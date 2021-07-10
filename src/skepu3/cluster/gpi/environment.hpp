#ifndef ENVIRONMENT_HPP
#define ENVIRONMENT_HPP

#include <GASPI.h>
#include <cassert>
#include <random>
#include <mutex>
#include <thread>
#include <omp.h>
#include <climits>
/*
* This singleton scheme allows for better control of global state.
* In particular it calls gaspi_init and terminate correctly.
*/
class Environment
{
private:

  ~Environment(){
    delete generator;

    // This barrier is only needed by mapreduce without any container
    gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);

    gaspi_proc_term(GASPI_BLOCK);
  };

  Environment() :
  init_called{false},
  mapred_seg_created{false},
  mapred_swapspace_ptr{nullptr}
  {};

  bool init_called;
  bool mapred_seg_created;
  std::mutex** rank_locks;


  gaspi_rank_t rank;
  gaspi_rank_t nr_nodes;
  inline static gaspi_pointer_t vclock_void;

  gaspi_pointer_t mapred_swapspace_ptr;
  gaspi_queue_id_t env_queue;

  inline static std::mt19937* generator;




  // Reads the remote vclock and updates our own
  // Take in the vclock offsets since we might read from a container with
  // a different type or size (and hence offset) than ourselves.
  void get_vclock(int dest_rank){
      if(dest_rank == get_rank()){
        // Reading our own vclock does nothing
        return;
      }

    auto res = gaspi_read_notify(
      0, // local seg
      sizeof(unsigned long) * nr_nodes, // local offset
      dest_rank,
      0, // dest seg
      0, // remote offset
      sizeof(unsigned long) * nr_nodes,
      100, // notif id
      env_queue,
      GASPI_BLOCK
    );

    if(res == GASPI_QUEUE_FULL){
      gaspi_wait(env_queue, GASPI_BLOCK);

      gaspi_read_notify(
        0, // local seg
        sizeof(unsigned long) * nr_nodes, // local offset
        dest_rank,
        0, // dest seg
        0, // remote offsett
        sizeof(unsigned long) * nr_nodes,
        100, // notif id
        env_queue,
        GASPI_BLOCK
      );
    }


    gaspi_notification_id_t notify_id;
    gaspi_notification_t notify_val;


    gaspi_notify_waitsome(
      0, // seg id
      100, //notif id
      1,
      &notify_id,
      GASPI_BLOCK
      );

    gaspi_notify_reset(0, notify_id, &notify_val);


    for(int i = 0; i < nr_nodes; i++){
      if(i == rank){
        vclock[i] = Environment::op_nr;
      }
      else{
        vclock[i] = std::max(vclock[i + nr_nodes], vclock[i]);
      }
    }
  };







public:

  static Environment& get_instance(){
    static Environment ins;
    return ins;
  }


  inline static unsigned long* vclock;

  inline static long unsigned op_nr = 0;

  void init(){
    if(!init_called){
      gaspi_proc_init(GASPI_BLOCK);
      init_called = true;

      gaspi_proc_num(&nr_nodes);
      gaspi_proc_rank(&rank);
      rank_locks = (std::mutex**) malloc(sizeof(std::mutex*) * nr_nodes);

      for(int i  = 0; i < nr_nodes; i++){
        rank_locks[i] = new std::mutex{};
      }

      assert(gaspi_segment_create(
        0,
        gaspi_size_t{2 * nr_nodes * sizeof(unsigned long)},
        GASPI_GROUP_ALL,
        GASPI_BLOCK,
        GASPI_MEM_INITIALIZED
      ) == GASPI_SUCCESS);


      gaspi_segment_ptr(0, &(Environment::vclock_void));
      Environment::vclock = (unsigned long*) Environment::vclock_void;

    }
   };

   static std::mt19937& get_generator(){
     if(generator == nullptr){
       generator = new std::mt19937(std::random_device{}() );
     }
     return *generator;
   }

   static int get_rank(){
     Environment& env = Environment::get_instance();
     if(!env.init_called){
       env.init();
     }

     return env.rank;
   }

   static int get_nr_nodes(){
     Environment& env = Environment::get_instance();
     if(!env.init_called){
       env.init();
     }

     return env.nr_nodes;
   }

   static long unsigned* get_vclock(){
     Environment& env = Environment::get_instance();
     if(!env.init_called){
       env.init();
     }

     return (long unsigned*) env.vclock_void;

   }


   /* This is a thread safe function which waits until dest_rank has reached
   * atleast operation number wait_for_op. The remaining arguments are used to
   * provide a GASPI segment point where the remote vclock can be stored which
   * allows for multiple concurrent reads.
   *
   * Note:
   * 1: This function prevents multiple vclock reads to the same rank. Instead
   *    one thread will do it while the rest wait for the result.
   *
   * 2: This function assumes that the read and writes operation to the long
   *    integer type are atomic. This assumption means that we do not need a
   *    read lock.
   *
   * 3: The reason that this function exists here and not in a container is that
   *    the vclock is shared between all containers and hence must be synchronized
   *    between all of them.
   */
   void wait_for_remote(
     int dest_rank,
     long unsigned wait_for_op,
     long unsigned* buffer,
     long unsigned seg_id,
     long unsigned loc_offset,
     gaspi_queue_id_t& queue
     ){

    long unsigned* vclock = (long unsigned*) vclock_void;

     bool is_done = vclock[dest_rank] >= wait_for_op;

     if(is_done)
      return;

     // Lock to prevent other threads from reading the same remote rank's vclock
     rank_locks[dest_rank]->lock();
     is_done = vclock[dest_rank] >= wait_for_op;

     // Check if the work has been done while we were sleeping
     if(is_done)
      return;

      // Fetch the remote vclock and store it in the given place
      auto res = gaspi_read_notify(
        seg_id, // local seg
        loc_offset, // local offset
        dest_rank,
        0, // dest seg
        0, // remote offset
        sizeof(unsigned long) * nr_nodes,
        omp_get_thread_num(), // notif id
        queue,
        GASPI_BLOCK
      );

      if(res == GASPI_QUEUE_FULL){
        gaspi_wait(queue, GASPI_BLOCK);

        gaspi_read_notify(
          seg_id, // local seg
          loc_offset, // local offset
          dest_rank,
          0, // dest seg
          0, // remote offset
          sizeof(unsigned long) * nr_nodes,
          omp_get_thread_num(), // notif id
          queue,
          GASPI_BLOCK
        );
      }

      gaspi_notification_id_t notify_id;
      gaspi_notification_t notify_val;

      gaspi_notify_waitsome(
        seg_id, // seg id
        omp_get_thread_num(),  //notif id
        1,
        &notify_id,
        GASPI_BLOCK
        );

      gaspi_notify_reset(seg_id, notify_id, &notify_val);

      // Lock to prevent concurrent writes to the vclock
      rank_locks[rank]->lock();

      for(int i = 0; i < nr_nodes; i++){
        vclock[i] = std::max(vclock[i], buffer[i]);

        if(i == dest_rank)
           is_done = vclock[i] >= wait_for_op;
      }

      rank_locks[rank]->unlock();
      rank_locks[dest_rank]->unlock();

      if(is_done){
        return;
      }
      else{
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        wait_for_remote(dest_rank, wait_for_op, buffer,  seg_id,  loc_offset, queue);
      }
   }


   // MapReduce may be called without a container existing hence the data must
   // be stored globally somewhere.
   template<typename T>
   static void create_mapred_swapspace(T){
     Environment& env = Environment::get_instance();


     if(!env.init_called){
       env.init();
     }

     if(!env.mapred_seg_created){

       assert(gaspi_segment_create(
         UCHAR_MAX - 1, // Highest allowed segment ID
         gaspi_size_t{sizeof(T) * env.nr_nodes},
         GASPI_GROUP_ALL,
         GASPI_BLOCK,
         GASPI_MEM_INITIALIZED
       ) == GASPI_SUCCESS);

       gaspi_segment_ptr(UCHAR_MAX - 1, &(env.mapred_swapspace_ptr));
       gaspi_queue_create(&env.env_queue, GASPI_BLOCK);
       env.mapred_seg_created = true;

     }

   }

   static gaspi_pointer_t get_mapred_swapspace(){
     Environment& env = Environment::get_instance();
     return env.mapred_swapspace_ptr;
   }

   static gaspi_queue_id_t get_queue(){
     Environment& env = Environment::get_instance();
     return env.env_queue;
   }


   static void wait_for_vclocks(unsigned long wait_val, int dest_rank){
     Environment& env = Environment::get_instance();

       while(true){
         env.get_vclock(dest_rank);
         if(env.vclock[dest_rank] >= wait_val){
           break;
         }
         else{
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
         }
       }
     }
};


#endif//ENVIRONMENT_HPP
