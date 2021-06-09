#ifndef ENVIRONMENT_HPP
#define ENVIRONMENT_HPP

#include <GASPI.h>
#include <cassert>
#include <random>
#include <mutex>
#include <thread>
#include <omp.h>
/*
* This singleton scheme allows for better control of global state.
* In particular it calls gaspi_init and terminate correctly.
*/
class Environment
{
private:

  ~Environment(){
    delete generator;
    gaspi_proc_term(GASPI_BLOCK);
  };

  Environment() : init_called{false}{};

  bool init_called;
  std::mutex* rank_locks;


  gaspi_rank_t rank;
  gaspi_rank_t nr_nodes;

  inline static gaspi_pointer_t vclock_void;


  inline static std::mt19937* generator;

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
      rank_locks = new std::mutex[nr_nodes];


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
     gaspi_queue_id_t queue
     ){

    long unsigned* vclock = (long unsigned*) vclock_void;

     bool is_done = vclock[dest_rank] >= wait_for_op;

     if(is_done)
      return;

     // Lock to prevent other threads from reading the same remote rank's vclock
     rank_locks[dest_rank].lock();
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
      rank_locks[rank].lock();

      for(int i = 0; i < nr_nodes; i++){
        vclock[i] = std::max(vclock[i], buffer[i]);

        if(i == dest_rank)
           is_done = vclock[i] >= wait_for_op;
      }

      rank_locks[rank].unlock();
      rank_locks[dest_rank].unlock();

      if(is_done){
        return;
      }
      else{
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        wait_for_remote(dest_rank, wait_for_op, buffer,  seg_id,  loc_offset, queue);
      }
   }


};



#endif//ENVIRONMENT_HPP
