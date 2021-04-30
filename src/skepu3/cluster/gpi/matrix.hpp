#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cassert>
#include <vector>
#include <array>
#include <type_traits>
#include <chrono>
#include <thread>
#include <cmath>
#include <cstring>
#include <atomic>
#include <climits>
#include <omp.h>

#include <GASPI.h>

#include "utils.hpp"
#include "container.hpp"

namespace skepu{
  namespace _gpi{
    struct build_tup_util;


    template<int, typename Tup, typename T>
    struct build_buff_helper;
  }

  template<typename T>
  class Matrix : public skepu::_gpi::Container{


    // All skeletons must be added as friend classes here
    // and in all other containers. This is an ad hoc solution due to being a demo
    template<typename TT>
    friend class Reduce1D;
    template<int, typename Ret, typename... Func_args>
    friend class Map1D;
    template<typename TT>
    friend class FilterClass;
    template<typename TT>
    friend class Matrix;


    template<typename TT>
    friend class Vec;

    friend class _gpi::build_tup_util;


    template<int, typename TTup, typename TT>
    friend class _gpi::build_buff_helper;

  private:
    // Used to store two versions of the data so that Map may use a random access
    // type within its function.
    T* local_buffer;

    int local_size;
    long global_size;

    unsigned long comm_offset;
    unsigned long comm_size;

    // Global information about the container
    int last_partition_size;
    long last_partition_comm_offset;

    int norm_partition_size;
    long norm_partition_comm_offset;


    unsigned long last_mod_op;
    unsigned long* last_flush;


    /* TODO Remove!
    // Global information regarding offset
    long last_partition_vclock_offset;
    long norm_vclock_offset;

    // This object's offset
    long vclock_offset;
    */
    // The OP number of the remote when its container was fetched
    std::atomic_ulong* comm_buffer_state;
    std::mutex* comm_buffer_locks;

    // Indiciates which indeces this partition handles
    long start_i;
    long end_i;

    int get_owner(unsigned long index){
      if(index >= global_size){
        return -1;
      }
      return std::min((int) std::floor(index / (double) norm_partition_size), nr_nodes - 1);
    }



    // Puts all of dest_cont's elements within the range in our
    //communication buffer starting at offset local_offset.
    // Inclusive range, [start, end]
    void read_range(
      int start,
      int end,
      unsigned long local_offset,
      Matrix& dest_cont,
      bool no_wait = false
      ){
        int lowest_rank = dest_cont.get_owner(start);
        int highest_rank = dest_cont.get_owner(end);

        int curr_seg_id;
        int ranks_last_elem;
        int ranks_first_elem;

        int nr_elems_to_send;
        int sent_elems = 0;

        for(int i = lowest_rank; i <= highest_rank; i++){
          curr_seg_id = dest_cont.segment_id - dest_cont.rank + i;

          ranks_last_elem = i == nr_nodes - 1 ?
            dest_cont.global_size - 1
            : (i + 1) * dest_cont.norm_partition_size - 1;

          ranks_last_elem = std::min(ranks_last_elem, end);
          ranks_first_elem = std::max(i * dest_cont.norm_partition_size, start);

          nr_elems_to_send = ranks_last_elem - ranks_first_elem + 1;

          gaspi_read(
            segment_id,
            comm_offset + local_offset + sizeof(T) * sent_elems,
            i,
            dest_cont.segment_id - rank + i,
            sizeof(T) * (ranks_first_elem % dest_cont.norm_partition_size), // remote offset
            sizeof(T) * nr_elems_to_send, // size
            dest_cont.queue,
            GASPI_BLOCK
          );
          sent_elems += nr_elems_to_send;
        }
        if(!no_wait){
          gaspi_wait(dest_cont.queue, GASPI_BLOCK);
        }
    }


    // Reads the remote vclock and updates our own
    // Take in the vclock offsets since we might read from a container with
    // a different type or size (and hence offset) than ourselves.
    void get_vclock(int dest_rank){

        if(dest_rank == rank){
          // Reading our own vclock does nothing
          return;
        }


      auto res = gaspi_read_notify(
        0, // local seg
        sizeof(unsigned long) * nr_nodes, // local offset
        dest_rank,
        0, // dest seg
        0, // remote offsett
        sizeof(unsigned long) * nr_nodes,
        100, // notif id
        queue,
        GASPI_BLOCK
      );

      if(res == GASPI_QUEUE_FULL){
        gaspi_wait(queue, GASPI_BLOCK);

        gaspi_read_notify(
          0, // local seg
          sizeof(unsigned long) * nr_nodes, // local offset
          dest_rank,
          0, // dest seg
          0, // remote offsett
          sizeof(unsigned long) * nr_nodes,
          100, // notif id
          queue,
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
          vclock[i] = op_nr;
        }
        else{
          vclock[i] = std::max(vclock[i + nr_nodes], vclock[i]);
        }
      }
    };


    void wait_for_vclocks(unsigned long wait_val){
      int curr_rank;
      bool done;

      for(int i = 0; i < wait_ranks.size(); i++){
        curr_rank = wait_ranks[i];

        if(curr_rank == rank || vclock[curr_rank] >= wait_val){
          continue;
        }


        while(true){
          get_vclock(curr_rank);
          done = vclock[curr_rank] >= wait_val;
          if(done){
            break;
          }
          else{
            // TODO Remove
            // std::cout << "Rank " << rank << " waiting on ranks: ";
            //
            // for(int i = 0; i < wait_ranks.size(); i++){
            //   std::cout << wait_ranks[i] << ", ";
            // }
            // std::cout << "  and value " << wait_val << std::endl;


            // Sleep TODO change val
            std::this_thread::sleep_for(std::chrono::milliseconds(5));

          }
        }
      }
      wait_ranks.clear();
    };


    /*
    * Gets the harshers constraints in the arg list and stores it in dest.
    * Clears the constraint of non dest members.
    */
    template<typename Curr, typename... Rest>
    void get_constraints(Curr& curr,
      Rest&... rest){
        get_constraints_helper(curr, curr == *this);
        get_constraints(rest...);
      }

    // Sink
    void get_constraints(){
      }

    template<typename Cont>
    void get_constraints_helper(Cont& curr, bool is_same){

      for(auto& c : curr.constraints){

        if(c.second > constraints[c.first]){
          constraints[c.first] = c.second;
        }
      }
      if(!is_same){
        curr.constraints.clear();
      }
    }


    void wait_for_constraints(){

      //if(op_nr == 2)
        //printf("Rank %d, op %lu waiting for: ", rank, op_nr);


      for(auto& c : constraints){

        //if(op_nr == 2){
          //printf("(%d, %lu ), ", c.first, c.second);
        //}

        wait_ranks.push_back(c.first);
        wait_for_vclocks(c.second);

      }

      //if(op_nr == 2)
        //printf("\n");

      constraints.clear();
    }



    // Flushes every container which has not flushed since the last modifying
    // operation was performed.
    template<typename Curr, typename... Rest>
    static void flush_rest(Curr& curr, Rest&... rest){

      curr.conditional_flush();
      flush_rest(rest...);
    }

    static void flush_rest(){}


    // Only flush if there are new changes
    void conditional_flush(){

      //printf("Attempting to flush: op_nr: %lu, last_flush %lu, mod_op %lu\n",
      //op_nr, last_flush[rank], last_mod_op);

      if(last_flush[rank] <= last_mod_op){
        flush();
      }
    }

    void flush(){

      bool empty_buff = last_flush[rank] == 0 && last_mod_op == 0;

      if(empty_buff || last_flush[rank] == op_nr){
        return;
      }

      //printf("Did flush: op_nr: %lu, last_flush %lu, mod_op %lu\n",
      //op_nr, last_flush[rank], last_mod_op);

      //printf("S %d, E %d, LS %d\n", start_i, end_i, local_size);

      //for(int i = 0; i == 10000; i++){
      //}

      /*
      std::memcpy(
        cont_seg_ptr,
        local_buffer,
        sizeof(T) * local_size);
        return;
        */

      for(int i = 0; i < nr_nodes; i++){
        last_flush[i] = op_nr;
      }




      #pragma omp parallel
      {
        size_t t_step = local_size / omp_get_num_threads();

        size_t size = omp_get_thread_num() == omp_get_num_threads() -1 ?
          t_step + local_size % omp_get_num_threads() :
          t_step;

        std::memcpy(
          (value_type*) cont_seg_ptr + omp_get_thread_num() * t_step,
          (value_type*) local_buffer + omp_get_thread_num() * t_step,
           sizeof(T) * size);
      }
    }



    /* Fetches all indeces between our own start_i to end_i from the remote
    * container cont. The values are put cont's communication buffer.
    *
    * The no_wait argument indicates whether we wait for the read to finish
    * or not. We may still have to wait for the remote ranks to reach our
    * current state.
    */
    template<typename Cont>
    void build_buffer_helper(bool no_wait, Cont& cont){
      int transfered_obj = 0;

      unsigned swap_offset = sizeof(T) * cont.rank * cont.norm_partition_size;

      if(cont.start_i > start_i){
        // transfer up to our start_i
        transfered_obj = cont.start_i - start_i;

        int first_owner = cont.get_owner(start_i);
        wait_ranks.clear();

        for(int i = first_owner; i < cont.rank; i++){
          wait_ranks.push_back(i);
        }
        wait_for_vclocks(op_nr);

        // read range is inclusive
        cont.read_range(start_i, cont.start_i - 1, swap_offset, cont, no_wait);
      }

      if(cont.end_i < end_i){
        // transfer up to our end_i
        int last_owner = cont.get_owner(end_i);
        wait_ranks.clear();

        for(int i = cont.rank + 1; i < last_owner; i++){
          wait_ranks.push_back(i);
        }
        wait_for_vclocks(op_nr);

        cont.read_range(cont.end_i + 1, end_i, swap_offset +
            transfered_obj * sizeof(T), cont, no_wait);

      }
    }

    template<typename First, typename ... Rest>
    static long smallest(First& first, Rest&... rest){
      return smallest(first.global_size, rest...);
    }

    template<typename First, typename ... Rest>
    static long smallest(long i, First& first, Rest&... rest){
      return smallest(std::min(first.global_size, i), rest...);
    }

    template<typename First>
    static long smallest(long i, First& first){
      return std::min(first.global_size, i);
    }

    template<typename First, typename ... Rest>
    static long lowest_shared_i(First& first, Rest&... rest){
      return lowest_shared_i(first.start_i, rest...);
    }

    template<typename First, typename ... Rest>
    static long lowest_shared_i(long i, First& first, Rest&... rest){
      return lowest_shared_i(std::max(first.start_i, i), rest...);
    }

    static long lowest_shared_i(long i){
      return i;
    }


    template<typename First, typename ... Rest>
    static long highest_shared_i(First& first, Rest&... rest){
      return highest_shared_i(first.end_i, rest...);
    }

    template<typename First, typename ... Rest>
    static long highest_shared_i(long i, First& first, Rest&... rest){
      return highest_shared_i(std::min(first.end_i, i), rest...);
    }

    static long highest_shared_i(long i){
      return i;
    }




    /* TODO Remove these functions!

    // Given a list of Matrices and other types returns the highest OP number
    // of the matrices.
    template<typename ... Rest>
    static unsigned long max_op(Rest&... rest){
      return max_op(int{0}, 0, rest...);
    }


    template<typename First, typename ... Rest>
    static unsigned long max_op(int SFINAE, unsigned long acc, Matrix<First>& first,
      Rest&... rest){
      return max_op(int{}, std::max(acc, first.op_nr), rest...);
    }

    template<typename First, typename ... Rest>
    static unsigned long max_op(double SFINAE, unsigned long acc, First&& first,
      Rest&... rest){
      return max_op(int{}, acc, rest...);
    }


    template<typename Last>
    static unsigned long max_op(int SFINAE, unsigned long acc, Matrix<Last>& last){
      return std::max(acc, last.op_nr);
    }

    template<typename Last>
    static unsigned long max_op(double SFINAE, unsigned long acc, Last&& last){
      return acc;
    }



  template<typename ... Rest>
  static void set_op_nr(unsigned long val, Rest&... rest){
    set_op_nr_helper(int{}, val, rest...);
  }

  template<typename First, typename ... Rest>
  static void set_op_nr_helper(int SFINAE, unsigned long val, Matrix<First>& first,
     Rest&... rest){
    first.op_nr = val;
    first.vclock[first.rank] = val;
    set_op_nr_helper(int{}, val, rest...);
  }


  template<typename First, typename ... Rest>
  static void set_op_nr_helper(double SFINAE, unsigned long val,
    First&& first, Rest&... rest){
    set_op_nr_helper(int{}, val, rest...);
  }

  template<typename Last>
  static void set_op_nr_helper(int SFINAE, unsigned long val,
      Matrix<Last>& last){
    last.op_nr = val;
    last.vclock[last.rank] = val;
  }

  template<typename Last>
  static void set_op_nr_helper(double SFINAE, unsigned long val,
      Last& last){
  }

  */

  long get_comm_offset(size_t inc_rank){
    if(inc_rank != nr_nodes -1){
      return norm_partition_comm_offset;
    }
    return last_partition_comm_offset;
  }


  size_t get_end(size_t inc_rank){
    if(inc_rank == nr_nodes - 1){
      return global_size - 1;
    }
    else{
      return (inc_rank + 1) * norm_partition_size - 1;
    }
  }


  size_t get_start(size_t inc_rank){
    return inc_rank * norm_partition_size;
  }


  /* Fetches a potentially remote value. If needed it fills the communication
  * buffer with a remote nodes values. The function is thread safe and assumes
  * a not changing op number during the multithreaded runs of it.
  */
  T proxy_get(const size_t i){
    T* comm_buffer = ((T*) comm_seg_ptr );

    if(i >= start_i && i <= end_i){



      // If we have flushed or if the container has never been modified before
      // we read from the container rather than buffer
      if(last_flush[rank] >= last_mod_op || last_flush[rank] == 0){
        return ((T*) cont_seg_ptr)[i - start_i];
      }
      else{
        return ((T*) local_buffer)[i - start_i];
      }
    }

    int remote_rank = get_owner(i);
    unsigned long remote_size = remote_rank == nr_nodes - 1 ?
      last_partition_size : norm_partition_size;


    if(comm_buffer_state[remote_rank] == op_nr){
      return comm_buffer[i];
    }


    comm_buffer_locks[remote_rank].lock();

    // Check if the work has been done while we were waiting on the lock
    if(comm_buffer_state[remote_rank] == op_nr){
      comm_buffer_locks[remote_rank].unlock();
      //return -100;
      return comm_buffer[i];
    }

    // We are the first to get the lock and need to fetch the remote values
    else{

      // TODO Wait for the remote rank to reach our current OP number
      comm_buffer_locks[rank].lock();
      wait_ranks.push_back(remote_rank);

      // if(last_flush[remote_rank] > last_mod_op){
      //   // Remote has flushed before this operation and we have a weaker restriction
      //
      //   //wait_for_vclocks(op_nr);
      //
      //   if(op_nr == 2)
      //   printf("Weak wait: op_nr: %lu, last_flush %lu, mod_op %lu\n",
      //   op_nr, last_flush[rank], last_mod_op);
      //
      //
      //   wait_for_vclocks(last_flush[remote_rank]);
      //
      // }
      // else{
      //   if(op_nr == 2)
      //   printf("Strong wait: op_nr: %lu, last_flush %lu, mod_op %lu\n",
      //   op_nr, last_flush[rank], last_mod_op);
      //
      //   // The remote has not flushed early and will have to flush during this op
      //   wait_for_vclocks(op_nr);
      // }

        if(op_nr == 34)
        printf("Wait: op_nr: %lu, last_flush %lu, mod_op %lu               ",
          op_nr, last_flush[rank], last_mod_op);

      wait_for_vclocks(last_flush[remote_rank]);

      comm_buffer_locks[rank].unlock();

      if(op_nr == 34)
        print_vclock();


      unsigned long read_size = remote_rank != nr_nodes - 1 ?
        (norm_partition_size ) * sizeof(T) :
        (last_partition_size ) * sizeof(T);


        if(op_nr == 244){
          printf("Seg %d reading seg %d\n", segment_id, segment_id + remote_rank - rank);
        }


      // A read notify might be more fitting here
      auto res = gaspi_read(
        segment_id,
        comm_offset + sizeof(T) * remote_rank * norm_partition_size, // local offset
        remote_rank,
        segment_id + remote_rank - rank, // dest_seg
        0, // remote offset
        read_size,
        queue,
        GASPI_BLOCK
      );


      // printf("Rank %d reading from rank %d for %lu elems. Local_size = %d"
      // " Last Size = %d\n", rank, remote_rank, read_size / sizeof(T),
      //   local_size, last_partition_size);


      if(res == GASPI_QUEUE_FULL){
        gaspi_wait(queue, GASPI_BLOCK);
        gaspi_read(
          segment_id,
          comm_offset + sizeof(T) * remote_rank * norm_partition_size, // local offset
          remote_rank,
          segment_id + remote_rank - rank, // dest_seg
          0, // remote offset
          read_size,
          queue,
          GASPI_BLOCK
        );
      }

      gaspi_wait(queue, GASPI_BLOCK);
      comm_buffer_state[remote_rank] = op_nr;

      comm_buffer_locks[remote_rank].unlock();


      bool b = remote_rank == nr_nodes - 1;
      int it = b ? last_partition_size : norm_partition_size;
      //printf("it = %d\n", it);

      if(op_nr == 244 && rank == 2){

        int s = norm_partition_size * remote_rank;
        int e = norm_partition_size * remote_rank + it -1;

        printf("index (%d, %d) remote: %d, flush %lu: ", s, e, remote_rank,
        last_flush[remote_rank]);
      for(int j = 0; j < it; j++){


        std::cout << comm_buffer[norm_partition_size * remote_rank + j] << ", ";
      }
      std::cout << "\n";
      print_vclock();
    }


      //comm_seg_ptr = ((T*) cont_seg_ptr) + local_size;
      //T* comm_buffer = ((T*) comm_seg_ptr );

      if(op_nr == 24)
        printf("At index %lu (pos 0) is val %d\n", i, comm_buffer[i]);

      //return 100;
      return comm_buffer[i];

    }
  }


  public:

    using value_type = T;
    using is_skepu_container = std::true_type;

    bool operator==(Matrix<T>& that){
      return this == &that;
    }


    Matrix(const Matrix&) = delete;
    Matrix& operator=(const Matrix&) = delete;


    Matrix(int rows, int cols){
      // Partition the matrix so that each rank receivs an even
      // amount of elements
      int step = (rows * cols) / nr_nodes;
      // The part which can not be evenly split
      int residual = (rows * cols) % nr_nodes;

      if(rank != nr_nodes - 1){
        start_i = rank * step;
        end_i = (rank + 1) * step - 1;
      }
      else{
        start_i = rank * step;
        end_i = (rank + 1) * step + residual - 1;

      }

      last_partition_size = step + residual;
      norm_partition_size = step;

      local_size = end_i - start_i + 1;
      global_size = rows * cols;


      comm_offset = sizeof(T) * local_size;
      comm_size = sizeof(T) * global_size;

      last_partition_comm_offset = sizeof(T) * last_partition_size;
      norm_partition_comm_offset = sizeof(T) * norm_partition_size;


      // Allocate enough memory to:
      // 1 - Store the local part of the container
      // 2 - Store remote values which we might read
      // 2.5 - Store remote values which we prefetch
      assert(gaspi_segment_create(
        segment_id,
        gaspi_size_t{sizeof(T) * (local_size + global_size)},
        GASPI_GROUP_ALL,
        GASPI_BLOCK,
        GASPI_MEM_INITIALIZED
      ) == GASPI_SUCCESS);


      gaspi_segment_ptr(segment_id, &cont_seg_ptr);
      comm_seg_ptr = ((T*) cont_seg_ptr) + local_size;

      local_buffer = (T*) malloc(local_size * sizeof(T));

      comm_buffer_state = (std::atomic_ulong*) malloc(sizeof(std::atomic_ulong) * nr_nodes);
      for(int i  = 0; i < nr_nodes; i++){
        comm_buffer_state[i] = ULONG_MAX;
      }

      comm_buffer_locks = new std::mutex[nr_nodes];


      last_mod_op = 0;
      last_flush = (unsigned long*) calloc(nr_nodes, sizeof(unsigned long));

    };


    Matrix(int rows, int cols, T init) : Matrix(rows, cols){
      set(init);
    }

    ~Matrix(){
      delete local_buffer;
      delete comm_buffer_state;
      delete comm_buffer_locks;
    }

    void set(T scalar){
      for(int i = 0; i < local_size; i ++){
        ((T*) cont_seg_ptr)[i] = scalar;
      }
    }


    // Randomly initializes an int or long matrix
    void rand(int from, int to){
      if(std::is_same<T, int>::value || std::is_same<T, long>::value){
        for(int i = 0; i < local_size; i ++){
          ((T*) cont_seg_ptr)[i] = from + std::rand() % (to - from);
        }

      }
      else{
        std::cout << "Rand only supports matrices of int or long\n";
      }
    }



    void set(int index, T value){
      if(index >= start_i && index <= end_i){
        ((T*) cont_seg_ptr)[index - start_i] = value;
      }
      //vclock[rank] = ++op_nr;
    }


    size_t size(){
      return global_size;
    }


    T operator[](const size_t index){

      if(index >= start_i && index <= end_i){
        // The value is local, wait until we have distributed the value

        wait_for_constraints();
        conditional_flush();
        vclock[rank] = ++op_nr;

        for(int i = 0; i < nr_nodes; i++){
          if(i == rank){
            continue;
          }
          constraints[i] = op_nr;
        }

        return ((T*) cont_seg_ptr)[index - start_i];
      }
      else{
        vclock[rank] = ++op_nr;

        // The rank holding the value
        int dest_rank = get_owner(index);
        wait_ranks.push_back(dest_rank);

        //long unsigned foo = std::min(last_flush[dest_rank])

        wait_for_vclocks(last_flush[dest_rank]);


        // Wait until the rank is on the current OP

        auto res = gaspi_read(
          segment_id,
          comm_offset,
          dest_rank,
          segment_id + dest_rank - rank, // remote segment id
          sizeof(T) * (index - norm_partition_size * dest_rank), // Remote offset
          sizeof(T),
          queue,
          GASPI_BLOCK
        );

        if(res == GASPI_QUEUE_FULL){
          gaspi_wait(queue, GASPI_BLOCK);

          gaspi_read(
            segment_id,
            comm_offset,
            dest_rank,
            segment_id + dest_rank - rank, // remote segment id
            sizeof(T) * (index - norm_partition_size * dest_rank), // Remote offset
            sizeof(T),
            queue,
            GASPI_BLOCK
          );
        }

        // Wait for read to finish
        gaspi_wait(queue, GASPI_BLOCK);

        // Notify the remote that we have read the value
        gaspi_notify(
          segment_id + dest_rank - rank, // remote segment
          dest_rank,
          notif_ctr * nr_nodes + rank, // notif id
          13, // notif val, not used
          queue,
          GASPI_BLOCK
        );

        ++notif_ctr;
        vclock[rank] = ++op_nr;
        return ((T*) comm_seg_ptr)[0];
      }
    }

    void print_vclock(){
      std::cout << "Rank " << rank << " has vclock: ";

      for(int i = 0; i < nr_nodes; i++){
        std::cout << vclock[i] << ", ";
      }

      std::cout << std::endl;
    }


    // WARNING Only use this for debugging, it has very poor performance
    void print(){
      gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);

      op_nr++;

      wait_for_constraints();
      //print_vclock();
      flush();

      for(int i = 0; i < nr_nodes; i++){
        if(i == rank){
          for(int j = 0; j < local_size; j++){
            std::cout << ((T*) cont_seg_ptr)[j] << ", ";
          }
          std::cout << std::endl;
        }
        gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);
      }

      // The sleep makes multiple prints prettier
      std::this_thread::sleep_for(std::chrono::milliseconds(500));
      if(rank == 0){
        printf("\n");
      }
    }

    void print_buff(){
      //gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);

      op_nr++;
      //flush();

      for(int i = 0; i < nr_nodes; i++){
        if(i == rank){
          for(int j = 0; j < local_size; j++){
            std::cout << ((T*) local_buffer)[j] << ", ";
          }
          std::cout << std::endl;
        }
        gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);
      }

      // The sleep makes multiple prints prettier
      std::this_thread::sleep_for(std::chrono::milliseconds(500));
      if(rank == 0){
        printf("\n");
      }
    }


    // TODO Remove
    void foobar(){
      //flush();

      T* ptr = (T*) cont_seg_ptr;

      for(int i = 0; i < local_size; i++){

        if(ptr[i] != local_buffer[i]){
          printf("cont: %d, buff %d\n", ptr[i], local_buffer[i]);
        }

      }

    }



    template<typename Lambda>
    void print(Lambda lamb){
      for(int i = 0; i < nr_nodes; i++){
        if(i == rank){
          for(int j = 0; j < local_size; j++){
            std::cout << lamb(((T*) cont_seg_ptr)[j]) << "\n";
          }
          std::cout << std::endl;
        }
        gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);
      }

      // The sleep makes multiple prints prettier
      std::this_thread::sleep_for(std::chrono::milliseconds(500));
      if(rank == 0){
        printf("\n");
      }
    }


};

  template<typename T>
  struct is_skepu_container<skepu::Matrix<T>> : std::true_type {};
}

#endif // MATRIX_HPP
