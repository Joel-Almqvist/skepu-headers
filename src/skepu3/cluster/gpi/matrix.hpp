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
#include <random>

#include <GASPI.h>

#include "utils.hpp"
#include "container.hpp"

// NOTE This included is only used for get_avg_time() and should
// eventually be removed.
#include "reduce.hpp"


namespace skepu{
  namespace _gpi{
    struct build_tup_util;


    template<int, typename Tup, typename T>
    struct build_buff_helper;

    template<int nr_args, int ctr, typename Tup>
    struct build_buffer_util;

  }

  template<typename T>
  class Matrix : public skepu::_gpi::Container{


    // All skeletons must be added as friend classes here
    // and in all other containers. This is an ad hoc solution due to being a demo
    template<typename ret_t, typename arg_t>
    friend class Reduce1D;

    template<int, typename Ret, typename... Func_args>
    friend class Map1D;

    template<int Arity, typename red_f_t, typename map_ret_t, typename... map_args>
    friend class MapReduce1D;

    template<typename TT>
    friend class FilterClass;

    template<typename TT>
    friend class Matrix;

    template<typename TT>
    friend class Vec;

    template<typename TT>
    friend class Mat;

    friend class _gpi::build_tup_util;

    template<int nr_args, int ctr, typename Tup>
    friend class _gpi::build_buffer_util;


  private:
    // Used to store two versions of the data so that Map may use a random access
    // type within its function.
    T* local_buffer;

    unsigned long local_size;
    unsigned long global_size;

    unsigned long comm_offset;
    unsigned long comm_size;

    // Global information about the container
    unsigned long last_partition_size;
    unsigned long last_partition_comm_offset;

    unsigned long norm_partition_size;
    unsigned long norm_partition_comm_offset;

    unsigned long* last_mod_op;
    unsigned long* last_flush;

    // The OP number of the remote when its container was fetched
    std::atomic_ulong* comm_buffer_state;
    std::mutex* comm_buffer_locks;

    // Indiciates which indeces this partition handles
    unsigned long start_i;
    unsigned long end_i;

    size_t row;
    size_t col;

    int get_owner(unsigned long index){
      if(index >= global_size){
        return -1;
      }
      return std::min((int) std::floor(index / (double) norm_partition_size), nr_nodes - 1);
    }

    size_t get_col(unsigned long index){
      return index % row;
    }

    size_t get_row(unsigned long index){
      return index / row;
    }


    void set(int index, T value){
      if(index >= start_i && index <= end_i){
        ((T*) cont_seg_ptr)[index - start_i] = value;
      }
    }


    // Puts all of dest_cont's elements within the range in our
    //communication buffer starting at offset local_offset.
    // Inclusive range, [start, end]
    void read_range(
      unsigned long start,
      unsigned long end,
      unsigned long local_offset,
      Matrix& dest_cont,
      bool no_wait = false
      ){
        int lowest_rank = dest_cont.get_owner(start);
        int highest_rank = dest_cont.get_owner(end);

        int curr_seg_id;
        unsigned long ranks_last_elem;
        unsigned long ranks_first_elem;

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
    void get_constraints(int sfinae, Matrix<Curr>& curr,
      Rest&... rest){
        get_constraints_helper(curr, curr == *this);
        get_constraints(int{}, rest...);
      }

      template<typename Curr, typename... Rest>
      void get_constraints(long sfinae, Curr& curr,
        Rest&... rest){
          get_constraints(int{}, rest...);
        }

    // Sink
    void get_constraints(int){}

    void get_constraints(long){}

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

      for(auto& c : constraints){

        if(vclock[c.first] < c.second){
          wait_ranks.push_back(c.first);
          wait_for_vclocks(c.second);
        }
      }
      constraints.clear();
    }



    // Flushes every container which has not flushed since the last modifying
    // operation was performed.
    template<typename Curr, typename... Rest>
    static void flush_rest(int sfinae, Matrix<Curr>& curr, Rest&... rest){

      curr.conditional_flush();

      for(int i = 0; i < curr.nr_nodes; i++){
        if(i != curr.rank){
          curr.remc_flush(i);
        }
      }

      flush_rest(sfinae, rest...);
    }

    template<typename Curr, typename... Rest>
    static void flush_rest(long sfinae, Curr& curr, Rest&... rest){
      flush_rest(sfinae, rest...);
    }

    static void flush_rest(int){}

    static void flush_rest(long){}


    // Flushes any changes from the double buffer to the GASPI segment if
    // there exists any non-flushed changes.
    void conditional_flush(){

      bool empty_buff = last_flush[rank] == 0 && last_mod_op[rank] == 0;
      if(last_flush[rank] <= last_mod_op[rank] && !empty_buff){

          last_flush[rank] = op_nr + 1;
          wait_for_constraints();

          internal_flush();
      }
    }

    // Moves the data from the double buffer to the GASPI segment making it
    // visible for remote ranks.
    void internal_flush(){

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


    // To be called whenever a remote rank calls a conditional flush in order
    // to track the remote's flush status.
    void remc_flush(int remote){

      bool empty_buff = last_flush[remote] == 0 && last_mod_op[remote] == 0;

      if(last_flush[remote] <= last_mod_op[remote] && !empty_buff){
          last_flush[remote] = op_nr + 1;
        }
    }


    /* Takes in a rank, index, accumulator, decision flag and a variadic list
    * of matrices mixed with constants. Given this it either returns the highest
    * or lowest rank which contains this index in any of the containers.
    */
    template<typename Curr, typename... Rest>
    static int find_rank_overlap(
      int sfinae,
      int remote_rank,
      long unsigned i,
      int acc,
      bool find_min,
      Matrix<Curr>& curr,
      Rest&... rest){

        int new_acc;
        if(find_min){
          new_acc = std::min(acc, curr.get_owner(i));
        }
        else{
          new_acc = std::max(acc, curr.get_owner(i));
        }

        return find_rank_overlap(sfinae, remote_rank, i, new_acc, find_min, rest...);

    }


    template<typename Curr, typename... Rest>
    static int find_rank_overlap(
      long sfinae,
      int remote_rank,
      long unsigned i,
      int acc,
      bool find_min,
      Curr& curr,
      Rest&... rest){

         return find_rank_overlap(sfinae, remote_rank, i, acc, find_min, rest...);

    }

    // Sink
    static int find_rank_overlap(long sfinae, int remote_rank, long unsigned,
    int acc, bool){
      return acc;
    }
    static int find_rank_overlap(int sfinae, int remote_rank, long unsigned,
    int acc, bool){
      return acc;
    }


    /* Fetches all indeces between our own start_i to end_i from the remote
    * container cont. The values are put cont's communication buffer.
    *
    * The no_wait argument indicates whether we wait for the read to finish
    * or not. We may still have to wait for the remote ranks to reach our
    * current state.
    */
    template<typename Cont>
    void build_buffer_helper(Cont& cont){

      if(cont == *this)
        return;

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
        cont.read_range(start_i, cont.start_i - 1, swap_offset, cont, false);
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
            transfered_obj * sizeof(T), cont, false);

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


  // Multi-threaded fetch for a potentially remote value
  T proxy_get(const size_t i){
    T* comm_buffer = ((T*) comm_seg_ptr );

    if(i >= start_i && i <= end_i){

      if(last_flush[rank] >= last_mod_op[rank]){

        // If we have flushed or if the container has never been modified before
        // we read from the container rather than buffer
        return ((T*) cont_seg_ptr)[i - start_i];
      }
      else{

        return ((T*) local_buffer)[i - start_i];

      }
    }

    int remote_rank = get_owner(i);
    unsigned long remote_size = remote_rank == nr_nodes - 1 ?
      last_partition_size : norm_partition_size;

    // Check if the remote value is already cached
    if(
        comm_buffer_state[remote_rank] > last_mod_op[remote_rank]
        && comm_buffer_state[remote_rank] > last_flush[remote_rank]
        && comm_buffer_state[remote_rank] != ULONG_MAX
      ){
      return comm_buffer[i];
    }

    comm_buffer_locks[remote_rank].lock();

    // Check if the work has been done while we were waiting on the lock
    if(comm_buffer_state[remote_rank] == op_nr){
      comm_buffer_locks[remote_rank].unlock();
      return comm_buffer[i];
    }

    else{

      // Use the thread-safe wait provided by Environment
      Environment& env = Environment::get_instance();
      T* buffer = ((T*) comm_seg_ptr + norm_partition_size * remote_rank);

      env.wait_for_remote(
        remote_rank,
        last_flush[remote_rank],
        (long unsigned*) buffer,
        segment_id,
        sizeof(T) * (norm_partition_size * remote_rank + local_size),
        queue
      );

      unsigned long read_size = remote_rank != nr_nodes - 1 ?
        (norm_partition_size ) * sizeof(T) :
        (last_partition_size ) * sizeof(T);


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

      return comm_buffer[i];

    }
  }


  public:

    using value_type = T;
    using is_skepu_container = std::true_type;

    bool operator==(Matrix<T>& that){
      return this == &that;
    }

    template<typename TT>
    bool operator==(TT& that){
      return false;
    }


    Matrix(const Matrix&) = delete;
    Matrix& operator=(const Matrix&) = delete;


    Matrix(size_t rows, size_t cols) : row{rows}, col{cols}{
      // Partition the matrix so that each rank receivs an even
      // amount of elements
      size_t step = (rows * cols) / nr_nodes;
      // The part which can not be evenly split
      size_t residual = (rows * cols) % nr_nodes;

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


      last_mod_op = (unsigned long*) calloc(nr_nodes, sizeof(unsigned long));
      last_flush = (unsigned long*) calloc(nr_nodes, sizeof(unsigned long));

    };


    Matrix(size_t rows, size_t cols, T init) : Matrix(rows, cols){
      set(init);
    }

    // NOTE This constructor is added to mimic the non implemented Vector class.
    // When Vector is implemented this should be removed.
    Matrix(size_t size) : Matrix(size, 1){}

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


  private:

    // Helps deduce whether we can cast the container's value type T and as
    // such whether we can compile the randomize function or not.
    //
    // Only allows for randomizing of real and disccrete numbers
    template<typename TT, typename Cont, bool real, bool discrete>
    struct randomize_helper{
      template<typename Arg>
      static void exec(Cont&, Arg, Arg){}
    };


    template<typename TT, typename Cont>
    struct randomize_helper<TT, Cont, true, false>{
      template<typename Arg>
      static void exec(Cont& cont, Arg from, Arg to){

        cont.template randomize_helper_real<Arg, TT>(from, to);
      }
    };

    template<typename TT, typename Cont>
    struct randomize_helper<TT, Cont, false, true>{
      template<typename Arg>
      static void exec(Cont& cont, Arg from, Arg to){

        cont.template randomize_helper_discrete<Arg, TT>(from, to);
      }
    };



    // Perfoms the actual randomizing
    template <typename Arg, typename TT>
    void randomize_helper_real(Arg from, Arg to) {

      static_assert(std::is_same<T, TT>::value);

      wait_for_constraints();
      conditional_flush();
      last_mod_op[rank] = ++op_nr;

      std::uniform_real_distribution<TT> dis(from, to);

      std::mt19937& generator = Environment::get_generator();

      T* lb = (T*) local_buffer;
      for(int i = 0; i < local_size; i++){
        lb[i] = dis(generator);
      }

      op_nr++;
    }

    // Perfoms the actual randomizing
    template <typename Arg, typename TT>
    void randomize_helper_discrete(Arg from, Arg to) {

      static_assert(std::is_same<T, TT>::value);

      wait_for_constraints();
      conditional_flush();
      last_mod_op[rank] = ++op_nr;


      std::uniform_int_distribution<TT> dis(from, to);

      std::mt19937& generator = Environment::get_generator();

      T* lb = (T*) local_buffer;
      for(int i = 0; i < local_size; i++){
        lb[i] = dis(generator);
      }

      op_nr++;
    }


  public:

    // Wrapper for randomize_helper_func
    // Capture the arguments with template parameter as they may be real or
    // discrete.
    template<typename Arg>
    void randomize(Arg from, Arg to){

      const bool is_real =
        std::is_same<T, float>::value ||
        std::is_same<T, double>::value ||
        std::is_same<T, long double>::value;

      const bool is_discrete =
        std::is_same<T, short>::value ||
        std::is_same<T, int>::value ||
        std::is_same<T, long >::value ||
        std::is_same<T, long long>::value ||
        std::is_same<T, unsigned short>::value ||
        std::is_same<T, unsigned int>::value ||
        std::is_same<T, unsigned long>::value ||
        std::is_same<T, unsigned long long>::value;

      randomize_helper<T, decltype(*this), is_real, is_discrete>
      ::exec(*this, from, to);
    }


    size_t size(){
      return global_size;
    }

    size_t size_i(){
      return row;
    }

    size_t size_j(){
      return col;
    }


    T operator[](const size_t index){

      // Local value
      if(index >= start_i && index <= end_i){

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

      // The value is remote

      int remote_rank = get_owner(index);
      remc_flush(remote_rank);

      T* comm_buffer = ((T*) comm_seg_ptr );
      vclock[rank] = ++op_nr;

      // Check if we have the value cached
      if(comm_buffer_state[remote_rank] <= last_mod_op[rank]
        || comm_buffer_state[remote_rank] == ULONG_MAX){

          wait_ranks.push_back(remote_rank);
          wait_for_vclocks(last_flush[remote_rank]);


          unsigned long remote_size = remote_rank == nr_nodes - 1 ?
          last_partition_size * sizeof(T): norm_partition_size * sizeof(T);

          auto res = gaspi_read(
            segment_id,
            comm_offset + norm_partition_size * remote_rank * sizeof(T),
            remote_rank,
            segment_id + remote_rank - rank, // remote segment id
            0, // Remote offset
            remote_size,
            queue,
            GASPI_BLOCK
          );

          if(res == GASPI_QUEUE_FULL){
            gaspi_wait(queue, GASPI_BLOCK);
            gaspi_read(
              segment_id,
              comm_offset + norm_partition_size * remote_rank * sizeof(T),
              remote_rank,
              segment_id + remote_rank - rank, // remote segment id
              0, // Remote offset
              remote_size,
              queue,
              GASPI_BLOCK
            );
          }
          gaspi_wait(queue, GASPI_BLOCK);
          comm_buffer_state[remote_rank] = op_nr;
        }

        return comm_buffer[index];
    }

    void print_vclock(){
      std::cout << "Rank " << rank << " has vclock: ";

      for(int i = 0; i < nr_nodes; i++){
        std::cout << vclock[i] << ", ";
      }

      std::cout << std::endl;
    }

    void print_flush(){
      std::cout << "Rank " << rank << " has last_flush: ";

      for(int i = 0; i < nr_nodes; i++){
        std::cout << last_flush[i] << ", ";
      }

      std::cout << std::endl;
    }


    void print_modop(){
      std::cout << "Rank " << rank << " has mod_op: ";

      for(int i = 0; i < nr_nodes; i++){
        std::cout << last_mod_op[i] << ", ";
      }

      std::cout << std::endl;
    }



    // WARNING Only use this for debugging, it has very poor performance
    void print(){

      op_nr++;
      wait_for_constraints();
      conditional_flush();
      for(int i = 0; i < nr_nodes; i++){
        if(i != rank)
        remc_flush(i);
      }


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

      op_nr++;

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


    template<typename Lambda>
    void print(Lambda lamb){

      op_nr++;
      wait_for_constraints();
      conditional_flush();

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


    bool single(){
      return rank == 0;
    }


    int get_rank(){
      return rank;
    }

    int get_nodes(){
      return nr_nodes;
    }

    // NOTE This function is only used for benchmarking and is itself very slow
    double get_avg_time(double delta){

      skepu::Matrix<double> helper(nr_nodes);

      helper.set(rank, delta / nr_nodes);

      auto plus = skepu::Reduce([](double a, double b){
        return a + b;
      });

      return plus(helper);

    }


    // NOTE This function is only used for benchmarking and is itself very slow
    double get_slowest_node(double delta){

      skepu::Matrix<double> helper(nr_nodes);

      helper.set(rank, delta);

      auto max = skepu::Reduce([](double a, double b){
        return a > b ? a : b;
      });


      return max(helper);

    }


    // ******************************************************************
    // All the function below are dummy functions for functionality which
    // has not been implemented but needed to compile.
    // ******************************************************************

    T* getAddress(){
      internal_flush();
      return (T*) cont_seg_ptr;
    }

    T* begin(){
      return (T*) cont_seg_ptr;
    }

    T* end(){
      return (T*) cont_seg_ptr + local_size - 1;
    }

    // Dummy function, see internal_flush for the real implementation
    void flush(){};

  };

  // Dummy function
  template<typename T>
  std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix){
    return os;
  }

  template<typename T>
  struct is_skepu_container<skepu::Matrix<T>> : std::true_type {};
}

#endif // MATRIX_HPP
