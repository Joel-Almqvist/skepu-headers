#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cassert>
#include <vector>
#include <array>
#include <iostream>
#include <GASPI.h>
#include <type_traits>
#include <chrono>
#include <thread>
#include <cmath>

#include <utils.hpp>
#include <container.hpp>
#include <reduce.hpp>
// TODO remove include iostream (among more?)

namespace skepu{

  template<typename T>
  class Matrix : public skepu::_gpi::Container{


    // All skeletons must be added as friend classes here
    // and in all other containers. This is an ad hoc solution due to being a demo
    template<typename TT>
    friend class Reduce1D;
    template<typename TT, int>
    friend class Map1D;
    template<typename TT>
    friend class FilterClass;
  private:

    using is_skepu_container = std::true_type;

    static const int COMM_BUFFER_NR_ELEMS = 50;

    int local_size;
    long global_size;

    // Global information about the container
    int last_partition_size;
    long last_partition_comm_offset;

    int norm_partition_size;
    long norm_partition_comm_offset;

    // Global information regarding offset
    long last_partition_vclock_offset;
    long norm_vclock_offset;

    // This object's offset
    long vclock_offset;

    long unsigned comm_size;

    // Indiciates which indeces this partition handles
    long start_i;
    long end_i;

    // TODO remove this, norm_partition_size == step
    int step;


    int get_owner(int index){
      return std::min((int) std::floor((float) index / step), nr_nodes - 1);
    }



    // Puts all of dest_cont's elements within the range in our
    //communication buffer starting at offset local_offset.
    // Inclusive range, [start, end]
    void read_range(
      int start,
      int end,
      int local_offset,
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
            : (i + 1) * dest_cont.step - 1;

          ranks_last_elem = std::min(ranks_last_elem, end);
          ranks_first_elem = std::max(i * dest_cont.step, start);

          nr_elems_to_send = ranks_last_elem - ranks_first_elem + 1;

          gaspi_read(
            segment_id,
            comm_offset + local_offset + sizeof(T) * sent_elems,
            i,
            dest_cont.segment_id - rank + i,
            sizeof(T) * (ranks_first_elem % dest_cont.step), // remote offset
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



    // Puts all elements from start to end (these are global indeces) into
    // the given GASPI segment. Many to one communication pattern
    //
    // WARNING This function is deprecated and will be removed eventually.
    // It has been replaced by read_range since due to being more intuitive.
    void get_range(
      int start,
      int end,
      int dest_rank,
      int dest_seg_id,
      int offset,
      int notify_id
      ) {

      int first_elem = -1;
      int last_elem = -1;

      if(end >= end_i && start <= start_i){
        // Middle part of the range
        first_elem = start_i;
        last_elem = end_i;
      }
      else if(start >= start_i && start <= end_i){
        // The start of the range
        first_elem = start;
        last_elem = end <= end_i ? end : end_i;
      }
      else if(end <= end_i && end >= start_i){
        // The end of the range
        first_elem = start >= start_i ? start : start_i;
        last_elem = end;
      }

      if(last_elem != -1){

       gaspi_write_notify(segment_id,
         sizeof(T) * (first_elem % local_size),
         dest_rank,
         dest_seg_id,
         offset + sizeof(T) * (first_elem - start), // offset remote
         sizeof(T) * (last_elem + 1 - first_elem), // size
         notify_id, // notification id
         rank + 1, // notification value, not used atm
         queue,
         GASPI_BLOCK);

      }
    }


    // Reads the remote vclock and updates our own
    void get_vclock(
      int dest_rank,
      int dest_seg_id
      ){

        if(dest_rank == rank && dest_seg_id == segment_id){
          // Reading our own vclock does nothing
          return;
        }

      unsigned long remote_offset = dest_rank == nr_nodes - 1 ?
        last_partition_vclock_offset :
        norm_vclock_offset;

      gaspi_read(
        segment_id, // local seg
        vclock_offset + sizeof(unsigned long) * nr_nodes, // local offset
        dest_rank,
        dest_seg_id,
        remote_offset,
        sizeof(unsigned long) * nr_nodes,
        queue,
        GASPI_BLOCK
      );

      gaspi_wait(queue, GASPI_BLOCK);

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
      int curr_seg_id;
      const int min_seg_id = segment_id - rank;

      for(int i = 0; i < wait_ranks.size(); i++){
        curr_rank = wait_ranks[i];
        curr_seg_id = min_seg_id + curr_rank;

        if(curr_rank == rank || vclock[curr_rank] >= wait_val){
          continue;
        }

        while(true){
          get_vclock(curr_rank, curr_seg_id);
          if(vclock[curr_rank] >= wait_val){
            break;
          }
          else{
            // Sleep
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
          }
        }
      }
    };


    template<typename First, typename ... Rest>
    void build_buffer(First& first, Rest&... rest){
      build_buffer(first);
      build_buffer(rest...);
    }

    template<typename Cont>
    void build_buffer(Cont& cont){

      int transfered_obj = 0;

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
        cont.read_range(start_i, cont.start_i - 1, 0, cont, true);
      }

      if(cont.end_i < end_i){
        // transfer up to our end_i
        int last_owner = cont.get_owner(end_i);
        wait_ranks.clear();

        for(int i = cont.rank + 1; i < last_owner; i++){
          wait_ranks.push_back(i);
        }
        wait_for_vclocks(op_nr);


        cont.read_range(cont.end_i + 1, end_i, transfered_obj * sizeof(T), cont, true);
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


    template<typename First, typename ... Rest>
    static unsigned long max_op(First& first, Rest&... rest){
      return max_op(first.op_nr, rest...);
    }


    template<typename First, typename ... Rest>
    static unsigned long max_op(unsigned long acc, First& first, Rest&... rest){
      return max_op(std::max(acc, first.op_nr), rest...);
    }

    template<typename Last>
    static unsigned long max_op(unsigned long acc, Last& last){
      return std::max(acc, last.op_nr);
    }

template<typename First, typename ... Rest>
static void set_op_nr(unsigned long val, First& first, Rest&... rest){
  first.op_nr = val;
  first.vclock[first.rank] = val;
  set_op_nr(val, rest...);
}

template<typename Last>
static void set_op_nr(unsigned long val, Last& last){
  last.op_nr = val;
  last.vclock[last.rank] = val;
}

  public:

    using value_type = T;

    bool operator==(Matrix<T>& that){
      return this == &that;
    }




    Matrix(){
      std::cout << "Empty constructor called\n";
    }

    Matrix(int rows, int cols){
      // Partition the matrix so that each rank receivs an even
      // amount of elements
      step = (rows * cols) / nr_nodes;
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
      last_partition_comm_offset = sizeof(T) * last_partition_size;
      norm_partition_comm_offset = sizeof(T) * norm_partition_size;


      // Guarantee that the comm buffer has size enough for Reduce to work
      assert(COMM_BUFFER_NR_ELEMS >= ((int) std::ceil(std::log2(nr_nodes))) + 1);

      comm_size = sizeof(T) * COMM_BUFFER_NR_ELEMS;


      norm_vclock_offset = sizeof(T) * norm_partition_size + comm_size;
      last_partition_vclock_offset = sizeof(T) * last_partition_size + comm_size;


      vclock_offset = comm_offset + comm_size;

      // 2 * nr_nodes * sizeof(unsigned long) is the size of the vector clock
      assert(gaspi_segment_create(
        segment_id,
        gaspi_size_t{sizeof(T) * local_size + comm_size
          + 2 * nr_nodes * sizeof(unsigned long)},
        GASPI_GROUP_ALL,
        GASPI_BLOCK,
        GASPI_ALLOC_DEFAULT
      ) == GASPI_SUCCESS);


      gaspi_segment_ptr(segment_id, &cont_seg_ptr);
      comm_seg_ptr = ((T*) cont_seg_ptr) + local_size;

      // Point vclock to the memory after communication segment
      vclock = (unsigned long*) (((T*) comm_seg_ptr) + COMM_BUFFER_NR_ELEMS);

      // TODO Ask Bernd, is this really necessary? Likely not
      for(int i = 0; i < 2*nr_nodes; i++){
        vclock[i] = op_nr;
      }

      gaspi_queue_create(&queue, GASPI_BLOCK);
    };


    Matrix(int rows, int cols, T init) : Matrix(rows, cols){
      set(init);
    }


    void set(T scalar){
      for(int i = 0; i < local_size; i ++){
        ((T*) cont_seg_ptr)[i] = scalar;
      }
      vclock[rank] = ++op_nr;
    }


    // Randomly initializes an int or long matrix
    void rand(int from, int to){
      if(std::is_same<T, int>::value || std::is_same<T, long>::value){
        for(int i = 0; i < local_size; i ++){
          ((T*) cont_seg_ptr)[i] = from + std::rand() % (to - from);
        }

        vclock[rank] = ++op_nr;
      }
      else{
        std::cout << "Rand only supports matrices of int or long\n";
      }
    }



    void set(int index, T value){
      if(index >= start_i && index <= end_i){
        ((T*) cont_seg_ptr)[index - start_i] = value;
      }
      vclock[rank] = ++op_nr;
    }


    // Gets a specific value from the Matrix. Either fetch the remote value
    // or wait untill all remotes have read our local value, then return.
    //
    // TODO A distributed distribution scheme of this value should be implemented.
    // currently we are likely to overload a single node
    T get(int index){

      if(index >= start_i && index <= end_i){
        // The value is local, wait until we have distributed the value

        gaspi_notification_id_t notify_id;
        gaspi_notification_t notify_val;

        // Wait for notification from all nodes
        for(int i = 0; i < nr_nodes - 1; ++i){

          // Note that we wait for notifs from start to start + #notifs - 1
          gaspi_notify_waitsome(
            segment_id,
            notif_ctr * nr_nodes, // notif start
            nr_nodes, // number of notifs
            &notify_id,
            GASPI_BLOCK
          );
          gaspi_notify_reset(segment_id, notify_id, &notify_val);

        }

        ++notif_ctr;
        return ((T*) cont_seg_ptr)[index - start_i];
      }
      else{
        // The rank holding the value
        int dest_rank = get_owner(index);

        // Wait until the rank is on the current OP
        wait_ranks.clear();
        wait_ranks.push_back(dest_rank);
        wait_for_vclocks(op_nr);

        gaspi_read(
          segment_id,
          comm_offset,
          dest_rank,
          segment_id + dest_rank - rank, // remote segment id
          sizeof(T) * (index - step * dest_rank), // Remote offset
          sizeof(T),
          queue,
          GASPI_BLOCK
        );

        gaspi_notify(
          segment_id + dest_rank - rank, // remote segment
          dest_rank,
          notif_ctr * nr_nodes + rank, // notif id
          13, // notif val, not used
          queue,
          GASPI_BLOCK
        );

        ++notif_ctr;
        gaspi_wait(queue, GASPI_BLOCK);

        return ((T*) comm_seg_ptr)[0];
      }
      vclock[rank] = ++op_nr;
    }



    // WARNING Only use this for debugging, it has very poor performance
    void print(){
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


    template<typename Lambda>
    void print(Lambda lamb){
      for(int i = 0; i < nr_nodes; i++){
        if(i == rank){
          for(int j = 0; j < local_size; j++){
            std::cout << lamb(((T*) cont_seg_ptr)[j]) << ",  ";
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
