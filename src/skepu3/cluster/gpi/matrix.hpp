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

    using is_skepu_container = decltype(true);

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
    int start_i;
    int end_i;

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
      Matrix& dest_cont
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
            queue,
            GASPI_BLOCK
          );
          sent_elems += nr_elems_to_send;
        }
        gaspi_wait(queue, GASPI_BLOCK);
    }






    // Puts all elements from start to end (these are global indeces) into
    // the given GASPI segment. Many to one communication pattern
    //
    // TODO Change this function to be read based, it is currently very
    // counter intuitive to use
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

      unsigned long remote_offset = dest_rank == nr_nodes - 1 ?
        last_partition_vclock_offset :
        norm_vclock_offset;

      //
      gaspi_read_notify(
        segment_id,
        vclock_offset + sizeof(unsigned long) * nr_nodes, // local offset
        dest_rank,
        dest_seg_id,
        remote_offset,
        sizeof(unsigned long) * nr_nodes,
        2 * nr_nodes + 1, // notif id
        queue,
        GASPI_BLOCK
      );

      gaspi_notification_t notify_val = 0;
      gaspi_notification_id_t first_id;

      gaspi_notify_waitsome(
        segment_id,
        2 * nr_nodes + 1, // notif begin
        1, // number of notif
        &first_id,
        GASPI_BLOCK
      );

      gaspi_notify_reset(segment_id, first_id, &notify_val);

      for(int i = 0; i < nr_nodes; i++){
        if(i == rank){
          vclock[i] = op_nr;
        }
        else{
          vclock[i] = std::max(vclock[i + nr_nodes], vclock[i]);
        }
      }

      if(false && rank == 1){
        std::cout << "Rank " << rank <<  " seg_id "<< (int) segment_id << " new vlock: ";
        for(int i = 0; i < nr_nodes; i++){
          std::cout << vclock[i] << ", ";
        }
        std::cout << std::endl;
      }
    };




  public:

    using value_type = T;


    bool operator==(Matrix<T>& that){
      return this == &that;
    }


        //TODO Make this private possibly?
        void wait_for_vclocks(int wait_val){
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


    T get(int index){

      if(index >= start_i && index <= end_i){
        // The value is local
        return ((T*) cont_seg_ptr)[index - start_i];
      }
      else{
        // The rank holding the value
        int dest_rank = std::floor(((double) index) / step);
        if(dest_rank >= nr_nodes){
          dest_rank = nr_nodes - 1;
        }

        wait_ranks.clear();
        wait_ranks.push_back(dest_rank);
        wait_for_vclocks(op_nr);


        gaspi_read_notify(
          segment_id,
          0,
          dest_rank,
          segment_id + dest_rank - rank, // remote segment id
          sizeof(T) * (index - step * dest_rank), // Remote offset
          sizeof(T),
          rank, // Notification id
          queue,
          GASPI_BLOCK
        );


        gaspi_notification_id_t notify_id;
        gaspi_notification_t notify_val = 0;
        gaspi_notify_waitsome(
          segment_id,
          rank,
          1,
          &notify_id,
          GASPI_BLOCK
        );
        gaspi_notify_reset(segment_id, notify_id, &notify_val);

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
    }


};

  template<typename T>
  struct is_skepu_container<skepu::Matrix<T>> : std::true_type {};
}

#endif // MATRIX_HPP
