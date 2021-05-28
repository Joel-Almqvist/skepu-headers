#ifndef REDUCE_HPP
#define REDUCE_HPP

#include <type_traits>
#include <numeric>
#include <cmath>

#include <GASPI.h>

#include "matrix.hpp"

namespace skepu{

  template<typename ReduceFunc>
  class Reduce1D{
  private:
    ReduceFunc func;
    bool always_use_buffer;

  protected:

    // Constructor used by MapReduce
    Reduce1D(ReduceFunc func, bool flag) : func{func}, always_use_buffer{flag} {};

  public:

    Reduce1D(ReduceFunc func) : func{func}, always_use_buffer{false} {};

    // Dummy
    template<typename T>
    void setBackend(T){};

    // If no initial value is given we use a default value
    template<typename Container>
    typename Container::value_type operator()(Container& cont){
      using T = typename Container::value_type;
      return operator()(cont, T{});
    }


    template<typename Container, typename Val>
    typename Container::value_type operator()(Container& cont, Val init){
      using T = typename Container::value_type;

      static_assert(std::is_same<T, Val>::value);

      cont.op_nr++;
      cont.vclock[cont.rank] = cont.op_nr;

      // We either read from the double buffer or the actual segment depending
      // on flush status. We do not flush.
      T* from;

      bool has_flushed = cont.last_flush[cont.rank] > cont.last_mod_op
        || cont.last_mod_op == 0;

      if(has_flushed && !always_use_buffer){
        from = (T*) cont.cont_seg_ptr;
      }
      else{
        from = (T*) cont.local_buffer;
      }

      // Calculate the local value
      #pragma omp parallel
      {

        T& to = *((T*) cont.comm_seg_ptr + cont.norm_partition_size * cont.rank
          + omp_get_thread_num());

        bool first_time = false;

        #pragma omp for schedule(static)
        for(size_t i = 0; i <= cont.end_i - cont.start_i; ++i){

          if(first_time){
            first_time = false;
            to = from[i];
          }
          else{
            to = func(to, from[i]);
          }
        }

        if(cont.rank == 0){
          #pragma omp single
          {
            to = func(to, init);
          }
        }

        #pragma omp barrier

        if(omp_get_thread_num() == 0){

          T* other_res = (T*) cont.comm_seg_ptr + cont.norm_partition_size *
            cont.rank;

          for(int i = 1; i < omp_get_num_threads(); ++i){
            to = func(to, other_res[i]);
          }
        }
      }

      // Indicate to other ranks that our value is ready to be used
      cont.vclock[cont.rank] = ++cont.op_nr;

      int iterations = std::ceil(std::log2(cont.nr_nodes));
      int dest_rank;
      int step;
      int prev_step;
      unsigned remote_offset;

      T& loc_val = *((T*) cont.comm_seg_ptr +
        cont.norm_partition_size * cont.rank);

      T& rec_val = *((T*) cont.comm_seg_ptr +
        cont.norm_partition_size * cont.rank + 1);

      // A distributed gather which in the end puts the data to node 0
      for(int i = 0; i < iterations; i++){

        step = pow(2, i + 1);
        prev_step = pow(2, i);

          if(cont.rank % step == 0){

            dest_rank = cont.rank + prev_step;

            if(dest_rank < cont.nr_nodes){
              cont.wait_ranks.push_back(dest_rank);
              cont.wait_for_vclocks(cont.op_nr);

              if(dest_rank == cont.nr_nodes - 1){
                remote_offset = cont.last_partition_comm_offset +
                  sizeof(T) * cont.norm_partition_size * dest_rank;
              }
              else{
                remote_offset = cont.norm_partition_comm_offset +
                  sizeof(T) * cont.norm_partition_size * dest_rank;
              }

              gaspi_read(
                cont.segment_id,
                cont.comm_offset + (cont.start_i + 1) * sizeof(T), // local offset
                dest_rank,
                cont.segment_id - cont.rank + dest_rank, // rem seg id
                remote_offset,
                sizeof(T), //size
                cont.queue,
                GASPI_BLOCK
              );

              gaspi_wait(cont.queue, GASPI_BLOCK);

              loc_val = func(loc_val, rec_val);

            }
          }
          cont.vclock[cont.rank] = ++cont.op_nr;
      }


      // A distributed broadcast from node 0 to all the other.
      for(int i = 0; i < iterations; i++){

        step = pow(2, i + 1);
        prev_step = pow(2, i);

        if(cont.rank >= prev_step && cont.rank < step){

          dest_rank = cont.rank - prev_step;

          cont.wait_ranks.push_back(dest_rank);
          cont.wait_for_vclocks(cont.op_nr);


          if(dest_rank == cont.nr_nodes - 1){
            remote_offset = cont.last_partition_comm_offset +
              sizeof(T) * cont.norm_partition_size * dest_rank;
          }
          else{
            remote_offset = cont.norm_partition_comm_offset +
              sizeof(T) * cont.norm_partition_size * dest_rank;
          }

          gaspi_read(
            cont.segment_id,
            cont.comm_offset + cont.start_i * sizeof(T), // local offset
            dest_rank,
            cont.segment_id - cont.rank + dest_rank, // rem seg id
            remote_offset,
            sizeof(T), //size
            cont.queue,
            GASPI_BLOCK
          );

          gaspi_wait(cont.queue, GASPI_BLOCK);

        }
        cont.vclock[cont.rank] = ++cont.op_nr;
      }

      return loc_val;
    }

     // Should take in a backend type
     void setBackend(){}

     // Need to be implemented
     void setReduceMode(){};
  };


  // Template deduction for classes are not allowed in c++11
  // This solves this problem
  template<typename ReduceFunc>
  Reduce1D<ReduceFunc> Reduce(ReduceFunc func){
    return Reduce1D<ReduceFunc>{func};
  }

} // end of namespace skepu
#endif // REDUCE_HPP
