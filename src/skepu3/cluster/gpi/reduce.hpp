#ifndef REDUCE_HPP
#define REDUCE_HPP

#include <type_traits>
#include <numeric>
#include <cmath>
#include <functional>

#include <GASPI.h>

#include "matrix.hpp"

namespace skepu{

  template<typename ret_t, typename arg_t>
  class Reduce1D{
  private:
    std::function<ret_t(arg_t, arg_t)> func;
    ret_t start_value;

  public:

    Reduce1D(std::function<ret_t(arg_t, arg_t)> func) : func{func},
      start_value{} {

    };

    void setStartValue(ret_t val){
      start_value = val;
    }


    // Dummy
    template<typename T>
    void setBackend(T){};



    template<typename Container>
    typename Container::value_type operator()(Container& cont){
      using T = typename Container::value_type;

      static_assert(std::is_same<T, ret_t>::value);

      cont.op_nr++;
      cont.vclock[cont.rank] = cont.op_nr;

      // We either read from the double buffer or the actual segment depending
      // on flush status. We do not flush.
      T* from;

      bool has_flushed = cont.last_flush[cont.rank] > cont.last_mod_op[cont.rank]
        || cont.last_mod_op[cont.rank] == 0;

      if(has_flushed){
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
            to = func(to, start_value);
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
        // This rank is reading during this iteration

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

        else if(cont.rank < prev_step && cont.rank < step){
          // This rank is the target of a read during this iteration

          dest_rank = cont.rank + prev_step;
          cont.wait_ranks.push_back(dest_rank);
        }

        cont.vclock[cont.rank] = ++cont.op_nr;
      }

      // To ensure that the remote partial sum is not overwritten we need to
      // wait for all remotes to get it.
      cont.wait_for_vclocks(cont.op_nr);

      return loc_val;
    }

     // Should take in a backend type
     void setBackend(){}

     // Need to be implemented
     void setReduceMode(){};
  };


  // The constructor wrapper for function pointers
  template <typename Ret, typename arg_t>
  Reduce1D<Ret, arg_t> Reduce(Ret(*arg)(arg_t, arg_t)){
      return Reduce1D<Ret, arg_t>((std::function<Ret(arg_t, arg_t)>)arg);
  }

  template <typename Ret, typename arg_t>
  Reduce1D<Ret, arg_t> _red_helper(std::function<Ret(arg_t, arg_t)> func){
    return Reduce1D<Ret, arg_t>(func);
  }


  template <typename Lambda>
  auto Reduce(Lambda&& lambda) -> decltype(_red_helper(lambda_cast(lambda))){

    return _red_helper(lambda_cast(lambda));
  }


} // end of namespace skepu
#endif // REDUCE_HPP
