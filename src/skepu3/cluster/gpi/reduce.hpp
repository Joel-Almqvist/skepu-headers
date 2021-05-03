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
  public:

    Reduce1D(ReduceFunc func) : func{func} {};


    // If no initial value is given we use a default value
    template<typename Container>
    typename Container::value_type apply(Container& cont){
      using T = typename Container::value_type;
      return apply(cont, T{});
    }



    template<typename Container, typename Val>
    typename Container::value_type apply(Container& cont, Val init){
      using T = typename Container::value_type;

      static_assert(std::is_same<T, Val>::value);


      cont.wait_for_constraints();
      //cont.conditional_flush();

      cont.op_nr++;
      cont.vclock[cont.rank] = cont.op_nr;

      // Calculate the local value
      #pragma omp parallel
      {

        // TODO this needs to change depending on flush status
        T* from = (T*) cont.local_buffer;

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

        #pragma omp single
        {
          to = func(to, init);
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
      cont.op_nr++;
      cont.vclock[cont.rank] = cont.op_nr;

      int i = cont.rank == cont.nr_nodes - 1 ? 0 : cont.rank + 1;
      bool has_looped = false;
      while(true){


        cont.wait_ranks.push_back(i);
        cont.wait_for_vclocks(cont.op_nr);

        // TODO Read from rank

        // TODO apply the function

        // TODO create a distributed version of this insteads


        ++i;
        if(i == cont.rank){
         ++i;
       }
        if(i == cont.nr_nodes -1 && !has_looped){
          i = 0;
        }
        else if(i >= cont.nr_nodes){
          break;
        }
      }

      return 4;


    }


     template<typename Container>
     typename Container::value_type operator()(Container& cont){
       using T = typename Container::value_type;

       if(is_skepu_container<Container>::value){

         cont.vclock[cont.rank] = ++cont.op_nr;
         gaspi_notification_id_t notify_id;

         T local_sum = func(((T*) cont.cont_seg_ptr)[0], ((T*) cont.cont_seg_ptr)[1]);

         for(int i = 2; i < cont.local_size; i++){
           local_sum = func(local_sum, ((T*) cont.cont_seg_ptr)[i]);
         }

         ((T*) cont.comm_seg_ptr)[0] = local_sum;


         int iterations = std::ceil(std::log2(cont.nr_nodes));

         bool received = true;
         int step;
         int remote_comm_offset;
         gaspi_notification_t notify_val = 0;


         // TODO add the following:
         // 1 - A proper wait which only does so before writing
         // 2 - Multithreaded
         // 3 - Better check for queue overflowing
         gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);
         gaspi_wait(cont.queue, GASPI_BLOCK);

         for(int i = 0; i < iterations; i++){

           step = pow(2, i);

           // Do work only if we received a value or if it is the first iteration
           if(received || i == 0){

             if(cont.rank % (step * 2) == step - 1){
               // Send

               // The last rank has a different offset
               remote_comm_offset = cont.rank + step == cont.nr_nodes - 1 ?
                  cont.last_partition_comm_offset :
                  cont.norm_partition_comm_offset;

                  // Make sure we do not overwrite the remote comm buffer
                  cont.wait_ranks.clear();
                  cont.wait_ranks.push_back(cont.rank + step);
                  cont.wait_for_vclocks(cont.op_nr);

                gaspi_write_notify(cont.segment_id, // local seg
                    cont.comm_offset, // local offset
                    cont.rank + step, // dest rank
                    cont.segment_id + step,
                    remote_comm_offset + (i + 1) * sizeof(T), // remote offset
                    sizeof(T),
                    i + 1, // notif ID
                    123,
                    cont.queue,
                    GASPI_BLOCK);

               received = false;
             }
             else if(cont.rank % (step * 2) == (step * 2) - 1){
               // Receive

                gaspi_notify_waitsome(
                  cont.segment_id,
                  i + 1,
                  1,
                  &notify_id,
                  GASPI_BLOCK);

                  gaspi_notify_reset(cont.segment_id, notify_id, &notify_val);

                ((T*) cont.comm_seg_ptr)[0] = func(((T*) cont.comm_seg_ptr)[0],
                      ((T*) cont.comm_seg_ptr)[i + 1]);
             }
             else{
               // Do nothing
             }
           }
         }


         // Distribute the reduces value
         for(int i = 0; i < iterations; i++){

           step = pow(2, i);

           if(cont.rank > (cont.nr_nodes - 1) - step){

             if(cont.rank - step >= 0){

               remote_comm_offset = cont.rank - step == cont.nr_nodes - 1 ?
                  cont.last_partition_comm_offset :
                  cont.norm_partition_comm_offset;

              gaspi_write_notify(cont.segment_id,
                  cont.comm_offset,
                  cont.rank - step,
                  cont.segment_id - step, // dest rank
                  remote_comm_offset, // remote offset
                  sizeof(T),
                  i + iterations,
                  123,
                  cont.queue,
                  GASPI_BLOCK);

            }
           }

           else if(cont.rank > (cont.nr_nodes - 1) - 2 * step){
               // receive
               gaspi_notify_waitsome(
                 cont.segment_id,
                 i + iterations,
                 1,
                 &notify_id,
                 GASPI_BLOCK);

                 gaspi_notify_reset(cont.segment_id, notify_id, &notify_val);

           }
         }


         return ((T*) cont.comm_seg_ptr)[0];
       }
       else{
         std::cout << "ERROR Non Skepu container\n";
         return typename Container::value_type{};
       }
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
