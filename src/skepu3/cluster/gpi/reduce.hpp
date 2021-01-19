#ifndef REDUCE_HPP
#define REDUCE_HPP

#include <type_traits>
#include <numeric>
#include <cmath>

#include <GASPI.h>
#include <matrix.hpp>
#include <skeleton_base.hpp>

namespace skepu{

  template<typename ReduceFunc>
  class Reduce1D : public _gpi::skeleton_base{
  private:
    ReduceFunc func;
  public:

    Reduce1D(ReduceFunc func) : func{func} {};

     template<typename Container>
     typename Container::value_type operator()(Container& cont){
       using T = typename Container::value_type;

       if(is_skepu_container<Container>::value){

         // TODO Find a better solution than a barrier.
         // Currently it prevents multiple operations from modifying
         // the communication and container segments and needs to be at
         // the start of all functions which use these.
         gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);

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
