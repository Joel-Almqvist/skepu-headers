#ifndef MAP_HPP
#define AMP_HPP

#include <type_traits>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <deque>
#include <utility>
#include <mutex>
#include <cstring>

#include <omp.h>
#include <GASPI.h>

#include "matrix.hpp"
#include "utils.hpp"
#include "proxy.hpp"
#include "build_tup.hpp"


namespace skepu{



  template<typename Function, int nr_args>
  class Map1D{
  private:
    Function func;
    using arg_tup_t = typename _gpi::get_tup_t<Function, nr_args - 1>::type;



    /* This is a help function for Map() with two arguments
    *
    * Fetches the values whithin the range [dest_cont.start_i, dest_cont.end_i]
    * which the specified rank owns.
    *
    * Avoids doing work which is known to have been done previously.
    *
    * A pair of values which are from different ranks are called orphans untill
    * their partner has been found. Any value read is compared with the existing
    * orphans in the orphans deque and either matched or added to it.
    *
    * The values are handled in one of the 3 following ways:
    *
    * 1 - Apply func() to c1_val and c2_val and store it in dest_cont. This
          occurs if both values happens to be read simultaneously in the buffer.
    *
    * 2 - The value is matched with a orphaned, func is applied and the result
          stored in dest_cont
    *
    * 3 - The value is orphaned and added to orphans deque
    */
    template<typename T>
    void apply_rank_unique(Matrix<T>& dest_cont, Matrix<T>& cont1, Matrix<T>& cont2,
      int rank, std::deque<std::tuple<unsigned int, unsigned int, T>>& orphans,
      bool first_time = true){

      int c1_last_elem_global_index = rank != cont1.nr_nodes - 1 ?
       (rank + 1) * cont1.step -1 :
       cont1.global_size - 1;

      int c1_first_elem_global_index = rank * cont1.step;

      // Are the reads limited by dest_cont or do we want the whole of c1?
      int c1_read_from = std::max(dest_cont.start_i, c1_first_elem_global_index);

      int c1_to = std::min(dest_cont.end_i, c1_last_elem_global_index);


      int c2_last_elem_global_index = rank != cont2.nr_nodes - 1 ?
       (rank + 1) * cont2.step - 1 :
       cont2.global_size - 1;

      int c2_first_elem_global_index = rank * cont2.step;

      int c2_read_from = std::max(dest_cont.start_i, c2_first_elem_global_index);

      int c2_to = std::min(dest_cont.end_i, c2_last_elem_global_index);


      int overlap_lower_bound = std::max({dest_cont.start_i, cont1.start_i,
        cont2.start_i});

      int overlap_upper_bound = std::min({dest_cont.end_i, cont1.end_i,
        cont2.end_i});

      if(overlap_upper_bound < overlap_lower_bound){
        // There is no overlap, no need to change our boundaries
      }
      else if((overlap_lower_bound < c1_read_from && overlap_upper_bound > c1_to)
      ||  (overlap_lower_bound < c2_read_from && overlap_upper_bound > c2_to)){
        // We reading from a non contiguous range

        if(first_time){
          // Handle the lower part of the non contiguous range
          apply_rank_unique(dest_cont, cont1, cont2, rank, orphans, false);

          // TODO Reevaluate if this case is handled correctly
          c1_to = overlap_lower_bound - 1;
          c2_to = overlap_lower_bound - 1;
        }
        else{
          // Handle upper part
          c1_read_from = overlap_upper_bound + 1;
          c2_read_from = overlap_upper_bound + 1;
        }

      }

      else if(overlap_lower_bound > c1_read_from ||
        overlap_lower_bound  > c2_read_from){
          // We are below the overlap
          if(c1_to >= c1_read_from)
            c1_to = overlap_lower_bound - 1;

          if(c2_to >= c2_read_from)
            c2_to = overlap_lower_bound - 1;

      }
      else if(overlap_upper_bound < c1_to || overlap_upper_bound < c2_to){
        // We are above the overlap
        c1_read_from = overlap_upper_bound + 1;
        c2_read_from = overlap_upper_bound + 1;
      }
      else{
        // TODO Reevaluate this case
        // There is a full overlap and the work has already been done
        return;
      }


      int transfer_size = dest_cont.COMM_BUFFER_NR_ELEMS / 2;

      // We may not want any values one of the two containers of this rank
      if(c1_to < dest_cont.start_i || c1_to < c1_read_from){
        c1_to = -2;
      }

      if(c2_to < dest_cont.start_i || c2_to < c2_read_from){
        c2_to = -2;
      }

      // Start at -1 so we can increment at start to make the while more readable
      int t = -1;

      int c1_send_up_to;
      int c2_send_up_to;

      int lowest_shared_i = std::max(c1_read_from, c2_read_from);

      // Globally indexed counters for which elements to send in every
      // transmission iteration
      int c1_rec_from = -1;
      int c1_rec_to;

      // Globally indexed
      int c2_rec_from = -1;
      int c2_rec_to;


      bool c1_has_unread = true;
      bool c2_has_unread = true;

      int c1_sent_elems = 0;
      int c2_sent_elems = 0;

      // Stores values which have not yet been paired up.
      // Tuple scheme: (global_index, 0/1 indicates removal, value)
      std::vector<std::tuple<unsigned int, unsigned int, T>> new_orphans{};
      std::mutex vlock;

      // Loop while there are elements to transfer from c1 or c2
      while(c1_has_unread || c2_has_unread){

        ++t;
        c1_has_unread = c1_sent_elems <= c1_to - c1_read_from;
        c2_has_unread = c2_sent_elems <= c2_to - c2_read_from;

        if(c1_has_unread){
          c1_rec_from = t * transfer_size + c1_read_from;
          c1_rec_to = std::min(c1_read_from + transfer_size * (t + 1) - 1, c1_to);

          c1_sent_elems += c1_rec_to - c1_rec_from + 1;

          dest_cont.read_range(c1_rec_from, c1_rec_to, 0, cont1);
        }

        if(c2_has_unread){
          c2_rec_from = t * transfer_size + c2_read_from;
          c2_rec_to = std::min(c2_read_from + transfer_size * (t + 1) - 1, c2_to);

          c2_sent_elems += c2_rec_to - c2_rec_from + 1;

          dest_cont.read_range(c2_rec_from, c2_rec_to, transfer_size * sizeof(T), cont2);
        }

        if(!c1_has_unread && !c2_has_unread){
          // Early break, we won't be transfering any more values
          break;
        }


        T* store_at = (T*) dest_cont.cont_seg_ptr;
        T* comm_ptr = (T*) dest_cont.comm_seg_ptr;

        #pragma omp parallel
        {
          for(int i = omp_get_thread_num(); i < transfer_size;
          i = i + omp_get_num_threads()){

            // c1_to < 0 indicates that nothing is being transfered
            if(c1_rec_from + i > c1_rec_to || c1_to < 0){
              // Base case
              break;
            }

            // if the matching value exists in the buffer
            if(c2_rec_from >= 0 && c1_rec_from + i >= c2_rec_from &&
                c1_rec_from + i <= c2_rec_to){
              int pair_offset = transfer_size + c1_rec_from + i - c2_rec_from;

              store_at[c1_rec_from - dest_cont.start_i + i]
              = func(comm_ptr[i], comm_ptr[pair_offset]);
            }

            else{
              vlock.lock();
              new_orphans.push_back(std::tuple<unsigned int, unsigned int, T>
                {c1_rec_from + i, 0, comm_ptr[i]});
              vlock.unlock();
            }
          }


          for(int i = omp_get_thread_num(); i < transfer_size;
          i = i + omp_get_num_threads()){

            // c2_to < 0 indicates that nothing is being transfered
            if(c2_rec_from + i > c2_rec_to || c2_to < 0){
              // Base case
              break;
            }

            // cx_rec_from may be negative to indicate that nothing is
            // transfered from the container
            if(c1_rec_from >= 0 && c2_rec_from + i >= c1_rec_from &&
                c2_rec_from + i <= c1_rec_to){
              // We have a matching pair in the buffer but the previous loop
              // already managed this case
              continue;
            }
            else{
              vlock.lock();
              new_orphans.push_back(std::tuple<unsigned int, unsigned int, T>
                {c2_rec_from + i, 0, comm_ptr[transfer_size + i]});
              vlock.unlock();
            }
          }

          #pragma omp barrier

          // Match the new orphans with the previous ones
          for(int i = omp_get_thread_num(); i < new_orphans.size();
          i = i + omp_get_num_threads()){

            auto it = std::find_if(orphans.begin(), orphans.end(), [&new_orphans, i]
            (typename std::tuple<unsigned int, unsigned int, T> tup){
              return std::get<0>(tup) == std::get<0>(new_orphans[i]);
              }
            );

            if(it != orphans.end()){
              int index = std::get<0>(*it);
              store_at[index - dest_cont.start_i] = func(std::get<2>(*it),
              std::get<2>(new_orphans[i]));

              // Mark the two ex orphans for removal
              std::get<1>(new_orphans[i]) = 1;
              std::get<1>(*it) = 1;
            }
            else{
            }
          }
        } // End of parallel region


        orphans.erase(
          std::remove_if(orphans.begin(), orphans.end(), []
          (typename std::tuple<unsigned int, unsigned int, T> tup){
            return std::get<1>(tup) == 1;
          }),
          orphans.end()
        );

        std::copy_if(
          new_orphans.begin(), new_orphans.end(), std::back_inserter(orphans),
          [](typename std::tuple<unsigned int, unsigned int, T> tup){
            return std::get<1>(tup) != 1;
          }
        );

        new_orphans.clear();

      } // end of while()
    } // end of apply_rank_unique()



    template<int ctr, typename Dest, typename Tup, typename Curr, typename... Rest>
    void build_tuple(int i, bool local_only, Dest& dest,
       Tup& tup, Curr& curr, Rest&... rest){

      if(local_only){
        build_tuple_local<ctr>(i, dest, tup, curr);

        // Early break if the tuple need a remote value
        //if(tup_flag != -1){
          build_tuple<ctr + 1>(i, local_only, dest, tup, rest...);
        //}
      }

      else{
        build_tuple_remote<ctr>(i, dest, tup, curr);

        // No early break possible for remote-valued tuples
        build_tuple<ctr + 1>(i, local_only, dest, tup, rest...);
      }
    }

    template<int ctr, typename Dest, typename Tup, typename Curr>
    void build_tuple(int i, bool local_only, Dest& dest,
       Tup& tup, Curr& curr){
         if(local_only){
           build_tuple_local<ctr>(i, dest, tup, curr);
         }
         else{
           build_tuple_remote<ctr>(i, dest, tup, curr);
         }
    }

    // Adds a local value to the tuple or indicate error
    template<int ctr, typename Dest, typename Tup, typename Curr>
    void build_tuple_local(int i, Dest& dest, Tup& tup, Curr& curr){
      using T = typename Curr::value_type;
      T* val_ptr = (T*) curr.cont_seg_ptr;

      if(i >= curr.start_i && i <= curr.end_i){
        // The data is local
        std::get<ctr>(tup) = val_ptr[i - curr.start_i];
      }
    }



    template<int ctr, typename Dest, typename Tup, typename Curr>
    void build_tuple_remote(int i, Dest& dest, Tup& tup, Curr& curr){
      using T = typename Curr::value_type;
      T* val_ptr = (T*) curr.cont_seg_ptr;
      T* comm_ptr = (T*) curr.comm_seg_ptr;

      if(i >= curr.start_i && i <= curr.end_i){
        // The data is local, parts of a tuple is allowed to be local
        std::get<ctr>(tup) = val_ptr[i - curr.start_i];
      }
      else if(i < curr.start_i){
        int offset = i - dest.start_i;

        gaspi_wait(curr.queue, GASPI_BLOCK);
        std::get<ctr>(tup) = comm_ptr[offset];
      }
      else{
        // i > curr.end_i
        int offset = std::max(curr.start_i - dest.start_i, long{0});

        // If curr.end_i < dest.end_i then we have transfered elements
        offset += std::max(dest.end_i - curr.end_i, long{0});
        gaspi_wait(curr.queue, GASPI_BLOCK);
        std::get<ctr>(tup) = comm_ptr[offset];
      }

    }



  public:


    Map1D(Function func) : func{func} {};


    /* Performs Map with 1 argument the following way:
    * 1 - Apply map to locally owned elements
    * 2 - Wait for all ranks which have elements we need to access remotely
    * 3 - Gather remote elements and apply Map on them
    * 4 - Wait until all remote ranks are done reading from us
    */
    template<typename DestCont, typename From>
     auto old1(DestCont& dest_cont, From& from) ->
     decltype(
       std::declval<typename DestCont::is_skepu_container>(),
       std::declval<typename From::is_skepu_container>(),

       // Check that the lambda takes one argument
       std::declval<Function>()(std::declval<typename DestCont::value_type>()),

       std::declval<void>()) {
         using T = typename DestCont::value_type;

         int lowest_local_i = std::max({from.start_i, dest_cont.start_i});
         int highest_local_i = std::min({from.end_i, dest_cont.end_i});

         T* dest_ptr = (T*) dest_cont.cont_seg_ptr;
         T* from_ptr = (T*) from.cont_seg_ptr;


         // Do the local work
         #pragma omp parallel shared(dest_cont, dest_ptr, from_ptr, \
           from, lowest_local_i, highest_local_i)
         {

           for(int i = lowest_local_i + omp_get_thread_num();
             i <= highest_local_i;
             i = i + omp_get_num_threads()){
               dest_ptr[i - dest_cont.start_i] = func(from_ptr[i - from.start_i]);
           }
         }


         int lowest_rank_i_depend_on = from.get_owner(dest_cont.start_i);
         int highest_rank_i_depend_on = from.get_owner(dest_cont.end_i);

         dest_cont.wait_ranks.clear();
         for(int i = lowest_rank_i_depend_on; i <= highest_rank_i_depend_on; i++){
           if (i != dest_cont.rank)
             dest_cont.wait_ranks.push_back(i);
         }


         dest_cont.wait_for_vclocks(dest_cont.op_nr);

         int transfer_size = dest_cont.COMM_BUFFER_NR_ELEMS;

         int range1_start = dest_cont.start_i;
         int range1_end = from.start_i;

         int range2_start = from.end_i;
         int range2_end = dest_cont.end_i;

         int start;
         int end;
         int t = -1;

         T* dest_cptr = (T*) dest_cont.comm_seg_ptr;
         if(range1_end > range1_start){
           start = range1_start;
           end = start;

           while(true){
             ++t;
             start = start + t * transfer_size;

             if(start >= range1_end){
               break;
             }

             // Minues one because read_range is inclusive
             end = std::min(end + (t + 1) * transfer_size, range1_end -1);
             dest_cont.read_range(start, end, 0, from);


             for(int i = 0; i <= end - start; i++){
               dest_ptr[start + i - dest_cont.start_i] = func(dest_cptr[i]);
             }
           }
         }


         if(range2_end > range2_start){

           start = range2_start;
           end = start;
           t = -1;

           while(true){
             ++t;
             start = start + t * transfer_size;

             if(start >= range2_end){
               break;
             }

             // Minues one because read_range is inclusive
             end = std::min(end + (t + 1) * transfer_size, range2_end -1);
             dest_cont.read_range(start, end, 0, from);

             for(int i = 0; i <= end - start; i++){
               dest_ptr[start + i - dest_cont.start_i] = func(dest_cptr[i]);
             }
           }
         }

         // Indicate that the first phase of the operation is done
         dest_cont.vclock[dest_cont.rank] = ++dest_cont.op_nr;


         // We must wait for other ranks which may want to read from us
         int lowest_rank_depending_on_me = dest_cont.get_owner(from.start_i);
         int highest_rank_depending_on_me = dest_cont.get_owner(from.end_i);

         dest_cont.wait_ranks.clear();

         for(int i = lowest_rank_depending_on_me; i <= highest_rank_depending_on_me; i++){
           if (i != dest_cont.rank)
             dest_cont.wait_ranks.push_back(i);
         }

         dest_cont.wait_for_vclocks(dest_cont.op_nr);

         dest_cont.vclock[dest_cont.rank] = ++dest_cont.op_nr;

       }


     /*
     * Applies Map with 2 arguments on 3 containers, which or may not be
     * unique.
     */
     template<typename DestCont, typename Cont1, typename Cont2>
      auto old2(DestCont& dest_cont, Cont1& cont1, Cont2& cont2) ->
      decltype(
        std::declval<typename DestCont::is_skepu_container>(),
        std::declval<typename Cont1::is_skepu_container>(),
        std::declval<typename Cont2::is_skepu_container>(),

        // Check that the lambda takes two arguments
        std::declval<Function>()(std::declval<typename DestCont::value_type>(),
        std::declval<typename DestCont::value_type>()),
        std::declval<void>()){

          using T = typename DestCont::value_type;

          int lowest_rank_i_depend_on = std::min(cont1.get_owner(dest_cont.start_i),
            cont2.get_owner(dest_cont.start_i));

          int highest_rank_i_depend_on = std::max(cont1.get_owner(dest_cont.end_i),
            cont2.get_owner(dest_cont.end_i));


          int lowest_local_i = std::max({cont1.start_i, cont2.start_i, dest_cont.start_i});
          int highest_local_i = std::min({cont1.end_i, cont2.end_i, dest_cont.end_i});


          // Do the work which does not need remote communication
          #pragma omp parallel shared(dest_cont, cont2, cont1, lowest_local_i, highest_local_i)
          {

            for(int i = lowest_local_i + omp_get_thread_num();
              i <= highest_local_i;
              i = i + omp_get_num_threads()){

              ((T*) dest_cont.cont_seg_ptr)[i - dest_cont.start_i] =
              func(((T*) cont1.cont_seg_ptr)[i - cont1.start_i],
                ((T*) cont2.cont_seg_ptr)[i - cont2.start_i]);
            }
          }

          // Manually update the vclock once since there can't be any ranks
          // which are ready at this point.
          dest_cont.get_vclock(lowest_rank_i_depend_on, dest_cont.segment_id - dest_cont.rank
           + lowest_rank_i_depend_on);

         bool got_ranks_part [highest_rank_i_depend_on - lowest_rank_i_depend_on
            + 1] = {};

         bool work_remaining = !(lowest_rank_i_depend_on == dest_cont.rank
         && highest_rank_i_depend_on == dest_cont.rank);


         // Orphans are unique for each destination rank, but shared for multiple
         // calls to the same rank.
         std::deque<std::tuple<unsigned int, unsigned int, T>> orphans{};
         int j;
         while(work_remaining){

           work_remaining = false;
           j = -1;

           for(int i = lowest_rank_i_depend_on; i <= highest_rank_i_depend_on; i++){
             j++;

             if(got_ranks_part[j] == true){
               continue;
             }

             else if(dest_cont.vclock[i] < dest_cont.op_nr){

               dest_cont.get_vclock(i, dest_cont.segment_id - dest_cont.rank + i);

               if(dest_cont.vclock[i] < dest_cont.op_nr){
                 // The updated vclock is still not ready
                 work_remaining = true;
                 continue;
               }
             }

             // We may or may not have updated the vclock but there is work
             // to be done at this point
             got_ranks_part[j] = true;
             apply_rank_unique(dest_cont, cont1, cont2, i, orphans);

           } // End of for

           std::this_thread::sleep_for(std::chrono::milliseconds(5));

         }

         assert(orphans.size() == 0);

         int lowest_rank_depending_on_me = std::min(dest_cont.get_owner(cont1.start_i),
           dest_cont.get_owner(cont2.start_i));

         int highest_rank_depending_on_me = std::max(dest_cont.get_owner(cont1.end_i),
           dest_cont.get_owner(cont2.end_i));

         dest_cont.wait_ranks.clear();
         dest_cont.vclock[dest_cont.rank] = ++dest_cont.op_nr;


         for(int i = lowest_rank_depending_on_me; i <= highest_rank_depending_on_me; i++){
           if (i != dest_cont.rank)
             dest_cont.wait_ranks.push_back(i);
         }

         // TODO This wait should be done at the start of Map() instead
         // since we are waiting unnecesarily currently. Add this after variadic
         dest_cont.wait_for_vclocks(dest_cont.op_nr);
         dest_cont.vclock[dest_cont.rank] = ++dest_cont.op_nr;

     } // end of apply_on_unique_conts



     // Should take in a backend type
     void setBackend(){}

     // Need to be implemented
     void setReduceMode(){};


     template<typename DestCont, typename ... Conts>
      auto apply_old(DestCont& dest_cont, Conts&... conts) ->
      decltype(
        std::declval<typename DestCont::is_skepu_container>(),
        std::declval<void>()){

          using T = typename DestCont::value_type;
          using tup_type = _gpi::tuple_of<nr_args, T>;

          static_assert(nr_args == sizeof...(Conts), "Missmatching number of arguments");

          if(DestCont::smallest(dest_cont, conts...) < dest_cont.global_size){
            throw std::logic_error("Can not map a smaller container into a larger one");
          }

          unsigned long this_op_nr = DestCont::max_op(dest_cont, conts...) + 1;
          DestCont::set_op_nr(this_op_nr, dest_cont, conts...);

          T* const dest_ptr = (T*) dest_cont.cont_seg_ptr;

          long lowest = DestCont::lowest_shared_i(conts...);
          long highest = DestCont::highest_shared_i(conts...);

          if(highest < lowest){
            // No overlap
            highest = -1;
            lowest = -1;
          }

          #pragma omp parallel
          {

            int glob_i;
            tup_type tup{};

            // Handle pure local indeces only if they exist
            if(highest != -1 && lowest != -1){

            // Dynamic scheduling as to not overload the thread building the
            // buffer. Another solution is to not give the buffer building
            // thread any work whatsoever and use static scheduling.
              #pragma omp for schedule(static)
              for(int glob_i = lowest; glob_i <= highest; glob_i++){
                int i = glob_i - dest_cont.start_i;

                // The template arg is used to traverse the tupel starting at 0
                build_tuple<0>(glob_i, true, dest_cont, tup, conts...);

                dest_ptr[i] = _gpi::dummy<0, true>::exec(func, tup);
              }
            }

            #pragma omp single
            {
              // TODO build buffer has been modified and needs to be
              // repurposed for this usage
              //dest_cont.build_buffer(true, double{}, conts...);
            }
              // After the above barrier we guarantee that:
            // 1 -  The pure local operations are done
            // 2 - Our read requests have been sent (but not neccesarily finished)
            // 3 - All the ranks we will read are on the same operation as us.


            // Handle non-purely-local indeces only if they exist
            if(!(lowest == dest_cont.start_i && highest == dest_cont.end_i)){


              // The remaining work all use remote information and may require
              // waiting, but static scheduling still seem to have superior
              // performance since the drift between nodes is generally very low.
              #pragma omp for schedule(static)
              for(int glob_i = dest_cont.start_i; glob_i < lowest; glob_i++){
                int i = glob_i - dest_cont.start_i;

                build_tuple<0>(glob_i, false, dest_cont, tup, conts...);

                dest_ptr[i] = _gpi::dummy<0, true>::exec(func, tup);
              }

              #pragma omp for schedule(static)
              for(int glob_i = highest + 1; glob_i <= dest_cont.end_i; glob_i++){
                int i = glob_i - dest_cont.start_i;

                build_tuple<0>(glob_i, false, dest_cont, tup, conts...);

                dest_ptr[i] = _gpi::dummy<0, true>::exec(func, tup);
              }
          }
        } // end of parallel region
      }


      template<typename DestCont, typename ... Args>
      void operator()(DestCont& dest, Args&&... args)
      {

        using arg_0_t = typename std::remove_reference<
          decltype(std::get<0>(std::declval<arg_tup_t>()))>::type;

        const bool case1 = std::is_same<arg_0_t, Index1D>::value &&
          sizeof...(args) == nr_args - 1;

        const bool case2 = !std::is_same<arg_0_t, Index1D>::value &&
          sizeof...(args) == nr_args;

        static_assert(case1 || case2);

        using T = typename DestCont::value_type;


        // We assume that no lambda-call is pure-local, we can deduce such
        // cases by checking if the lambda does not contain any random access
        // iterator. The old variadic solution already have performance
        // improvement regarding this and hence the easiest solution is to
        // update it to allow for indeces and let it handle all calls with no
        // random access iterator.
        //
        // Summarized: Any calls without a random acess interator should be
        // forwarded to the old solution in the future. Here we assume that
        // there exists a random access iterator argument.

        unsigned long this_op_nr = DestCont::max_op(dest, args...) + 1;

        DestCont::set_op_nr(this_op_nr, dest, args...);
        // Reset all proxy caches since they may be invalid
        DestCont::init_proxy_cache(args...);

        // Pre fetch all remote values we know that we want
        _gpi::build_buff_helper<0, arg_tup_t, typename _gpi::is_skepu_proxy_type<arg_0_t>::type>
        ::build(false, double{}, dest, args...);

        #pragma omp parallel
        {
          arg_tup_t tup{};



          #pragma omp for schedule(static)
          for(size_t i = dest.start_i; i <= dest.end_i; ++i){

            _gpi::helper<0, arg_0_t>::build_tuple( i, tup, dest, args...);

            _gpi::apply_helper<T, 0, true>::exec(
              dest.local_buffer + i - dest.start_i, func, tup);
          }


        #pragma omp single
        {
          // Due to the random access iterator we need a barrier,
          // for more discussion read the long comment above
          gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);
        }

        size_t t_step = dest.local_size / omp_get_num_threads();

        size_t size = omp_get_thread_num() == omp_get_num_threads() -1 ?
          t_step + dest.local_size % omp_get_num_threads() :
          t_step;


        std::memcpy(
          (T*) dest.cont_seg_ptr + omp_get_thread_num() * t_step,
          (T*) dest.local_buffer + omp_get_thread_num() * t_step,
           sizeof(T) * size);

        }

        DestCont::set_op_nr(this_op_nr + 1, dest, args...);
      }
  };



  // Template deduction for classes are not allowed in c++11
  // This solves this problem
  template<int nr_args, typename Function>
  auto Map(Function func) ->
    decltype(
      std::declval<
        typename std::remove_reference<Map1D<Function, nr_args>>::type
      >()
    )
  {
    return std::move(Map1D<Function, nr_args>{func});
  }




} // end of namespace skepu
#endif // MAP_HPP
