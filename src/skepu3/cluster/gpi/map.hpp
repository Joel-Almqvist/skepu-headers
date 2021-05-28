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
#include <functional>

#include <omp.h>
#include <GASPI.h>

#include "matrix.hpp"
#include "utils.hpp"
#include "proxy_vec.hpp"
#include "proxy_mat.hpp"
#include "build_tup.hpp"


namespace skepu{



  template<int Arity, typename Ret, typename... Func_args>
  class Map1D{
  private:
    std::function<Ret(Func_args...)> func;

    static const int nr_args = sizeof...(Func_args);

    const bool uses_random_access;

    using arg_tup_t = typename std::tuple<Func_args...>;


  public:
    Map1D(std::function<Ret(Func_args...)> func) : func{func},
    uses_random_access{has_random_access<Func_args...>()} {};


    // Dummy
    template<typename T>
    void setBackend(T){};


  private:
    template<int ctr, typename Dest, typename Tup, typename Curr, typename... Rest>
    void build_tuple(int i, bool local_only, Dest& dest,
       Tup& tup, Curr& curr, Rest&... rest){

      if(local_only){
        build_tuple_local<ctr>(i, dest, tup, curr);

        // Early break if the tuple need a remote value
          build_tuple<ctr + 1>(i, local_only, dest, tup, rest...);
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



    template<typename Curr, typename... Rest>
    static bool has_random_access(){
      return has_random_access(int{}, std::tuple<Curr, Rest...>{});
    }


    template<typename Curr, typename... Rest>
    static auto has_random_access(int sfinae, std::tuple<Curr, Rest...>)
      -> decltype((typename Curr::is_proxy_type){}, true)
    {
      return true;
    }


    template<typename Curr, typename... Rest>
    static bool has_random_access(long sfinae, std::tuple<Curr, Rest...>) {

      return has_random_access(int{}, std::tuple<Rest...>{});
    }


    static bool has_random_access(int sfinae, std::tuple<>) {
      return false;
    }



  public:

     // Should take in a backend type
     void setBackend(){}

     // Need to be implemented
     void setReduceMode(){};



      template<typename DestCont, typename ... Args>
      void operator()(DestCont& dest, Args&&... args)
      {

        using arg_0_t = typename std::remove_reference<
          decltype(std::get<0>(std::declval<arg_tup_t>()))>::type;

        const bool case1 = std::is_same<arg_0_t, Index1D>::value &&
          sizeof...(args) == nr_args - 1;

        const bool case2 = std::is_same<arg_0_t, Index2D>::value &&
          sizeof...(args) == nr_args - 1;


        const bool case3 = !std::is_same<arg_0_t, Index1D>::value &&
          !std::is_same<arg_0_t, Index2D>::value &&
          sizeof...(args) == nr_args;

        static_assert(case1 || case2 || case3);

        using T = typename DestCont::value_type;


        dest.get_constraints(int{}, args...);
        dest.wait_for_constraints();



        DestCont::flush_rest(int{}, dest, args...);

        // Increment after the flush to guarantee the reads to be safe
        dest.op_nr++;

        // Remember that these variables are globally shared among containers
        dest.vclock[dest.rank] = dest.op_nr;


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
        }

        dest.last_mod_op = dest.op_nr;
        dest.vclock[dest.rank] = ++dest.op_nr;


        // Add constraints
        if(uses_random_access){
          // Due to the random access we have constraints to all nodes
          for(int i = 0; i < dest.nr_nodes; i++){
            if(i == dest.rank){
              continue;
            }
            dest.constraints[i] = dest.op_nr;
          }
        }

        else{

          long unsigned lowest_i;
          long unsigned highest_i;

          int lowest_rank;
          int highest_rank;

          // Calculate which nodes every remote is going to contact and create
          // constraints to any node which will read from us.
          for(int i = 0; i < dest.nr_nodes; i++){

            if(i == dest.rank)
              continue;

            lowest_i = i * dest.norm_partition_size;

            if(i != dest.nr_nodes - 1 ){
              highest_i = (i + 1) * dest.norm_partition_size;
            }
            else{
              highest_i = dest.global_size - 1;
            }

            lowest_rank = DestCont::find_rank_overlap(
              int{},
              i,
              lowest_i,
              dest.nr_nodes, // start the accumulator at impossible value
              true,
              dest, args...);

            highest_rank = DestCont::find_rank_overlap(
              int{},
              i,
              lowest_i,
              -1, // start the accumulator at impossible value
              false,
              dest, args...);

            if(lowest_rank <= dest.rank && dest.rank <= highest_rank){
              dest.constraints[i] = dest.op_nr;
            }
          }
        }
      }
  };


  // The constructor wrapper for function pointers
  template <int Arity = -1, typename Ret, typename... Args>
  Map1D<Arity, Ret, Args...> Map(Ret(*map_arg)(Args...)){
      return Map1D<Arity, Ret, Args...>((std::function<Ret(Args...)>)map_arg);
  }


  // Helper for the lambda and functor constructor wrapper
  template <int Arity = -1, typename Ret, typename... Args>
  Map1D<Arity, Ret, Args...> _map_helper(std::function<Ret(Args...)> map_arg){
    return Map1D<Arity, Ret, Args...>(map_arg);
  }


  // The constructor wrapper for lambda and functors
  template <int Arity = -1, typename Lambda>
  auto Map(Lambda&& lambda) -> decltype(_map_helper<Arity>(lambda_cast(lambda))){

    return _map_helper<Arity>(lambda_cast(lambda));
  }



} // end of namespace skepu
#endif // MAP_HPP
