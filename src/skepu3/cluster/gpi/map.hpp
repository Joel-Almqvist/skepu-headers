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
#include "proxy.hpp"
#include "build_tup.hpp"


namespace skepu{



  template<int Arity, typename Ret, typename... Func_args>
  class Map1D{
  private:
    std::function<Ret(Func_args...)> func;

    static const int nr_args = sizeof...(Func_args);

    const bool uses_random_access;

    //using arg_tup_t = typename _gpi::get_tup_t<Function, nr_args - 1>::type;

    using arg_tup_t = typename std::tuple<Func_args...>;


  public:
    Map1D(std::function<Ret(Func_args...)> func) : func{func},
    uses_random_access{has_random_access<Func_args...>()} {};


  private:
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

          // These functions have been removed
          //unsigned long this_op_nr = DestCont::max_op(dest_cont, conts...) + 1;
          //DestCont::set_op_nr(this_op_nr, dest_cont, conts...);

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


      // Work in progress
      template<typename DestCont, typename ... Args>
      void apply(DestCont& dest, Args&&... args)
      {
        if(uses_random_access){
          std::cout << "Applying assuming atleast one random access iterator\n";
        }
        else{
          std::cout << "Applying without random access iterators\n";
        }
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


        dest.get_constraints(args...);
        // TODO here we must flush if needed
        dest.wait_for_constraints();

        dest.op_nr++;

        // Do not flush if this is the first operation
        if(dest.state != 0){
          dest.flush();
        }

        DestCont::flush_rest(args...);

        dest.vclock[dest.rank] = dest.op_nr;

        dest.state = dest.op_nr;

        /*
        max_op_nr = ...

        set_all_containers_state(max_op_nr);
        */

        /* Det enda värdet i c1 som är har en betydelse i c2 är op_nr för just
        * denna operation. Nej eftersom alla noder gör arbete på alla operationer


        TODO Ändra op_nr till en global variabel i container klassen.
        Nu har vi en global ordning av alla operationer, dvs alla vclocks kan
        jämföras med varandra.

        Ändra vclock till en global array som alla kontainrar delar.


        SEN så kan vi göra:

        1 - Hitta de starkaste beronden i alla requests

        2 - Jämför med vclock och vänta ifall det behövs

        3 - Flush

        4 - Sätt vclock[rank] till nuvarande op_nr

        5 - Gör operationen... vänta på remotes vi vill läsa ifrån, spara
        resultatet i en temp buff.

        6 - Lägg till beroenden, dvs alla remotes som vill eller kanske vill
        läsa ifrån oss läggs till i dest kontainern ENBART.

        Så i icke RAO fallet är det våra "grannar", i RAO fallet är det alla

        ----

        I praktiken så väntar vi då i början på att alla remotes ska ha läst av
        oss ifall det behövs.

        Under OP väntar vi på att remote ska göra förändringar som vi vill läsa.

         Och i slutet så lägger vi till vilka remotes som vill läsa av oss.


         Skillnaden blir först bugfix för op_nr och sen att vi inte måste vänta
         till att operationen ska nå JUST denna op_nr utan vi väntar enbart på
         att beroenderna ska lösas.

        */



        /*
        op_nr++;
        update_state(op_nr);

        flush();
        state = op_nr;
        vclock[rank] = state;

        // NOTE We assume random access operator

        for(int i = 0; i < nr_nodes; i++){
          if(i != rank){
            state_tracker[i] = op_nr;
          }
        }
        */

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


        // TODO write over changes if they exist before incrementing op_nr

        // TODO In non RAO we must wait for nodes which we know read from us

        // Idea 1 - Save which op_nr was the last random access operation and
        // wait for every node to reach atleast that state.


        /* Idea 2:

        * create state_nr which indicates which state the container chunk is in
        * states are indicated by op_nr
        * Ie (op_nr, state_nr), (12, 9)

        * wait_for_vclock and similar functions now reads the state_nr

        * Increment op_nr at the end of every operation

        * At the start of every operation where a remote reads our local value:
        *
        * 1 - Check the requests, read the op_nr of all requests and wait if
        * needed to guarantee that they have read our local value. We now know
        * that no ranks wants what is inside the gaspi_segment.
        *
        * 2 - Write over the changes, update state_nr


        * At the start of every operation WHICH MODIFIES THE LOCAL BUFFER
        *
        * 1 - Check the requests, read the op_nr of all requests and wait if
        * needed to guarantee that they have read our local value. We now know
        * that no ranks wants what is inside the gaspi_segment.
        *
        * 2 - Write over the changes, update state_nr
        *
        * 3 - Add which ever remote ranks might want to read from us (for example
        *  in map with a random access index that would be all ranks).
        *
        * 4 - Perform work and put it in the buffer, update op_nr.


         * At every operation calculate which ranks want to read my data
         *
        */
        // [(rank, op_nr), (rank, op_nr) ...]




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

        // Due to the random access we have constraints to all nodes
        for(int i = 0; i < dest.nr_nodes; i++){
          if(i == dest.rank){
            continue;
          }
          dest.constraints[i] = dest.op_nr;
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
