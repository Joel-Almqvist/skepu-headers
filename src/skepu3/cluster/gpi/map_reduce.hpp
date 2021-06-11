#ifndef MAP_RED_HPP
#define MAP_RED_HPP
#include "map.hpp"
#include "reduce.hpp"

namespace skepu{


  // MapReduce inherits from Map and Reduce their respective functions. To
  // signal slightly different behavior MapReduce uses a special constructor
  // for its super classes.
  template<int Arity, typename red_ret_t, typename map_ret_t, typename... map_args>
  class MapReduce1D {
    private:

    std::function<map_ret_t(map_args...)> map_func;
    using arg_tup_t = typename std::tuple<map_args...>;
    static const int nr_args = sizeof...(map_args);


    std::function<red_ret_t(map_ret_t, map_ret_t)> red_func;
    red_ret_t start_value;

  public:
    MapReduce1D(
      std::function<map_ret_t(map_args...)> m_func,
      std::function<red_ret_t(map_ret_t, map_ret_t)> r_func
      ) :
      map_func{m_func},
      red_func{r_func},
      start_value{}
      {};

      void setStartValue(red_ret_t val){
        start_value = val;
      }



     template<typename FirstCont, typename ... Args>
     typename FirstCont::value_type operator()(
       FirstCont& first, Args&&... args) {


       using arg_0_t = typename std::remove_reference<
         decltype(std::get<0>(std::declval<arg_tup_t>()))>::type;

       const bool case1 = std::is_same<arg_0_t, Index1D>::value &&
         (1 + sizeof...(args)) == nr_args - 1;

       const bool case2 = std::is_same<arg_0_t, Index2D>::value &&
         (1 + sizeof...(args)) == nr_args - 1;


       const bool case3 = !std::is_same<arg_0_t, Index1D>::value &&
         !std::is_same<arg_0_t, Index2D>::value &&
         (1 + sizeof...(args)) == nr_args;

       static_assert(case1 || case2 || case3);


       // NOTE: Here we must flush old changes and wait if we allow random access
       // since we might access any value from any container.
       first.get_constraints(int{}, args...);
       first.wait_for_constraints();
       FirstCont::flush_rest(int{}, first, args...);


       first.vclock[first.rank] = ++first.op_nr;

       // Pre fetch all remote values we know that we want
       _gpi::build_buffer_util<nr_args, 0, arg_tup_t>
       ::build(
         int{}, int{}, int{}, int{}, int{},
         Index1D{}, first, first, args...);


       map_ret_t* all_accs;
       size_t nthreads;

       #pragma omp parallel
       {

         #pragma omp single
         {

           nthreads = omp_get_num_threads();
           all_accs = (map_ret_t*) malloc(sizeof(map_ret_t) * nthreads);

         }

         map_ret_t& acc = all_accs[omp_get_thread_num()];
         map_ret_t temp;
         size_t trank = omp_get_thread_num();

         bool first_run = true;

         arg_tup_t tup{};

         #pragma omp for schedule(static)
         for(size_t i = first.start_i; i <= first.end_i; ++i){

           _gpi::helper<0, arg_0_t>::build_tuple( i, tup, first, first, args...);


           if(first_run){
             _gpi::apply_helper<red_ret_t, 0, true>::exec(
               &acc, map_func, tup);
             first_run = false;


             if(first.rank == 0 && omp_get_thread_num() == 0){
               acc = red_func(acc, start_value);
             }

           }
           else{
             _gpi::apply_helper<red_ret_t, 0, true>::exec(
               &temp, map_func, tup);

              acc = red_func(acc, temp);

           }
         }
       } // end of parallel region

       red_ret_t local_sum = all_accs[0];

       for(int i = 1; i < nthreads; i++){
         local_sum = red_func(local_sum, all_accs[i]);
       }

       // Indicate to other ranks that our value is ready to be used
       first.vclock[first.rank] = ++first.op_nr;

       int iterations = std::ceil(std::log2(first.nr_nodes));
       int dest_rank;
       int step;
       int prev_step;
       unsigned remote_offset;

       using T = map_ret_t;

       T& loc_val = *((T*) first.comm_seg_ptr +
         first.norm_partition_size * first.rank);

       loc_val = local_sum;

       T& rec_val = *((T*) first.comm_seg_ptr +
         first.norm_partition_size * first.rank + 1);


         // A distributed gather which in the end puts the data to node 0
         for(int i = 0; i < iterations; i++){

           step = pow(2, i + 1);
           prev_step = pow(2, i);

             if(first.rank % step == 0){

               dest_rank = first.rank + prev_step;

               if(dest_rank < first.nr_nodes){
                 first.wait_ranks.push_back(dest_rank);
                 first.wait_for_vclocks(first.op_nr);

                 if(dest_rank == first.nr_nodes - 1){
                   remote_offset = first.last_partition_comm_offset +
                     sizeof(T) * first.norm_partition_size * dest_rank;
                 }
                 else{
                   remote_offset = first.norm_partition_comm_offset +
                     sizeof(T) * first.norm_partition_size * dest_rank;
                 }

                 gaspi_read(
                   first.segment_id,
                   first.comm_offset + (first.start_i + 1) * sizeof(T), // local offset
                   dest_rank,
                   first.segment_id - first.rank + dest_rank, // rem seg id
                   remote_offset,
                   sizeof(T), //size
                   first.queue,
                   GASPI_BLOCK
                 );

                 gaspi_wait(first.queue, GASPI_BLOCK);

                 loc_val = red_func(loc_val, rec_val);

               }
             }
             first.vclock[first.rank] = ++first.op_nr;
         }

         // A distributed broadcast from node 0 to all the other.
         for(int i = 0; i < iterations; i++){

           step = pow(2, i + 1);
           prev_step = pow(2, i);

           if(first.rank >= prev_step && first.rank < step){

             dest_rank = first.rank - prev_step;

             first.wait_ranks.push_back(dest_rank);
             first.wait_for_vclocks(first.op_nr);


             if(dest_rank == first.nr_nodes - 1){
               remote_offset = first.last_partition_comm_offset +
                 sizeof(T) * first.norm_partition_size * dest_rank;
             }
             else{
               remote_offset = first.norm_partition_comm_offset +
                 sizeof(T) * first.norm_partition_size * dest_rank;
             }

             gaspi_read(
               first.segment_id,
               first.comm_offset + first.start_i * sizeof(T), // local offset
               dest_rank,
               first.segment_id - first.rank + dest_rank, // rem seg id
               remote_offset,
               sizeof(T), //size
               first.queue,
               GASPI_BLOCK
             );

             gaspi_wait(first.queue, GASPI_BLOCK);

           }
           first.vclock[first.rank] = ++first.op_nr;
         }

       return loc_val;
     }
  };


  // The 6 functions below are used to instantiate a MapReduce object. First
  // the map argument is converted, secondly the reduce argument.
  //
  // Reduce arguments conversion:
  // Second step of reduce lambda handling
  template<int Arity = -1, typename red_ret_t, typename map_ret_t, typename... map_args_t>
  MapReduce1D<Arity, red_ret_t, map_ret_t, map_args_t...> _mr_red_deducer(
    std::function<map_ret_t(map_args_t...)> m_func,
    std::function<red_ret_t(map_ret_t, map_ret_t)> r_func
   ){

     return MapReduce1D<Arity, red_ret_t, map_ret_t, map_args_t...>{m_func, r_func};

  }

  // Handles reduce function pointers
  template<int Arity = -1, typename red_ret_t, typename map_ret_t, typename... map_args_t>
  MapReduce1D<Arity, red_ret_t, map_ret_t, map_args_t...> _mr_red_deducer(
    std::function<map_ret_t(map_args_t...)> m_func,
    red_ret_t(*arg)(map_ret_t, map_ret_t)
   ){

     return MapReduce1D<Arity, red_ret_t, map_ret_t, map_args_t...>{
        m_func, (std::function<red_ret_t(map_ret_t, map_ret_t)>)arg};

  }


  // First step of reduce lambda handling
  template<int Arity = -1, typename lambda_t, typename map_ret_t, typename... map_args_t>
  auto _mr_red_deducer(
    std::function<map_ret_t(map_args_t...)> m_func,
    lambda_t lambda
  ) -> decltype(_mr_red_deducer(m_func, lambda_cast(lambda))){

      return _mr_red_deducer(m_func, lambda_cast(lambda));
  }


  // Map argument conversions
  //
  // Handles Map function pointers
  template<int Arity = -1, typename red_f_t, typename map_ret_t, typename... map_args_t>
  auto MapReduce(
    map_ret_t(*map_args)(map_args_t...),
    red_f_t red_func
  ) -> decltype(_mr_red_deducer(
      (std::function<map_ret_t(map_args_t...)>) map_args, red_func)
      ){

    return _mr_red_deducer(
      (std::function<map_ret_t(map_args_t...)>) map_args,
      red_func);
  }


  // Second step of handling Map lambdas and functors
  template<int Arity, typename red_f_t, typename map_ret_t, typename... map_args_t>
  auto _mr_lambda_deducer(
    std::function<map_ret_t(map_args_t...)> m_func,
    red_f_t red_func
  ) -> decltype(_mr_red_deducer(m_func, red_func)) {

    return _mr_red_deducer(m_func, red_func);

  }

  // First step of handling Map lambdas and functors
  template<int Arity = -1, typename red_f_t, typename lambda_t>
  auto MapReduce(lambda_t&& lambda, red_f_t red_func) ->
    decltype(_mr_lambda_deducer<Arity>(lambda_cast(lambda), red_func)){

    return _mr_lambda_deducer<Arity>(lambda_cast(lambda), red_func);
  }

}

#endif //MAP_RED_HPP
