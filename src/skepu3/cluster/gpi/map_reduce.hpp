#ifndef MAP_RED_HPP
#define MAP_RED_HPP
#include "map.hpp"
#include "reduce.hpp"

namespace skepu{


  // MapReduce inherits from Map and Reduce their respective functions. To
  // signal slightly different behavior MapReduce uses a special constructor
  // for its super classes.
  template<int Arity, typename red_ret_t, typename map_ret_t, typename... map_args>
  class MapReduce1D : public Map1D<Arity, map_ret_t, map_args...>,
    public Reduce1D<red_ret_t, map_ret_t>
    {
    private:
    using map_t = Map1D<Arity, map_ret_t, map_args...>;
    using red_t = Reduce1D<red_ret_t, map_ret_t>;

  public:
    MapReduce1D(
      std::function<map_ret_t(map_args...)> m_func,
      std::function<red_ret_t(map_ret_t, map_ret_t)> r_func
      ) :
      map_t{m_func, true}, red_t{r_func, true} {};


     template<typename DestCont, typename ... Args>
     typename DestCont::value_type operator()(
       DestCont& dest, Args&&... args) {

      dest.wait_for_constraints();
      dest.internal_flush();

       map_t::operator()(dest, dest, args...);
       return red_t::operator()(dest);
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
