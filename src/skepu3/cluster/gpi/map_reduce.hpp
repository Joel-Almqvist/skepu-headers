#ifndef MAP_RED_HPP
#define MAP_RED_HPP
#include "map.hpp"
#include "reduce.hpp"

namespace skepu{


  // MapReduce inherits from Map and Reduce their respective functions. To
  // signal slightly different behavior MapReduce uses a special constructor
  // for its super classes.
  template<int Arity, typename red_f_t, typename map_ret_t, typename... map_args>
  class MapReduce1D : public Map1D<Arity, map_ret_t, map_args...>,
    public Reduce1D<red_f_t>
    {
    private:
    using map_t = Map1D<Arity, map_ret_t, map_args...>;
    using red_t = Reduce1D<red_f_t>;

  public:
    MapReduce1D(std::function<map_ret_t(map_args...)> m_func, red_f_t red_func) :
    map_t{m_func, true},
    Reduce1D<red_f_t>{red_func, true}
     {};


     template<typename DestCont, typename ... Args>
     typename DestCont::value_type operator()(
       DestCont& dest, Args&&... args) {

      dest.wait_for_constraints();
      dest.internal_flush();



      // TODO Add support for initial values in reduce by
      // instance.setStartValue(value )
      // NOTE: Need to add more type deduction in reduce constructor.

       map_t::operator()(dest, dest, args...);
       return red_t::operator()(dest);
     }
  };


  // Handle function pointers
  template<int Arity = -1, typename red_f_t, typename map_ret_t, typename... map_args_t>
  MapReduce1D<Arity, red_f_t, map_ret_t, map_args_t...>
  MapReduce(map_ret_t(*map_args)(map_args_t...), red_f_t red_func){


    return MapReduce1D<Arity, red_f_t, map_ret_t, map_args_t...>{
      (std::function<map_ret_t(map_args_t...)>)map_args, red_func };
  }


  // Handle lambdas and functors
  template<int Arity, typename red_f_t, typename map_ret_t, typename... map_args_t>
  MapReduce1D<Arity, red_f_t, map_ret_t, map_args_t...>
  _map_reduce_helper(std::function<map_ret_t(map_args_t...)> m_func, red_f_t red_func){

    return MapReduce1D<Arity, red_f_t, map_ret_t, map_args_t...>{
      m_func, red_func };
  }

  // Handle lambdas and functors
  template<int Arity = -1, typename red_f_t, typename lambda_t>
  auto MapReduce(lambda_t&& lambda, red_f_t red_func) ->
    decltype(_map_reduce_helper<Arity>(lambda_cast(lambda), red_func)){

    return _map_reduce_helper<Arity>(lambda_cast(lambda), red_func);
  }

}

#endif //MAP_RED_HPP
