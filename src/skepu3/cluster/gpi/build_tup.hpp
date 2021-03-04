#ifndef BUILD_TUP_HPP
#define BUILD_TUP_HPP

#include <type_traits>
#include "matrix.hpp"


/* This file contains help functions for Map.
*/

struct Particle;

namespace skepu{
  namespace _gpi{


    struct build_tup_util{

      // Random access iterator type case
      template<int tup_arg_ctr, typename Tup, typename Dest, typename Curr, typename... Rest>
      static auto build_tuple_helper(double sfinae_param, size_t i, Tup& tup,
        Dest& dest, Matrix<Curr>& curr, Rest&... rest)
        -> decltype(
          std::declval<
            typename std::remove_reference<decltype(std::get<tup_arg_ctr>(tup))>::type::is_proxy_type>(),
        std::declval<void>())
      {

        std::get<tup_arg_ctr>(tup) = std::move(Vec<Curr>(curr));
      }

      // Index argument case
      template<int tup_arg_ctr, typename Tup, typename Dest, typename Curr, typename... Rest>
      static auto build_tuple_helper(double sfinae_param, size_t i, Tup& tup,
        Dest& dest, Matrix<Curr>& curr, Rest&... rest)
        -> decltype(
          std::declval<
            typename std::remove_reference<decltype(std::get<tup_arg_ctr>(tup))>::type::is_skepu_index>(),
        std::declval<void>())
      {

        std::get<tup_arg_ctr>(tup) = Index1D{i};

      }


      // Scalar value from container case
      template<int tup_arg_ctr, typename Tup, typename Dest, typename Curr, typename... Rest>
      static auto build_tuple_helper(int sfinae_param, size_t i, Tup& tup,
        Dest& dest, Matrix<Curr>& curr, Rest&... rest) -> decltype(std::declval<void>())
      {
        using T = typename Matrix<Curr>::value_type;

        T* val_ptr = (T*) curr.cont_seg_ptr;
        T* comm_ptr = (T*) curr.comm_seg_ptr;


        if(curr.start_i > i){
          std::get<tup_arg_ctr>(tup) = comm_ptr[i - dest.start_i];
        }

        else if (i > curr.end_i){
          // Every added elem increases the offset by one, and we want offset 0
          // to be possible
          int offset = -1;

          if(curr.start_i > dest.start_i){
            // We have sent lower missing chunks
            offset += curr.start_i - dest.start_i;
          }

          offset += i - curr.end_i;
          std::get<tup_arg_ctr>(tup) = comm_ptr[offset];
        }

        else{
          std::get<tup_arg_ctr>(tup) = val_ptr[i - curr.start_i];
        }
      }

      // Constant value case
      template<int tup_arg_ctr, typename Tup, typename Dest, typename Curr_scalar>
      static void build_tuple_helper(int sfinae_param, size_t i, Tup& tup, Dest& dest, Curr_scalar& curr)
      {
        std::get<tup_arg_ctr>(tup) = curr;
      }



    };

    template<int tup_arg_ctr, typename T>
    struct helper{

      // Traverses through the argument list and calls a helper for function on
      // every argument.
      template<typename Tup, typename Dest, typename Curr, typename... Rest>
      static void build_tuple( size_t i, Tup& tup, Dest& dest,
        Curr& curr, Rest&... rest)
      {
        build_tup_util::build_tuple_helper<tup_arg_ctr>(double{}, i, tup, dest, curr);
        helper<tup_arg_ctr + 1, T>::build_tuple( i, tup, dest, rest...);
      }

      template<typename Tup, typename Dest, typename Curr>
      static void build_tuple( size_t i, Tup& tup, Dest& dest, Curr& curr)
      {
        build_tup_util::build_tuple_helper<tup_arg_ctr>(double{}, i, tup, dest, curr);
      }


    };

    template<>
    struct helper<0, Index1D>{


      // Traverses through the argument list and calls a helper for function on
      // every argument.
      template<typename Tup, typename Dest, typename Curr, typename... Rest>
      static void build_tuple( size_t i, Tup& tup, Dest& dest,
        Curr& curr, Rest&... rest)
      {
        // We won't use the curr argument in this call but we want another
        // skepu::Matrix for overload resoultion so send in dest twice.
        build_tup_util::build_tuple_helper<0>(double{}, i, tup, dest, dest);

        helper<1, int>::build_tuple( i, tup, dest, curr, rest...);
      }


      template<typename Tup, typename Dest, typename... Rest>
      static void build_tuple( size_t i, Tup& tup, Dest& dest)
      {
        // Throw in dest as argument twice to help with SFINAE
        build_tup_util::build_tuple_helper<0>(double{}, i, tup, dest, dest);
      }


    };


    template <typename T, int ctr, bool done>
    struct apply_helper;

    template <typename T, int ctr>
    struct apply_helper<T, ctr, true> {

        template <typename Func, typename...Args, typename... Exp>
        static void exec(T* t, Func func, std::tuple<Args...>& tup, Exp&&... exp)
         {

          const bool not_done = ctr < sizeof...(Args) - 1;
          apply_helper<T, ctr + 1, not_done>::exec(t, func, tup, exp..., std::get<ctr>(tup));
        }
    };

    template <typename T, int ctr>
    struct apply_helper<T, ctr, false> {

        template<typename Func, typename Tup, typename...Exp>
        static void exec(T* t, Func func, Tup& tup, Exp&&... exp)
        {
          *t = func(exp...);
          }
    };


    // The current argument type is a proxy and we never want to call build_buffer
    template<int ctr, typename Tup, typename T>
    struct build_buff_helper{

      template<typename Dest, typename First, typename ... Rest>
      static void build(bool no_wait, int SFINAE_param, Dest& dest,
          First& first, Rest&... rest){

          using NextT = decltype(std::get<ctr + 1>(std::declval<Tup>()));

          build_buff_helper<ctr + 1, Tup, typename is_skepu_proxy_type<NextT>::type>
          ::build(no_wait, double{}, dest, rest...);
      }

      // Sink
      template<typename Dest, typename First>
      static void build(bool no_wait, int SFINAE_param, Dest& dest,
          Matrix<First>& first){}


      // Sink
      template<typename Dest>
      static void build(bool no_wait, double SFINAE_param, Dest& dest){}


    };

    // The current lambda argument is not a proxy type, but we need to make sure
    // that it also is not a constant.
    template<int ctr, typename Tup>
    struct build_buff_helper<ctr, Tup, std::false_type>{

      template<typename Dest, typename First, typename ... Rest>
      static void build(bool no_wait, double SFINAE_param, Dest& dest,
          Matrix<First>& first, Rest&... rest){

        dest.build_buffer_helper(no_wait, first);

        using NextT = decltype(std::get<ctr + 1>(std::declval<Tup>()));

        build_buff_helper<ctr + 1, Tup, typename is_skepu_proxy_type<NextT>::type>
        ::build(no_wait, double{}, dest, rest...);
      }


      // Not a proxy and not a matrix -> constant
      template<typename Dest, typename First, typename ... Rest>
      static void build(bool no_wait, double SFINAE_param, Dest& dest,
          First&& first, Rest&... rest){

        using NextT = decltype(std::get<ctr + 1>(std::declval<Tup>()));
        build_buff_helper<ctr + 1, Tup, typename is_skepu_proxy_type<NextT>::type>
        ::build(no_wait, double{}, dest, rest...);
      }

      template<typename Dest, typename First>
      static void build(bool no_wait, double SFINAE_param, Dest& dest,
          Matrix<First>& first){
        dest.build_buffer_helper(no_wait, first);
      }

      // Not a proxy and not a matrix -> constant
      template<typename Dest, typename First>
      static void build(bool no_wait, double SFINAE_param, Dest& dest,
          First&& first){}

      // Sink
      template<typename Dest>
      static void build(bool no_wait, double SFINAE_param, Dest& dest){}



    };
  } // end of namespace _gpi
}

#endif // BUILD_TUP_HPP
