#ifndef BUILD_TUP_HPP
#define BUILD_TUP_HPP

#include <type_traits>
#include "matrix.hpp"



namespace skepu{
  namespace _gpi{

    struct build_tup_util{

      // Random access iterator type case
      template<int tup_arg_ctr, typename Tup, typename Dest, typename Curr>
      static auto build_tuple_helper(Vec<Curr> sfinae, size_t i, Tup& tup,
        Dest& dest, Matrix<Curr>& curr)
        -> decltype(
        std::declval<void>())
      {

        std::get<tup_arg_ctr>(tup) = std::move(Vec<Curr>(curr));
      }


      // Random access iterator type case
      template<int tup_arg_ctr, typename Tup, typename Dest, typename Curr>
      static auto build_tuple_helper(Mat<Curr> sfinae, size_t i, Tup& tup,
        Dest& dest, Matrix<Curr>& curr)
        -> decltype(
        std::declval<void>())
      {

        std::get<tup_arg_ctr>(tup) = std::move(Mat<Curr>(curr));
      }




      // Index argument case
      template<int tup_arg_ctr, typename Tup, typename Dest, typename Curr>
      static auto build_tuple_helper(Index1D sfinae_param, size_t i, Tup& tup,
        Dest& dest, Curr& curr)
        -> decltype(
        std::declval<void>())
      {

        std::get<tup_arg_ctr>(tup) = Index1D{i};

      }


      // Index argument case
      template<int tup_arg_ctr, typename Tup, typename Dest, typename Curr>
      static auto build_tuple_helper(Index2D sfinae_param, size_t i, Tup& tup,
        Dest& dest, Curr& curr)
        -> decltype(
        std::declval<void>())
      {

        std::get<tup_arg_ctr>(tup) = Index2D{dest.get_row(i), dest.get_col(i)};

      }



      // Scalar value from container case
      template<int tup_arg_ctr, typename Tup, typename Dest, typename Curr>
      static auto build_tuple_helper(typename Matrix<Curr>::value_type sfinae_param, size_t i, Tup& tup,
        Dest& dest, Matrix<Curr>& curr) -> decltype(std::declval<void>())
      {
        using T = typename Matrix<Curr>::value_type;

        T* val_ptr = (T*) curr.cont_seg_ptr;
        T* swap_ptr = ((T*) curr.comm_seg_ptr) + curr.rank * curr.norm_partition_size;

        if(curr.start_i > i){
          std::get<tup_arg_ctr>(tup) = swap_ptr[i - dest.start_i];
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
          std::get<tup_arg_ctr>(tup) = swap_ptr[offset];
        }

        else{
          std::get<tup_arg_ctr>(tup) = val_ptr[i - curr.start_i];
        }
      }

      // Constant value case
      template<int tup_arg_ctr, typename Tup, typename Dest, typename Curr_scalar, typename SFINAE>
      static void build_tuple_helper(SFINAE sfinae_param, size_t i, Tup& tup, Dest& dest, Curr_scalar& curr)
      {
        std::get<tup_arg_ctr>(tup) = curr;
      }



    };


    // This struct iterates two ways, once through the input arguments to map
    // and once through the argument tuple to the embedded function itself.
    //
    // The specialization are needed since during an index arg we step once in
    // the embedded function chain but not in the map arg chain.
    template<int tup_arg_ctr, typename T>
    struct helper{

      // Traverses through the argument list and calls a helper for function on
      // every argument.
      template<typename Tup, typename Dest, typename Curr, typename... Rest>
      static void build_tuple( size_t i, Tup& tup, Dest& dest,
        Curr& curr, Rest&... rest)
      {

        using curr_arg_t = typename std::remove_reference<decltype(std::get<tup_arg_ctr>(tup))>::type;


        build_tup_util::build_tuple_helper<tup_arg_ctr>(curr_arg_t{}, i, tup, dest, curr);


        helper<tup_arg_ctr + 1, T>::build_tuple( i, tup, dest, rest...);
      }

      template<typename Tup, typename Dest, typename Curr>
      static void build_tuple( size_t i, Tup& tup, Dest& dest, Curr& curr)
      {

        using curr_arg_t = typename std::remove_reference<decltype(std::get<tup_arg_ctr>(tup))>::type;

        build_tup_util::build_tuple_helper<tup_arg_ctr>(curr_arg_t{}, i, tup, dest, curr);
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

        using curr_arg_t = typename std::remove_reference<decltype(std::get<0>(tup))>::type;


        // We won't use the curr argument in this call but we want another
        // skepu::Matrix for overload resoultion so send in dest twice.
        build_tup_util::build_tuple_helper<0>(curr_arg_t{}, i, tup, dest, dest);

        // Reuse curr
        helper<1, int>::build_tuple( i, tup, dest, curr, rest...);
      }


      template<typename Tup, typename Dest, typename... Rest>
      static void build_tuple( size_t i, Tup& tup, Dest& dest)
      {
        using curr_arg_t = typename std::remove_reference<decltype(std::get<0>(tup))>::type;

        // Throw in dest as argument twice to help with SFINAE
        build_tup_util::build_tuple_helper<0>(curr_arg_t{}, i, tup, dest, dest);
      }
    };


    // This is a copy of the case above, ugly solution but it works
    template<>
    struct helper<0, Index2D>{


      // Traverses through the argument list and calls a helper for function on
      // every argument.
      template<typename Tup, typename Dest, typename Curr, typename... Rest>
      static void build_tuple( size_t i, Tup& tup, Dest& dest,
        Curr& curr, Rest&... rest)
      {

        using curr_arg_t = typename std::remove_reference<decltype(std::get<0>(tup))>::type;

        // We won't use the curr argument in this call but we want another
        // skepu::Matrix for overload resoultion so send in dest twice.
        build_tup_util::build_tuple_helper<0>(curr_arg_t{}, i, tup, dest, dest);

        // Reuse
        helper<1, int>::build_tuple( i, tup, dest, curr, rest...);
      }


      template<typename Tup, typename Dest>
      static void build_tuple( size_t i, Tup& tup, Dest& dest)
      {

        using curr_arg_t = typename std::remove_reference<decltype(std::get<0>(tup))>::type;

        // Throw in dest as argument twice to help with SFINAE
        build_tup_util::build_tuple_helper<0>(curr_arg_t{}, i, tup, dest, dest);
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


    // Helps with compilation by adding a out of bounds check for the tuple type.
    template<int ctr, typename Tup, bool>
    struct next_t_helper{};

    template<int ctr, typename Tup>
    struct next_t_helper<ctr, Tup, true>{
      using type = typename std::remove_reference
      <
      decltype(std::get<ctr>(std::declval<Tup>()))
      >::type;

    };

    template<int ctr, typename Tup>
    struct next_t_helper<ctr, Tup, false>{
      using type = std::false_type;

    };


    template<int nr_args, int ctr, typename Tup>
    struct build_buffer_util{

      /*
      * This helper struct is used to decide which build to call. There are 6
      * cases which are sorted from highest SFINAE priority to lowest. The
      * functions attempts to deduce whether a argument which may or may not be
      * a container is using the random access pattern or not. If not the
      * function calls build buffer on the container.
      *
      * Note that the only case which in fact calls build buffer is case 5.
      *
      * Note that an index argument makes us traverse one step in the argument
      * types chain but no step in the variadic types chain.
      *
      */


      // Case 1 A
      template<typename Dest, typename First, typename... Rest>
      static void build(
        int, int, int, int, int,
        Index1D, Dest& dest, First& first, Rest&... rest){
        // Reuse first for recursive call


        using next_t = typename next_t_helper<ctr, Tup, (nr_args - 1) == ctr>::type;


        build_buffer_util<nr_args, ctr + 1, Tup>::build(
          int{}, int{}, int{}, int{}, int{},
          next_t{}, dest, first, rest...
        );
      }

      // Case 1 B
      template<typename Dest, typename First>
      static void build(
        int, int, int, int, int,
        Index1D, Dest& dest, First& first){
      }


      // Case 2 A
      template<typename Dest, typename First, typename... Rest>
      static void build(
        int, int, int, int, long,
        Index2D, Dest& dest, First& first, Rest&... rest){
        // Reuse first for recursive call

        using next_t = typename next_t_helper<ctr, Tup, (nr_args - 1) == ctr>::type;

        build_buffer_util<nr_args, ctr + 1, Tup>::build(
          int{}, int{}, int{}, int{}, int{},
          next_t{}, dest, first, rest...
        );

      }

      // Case 2 B
      template<typename Dest, typename First>
      static void build(
        int, int, int, int, long,
        Index2D, Dest& dest, First& first){
      }


      // Case 3 A
      template<typename T, typename Dest, typename First, typename... Rest>
      static void build(
        int, int, int, long, long,
        Vec<T>, Dest& dest, Matrix<First>& first, Rest&... rest){
        // Only recursive call

        using next_t = typename next_t_helper<ctr, Tup, (nr_args - 1) == ctr>::type;

        build_buffer_util<nr_args, ctr + 1, Tup>::build(
          int{}, int{}, int{}, int{}, int{},
          next_t{}, dest, rest...
        );

      }

      // Case 3 B
      template<typename T, typename Dest, typename First>
      static void build(
        int, int, int, long, long,
        Vec<T>, Dest& dest, Matrix<First>& first){
      }

      // Case 4 A
      template<typename T, typename Dest, typename First, typename... Rest>
      static void build(
        int, int, long, long, long,
        Mat<T>, Dest& dest, Matrix<First>& first, Rest&... rest){
        // Only recursive call

        using next_t = typename next_t_helper<ctr, Tup, (nr_args - 1) == ctr>::type;

        build_buffer_util<nr_args, ctr + 1, Tup>::build(
          int{}, int{}, int{}, int{}, int{},
          next_t{}, dest, rest...
        );

      }

      // Case 4 B
      template<typename T, typename Dest, typename First>
      static void build(
        int, int, long, long, long,
        Mat<T>, Dest& dest, Matrix<First>& first){
      }


      // Case 5 A
      // This case is the only one which actually calls build buffer
      template<typename T, typename Dest, typename First, typename... Rest>
      static void build(
        int, long, long, long, long,
        T, Dest& dest, Matrix<First>& first, Rest&... rest){

        using next_t = typename next_t_helper<ctr, Tup, (nr_args - 1) == ctr>::type;

        build_buffer_util<nr_args, ctr + 1, Tup>::build(
          int{}, int{}, int{}, int{}, int{},
          next_t{}, dest, rest...
        );

        dest.build_buffer_helper(first);
      }

      // Case 5 B
      // This case is the only one which actually calls build buffer
      template<typename T, typename Dest, typename First>
      static void build(
        int, long, long, long, long,
        T, Dest& dest, Matrix<First>& first){

        dest.build_buffer_helper(first);

      }


      // Case 6 A
      template<typename T, typename Dest, typename First, typename... Rest>
      static void build(
        long, long, long, long, long,
        T, Dest& dest, First& first, Rest&... rest){
        // Only recursive call

        using next_t = typename next_t_helper<ctr, Tup, (nr_args - 1) == ctr>::type;

        build_buffer_util<nr_args, ctr + 1, Tup>::build(
          int{}, int{}, int{}, int{}, int{},
          next_t{}, dest, rest...
        );

      }

      // Case 6 B
      template<typename T, typename Dest, typename First>
      static void build(
        long, long, long, long, long,
        T, Dest& dest, First& first){
        // Only recursive call
      }
    };





  } // end of namespace _gpi
}

#endif // BUILD_TUP_HPP
