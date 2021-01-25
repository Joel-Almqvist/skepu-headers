#ifndef SKELETON_BASE_HPP
#define SKELETON_BASE_HPP

namespace skepu{

  namespace _gpi{

    // TODO Move this functionality to Matrix or remove

    // Globally unique ID for each operation
    long unsigned op_nr = 0;

    class skeleton_base{
    protected:
      void increment_op_nr(){
        op_nr++;
      }

      long unsigned get_op_nr(){
        return op_nr;
      }


    };
  }
}


#endif //SKELETON_BASE_HPP
