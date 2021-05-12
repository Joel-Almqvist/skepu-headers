#ifndef GPI_DUMMIES_HPP
#define GPI_DUMMIES_HPP

/* This file contains dummy functions and or structs needed to compile Skepu
*  code but which is not yet implemented in the GPI backend.
*/


namespace skepu{

  // WARNING This might cause many issues
  #define Vector Matrix

  namespace cluster{
    bool mpi_rank(){
      return true;
    }
  }

  namespace Backend{

    template<typename... T>
    int typeFromString(T...){
      return 0;
    }
  }


  struct BackendSpec{
    template <typename... T>
    BackendSpec(T...){};


    std::string type(){
      return "dummy";
    }

  };

  template<typename... T>
  void setGlobalBackendSpec(T...){}


}


#endif //GPI_DUMMIES_HPP
