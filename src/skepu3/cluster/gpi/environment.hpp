#ifndef ENVIRONMENT_HPP
#define ENVIRONMENT_HPP

#include <GASPI.h>
#include <cassert>
#include <random>
/*
* This singleton scheme allows for better control of global state.
* In particular it calls gaspi_init and terminate correctly.
*/
class Environment
{
private:

  ~Environment(){
    delete generator;
    gaspi_proc_term(GASPI_BLOCK);
  };

  Environment() : init_called{false}{};

  bool init_called;

  inline static gaspi_pointer_t vclock_void;

  inline static std::mt19937* generator;

public:
  static Environment& get_instance(){
    static Environment ins;
    return ins;
  }


  inline static unsigned long* vclock;


  inline static long unsigned op_nr = 0;

  void init(){
    if(!init_called){
      gaspi_proc_init(GASPI_BLOCK);
      init_called = true;

      //generator = new std::mt19937( time(NULL));



      gaspi_rank_t nr_nodes;
      gaspi_proc_num(&nr_nodes);

      assert(gaspi_segment_create(
        0,
        gaspi_size_t{2 * nr_nodes * sizeof(unsigned long)},
        GASPI_GROUP_ALL,
        GASPI_BLOCK,
        GASPI_MEM_INITIALIZED
      ) == GASPI_SUCCESS);


      gaspi_segment_ptr(0, &(Environment::vclock_void));
      Environment::vclock = (unsigned long*) Environment::vclock_void;

    }
   };

   static std::mt19937& get_generator(){
     if(generator == nullptr){
       generator = new std::mt19937(std::random_device{}() );
     }
     return *generator;
   }


};



#endif//ENVIRONMENT_HPP
