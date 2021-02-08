#ifndef ENVIRONMENT_HPP
#define ENVIRONMENT_HPP


/*
* This singleton scheme allows for better control of global state.
* In particular it calls gaspi_init and terminate correctly.
*/
class Environment
{
private:

  ~Environment(){
    gaspi_proc_term(GASPI_BLOCK);
  };

  Environment() : init_called{false}{};

  bool init_called;

public:
  static Environment& get_instance(){
    static Environment ins;
    return ins;
  }

  void init(){
    if(!init_called){
      gaspi_proc_init(GASPI_BLOCK);
      init_called = true;
    }
   };
};



#endif//ENVIRONMENT_HPP
