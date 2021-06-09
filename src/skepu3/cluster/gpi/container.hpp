#ifndef CONTAINER_HPP
#define CONTAINER_HPP
#include <vector>
#include <mutex>
#include <map>

#include <GASPI.h>

#include "environment.hpp"

namespace skepu{

  namespace _gpi{

   int curr_containers = 0;
   int created_containers = 0;

    class Container{

      // This is only a temporary, non scalable solution
      friend class Reduce1D;
      template<int, typename Ret, typename... Func_args>
      friend class Map1D;
      friend class FilterClass;
      friend class build_tup_util;
    private:

    protected:
      gaspi_rank_t rank;
      gaspi_rank_t nr_nodes;
      gaspi_queue_id_t queue;

      // Contains all ranks which we wait for before an operation is started
      std::vector<int> wait_ranks;
      gaspi_segment_id_t segment_id;

      gaspi_pointer_t cont_seg_ptr;
      gaspi_pointer_t comm_seg_ptr;

      unsigned long& op_nr;

      // Must be initialized by derived classes
      unsigned long* vclock;

      // A constraint follows the following scheme (rank, op_nr) and should be
      // interpreted as the lowest op_nr the remote rank must reach before we
      // are allowed to modify the containers data.
      std::map<int, unsigned long> constraints;

      Container() : wait_ranks{}, op_nr{Environment::op_nr}{

        Environment& env = Environment::get_instance();
        env.init();

        vclock = Environment::vclock;

        gaspi_proc_rank(&rank);
        gaspi_proc_num(&nr_nodes);
        gaspi_queue_create(&queue, GASPI_BLOCK);
        segment_id = created_containers * nr_nodes + rank + 1;
        //segment_id = created_containers;

        curr_containers++;
        created_containers++;
      }

    public:
      virtual ~Container(){
        curr_containers--;
        gaspi_wait(queue, GASPI_BLOCK);
        gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);
        gaspi_segment_delete(segment_id);
      }
    };
  }
}

#endif //CONTAINER_HPP
