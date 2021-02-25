#ifndef CONTAINER_HPP
#define CONTAINER_HPP
#include <GASPI.h>
#include <vector>

#include "environment.hpp"

namespace skepu{

  namespace _gpi{

   int curr_containers = 0;
   int created_containers = 0;

    class Container{

      // This is only a temporary, non scalable solution
      friend class Reduce1D;
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
      long comm_offset;

      unsigned int comm_buffer_free_slot;

      unsigned long op_nr;

      // Use a short due to gaspi_notification_id_t being of that type
      // and all calculations using this will be converted to that type.
      short unsigned notif_ctr;

      // Must be initialized by derived classes
      unsigned long* vclock;

      Container() : wait_ranks{}{

        Environment& env = Environment::get_instance();
        env.init();

        gaspi_proc_rank(&rank);
        gaspi_proc_num(&nr_nodes);
        segment_id = created_containers * nr_nodes + rank;
        //segment_id = created_containers;

        op_nr = 0;
        notif_ctr = 0;
        comm_buffer_free_slot = 0;
        curr_containers++;
        created_containers++;
      }

    public:
      virtual ~Container(){
        curr_containers--;
      }
    };
  }
}

#endif //CONTAINER_HPP
