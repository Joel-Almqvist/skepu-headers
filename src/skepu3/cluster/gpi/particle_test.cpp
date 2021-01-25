#include <iostream>
#include <GASPI.h>
#include <utility>

#include <matrix.hpp>
#include <map.hpp>
#include <omp.h>

struct Particle
{
	float x, y, z;
	float vx, vy, vz;
	float m;
};



// Helps with the prints for the demonstration
auto part_to_string = [](Particle p) -> std::string{
	std::string str{};
	str += std::to_string((int) p.x);

	str += " ";
	str += std::to_string((int) p.y);

	str += " ";
	str += std::to_string((int) p.z);

	return str;
};


static const int DELTA = 1;


int main(){

  skepu::Matrix<Particle> m{4,4};

	// Initialize the particles.
	// This function is a demo function and is not part of final the SkePU interface
  m.set(Particle{1, 2, 3, 0, 0, 0, 10});



	auto map_part_move = skepu::Map<1>([](Particle p){
		float new_x = p.x + p.vx * DELTA;
		float new_y = p.y + p.vy * DELTA;
		float new_z = p.z + p.vz * DELTA;

		return Particle{new_x, new_y, new_z, p.vx, p.vy, p.vz, p.m};
	});



  auto map_part_acc = skepu::Map<1>([](Particle p){
    return Particle{p.x, p.y, p.z, p.vx + 1, p.vy + 1, p.vz + 1, p.m};
  });


	map_part_move(m, m);
	m.print(part_to_string);


	map_part_acc(m, m);
	map_part_move(m, m);
	map_part_move(m, m);
	map_part_move(m, m);

	// Acceleration is now 1 and we have moved 3 times
	m.print(part_to_string);


	map_part_acc(m, m);
	map_part_acc(m, m);
	map_part_acc(m, m);
	map_part_acc(m, m);
	map_part_move(m, m);

	// Acceleration is now 5 and we have moved 1 more time
	m.print(part_to_string);

  return 0;
}
