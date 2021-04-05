#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <ctime>
#include <float.h>

#include <skepu>

// Particle data structure that is used as an element type.
struct Particle
{
	float x, y, z;
	float vx, vy, vz;
	float m;
};

struct Acc
{
	float ax, ay, az;
};



constexpr float G = 1;
constexpr float delta_t = 0.1;


// Generate user-function that is used for initializing particles array.
Particle init(skepu::Index1D index, size_t np)
{
	int s = index.i;
	int d = np / 2 + 1;
	int i = s % np;
	int j = ((s - i) / np) % np;
	int k = (((s - i) / np) - j) / np;


	Particle p;

	p.x = i - d + 1;
	p.y = j - d + 1;
	p.z = k - d + 1;

	p.vx = 0.0;
	p.vy = 0.0;
	p.vz = 0.0;

	p.m = 1;

	return p;
}



int main(int argc, char *argv[])
{

	int c;
	int r;

  size_t iterations;

  if(argc == 4){
    c = std::stoul(argv[1]);
		r = std::stoul(argv[2]);
    iterations = std::stoul(argv[3]);
  }
  else{
    c = 7;
		r = 7;
    iterations = 4;
  }

	size_t np = c * r;

  auto start = std::chrono::system_clock::now();

	// Particle vectors....
	skepu::Matrix<Particle> particles(c, r);
	skepu::Matrix<Acc> particles_delta(c, r);


	auto init_map = skepu::Map<2>(init);
	init_map(particles, np);

	// Capture these in the lambdas
	Particle p;
	size_t p_index;

	auto acc = skepu::Map<2>([&p, &p_index](skepu::Index1D index, Particle pi)
	{
		Acc a;

		// A particle does not accelerate itself
		if(index.i == p_index){
			a.ax = 0;
			a.ay = 0;
			a.az = 0;

			return a;
		}

		float rij = sqrt((pi.x - p.x) * (pi.x - p.x)
									 + (pi.y - p.y) * (pi.y - p.y)
									 + (pi.z - p.z) * (pi.z - p.z));


		float dum = G * pi.m * p.m / pow(rij, 3);

		a.ax = dum * (pi.x - p.x);
		a.ay = dum * (pi.y - p.y);
		a.az = dum * (pi.z - p.z);

		return a;
	});


	auto red = skepu::Reduce([](Acc a, Acc b){
		Acc res;
		res.ax = a.ax + b.ax;
		res.ay = a.ay + b.ay;
		res.az = a.az + b.az;
		return res;
	});


	Acc curr_acc;

	for(int j = 0; j < iterations; j++){
		for(int i = 0; i < particles.size(); i++){
			p = particles[i];
			p_index = i;

			acc(particles_delta, particles);

			curr_acc = red(particles_delta);


			p.x = p.x + delta_t * p.vx + delta_t * delta_t * curr_acc.ax / 2;
			p.y = p.y + delta_t * p.vy + delta_t * delta_t * curr_acc.ay / 2;
			p.z = p.z + delta_t * p.vz + delta_t * delta_t * curr_acc.az / 2;

			p.vx = p.vx + delta_t * curr_acc.ax;
			p.vy = p.vy + delta_t * curr_acc.ay;
			p.vz = p.vz + delta_t * curr_acc.az;


			particles.set(i, p);
		}

		/*
		if(j == iterations - 1){
			particles.print([](Particle pp){
				return "x = "+std::to_string(pp.x)+
				" y = " + std::to_string(pp.y)+
				" z = "+std::to_string(pp.z)+
				" vx = "+std::to_string(pp.vx);
			});
		}
		*/
	}


  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;

  std::cout << "Exec time = " << elapsed_seconds.count() << " c = " << c
  << ", r = " << r <<", iterations = " << iterations << std::endl;

	return 0;
}
