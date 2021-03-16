#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <ctime>


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


//constexpr float G [[skepu::userconstant]] = 1;
//constexpr float delta_t [[skepu::userconstant]] = 0.1;

constexpr float G = 1;
constexpr float delta_t = 0.1;



/*
 * Array user-function that is used for applying nbody computation,
 * All elements from parr and a single element (named 'pi') are accessible
 * to produce one output element of the same type.
 */
Particle move(skepu::Index1D index, Particle pi, skepu::Vec<Particle> parr)
{
	size_t i = index.i;

	float ax = 0.0, ay = 0.0, az = 0.0;
	size_t np = parr.size;

	for (size_t j = 0; j < np; ++j)
	{
		if (i != j)
		{
			Particle pj = parr[j];

			float rij = sqrt((pi.x - pj.x) * (pi.x - pj.x)
			               + (pi.y - pj.y) * (pi.y - pj.y)
			               + (pi.z - pj.z) * (pi.z - pj.z));

			float dum = G * pi.m * pj.m / pow(rij, 3);

			ax += dum * (pi.x - pj.x);
			ay += dum * (pi.y - pj.y);
			az += dum * (pi.z - pj.z);
		}
	}

//	std::cout << "i = " << i << ": ax = " << ax << ", ay = " << ay << ", az = " << az << "\n";

	Particle newp;
	newp.m = pi.m;

	newp.x = pi.x + delta_t * pi.vx + delta_t * delta_t / 2 * ax;
	newp.y = pi.y + delta_t * pi.vy + delta_t * delta_t / 2 * ay;
	newp.z = pi.z + delta_t * pi.vz + delta_t * delta_t / 2 * az;

	newp.vx = pi.vx + delta_t * ax;
	newp.vy = pi.vy + delta_t * ay;
	newp.vz = pi.vz + delta_t * az;

	return newp;
}


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
    c = 8;
		r = 8;
    iterations = 10;
  }

	size_t np = c * r;

	//auto spec = skepu::BackendSpec{skepu::Backend::typeFromString(argv[3])};

  auto start = std::chrono::system_clock::now();

	// Particle vectors....
	skepu::Matrix<Particle> particles(c, r);
	skepu::Matrix<Acc> particles_delta(c, r);



	auto init_map = skepu::Map<2>(init);
	init_map(particles, np);


	Particle p;
	auto acc = skepu::Map<1>([&p](Particle pi){

		Acc a;

		float rij = sqrt((pi.x - p.x) * (pi.x - p.x)
									 + (pi.y - p.y) * (pi.y - p.y)
									 + (pi.z - p.z) * (pi.z - p.z));

		float dum = G * pi.m * p.m / pow(rij, 3);

		a.ax = dum * (pi.x - p.x);
		a.ay += dum * (pi.y - p.y);
		a.az += dum * (pi.z - p.z);

		return a;
	});


	auto red = skepu::Reduce([](Acc a, Acc b){
		Acc res;
		res.ax = a.ax + b.ax;
		return res;
	});


	Acc curr_acc;

	for(int j = 0; j < iterations; j++){
		for(int i = 0; i < particles.size(); i++){
			p = particles[i];
			acc(particles_delta, particles);

			curr_acc = red(particles_delta);

			// Here we do more work than actually needed
			p.x = p.x + delta_t * p.vx + delta_t * delta_t * curr_acc.ax / 2;
			p.y = p.y + delta_t * p.vy + delta_t * delta_t * curr_acc.ay / 2;
			p.z = p.z + delta_t * p.vz + delta_t * delta_t * curr_acc.az / 2;

			particles.set(i, p);
		}
	}


  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;

  std::cout << "Exec time = " << elapsed_seconds.count() << " c = " << c
  << ", r = " << r <<", iterations = " << iterations << std::endl;

	return 0;
}
