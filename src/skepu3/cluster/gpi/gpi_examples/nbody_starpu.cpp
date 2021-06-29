#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>

#include <skepu>


// Particle data structure that is used as an element type.
struct Particle
{
	float x, y, z;
	float vx, vy, vz;
	float m;
};


constexpr float G [[skepu::userconstant]] = 1;
constexpr float delta_t [[skepu::userconstant]] = 0.1;


/*
 * Array user-function that is used for applying nbody computation,
 * All elements from parr and a single element (named 'pi') are accessible
 * to produce one output element of the same type.
 */
Particle move(skepu::Index1D index, Particle pi, const skepu::Vec<Particle> parr)
{
	size_t i = index.i;

	float ax = 0.0, ay = 0.0, az = 0.0;
	size_t np = parr.size;

	for (size_t j = 0; j < np; ++j)
	{
		if (i != j)
		{
			Particle pj = parr[j];

			if(pj.x == 0){
				//printf("x is 0 at %lu\n", i);
			}


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



// A helper function to write particle output values to standard output stream.
void save_step(skepu::Vector<Particle> &particles, std::ostream &os = std::cout)
{
	int i = 0;

	os
		<< std::setw(4) << "#" << "  "
		<< std::setw(15) << "x"
		<< std::setw(15) << "y"
		<< std::setw(15) << "z"
		<< std::setw(15) << "vx"
		<< std::setw(15) << "vy"
		<< std::setw(15) << "vz" << "\n"
		<< std::string(96,'=') << "\n";
	for (Particle &p : particles)
	{
		os << std::setw( 4) << i++ << ": "
			<< std::setw(15) << p.x
			<< std::setw(15) << p.y
			<< std::setw(15) << p.z
			<< std::setw(15) << p.vx
			<< std::setw(15) << p.vy
			<< std::setw(15) << p.vz << "\n";
	}
}

//! A helper function to write particle output values to a file.
void save_step(skepu::Vector<Particle> &particles, const std::string &filename)
{
	std::ofstream out(filename);

	if (out.is_open())
		save_step(particles, out);
	else
		std::cerr << "Error: cannot open this file: " << filename << "\n";
}


auto nbody_init = skepu::Map<0>(init);
auto nbody_simulate_step = skepu::Map(move);

void nbody(skepu::Vector<Particle> &particles, size_t iterations, skepu::BackendSpec *spec = nullptr)
{
	size_t np = particles.size();
	skepu::Vector<Particle> doublebuffer(particles.size());

	if (spec) skepu::setGlobalBackendSpec(*spec);

	// particle vectors initialization
	nbody_init(particles, np);





	for (size_t i = 0; i < iterations; i += 2)
	{
		nbody_simulate_step(doublebuffer, particles, particles);
		nbody_simulate_step(particles, doublebuffer, doublebuffer);
	}
}

int main(int argc, char *argv[])
{
	if (argc < 4)
	{
		if(!skepu::cluster::mpi_rank())
			std::cout << "Usage: " << argv[0] << " particles-per-dim iterations backend\n";
		exit(1);
	}

	const size_t np = std::stoul(argv[1]);
	const size_t iterations = std::stoul(argv[2]);
	auto spec = skepu::BackendSpec{skepu::Backend::typeFromString(argv[3])};

	auto start = std::chrono::system_clock::now();


	// Particle vectors....
	skepu::Vector<Particle> particles(np);

	nbody(particles, iterations, &spec);
	particles.flush();

	//starpu_task_wait_for_all();
	auto end = std::chrono::system_clock::now();
	double rtime = std::chrono::duration<double>{end - start}.count();



	std::cout << "My time = " << rtime <<
	", size = " << np <<
	", iterations = " << iterations << std::endl;;



	// NOTE: Code below are working prints to verify the result
	/*
	std::stringstream outfile2;
	outfile2 << "output" << spec.type() << ".txt";

	if(!skepu::cluster::mpi_rank()){
		save_step(particles, outfile2.str());
	}
	else{
		particles.print([](Particle p){

			std::string str{};

			str += "x = ";
			str += std::to_string(p.x).substr(0,8);

			str += ", y = ";
			str += std::to_string(p.y).substr(0,6);

			str += ", z = ";
			str += std::to_string(p.z).substr(0,6);

			str += ", vx = ";
			str += std::to_string(p.vx).substr(0,8);

			str += ", vy = ";
			str += std::to_string(p.vy).substr(0,6);

			str += ", vz = ";
			str += std::to_string(p.vz).substr(0,6);


			return str;
		});
	}
	*/

	return 0;
}
