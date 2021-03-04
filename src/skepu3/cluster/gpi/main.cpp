#include <iostream>
#include <skepu>
#include <vector>
#include <utility>


struct Particle
{
	float x, y, z;
	float vx, vy, vz;
	float m;

};

auto lambda = [](
	skepu::Index1D i, Particle p) -> Particle{
	return Particle{i.i,i.i,i.i,i.i,i.i,i.i,i.i};
};

int main(){


  skepu::Matrix<int> m1{7,7,1};
	 skepu::Matrix<int> m2{7,7,2};
   skepu::Matrix<Particle> m3{7,7};

	auto mult = skepu::Map<3>([](
		skepu::Index1D i, int b, skepu::Vec<Particle> p) int {
		return b +p[48 - i.i].x;
	});

	auto init = skepu::Map<2>(lambda);

	for(int i = 0; i < 49; i++){
		m2.set(i, i);
	}

	init(m3, m3);
	mult( m1, m2, m3);
	m1.print();

  return 0;
}
