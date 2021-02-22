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



auto f = ([](int a, int b) {
	return new char{};
});



int main(){

	// auto x = []( int, long, bool ){};
  // std::cout << typeid( function_argument_type< decltype(x), 0 >::type ).name() << '\n';
  // std::cout << typeid( function_argument_type< decltype(x), 1 >::type ).name() << '\n';
  // std::cout << typeid( function_argument_type< decltype(x), 2 >::type ).name() << '\n';
	//using t1 = function_argument_type< decltype(x), 0 >::type;
	//bool const b1 = std::is_same<t1, int>::value;


  skepu::Matrix<int> m1{4,4,1};
  skepu::Matrix<int> m2{4,4,2};
  skepu::Matrix<int> m3{5,5,3};

	auto mult = skepu::Map<4>([](
		skepu::Index1D i, skepu::Vec<int> a, int b, int c) int {

		return 5;
	});

	mult.apply( m1, m2, m3, m3, 5);


  return 0;
}
