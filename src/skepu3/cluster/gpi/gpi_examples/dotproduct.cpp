#include <iostream>
#include <utility>
#include <cfloat>

#include <skepu>

template<typename T>
T mult(T a, T b)
{
	return a * b;
}

template<typename T>
T add(T a, T b)
{
	return a + b;
}


int main(int argc, char *argv[])
{
	if (argc < 3)
	{
		std::cout << "Usage: " << argv[0] << " size backend\n";
		exit(1);
	}

	const size_t size = atoi(argv[1]);
	auto spec = skepu::BackendSpec{argv[2]};


	auto dotprod_map = skepu::Map(mult<float>);
	auto dotprod_red = skepu::Reduce(add<float>);


	dotprod_map.setBackend(spec);

	skepu::Vector<float> a(size), b(size);


	a.randomize(0, 3);
	b.randomize(0, 2);


	dotprod_map(a, a, b);
	float res = 0;

	res = dotprod_red(a);

	/*
	skepu::external(
		skepu::read(a,b),
		[&]{
			std::cout << a << "\n";
			std::cout << b << "\n";
		});
		*/

	//float res = dotprod(a, b);

	std::cout << "res: " << res << std::endl;

	/*
	skepu::external(
		[&]{ std::cout << res << "\n"; });
		*/
	return 0;
}
