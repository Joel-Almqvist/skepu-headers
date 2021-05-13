/*!
 *  PPMCC stands for "Pearson product-moment correlation coefficient"
 *  In statistics, the Pearson coefficient of correlation is a measure by the
 *  linear dependence between two variables X and Y. The mathematical
 *  expression of the Pearson coefficient of correlation is as follows:
 *   r = ( (n*sum(X.Y)-sum(X)*sum(Y))/((n*sum(X^2)-(sum(X))^2)*(n*sum(Y^2)-(sum(Y))^2)) )
 */

#include <iostream>
#include <cmath>

#include <skepu>

// Unary user-function used for mapping
template<typename T>
T square(T a)
{
	return a * a;
}

// Binary user-function used for mapping
template<typename T>
T mult(T a, T b)
{
	return a * b;
}

// User-function used for reduction
template<typename T>
T plus(T a, T b)
{
	return a + b;
}

using T = float;

T ppmcc(skepu::Vector<T> &x, skepu::Vector<T> &y)
{
	// Skeleton definitions
	auto sum = skepu::Reduce(plus<T>);

	//auto dotProd_map = skepu::MapReduce(mult<T>, plus<T>);
	//auto sumSquare = skepu::MapReduce(square<T>, plus<T>);

	auto mult_map = skepu::Map(mult<T>);
	auto square_map = skepu::Map(square<T>);

	size_t N = x.size();

	T sumX = sum(x);
	T sumY = sum(y);

	skepu::Vector<T> swap_space(N);

	mult_map(swap_space, x, y);
	T dot_prod = sum(swap_space);

	square_map(x, x);
	T sumSquare_x = sum(x);

	square_map(y, y);
	T sumSquare_y = sum(y);

	return (N * dot_prod - sumX * sumY)
		/ sqrt((N * sumSquare_x - pow(sumX, 2)) * (N * sumSquare_y - pow(sumY, 2)));

}

int main(int argc, char *argv[])
{
	if (argc < 3)
	{
		/*
		skepu::external(
			[&]{
				std::cout << "Usage: " << argv[0] << " input_size backend\n";
			});
		*/
		exit(1);
	}

	const size_t size = std::stoul(argv[1]);
	auto spec = skepu::BackendSpec{skepu::Backend::typeFromString(argv[2])};
	skepu::setGlobalBackendSpec(spec);

	// Vector operands
	skepu::Vector<T> x(size), y(size);
	x.randomize(1, 3);
	y.randomize(2, 4);

	/*
	skepu::external(
		skepu::read(x,y),
		[&]{
			std::cout << "X: " << x << "\n";
			std::cout << "Y: " << y << "\n";
		});
	*/

	T res = ppmcc(x, y);

	std::cout << "res: " << res << std::endl;

	/*
	skepu::external([&]{
		std::cout << "res: " << res << "\n";});
	*/
	return 0;
}
