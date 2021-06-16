#include <iostream>
#include <skepu>


template<typename T>
T mmmult_f(const skepu::MatRow<T> ar, const skepu::MatCol<T> bc)
{
	T res = 0;
	for (size_t k = 0; k < ar.cols; ++k)
		res += ar(k) * bc(k);
	return res;
}

// A helper function to calculate dense matrix-matrix product. Used to verify that the SkePU output is correct.
template<typename T>
void directMM(skepu::Matrix<T> &lhs, skepu::Matrix<T> &rhs, skepu::Matrix<T> &res)
{
	int rows  = res.size_i();
	int cols  = res.size_j();
	int inner = lhs.total_cols();

	for (int i = 0; i < rows; ++i)
		for (int j = 0; j < cols; ++j)
		{
			T sum{};
			for (int k = 0; k < inner; ++k)
				sum += lhs(i, k) * rhs(k, j);

			res(i, j) = sum;
		}
}

int main(int argc, char *argv[])
{
	if (argc < 4)
	{
		skepu::external([&]{
			std::cout << "Usage: " << argv[0] << " height width inner backend\n";
		});
		exit(1);
	}

	size_t height = atoi(argv[1]);
	size_t width = atoi(argv[2]);
	size_t inner = atoi(argv[3]);
	auto spec = skepu::BackendSpec{skepu::Backend::typeFromString(argv[4])};
	skepu::setGlobalBackendSpec(spec);

	auto start = std::chrono::system_clock::now();

	skepu::Matrix<float> lhs(height, inner), rhs(inner, width), res(height, width), res2(height, width);
	lhs.randomize(3, 9);
	rhs.randomize(0, 9);

	auto mmprod = skepu::Map(mmmult_f<float>);
	mmprod(res, lhs, rhs);

	auto end = std::chrono::system_clock::now();
	double rtime = std::chrono::duration<double>{end - start}.count();

	printf("My time = %f, width = %lu, height = %lu\n",
		rtime, width, height);

	return 0;
}
