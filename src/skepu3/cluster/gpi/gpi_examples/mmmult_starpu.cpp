#include <iostream>
#include <skepu>
#include <fstream>

template<typename T>
T mmmult_f(skepu::Index2D idx, const skepu::Mat<T> lhs, const skepu::Mat<T> rhs)
{
	T res = 0;
	for (size_t i = 0; i < lhs.cols; ++i)
		res += lhs.data[idx.row * lhs.cols + i] * rhs.data[i * rhs.cols + idx.col];
	return res;
}



int main(int argc, char *argv[])
{
	if (argc < 4)
	{
		if(!skepu::cluster::mpi_rank())
			std::cout << "Usage: " << argv[0] << " height width inner backend\n";
		exit(1);
	}

	size_t height = atoi(argv[1]);
	size_t width = atoi(argv[2]);
	size_t inner = atoi(argv[3]);
	auto spec = skepu::BackendSpec{skepu::Backend::typeFromString(argv[4])};
	skepu::setGlobalBackendSpec(spec);



	skepu::Matrix<float> lhs(height, inner), rhs(inner, width), res(height, width);
	lhs.randomize(3, 9);
	rhs.randomize(0, 9);

	lhs.flush();
	rhs.flush();


	auto mmprod = skepu::Map(mmmult_f<float>);
	auto start = std::chrono::system_clock::now();
	mmprod(res, lhs, rhs);

	res.flush();
	//starpu_task_wait_for_all();
	auto end = std::chrono::system_clock::now();
	double rtime = std::chrono::duration<double>{end - start}.count();


	std::cout << "My time = " << rtime <<
	", height = " << height <<
	", width = " << width <<
	", inner = " << inner << std::endl;


	return 0;
}
