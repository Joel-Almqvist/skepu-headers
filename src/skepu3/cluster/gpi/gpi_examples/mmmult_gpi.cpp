#include <iostream>
#include <skepu>
#include <fstream>

template<typename T>
T mmmult_f(skepu::Index2D idx, const skepu::Mat<T> lhs, const skepu::Mat<T> rhs)
{
	T res = 0;
	for (size_t i = 0; i < lhs.cols; ++i)
		res += lhs[idx.row * lhs.cols + i] * rhs[i * rhs.cols + idx.col];
	return res;
}



int main(int argc, char *argv[])
{
	if (argc < 5)
	{
		if(!skepu::cluster::mpi_rank())
			std::cout << "Usage: " << argv[0] << " height width inner path\n";
		exit(1);
	}

	size_t height = atoi(argv[1]);
	size_t width = atoi(argv[2]);
	size_t inner = atoi(argv[3]);
	std::string path = argv[4];
	auto spec = skepu::BackendSpec{skepu::Backend::typeFromString(argv[4])};
	skepu::setGlobalBackendSpec(spec);


	auto start = std::chrono::system_clock::now();

	skepu::Matrix<float> lhs(height, inner), rhs(inner, width), res(height, width);
	lhs.randomize(3, 9);
	rhs.randomize(0, 9);

	lhs.flush();
	rhs.flush();


	auto mmprod = skepu::Map(mmmult_f<float>);
	mmprod(res, lhs, rhs);

	res.flush();
	auto end = std::chrono::system_clock::now();
	double rtime = std::chrono::duration<double>{end - start}.count();

	double slowest_rtime = rhs.get_slowest_node(rtime);
	double avg_rtime = rhs.get_avg_time(rtime);

	std::ofstream ofile(path+"/mmmult_gpi.txt");

	ofile << "Slowest = " << slowest_rtime <<
	", Avg time = " << avg_rtime <<
	", my time = " << rtime <<
	", height = " << height <<
	", width = " << width <<
	", inner = " << inner;

	ofile.close();

	return 0;
}
