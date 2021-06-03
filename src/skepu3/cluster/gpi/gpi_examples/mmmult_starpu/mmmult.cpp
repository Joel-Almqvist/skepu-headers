#define SKEPU_PRECOMPILED 1
#define SKEPU_OPENMP 1
#define SKEPU_STARPU_MPI 1
#line 1 "/home/joel/Documents/exjobb/skepu/skepu_fork/skepu/skepu-headers/src/skepu3/cluster/gpi/gpi_examples/mmmult_starpu.cpp"
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




struct skepu_userfunction_skepu_skel_0mmprod_mmmult_f_float
{
using T = float;
constexpr static size_t totalArity = 3;
constexpr static size_t outArity = 1;
constexpr static bool indexed = 1;
using IndexType = skepu::Index2D;
using ElwiseArgs = std::tuple<>;
using ContainerArgs = std::tuple<const skepu::Mat<float>, const skepu::Mat<float>>;
using UniformArgs = std::tuple<>;
typedef std::tuple<skepu::ProxyTag::Default, skepu::ProxyTag::Default> ProxyTags;
constexpr static skepu::AccessMode anyAccessMode[] = {
skepu::AccessMode::Read, skepu::AccessMode::Read, };

using Ret = float;

constexpr static bool prefersMatrix = 1;

#define SKEPU_USING_BACKEND_OMP 1
#undef VARIANT_CPU
#undef VARIANT_OPENMP
#undef VARIANT_CUDA
#define VARIANT_CPU(block)
#define VARIANT_OPENMP(block) block
#define VARIANT_CUDA(block)
static inline SKEPU_ATTRIBUTE_FORCE_INLINE float OMP(skepu::Index2D idx, const skepu::Mat<float> lhs, const skepu::Mat<float> rhs)
{
	T res = 0;
	for (size_t i = 0; i < lhs.cols; ++i)
		res += lhs.data[idx.row * lhs.cols + i] * rhs.data[i * rhs.cols + idx.col];
	return res;
}
#undef SKEPU_USING_BACKEND_OMP

#define SKEPU_USING_BACKEND_CPU 1
#undef VARIANT_CPU
#undef VARIANT_OPENMP
#undef VARIANT_CUDA
#define VARIANT_CPU(block) block
#define VARIANT_OPENMP(block)
#define VARIANT_CUDA(block) block
static inline SKEPU_ATTRIBUTE_FORCE_INLINE float CPU(skepu::Index2D idx, const skepu::Mat<float> lhs, const skepu::Mat<float> rhs)
{
	T res = 0;
	for (size_t i = 0; i < lhs.cols; ++i)
		res += lhs.data[idx.row * lhs.cols + i] * rhs.data[i * rhs.cols + idx.col];
	return res;
}
#undef SKEPU_USING_BACKEND_CPU
};

#line 17 "/home/joel/Documents/exjobb/skepu/skepu_fork/skepu/skepu-headers/src/skepu3/cluster/gpi/gpi_examples/mmmult_starpu.cpp"
int main(int argc, char *argv[])
{
	if (argc < 5)
	{
		if(!skepu::cluster::mpi_rank())
			std::cout << "Usage: " << argv[0] << " height width inner backend path\n";
		exit(1);
	}

	size_t height = atoi(argv[1]);
	size_t width = atoi(argv[2]);
	size_t inner = atoi(argv[3]);
	auto spec = skepu::BackendSpec{skepu::Backend::typeFromString(argv[4])};
	std::string path = argv[5];
	skepu::setGlobalBackendSpec(spec);


	auto start = std::chrono::system_clock::now();

	skepu::Matrix<float> lhs(height, inner), rhs(inner, width), res(height, width);
	lhs.randomize(3, 9);
	rhs.randomize(0, 9);

	lhs.flush();
	rhs.flush();


	skepu::backend::Map<0, skepu_userfunction_skepu_skel_0mmprod_mmmult_f_float, bool, void> mmprod(false);
	mmprod(res, lhs, rhs);

	res.flush();

	auto end = std::chrono::system_clock::now();
	double rtime = std::chrono::duration<double>{end - start}.count();

	std::ofstream ofile(path+"/mmmult_starpu.txt");

	ofile << "My time = " << rtime <<
	", height = " << height <<
	", width = " << width <<
	", inner = " << inner;


	return 0;
}
