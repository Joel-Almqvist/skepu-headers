#define SKEPU_PRECOMPILED 1
#define SKEPU_OPENMP 1
#define SKEPU_STARPU_MPI 1
#line 1 "/home/joel/Documents/exjobb/skepu/skepu_fork/skepu/skepu-headers/src/skepu3/cluster/gpi/gpi_examples/mmmult_row_col.cpp"
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


struct skepu_userfunction_skepu_skel_0mmprod_mmmult_f_float
{
using T = float;
constexpr static size_t totalArity = 2;
constexpr static size_t outArity = 1;
constexpr static bool indexed = 0;
using IndexType = void;
using ElwiseArgs = std::tuple<>;
using ContainerArgs = std::tuple<const skepu::MatRow<float>, const skepu::MatCol<float>>;
using UniformArgs = std::tuple<>;
typedef std::tuple<skepu::ProxyTag::MatRow, skepu::ProxyTag::MatCol> ProxyTags;
constexpr static skepu::AccessMode anyAccessMode[] = {
skepu::AccessMode::Read, skepu::AccessMode::Read, };

using Ret = float;

constexpr static bool prefersMatrix = 0;

#define SKEPU_USING_BACKEND_OMP 1
#undef VARIANT_CPU
#undef VARIANT_OPENMP
#undef VARIANT_CUDA
#define VARIANT_CPU(block)
#define VARIANT_OPENMP(block) block
#define VARIANT_CUDA(block)
static inline SKEPU_ATTRIBUTE_FORCE_INLINE float OMP(const skepu::MatRow<float> ar, const skepu::MatCol<float> bc)
{
	T res = 0;
	for (size_t k = 0; k < ar.cols; ++k)
		res += ar(k) * bc(k);
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
static inline SKEPU_ATTRIBUTE_FORCE_INLINE float CPU(const skepu::MatRow<float> ar, const skepu::MatCol<float> bc)
{
	T res = 0;
	for (size_t k = 0; k < ar.cols; ++k)
		res += ar(k) * bc(k);
	return res;
}
#undef SKEPU_USING_BACKEND_CPU
};

#line 33 "/home/joel/Documents/exjobb/skepu/skepu_fork/skepu/skepu-headers/src/skepu3/cluster/gpi/gpi_examples/mmmult_row_col.cpp"
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

	skepu::backend::Map<0, skepu_userfunction_skepu_skel_0mmprod_mmmult_f_float, bool, void> mmprod(false);
	mmprod(res, lhs, rhs);

	auto end = std::chrono::system_clock::now();
	double rtime = std::chrono::duration<double>{end - start}.count();

	printf("My time = %f, width = %lu, height = %lu\n",
		rtime, width, height);

	return 0;
}
