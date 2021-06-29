#define SKEPU_PRECOMPILED 1
#define SKEPU_OPENMP 1
#define SKEPU_STARPU_MPI 1
#line 1 "/home/joel/Documents/exjobb/skepu/skepu_fork/skepu/skepu-headers/src/skepu3/cluster/gpi/gpi_examples/taylor_starpu.cpp"
/*
 * Taylor series calculation, natural log(1+x)  sum(1:N) (((-1)^(i+1))/i)*x^i
 */

#include <iostream>
#include <cmath>

#include <skepu>


float nth_term(skepu::Index1D index, float x)
{
	int k = index.i + 1;
	float temp_x = pow(x, k);
	int sign = (k % 2 == 0) ? -1 : 1;
	return sign * temp_x / k;
}

float plus(float a, float b)
{
	return a + b;
}


struct skepu_userfunction_skepu_skel_0taylor_nth_term
{
constexpr static size_t totalArity = 2;
constexpr static size_t outArity = 1;
constexpr static bool indexed = 1;
using IndexType = skepu::Index1D;
using ElwiseArgs = std::tuple<>;
using ContainerArgs = std::tuple<>;
using UniformArgs = std::tuple<float>;
typedef std::tuple<> ProxyTags;
constexpr static skepu::AccessMode anyAccessMode[] = {
};

using Ret = float;

constexpr static bool prefersMatrix = 0;

#define SKEPU_USING_BACKEND_OMP 1
#undef VARIANT_CPU
#undef VARIANT_OPENMP
#undef VARIANT_CUDA
#define VARIANT_CPU(block)
#define VARIANT_OPENMP(block) block
#define VARIANT_CUDA(block)
static inline SKEPU_ATTRIBUTE_FORCE_INLINE float OMP(skepu::Index1D index, float x)
{
	int k = index.i + 1;
	float temp_x = pow(x, k);
	int sign = (k % 2 == 0) ? -1 : 1;
	return sign * temp_x / k;
}
#undef SKEPU_USING_BACKEND_OMP

#define SKEPU_USING_BACKEND_CPU 1
#undef VARIANT_CPU
#undef VARIANT_OPENMP
#undef VARIANT_CUDA
#define VARIANT_CPU(block) block
#define VARIANT_OPENMP(block)
#define VARIANT_CUDA(block) block
static inline SKEPU_ATTRIBUTE_FORCE_INLINE float CPU(skepu::Index1D index, float x)
{
	int k = index.i + 1;
	float temp_x = pow(x, k);
	int sign = (k % 2 == 0) ? -1 : 1;
	return sign * temp_x / k;
}
#undef SKEPU_USING_BACKEND_CPU
};

#line 24 "/home/joel/Documents/exjobb/skepu/skepu_fork/skepu/skepu-headers/src/skepu3/cluster/gpi/gpi_examples/taylor_starpu.cpp"

struct skepu_userfunction_skepu_skel_0taylor_plus
{
constexpr static size_t totalArity = 2;
constexpr static size_t outArity = 1;
constexpr static bool indexed = 0;
using IndexType = void;
using ElwiseArgs = std::tuple<float, float>;
using ContainerArgs = std::tuple<>;
using UniformArgs = std::tuple<>;
typedef std::tuple<> ProxyTags;
constexpr static skepu::AccessMode anyAccessMode[] = {
};

using Ret = float;

constexpr static bool prefersMatrix = 0;

#define SKEPU_USING_BACKEND_OMP 1
#undef VARIANT_CPU
#undef VARIANT_OPENMP
#undef VARIANT_CUDA
#define VARIANT_CPU(block)
#define VARIANT_OPENMP(block) block
#define VARIANT_CUDA(block)
static inline SKEPU_ATTRIBUTE_FORCE_INLINE float OMP(float a, float b)
{
	return a + b;
}
#undef SKEPU_USING_BACKEND_OMP

#define SKEPU_USING_BACKEND_CPU 1
#undef VARIANT_CPU
#undef VARIANT_OPENMP
#undef VARIANT_CUDA
#define VARIANT_CPU(block) block
#define VARIANT_OPENMP(block)
#define VARIANT_CUDA(block) block
static inline SKEPU_ATTRIBUTE_FORCE_INLINE float CPU(float a, float b)
{
	return a + b;
}
#undef SKEPU_USING_BACKEND_CPU
};

#line 24 "/home/joel/Documents/exjobb/skepu/skepu_fork/skepu/skepu-headers/src/skepu3/cluster/gpi/gpi_examples/taylor_starpu.cpp"
int main(int argc, char *argv[])
{
	if (argc < 4)
	{
		std::cout << "Usage: " << argv[0] << " x-value number-of-terms backend\n";
		exit(1);
	}

	auto spec = skepu::BackendSpec{skepu::Backend::typeFromString(argv[3])};
	skepu::setGlobalBackendSpec(spec);
	float x = atof(argv[1]);
	size_t N = std::stoul(argv[2]);

	auto start = std::chrono::system_clock::now();

	skepu::backend::MapReduce<0, skepu_userfunction_skepu_skel_0taylor_nth_term, skepu_userfunction_skepu_skel_0taylor_plus, bool, bool, void> taylor(false, false);
	taylor.setDefaultSize(N);

	starpu_task_wait_for_all();
	auto end = std::chrono::system_clock::now();

	double rtime = std::chrono::duration<double>{end - start}.count();

	float res = taylor(x - 1);

	std::cout << "Res: ln(" << x << ") = " << res <<
	", my time = " << rtime <<
	", N = " << N <<  std::endl;


	//std::cout << "Result: ln(" << x << ") = " << taylor(x - 1) << "\n";

	return 0;

}
