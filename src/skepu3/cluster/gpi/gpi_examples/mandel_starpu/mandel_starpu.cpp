#define SKEPU_PRECOMPILED 1
#define SKEPU_OPENMP 1
#define SKEPU_STARPU_MPI 1
#line 1 "/home/joel/Documents/exjobb/skepu/skepu_fork/skepu/skepu-headers/src/skepu3/cluster/gpi/gpi_examples/mandelbrot_starpu.cpp"
/*!
Mandelbrot fractals. The Mandelbrot set
{
"B. B. Mandelbrot. Fractal aspects of the iteration of z → λz(1 − z) for complex λ and z.
Annals of the New York Academy of Sciences, 357:249–259, December 1980."
}
is a set of complex numbers which boundary draws a fractal in the complex numbers plane. A complex number c lies
within the Mandelbrot set, if the sequence
z_{i+1} = z_{i}2 + c
with i ∈ N, starting with z0 = 0 does not escape to infinity, otherwise c is not part of
the Mandelbrot set.

When computing a Mandelbrot fractal, the sequence in equation 3.1 is calculated
for every pixel of an image representing a section of the complex numbers plane. If a
given threshold is crossed, it is presumed that the sequence will escape to infinity and
that the pixel is not inside the Mandelbrot set. If the threshold is not crossed for a given
number of steps in the sequence, the pixel is taken as a member of the Mandelbrot
set. A pixel within the Mandelbrot set painted in black, other pixels are given a color
that corresponds to the number of sequence steps that have been calculated before
excluding the pixel from the Mandelbrot set. By setting the threshold and the number
of sequence steps accordingly, the calculation of the fractal can be a time-consuming
task. However, as all pixels are calculated independently, it is a common benchmark
application for data-parallel computations.
*/

#include <iostream>
#include <fstream>
#include <cstdlib>

#include <skepu>

#include "bitmap_image.hpp"

template<typename T>
void save_image(size_t width, size_t height, T *buf, T max)
{

	bitmap_image image(width, height);

	for (size_t y = 0; y < height; ++y)
	{
		for (size_t x = 0; x < width; ++x)
		{
			float val = buf[y*width + x];
			unsigned char shade = val / max * 255;
			rgb_store col = prism_colormap[shade];
			image.set_pixel(x, y, col.red, col.green, col.blue);
		}
	}

	image.save_image("generated_image.bmp");
}


[[skepu::userconstant]] constexpr float
	CENTER_X = -.5f,
	CENTER_Y = 0.f,
	SCALE = 2.5f;

[[skepu::userconstant]] constexpr size_t
	MAX_ITERS = 10000;

struct cplx
{
	float a, b;
};

cplx mult_c(cplx lhs, cplx rhs)
{
	cplx r;
	r.a = lhs.a * rhs.a - lhs.b * rhs.b;
	r.b = lhs.b * rhs.a + lhs.a * rhs.b;
	return r;
}

cplx add_c(cplx lhs, cplx rhs)
{
	cplx r;
	r.a = lhs.a + rhs.a;
	r.b = lhs.b + rhs.b;
	return r;
}

size_t mandelbrot_f(skepu::Index2D index, float height, float width)
{
	cplx a;
	a.a = SCALE / height * (index.col - width/2.f) + CENTER_X;
	a.b = SCALE / height * (index.row - width/2.f) + CENTER_Y;
	cplx c = a;

	for (size_t i = 0; i < MAX_ITERS; ++i)
	{
		a = add_c(mult_c(a, a), c);
		if ((a.a * a.a + a.b * a.b) > 9)
			return i;
	}
	return MAX_ITERS;
}


struct skepu_userfunction_skepu_skel_0mandelbroter_mult_c
{
constexpr static size_t totalArity = 2;
constexpr static size_t outArity = 1;
constexpr static bool indexed = 0;
using IndexType = void;
using ElwiseArgs = std::tuple<>;
using ContainerArgs = std::tuple<>;
using UniformArgs = std::tuple<cplx, cplx>;
typedef std::tuple<> ProxyTags;
constexpr static skepu::AccessMode anyAccessMode[] = {
};

using Ret = cplx;

constexpr static bool prefersMatrix = 0;

#define SKEPU_USING_BACKEND_OMP 1
#undef VARIANT_CPU
#undef VARIANT_OPENMP
#undef VARIANT_CUDA
#define VARIANT_CPU(block)
#define VARIANT_OPENMP(block) block
#define VARIANT_CUDA(block)
static inline SKEPU_ATTRIBUTE_FORCE_INLINE cplx OMP(cplx lhs, cplx rhs)
{
	cplx r;
	r.a = lhs.a * rhs.a - lhs.b * rhs.b;
	r.b = lhs.b * rhs.a + lhs.a * rhs.b;
	return r;
}
#undef SKEPU_USING_BACKEND_OMP

#define SKEPU_USING_BACKEND_CPU 1
#undef VARIANT_CPU
#undef VARIANT_OPENMP
#undef VARIANT_CUDA
#define VARIANT_CPU(block) block
#define VARIANT_OPENMP(block)
#define VARIANT_CUDA(block) block
static inline SKEPU_ATTRIBUTE_FORCE_INLINE cplx CPU(cplx lhs, cplx rhs)
{
	cplx r;
	r.a = lhs.a * rhs.a - lhs.b * rhs.b;
	r.b = lhs.b * rhs.a + lhs.a * rhs.b;
	return r;
}
#undef SKEPU_USING_BACKEND_CPU
};

#line 100 "/home/joel/Documents/exjobb/skepu/skepu_fork/skepu/skepu-headers/src/skepu3/cluster/gpi/gpi_examples/mandelbrot_starpu.cpp"

struct skepu_userfunction_skepu_skel_0mandelbroter_add_c
{
constexpr static size_t totalArity = 2;
constexpr static size_t outArity = 1;
constexpr static bool indexed = 0;
using IndexType = void;
using ElwiseArgs = std::tuple<>;
using ContainerArgs = std::tuple<>;
using UniformArgs = std::tuple<cplx, cplx>;
typedef std::tuple<> ProxyTags;
constexpr static skepu::AccessMode anyAccessMode[] = {
};

using Ret = cplx;

constexpr static bool prefersMatrix = 0;

#define SKEPU_USING_BACKEND_OMP 1
#undef VARIANT_CPU
#undef VARIANT_OPENMP
#undef VARIANT_CUDA
#define VARIANT_CPU(block)
#define VARIANT_OPENMP(block) block
#define VARIANT_CUDA(block)
static inline SKEPU_ATTRIBUTE_FORCE_INLINE cplx OMP(cplx lhs, cplx rhs)
{
	cplx r;
	r.a = lhs.a + rhs.a;
	r.b = lhs.b + rhs.b;
	return r;
}
#undef SKEPU_USING_BACKEND_OMP

#define SKEPU_USING_BACKEND_CPU 1
#undef VARIANT_CPU
#undef VARIANT_OPENMP
#undef VARIANT_CUDA
#define VARIANT_CPU(block) block
#define VARIANT_OPENMP(block)
#define VARIANT_CUDA(block) block
static inline SKEPU_ATTRIBUTE_FORCE_INLINE cplx CPU(cplx lhs, cplx rhs)
{
	cplx r;
	r.a = lhs.a + rhs.a;
	r.b = lhs.b + rhs.b;
	return r;
}
#undef SKEPU_USING_BACKEND_CPU
};

#line 100 "/home/joel/Documents/exjobb/skepu/skepu_fork/skepu/skepu-headers/src/skepu3/cluster/gpi/gpi_examples/mandelbrot_starpu.cpp"

struct skepu_userfunction_skepu_skel_0mandelbroter_mandelbrot_f
{
constexpr static size_t totalArity = 3;
constexpr static size_t outArity = 1;
constexpr static bool indexed = 1;
using IndexType = skepu::Index2D;
using ElwiseArgs = std::tuple<>;
using ContainerArgs = std::tuple<>;
using UniformArgs = std::tuple<float, float>;
typedef std::tuple<> ProxyTags;
constexpr static skepu::AccessMode anyAccessMode[] = {
};

using Ret = unsigned long;

constexpr static bool prefersMatrix = 1;

#define SKEPU_USING_BACKEND_OMP 1
#undef VARIANT_CPU
#undef VARIANT_OPENMP
#undef VARIANT_CUDA
#define VARIANT_CPU(block)
#define VARIANT_OPENMP(block) block
#define VARIANT_CUDA(block)
static inline SKEPU_ATTRIBUTE_FORCE_INLINE unsigned long OMP(skepu::Index2D index, float height, float width)
{
	cplx a;
	a.a = SCALE / height * (index.col - width/2.f) + CENTER_X;
	a.b = SCALE / height * (index.row - width/2.f) + CENTER_Y;
	cplx c = a;

	for (size_t i = 0; i < MAX_ITERS; ++i)
	{
		a = skepu_userfunction_skepu_skel_0mandelbroter_add_c::OMP(skepu_userfunction_skepu_skel_0mandelbroter_mult_c::OMP(a, a), c);
		if ((a.a * a.a + a.b * a.b) > 9)
			return i;
	}
	return MAX_ITERS;
}
#undef SKEPU_USING_BACKEND_OMP

#define SKEPU_USING_BACKEND_CPU 1
#undef VARIANT_CPU
#undef VARIANT_OPENMP
#undef VARIANT_CUDA
#define VARIANT_CPU(block) block
#define VARIANT_OPENMP(block)
#define VARIANT_CUDA(block) block
static inline SKEPU_ATTRIBUTE_FORCE_INLINE unsigned long CPU(skepu::Index2D index, float height, float width)
{
	cplx a;
	a.a = SCALE / height * (index.col - width/2.f) + CENTER_X;
	a.b = SCALE / height * (index.row - width/2.f) + CENTER_Y;
	cplx c = a;

	for (size_t i = 0; i < MAX_ITERS; ++i)
	{
		a = skepu_userfunction_skepu_skel_0mandelbroter_add_c::CPU(skepu_userfunction_skepu_skel_0mandelbroter_mult_c::CPU(a, a), c);
		if ((a.a * a.a + a.b * a.b) > 9)
			return i;
	}
	return MAX_ITERS;
}
#undef SKEPU_USING_BACKEND_CPU
};

#line 100 "/home/joel/Documents/exjobb/skepu/skepu_fork/skepu/skepu-headers/src/skepu3/cluster/gpi/gpi_examples/mandelbrot_starpu.cpp"
int main(int argc, char* argv[])
{
	if (argc < 4)
	{
		if(!skepu::cluster::mpi_rank())
			std::cout << "Usage: " << argv[0] << " width height backend\n";
		exit(1);
	}

	const size_t width = std::stoul(argv[1]);
	const size_t height = std::stoul(argv[2]);
	auto spec = skepu::BackendSpec{skepu::Backend::typeFromString(argv[3])};
	skepu::setGlobalBackendSpec(spec);

	auto start = std::chrono::system_clock::now();

	skepu::Matrix<size_t> iterations(height, width);

	skepu::backend::Map<0, skepu_userfunction_skepu_skel_0mandelbroter_mandelbrot_f, bool, void> mandelbroter(false);
	mandelbroter(iterations, height, width);
	iterations.flush();

	starpu_task_wait_for_all();
	auto end = std::chrono::system_clock::now();
	double rtime = std::chrono::duration<double>{end - start}.count();

	printf("My time = %f, width = %lu, height = %lu\n",
		rtime, width, height);



	//if(!skepu::cluster::mpi_rank())
	//if(iterations.single())
	//	save_image(width, height, iterations.getAddress(), MAX_ITERS);


}
