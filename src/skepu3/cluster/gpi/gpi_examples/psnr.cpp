/*!
* PSNR (Peak Signal to Noise Ratio): The PSNR represents a measure of the peak error between the compressed and the original image.
* It is clsely related to MSE which represents the cumulative squared error between the compressed and the original image.
*/

#include <iostream>
#include <cmath>

#include <skepu>

[[skepu::userconstant]] constexpr int
	MAX = 255,
	NOISE = 10;

float diff_squared(int a, int b)
{
	return (a - b) * (a - b);
}

template<typename T>
T sum(T a, T b)
{
	return a + b;
}

template<typename T>
T clamp_sum(T a, T b)
{
	T temp = a + b;
	return temp < 0 ? 0 : (temp > MAX ? MAX : temp);
}

float psnr(skepu::Matrix<int> &img, skepu::Matrix<int>& noise)
{
	const size_t rows = img.size_i();
	const size_t cols = img.size_j();

	auto clamped_sum = skepu::Map(clamp_sum<int>);
	//auto squared_diff_sum = skepu::MapReduce(diff_squared, sum<float>);
	skepu::Matrix<int> comp_img(rows, cols);

	skepu::Matrix<float> swap_space(rows, cols);

	auto diff_sqr = skepu::Map(diff_squared);
	auto sum_float = skepu::Reduce(sum<float>);

	// Add noise
	clamped_sum(comp_img, img, noise);

	diff_sqr(swap_space, img, noise);
	float squared_diff_sum = sum_float(swap_space);


//	std::cout << "Compressed image: " << comp_img << "\n";

	float mse = squared_diff_sum / (rows * cols);
	return 10 * log10((MAX * MAX) / mse);
}

int main(int argc, char *argv[])
{
	if (argc < 4)
	{
		/*
		skepu::external([&]{
			std::cout << "Usage: " << argv[0] << " rows cols backend\n";});
		*/
		exit(1);
	}

	const size_t rows = std::stoul(argv[1]);
	const size_t cols = std::stoul(argv[2]);
	auto spec = skepu::BackendSpec{skepu::Backend::typeFromString(argv[3])};
	skepu::setGlobalBackendSpec(spec);

	skepu::Matrix<int> img(rows, cols), noise(rows, cols);

	// Generate random image and random noise
	img.randomize(0, MAX);
	noise.randomize(-NOISE, NOISE);

	/*
	skepu::external(
		skepu::read(img),
		[&]{
			std::cout << "Actual image: " << img << "\n";
		});
		*/

	float psnrval = psnr(img, noise);

	std::cout << "psnr = " << psnrval << std::endl;

	/*
	skepu::external([&]{
		std::cout << "PSNR of two images: " << psnrval << "\n";});
	*/
	return 0;
}
