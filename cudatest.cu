//
// Created by af on 25.01.21.
//
#include "amathlib.h"
#include <iostream>

__global__ void kernel(double height, double width) {
	MatrixDouble4X4 m;
	//m.identity();
	//MatrixDouble4X4 m2({height, width, height, width}, {height, width, height, width},
	//				   {height, width, height, width}, { height, width, height, width });
	//MatrixDouble4X4 m3 = m2 * m;
}


int main() {

#if defined(__global__)
	std::cout << "cuda" << std::endl;

#endif
#if defined(__device__)
	std::cout << "cuda" << std::endl;

#endif

	// Run kernel
	dim3 blockDim(1, 1, 1);
	dim3 gridDim((1 + blockDim.x - 1) / blockDim.x, (1 + blockDim.y - 1) / blockDim.y, 1);
	kernel<<< gridDim, blockDim, 0 >>>(1, 1);

}
