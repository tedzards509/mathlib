#include <iostream>
#define USE_SSE42
#include "amathlib.h"


int main() {

#if defined(__NVCC__)
	std::cout << "USE CUDA" << std::endl;
#endif
#if defined(__EMSCRIPTEN__)
	std::cout << "WEB" << std::endl;
#endif

#if defined(DEBUG)
	std::cout << "debug" << std::endl;
#endif

	std::cout << "size : " << sizeof(VectorDouble4D) << " align : " << alignof(VectorDouble4D) << std::endl;
	std::cout << "size double : " << sizeof(double) << " size float : " << sizeof(float) << std::endl;
	VectorDouble4D vector;
	VectorDouble4D vector2(0.4, 21.0, 2.4, 2.5);
	vector += vector2;
	std::cout << vector[0] << " " << vector[1] << " " << vector[2] << " " << vector[3] << std::endl;
	vector.normalize();
	std::cout << vector[0] << " " << vector[1] << " " << vector[2] << " " << vector[3] << std::endl;
	vector.inverse();
	std::cout << vector[0] << " " << vector[1] << " " << vector[2] << " " << vector[3] << std::endl;
	VectorDouble4D vector3(2);
	vector *= vector3;
	vector.inverse();
	std::cout << vector[0] << " " << vector[1] << " " << vector[2] << " " << vector[3] << std::endl;

	vector.forEachSqrt();
	std::cout << vector[0] << " " << vector[1] << " " << vector[2] << " " << vector[3] << std::endl;

	vector.forEachSin()->add({0.0, 0.0, 0.0, 0.0});
	std::cout << vector[0] << " " << vector[1] << " " << vector[2] << " " << vector[3] << std::endl;


	double pixel = 100;
	VectorDouble4D *image = new VectorDouble4D[(int) pixel];
	for (int i = 0; i < (int) pixel; i++) {
		image[i] = VectorDouble4D(1.0 + i, 2.0 + i, 3.0 + i, 4.0 + i);
	}
	//get average per channel;
	VectorDouble4D sum;
	for (int i = 0; i < (int) pixel; i++) {
		sum += image[i];
	}
	sum *= (1.0 / pixel);
	std::cout << sum[0] << " " << sum[1] << " " << sum[2] << " " << sum[3] << std::endl;

	VectorDouble4D x(-20.0, 100.0, 80.0, 200.0);

	std::cout << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << std::endl;

	x.map(212.0, 32.0, 100.0, 0.0);// fahrenheit to celsius

	std::cout << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << std::endl;

	x = VectorDouble4D(6.0);

	std::cout << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << std::endl;

	x.map(6.56168, 3.28084, 200.0000064, 100.0000032);// feet to cm

	std::cout << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << std::endl;


	x = VectorDouble4D(0.0, 0.3333, 0.5, 1.0);

	std::cout << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << std::endl;

	x.map(1.0, 0.0, 255, 0.0);// NORM RGB -> 8-bit RGB

	std::cout << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << std::endl;


	VectorDouble8D vector8_1(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0);
	VectorDouble8D vector8_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);

	vector8_1 += vector8_2;

	std::cout << vector8_1[0] << " " << vector8_1[1] << " " << vector8_1[2] << " " << vector8_1[3] << " "
			  << vector8_1[4] << " " << vector8_1[5] << " " << vector8_1[6] << " " << vector8_1[7] << std::endl;

	vector8_1 += vector8_1;

	std::cout << vector8_1[0] << " " << vector8_1[1] << " " << vector8_1[2] << " " << vector8_1[3] << " "
			  << vector8_1[4] << " " << vector8_1[5] << " " << vector8_1[6] << " " << vector8_1[7] << std::endl;

	Complex64 cmp1(2, 4);
	Complex64 cmp2(4, 1);

	cmp2.divide(cmp1);

	std::cout << cmp2.real() << " + " << cmp2.imaginary() << "i" << std::endl;


	cmp2.multiply(cmp1);

	std::cout << cmp2.real() << " + " << cmp2.imaginary() << "i" << std::endl;

	cmp2.sqrt();
	std::cout << cmp2.real() << " + " << cmp2.imaginary() << "i" << std::endl;

	cmp2.multiply(cmp2);

	std::cout << cmp2.real() << " + " << cmp2.imaginary() << "i" << std::endl;


	MatrixDouble4X4 m;
	m.identity();

	VectorDouble4D matVec = m[3];

	std::cout << matVec[0] << " " << matVec[1] << " " << matVec[2] << " " << matVec[3] << std::endl;

	MatrixDouble4X4 m2 = {{1,  2,  3,  4},
						  {5,  6,  7,  8},
						  {9,  10, 11, 12},
						  {13, 14, 15, 16}};

	//MatrixDouble4X4 m2;

	//m2.identity();

	MatrixDouble4X4 m3 = m2 * m2;

	matVec = m3[3];

	std::cout << (float) matVec[0] << " " << (float) matVec[1] << " " << (float) matVec[2] << " " << (float) matVec[3]
			  << std::endl;

	m2.identity();

	m3 = m3 * m2;

	matVec = m3[3];

	std::cout << (float) matVec[0] << " " << (float) matVec[1] << " " << (float) matVec[2] << " " << (float) matVec[3]
			  << std::endl;

	VectorDouble4D matVec2 = {0, 1, 2, 3};

	matVec2 = m3 * matVec2;

	std::cout << (float) matVec2[0] << " " << (float) matVec2[1] << " " << (float) matVec2[2] << " "
			  << (float) matVec2[3] << std::endl;

	return 0;

}
