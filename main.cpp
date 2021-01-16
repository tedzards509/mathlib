#include <iostream>
#include "amathlib.h"


int main() {
	std::cout << "size : " << sizeof(VectorDouble4D) << " align : " << alignof(VectorDouble4D) << std::endl;
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

	return 0;
}
