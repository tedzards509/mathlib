#include <iostream>
#include <chrono>
#include "amathlib.h"
#include "aml_lua_binding.h"


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

	lua_State *L = luaL_newstate();
	luaL_openlibs(L);

	initAmlLua(L);
	int result = luaL_dofile(L, "../test.lua");
	if (result != LUA_OK) {
		std::cout << "ERROR" << std::endl;
	}

	lua_close(L);

	uint64_t beginTimer = std::chrono::steady_clock::now().time_since_epoch().count();

	for (int i = 0; i < 10; i++) {

		VectorDouble4D vecABench(1.0, 2.0, 3.0, 4.0);
		VectorDouble4D vecBBench(-20.0, 100.0, 80.0, 200.0);

		vecABench.add(vecBBench)->map(212.0, 32.0, 100.0, 0.0);

		volatile double a = vecABench[0];

	}


	uint64_t endTimer = std::chrono::steady_clock::now().time_since_epoch().count();
	uint64_t timeDelta = (endTimer - beginTimer);

	std::cout << timeDelta / 1000000000.0f << std::endl;

	return 0;

}

VectorDouble8D addV1(VectorDouble8D a, VectorDouble8D b) {
	VectorDouble8D ret;
	ret += a;
	ret += b;
	return ret;
}
