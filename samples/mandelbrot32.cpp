//
// Created by af on 02.02.21.

#define AML_USE_STD_COMPLEX

#include <stdint.h>
#include "../amathlib.h"
#include <iostream>
#include <chrono>

/*
 * This sample shows how to implement the mandelbrot set in C++. It is implemented in two ways a fast and a slow one.
 * It is recommended to know what the mandelbrot set is, not only because it makes it easy to follow, but because it is a
 * beautiful part of mathematics. Refer to https://en.wikipedia.org/wiki/Mandelbrot_set when needed.
 *
 * To compile on x86-32 choose clang or gcc and run
 * 		clang++ -O3 -ffast-math -march=native mandelbrot.cpp -o mandelbrot
 *	OR
 *		g++ -O3 -ffast-math -march=native mandelbrot.cpp -o mandelbrot
 *
 * "-O3" and "-ffast-math" are optimisation flags, you can try to run the code with just "-O0" or without -march=native
 *
 * On non x86-32 platforms remove -march=native
 */

uint64_t getTime() {
	return std::chrono::duration_cast<std::chrono::milliseconds>(
			std::chrono::system_clock::now().time_since_epoch()).count();
}

int width = 250;          // example values
int height = 80;          // example values
int accuracy = 100000;    // example values way to high, just to get useful timing
int main() {
	std::cout << "width recommended 250" << std::endl;
	std::string input;
	//std::getline(std::cin, input);
	//if (!input.empty()) {
	//	std::istringstream stream(input);
	//	stream >> width;
	//}
	std::cout << "height recommended 80" << std::endl;
	//std::getline(std::cin, input);
	//if (!input.empty()) {
	//	std::istringstream stream(input);
	//	stream >> height;
	//}
	std::cout << "accuracy recommended 100000" << std::endl;
	//std::getline(std::cin, input);
	//if (!input.empty()) {
	//	std::istringstream stream(input);
	//	stream >> accuracy;
	//}
	uint64_t begin;
	uint64_t end;

	begin = getTime();
	for (int x = 0; x < height; x++) {
		for (int y = 0; y < width; y++) {
			std::complex<float> c(((double) x) / ((double) height / 2.0f) - 1.5f,
								  ((double) y) / ((double) width / 2.0f) - 1.0f);
			std::complex<float> z = c;
			int i = 0;
			int result = accuracy;
			for (; i < accuracy; ++i) {
				z = z * z + c;
				if (norm(z) > 4) {
					result = i;
					break;
				}
			}
			if (result >= accuracy) {
				std::cout << "#";
			} else {
				std::cout << " ";
			}
		}
		std::cout << "\n";
	}
	end = getTime();

	std::cout << "Time 0 : " << ((double) (end - begin)) / 1000.0f << std::endl;

	begin = getTime();
	for (int x = 0; x < height; x++) {
		for (int y = 0; y < width; y++) {
			Complex32 c(((double) x) / ((double) height / 2.0f) - 1.5f, ((double) y) / ((double) width / 2.0f) - 1.0f);
			Complex32 z = c;
			int i = 0;
			int result = accuracy;
			for (; i < accuracy; ++i) {
				z = z * z + c;
				if (z.abs_gt(2)) {
					result = i;
					break;
				}
			}
			if (result >= accuracy) {
				std::cout << "#";
			} else {
				std::cout << " ";
			}
		}
		std::cout << "\n";
	}
	end = getTime();

	std::cout << "Time 1 : " << ((double) (end - begin)) / 1000.0f << std::endl;
	//

	begin = getTime();
	for (int x = 0; x < height; x++) {
		for (int y = 0; y < width / MIN_COMPLEX_32_SIZE; y++) {
			MIN_COMPLEX_32_TYPE complex32_C;
			MIN_COMPLEX_32_TYPE complex32_Z;
			for (int index = 0; index < MIN_COMPLEX_32_SIZE; index++) {
				complex32_C.set(index, Complex32(AML::mapLinear((float) x, 0.0f, (float) height, -1.5f, 0.5f),
												 AML::mapLinear((float) y * MIN_COMPLEX_32_SIZE + index, 0.0f,
																(float) width,
																-1.0f, 1.0f)));
			}
			complex32_Z = complex32_C;
			MIN_COMPLEX_32_VECTOR_TYPE result(accuracy);
			int i = 0;
			MIN_COMPLEX_32_MASK_TYPE finished;
			MIN_COMPLEX_32_MASK_TYPE nowFinished;
			MIN_COMPLEX_32_MASK_TYPE alreadyFinished;
			for (; i < accuracy; ++i) {
				complex32_Z = (complex32_Z * complex32_Z) + complex32_C;
				finished = complex32_Z.abs_gt(2);
				bool anyFinished = finished.anyTrue();
				if (anyFinished) {
					alreadyFinished = finished;
					result.set(i, finished);
					break;
				}
			}
			if (!(finished.allTrue())) {
				for (; i < accuracy; ++i) {
					complex32_Z.square(!finished)->add(complex32_C, !finished);
					//complex32_Z = (complex32_Z * complex32_Z) + complex32_C; Identical result but slower on my machine
					finished = complex32_Z.abs_gt(2);
					nowFinished = finished && !alreadyFinished;
					if (nowFinished.anyTrue()) {
						alreadyFinished = finished;
						result.set(i, nowFinished);
					}
					if (finished.allTrue()) {
						break;
					}
				}
			}
			for (int index = 0; index < MIN_COMPLEX_32_SIZE; index++) {
				if (result[index] >= accuracy) {
					std::cout << "#";
				} else {
					std::cout << " ";
				}
			}
		}
		std::cout << "\n";
	}
	end = getTime();

	std::cout << "Time 2 : " << ((double) (end - begin)) / 1000.0f << std::endl;
}