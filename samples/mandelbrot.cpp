//
// Created by af on 26.01.21.
//
#include <stdlib.h>
#include <stdint.h>
#include "../amathlib.h"
#include <iostream>
#include <chrono>

/*
 * This sample shows how to implement the mandelbrot set in C++. It is implemented in two ways a fast and a slow one.
 * It is recommended to know what the mandelbrot set is, not only because it makes it easy to follow, but because it is a
 * beautiful part of mathematics. Refer to https://en.wikipedia.org/wiki/Mandelbrot_set when needed.
 *
 * To compile on x86-64 choose clang or gcc and run
 * 		clang++ -O3 -ffast-math -march=native mandelbrot.cpp -o mandelbrot
 *	OR
 *		g++ -O3 -ffast-math -march=native mandelbrot.cpp -o mandelbrot
 *
 * "-O3" and "-ffast-math" are optimisation flags, you can try to run the code with just "-O0" or without -march=native
 *
 * On non x86-64 platforms remove -march=native
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
	std::cin >> width;
	std::cout << "height recommended 80" << std::endl;
	std::cin >> height;
	std::cout << "accuracy recommended 100000" << std::endl;
	std::cin >> accuracy;
	uint32_t begin = getTime();
	uint32_t end = 0;
	for (int x = 0; x < height; x++) {
		for (int y = 0; y < width; y++) {
			Complex64 c(((double) x) / ((double) height / 2.0f) - 1.5f, ((double) y) / ((double) width / 2.0f) - 1.0f);
			Complex64 z = c;
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
		for (int y = 0; y < width / IDEAL_COMPLEX_64_SIZE; y++) {
			IDEAL_COMPLEX_64_TYPE complex64_C;
			IDEAL_COMPLEX_64_TYPE complex64_Z;
			for (int index = 0; index < IDEAL_COMPLEX_64_SIZE; index++) {
				complex64_C.set(index, Complex64(((double) x) / ((double) height / 2.0f) - 1.5f,
												 (((double) y * IDEAL_COMPLEX_64_SIZE + index) /
												  ((double) width / 2.0f)) - 1.0f));
			}
			complex64_Z = complex64_C;
			IDEAL_COMPLEX_64_VECTOR_TYPE result(accuracy);
			int i = 0;
			IDEAL_COMPLEX_64_MASK_TYPE finished;
			IDEAL_COMPLEX_64_MASK_TYPE nowFinished;
			IDEAL_COMPLEX_64_MASK_TYPE alreadyFinished;
			for (; i < accuracy; ++i) {
				complex64_Z = (complex64_Z * complex64_Z) + complex64_C;
				finished = complex64_Z.abs_gt(2);
				bool anyFinished = finished.anyTrue();
				if (anyFinished) {
					alreadyFinished = finished;
					result.set(i, finished);
					break;
				}
			}
			if (!(finished.allTrue())) {
				for (; i < accuracy; ++i) {
					complex64_Z.multiply(complex64_Z, !finished)->add(complex64_C, !finished);
					//complex64_Z = (complex64_Z * complex64_Z) + complex64_C; Identical result but slower on my machine
					finished = complex64_Z.abs_gt(2);
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
			for (int index = 0; index < IDEAL_COMPLEX_64_SIZE; index++) {
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