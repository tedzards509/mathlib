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
 * It is recommended to know what the mandelbrot set is not only because it makes it easy to follow, but because it is a
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
	return std::chrono::steady_clock::now().time_since_epoch().count(); // this helps us to profile our code
}

int width = 250;           // example values
int height = 80;           // example values
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
	for (int x = 0; x < height; x++) {        // Imagine the console as an image with dimensions width X height.
		for (int y = 0; y < width; y++) {    // For each pixel we check weather it's in the set or not.
			Complex64 c(((double) x) / ((double) height / 2.0f) - 1.5f, ((double) y) / ((double) width / 2.0f) - 1.0f);
			Complex64 z = c; // Initializing the variables c and z by mapping pixel space onto the complex plane.
			int i = 0;
			int result = accuracy;
			for (; i < accuracy; ++i) {
				z = z * z + c;                // The central iteration / sequence
				if (z.abs_gt(2)) {            // check for diversion
					result = i;
					break;
				}
			}
			if (result >= accuracy) {    // check weather the sequence went through the whole loop
				std::cout << "#";    // c is part of the mandelbrot set
			} else {
				std::cout << " ";    // c is not
			}
		}
		std::cout << "\n";
	}
	end = getTime();

	std::cout << "Time 1 : " << ((double) (end - begin)) / 100000000.0f << std::endl; // display time
	//

	begin = getTime();
	for (int x = 0; x < height; x++) { // The same functionality as before, but now we group the pixel into groups of
		for (int y = 0; y < width / IDEAL_COMPLEX_64_SIZE; y++) {//IDEAL_COMPLEX_64_SIZE. This value is system dependent
			IDEAL_COMPLEX_64_TYPE complex64_C; // Again do initialising
			IDEAL_COMPLEX_64_TYPE complex64_Z; // Read through  the next lines until you see that it's basically the
			for (int index = 0; index < IDEAL_COMPLEX_64_SIZE; index++) { // same code
				complex64_C.set(index, Complex64(((double) x) / ((double) height / 2.0f) - 1.5f,
												 (((double) y * IDEAL_COMPLEX_64_SIZE + index) /
												  ((double) width / 2.0f)) - 1.0f));
			}
			complex64_Z = complex64_C;
			IDEAL_COMPLEX_64_VECTOR_TYPE result(accuracy); // equivalent to "int result = accuracy" from before
			int i = 0;
			IDEAL_COMPLEX_64_MASK_TYPE finished; // represents the finished indices of our pixels
			IDEAL_COMPLEX_64_MASK_TYPE nowFinished; // represents the newly finished indices of our pixels
			IDEAL_COMPLEX_64_MASK_TYPE alreadyFinished; // represents the already finished indices of our pixels
			for (; i < accuracy; ++i) { // iterate through the sequence until the first pixels is finished
				complex64_Z = (complex64_Z * complex64_Z) + complex64_C; // equivalent to "z = z * z + c;"
				finished = complex64_Z.abs_gt(2);
				bool anyFinished = finished.anyTrue();
				if (anyFinished) {
					alreadyFinished = finished;
					result.set(i, finished);
					break;
				}
			}
			if (!(finished.allTrue())) { // In the case not all pixel are finished
				for (; i < accuracy; ++i) {
					complex64_Z.multiply(complex64_Z, !finished)->add(complex64_C, !finished); // only iterate over the
					finished = complex64_Z.abs_gt(2); // pixels left (not finished)
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
			for (int index = 0; index < IDEAL_COMPLEX_64_SIZE; index++) { // Now that you have accomplished the hard
				if (result[index] >= accuracy) { // part, we just have to print the result
					std::cout << "#";
				} else {
					std::cout << " ";
				}
			}
		}
		std::cout << "\n";
	}
	end = getTime();

	std::cout << "Time 2 : " << ((double) (end - begin)) / 100000000.0f << std::endl;
}