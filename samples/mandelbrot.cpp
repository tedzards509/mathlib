//
// Created by af on 26.01.21.
//

#include <stdlib.h>
#include <stdint.h>
#include "../amathlib.h"
#include <iostream>
#include <chrono>


uint64_t getTime() {
	return std::chrono::steady_clock::now().time_since_epoch().count();
}


int main() {
	uint32_t totalTime = 0;
	uint32_t begin = getTime();
	uint32_t end = 0;
	for (int x = 0; x < 40; x++) {
		for (int y = 0; y < 200; y++) {
			Complex64 c(((double) x) / 20.0f - 1.5f, ((double) y) / 100.0f - 1.0f);
			Complex64 z = c;
			int i = 0;
			int result = 1000000;
			for (; i < 1000000; ++i) {
				z = z * z + c;
				if (z.abs_gt(2)) {
					result = i;
					break;
				}
			}
			if (result >= 1000000) {
				std::cout << "#";
			} else {
				std::cout << " ";
			}
		}
		std::cout << std::endl;
	}
	end = getTime();

	std::cout << "Time 1 : " << ((double) (end - begin)) / 100000000.0f << std::endl;
	//

	begin = getTime();
	for (int x = 0; x < 40; x++) {
		for (int y = 0; y < 200 / IDEAL_COMPLEX_64_SIZE; y++) {
			IDEAL_COMPLEX_64_TYPE complex64_C;
			IDEAL_COMPLEX_64_TYPE complex64_Z;
			for (int index = 0; index < IDEAL_COMPLEX_64_SIZE; index++) {
				complex64_C.set(index, Complex64(((double) x) / 20.0f - 1.5f,
												 ((double) y * IDEAL_COMPLEX_64_SIZE + index) / 100.0f - 1.0f));
			}
			complex64_Z = complex64_C;
			IDEAL_COMPLEX_64_VECTOR_TYPE result(1000000);
			int i = 0;
			IDEAL_COMPLEX_64_MASK_TYPE mask;
			IDEAL_COMPLEX_64_MASK_TYPE mask2;
			IDEAL_COMPLEX_64_MASK_TYPE oldMask;
			for (; i < 1000000; ++i) {
				complex64_Z = (complex64_Z * complex64_Z) + complex64_C;
				mask = complex64_Z.abs_gt(2);
				bool anyFinished = mask.anyTrue();
				if (anyFinished) {
					oldMask = mask;
					result.set(i, mask);
					break;
				}
			}
			if (!(mask.allTrue())) {
				for (; i < 1000000; ++i) {
					complex64_Z.multiply(complex64_Z, !mask)->add(complex64_C, !mask);
					mask = complex64_Z.abs_gt(2);
					mask2 = mask && !oldMask;
					if (mask2.anyTrue()) {
						oldMask = mask;
						result.set(i, mask2);

					}
					if (mask.allTrue()) {
						break;
					}
				}
			}
			for (int index = 0; index < IDEAL_COMPLEX_64_SIZE; index++) {
				if (result[index] >= 1000000) {
					std::cout << "#";
				} else {
					std::cout << " ";
				}
			}
		}
		std::cout << std::endl;
	}
	end = getTime();

	std::cout << "Time 1 : " << ((double) (end - begin)) / 100000000.0f << std::endl;
}