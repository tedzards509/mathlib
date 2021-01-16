//
// Created by af on 16.01.21.
//

#ifndef MATH_LIB_A_MATH_LIB_H
#define MATH_LIB_A_MATH_LIB_H


#ifdef X86_64
#define USE_AVX
#define USE_SSE

#if defined(USE_AVX)

#include <immintrin.h>

#define USE_SSE

#endif
#ifdef USE_SSE

#include <emmintrin.h>

#endif
#endif

#include <cmath>


union doublevec4 {
#ifdef USE_AVX
	__m256d avx;
#endif
#ifdef USE_SSE
	__m128d sse[2];
#endif
	double c[4];
};

union doublemat4x4 {
#ifdef USE_AVX
	__m256d avx[4];
#endif
#ifdef USE_SSE
	__m128d sse[8];
#endif
	double c[16];
};

union doublevec2 {
#ifdef USE_SSE
	__m128d sse;
#endif
	double c[2];
};


class VectorDouble4D {
private:

public:
	doublevec4 v{};

	inline double operator[](uint32_t position) {
		return v.c[position];//TODO maybe error handling
	}

	inline void operator+=(VectorDouble4D vec2) {
#if defined(USE_AVX)
		v.avx = _mm256_add_pd(v.avx, vec2.v.avx);
#elif defined(USE_SSE) // SSE2
		v.sse[0] = _mm_add_pd(v.sse[0], vec2.v.sse[0]);
		v.sse[1] = _mm_add_pd(v.sse[1], vec2.v.sse[1]);
#else
		v.c[0] += vec2[0];
		v.c[1] += vec2[1];
		v.c[2] += vec2[2];
		v.c[3] += vec2[3];
#endif


	}

	inline VectorDouble4D operator+(VectorDouble4D vec2) {
#if defined(USE_AVX)
		VectorDouble4D ret;
		ret.v.avx = _mm256_add_pd(v.avx, vec2.v.avx);
		return ret;
#elif defined(USE_SSE) // SSE2
		VectorDouble4D ret;
		ret.v.sse[0] = _mm_add_pd(v.sse[0], vec2.v.sse[0]);
		ret.v.sse[1] = _mm_add_pd(v.sse[1], vec2.v.sse[1]);
#else
		VectorDouble4D ret(v.c[0] + vec2.v.c[0], v.c[1] + vec2.v.c[1], v.c[2] + vec2.v.c[2], v.c[3] + vec2.v.c[3]);
		return ret;
#endif


	}

	inline VectorDouble4D operator+(double a) {
#if defined(USE_AVX)
		VectorDouble4D ret(a);
		ret.v.avx = _mm256_add_pd(v.avx, ret.v.avx);
		return ret;
#elif defined(USE_SSE) // SSE2
		VectorDouble4D ret(a);
		ret.v.sse[0] = _mm_add_pd(v.sse[0], ret.v.sse[0]);
		ret.v.sse[1] = _mm_add_pd(v.sse[1], ret.v.sse[1]);
#else
		VectorDouble4D ret(v.c[0] +a, v.c[1] + a, v.c[2] + a, v.c[3] + a);
		return ret;
#endif


	}

	inline VectorDouble4D *add(VectorDouble4D a) {
#if defined(USE_AVX)
		v.avx = _mm256_add_pd(v.avx, a.v.avx);
#elif defined(USE_SSE) // SSE2
		v.sse[0] = _mm_add_pd(v.sse[0], a.v.sse[0]);
		v.sse[1] = _mm_add_pd(v.sse[1], a.v.sse[1]);
#else
		v.c[0] += a.v.c[0]; v.c[1] += a.v.c[1]; v.c[2] += a.v.c[2]; v.c[3] += a.v.c[3];
#endif
		return this;
	}

	inline void inverse() {

#if defined(USE_AVX)
		doublevec4 a = {0.0f, 0.0f, 0.0f, 0.0f};
		v.avx = _mm256_sub_pd(a.avx, v.avx);
#elif defined(USE_SSE) // SSE2
		double a[2] = {0.0f, 0.0f};
		__m128d b = _mm_loadu_pd(a);
		v.sse[0] = _mm_sub_pd(b, v.sse[0]);
		v.sse[1] = _mm_sub_pd(b, v.sse[1]);
#else
		v.c[0] = 0 - v.c[0];
		v.c[1] = 0 - v.c[1];
		v.c[2] = 0 - v.c[2];
		v.c[3] = 0 - v.c[3];
#endif
	}

	inline void operator-=(VectorDouble4D vec2) {
#if defined(USE_AVX)
		v.avx = _mm256_sub_pd(v.avx, vec2.v.avx);
#else
		v.c[0] -= vec2[0];
		v.c[1] -= vec2[1];
		v.c[2] -= vec2[2];
		v.c[3] -= vec2[3];
#endif


	}


	inline VectorDouble4D operator-(VectorDouble4D vec2) {
#if defined(USE_AVX)
		VectorDouble4D ret;
		ret.v.avx = _mm256_sub_pd(v.avx, vec2.v.avx);
		return ret;
#else
		VectorDouble4D ret(v.c[0] - vec2.v.c[0], v.c[1] - vec2.v.c[1], v.c[2] - vec2.v.c[2], v.c[3] - vec2.v.c[3]);
		return ret;
#endif


	}

	inline VectorDouble4D operator-(double a) {
		VectorDouble4D ret(a);
#if defined(USE_AVX)
		ret.v.avx = _mm256_sub_pd(v.avx, ret.v.avx);
		return ret;
#else
		ret = VectorDouble4D(v.c[0] - a, v.c[1] - a, v.c[2] - a, v.c[3] - a);
		return ret;
#endif


	}

	inline double length() {
		return sqrt(v.c[0] * v.c[0] + v.c[1] * v.c[1] + v.c[2] * v.c[2] + v.c[3] * v.c[3]);
	}

	inline void normalize() {
		//TODO check length==0
		double vecLength = 1 / length();
		v.c[0] *= vecLength;
		v.c[1] *= vecLength;
		v.c[2] *= vecLength;
		v.c[3] *= vecLength;
	}

	inline VectorDouble4D *forEachSin() {
#if defined(USE_AVX) && defined(__INTEL_COMPILER)
		v.avx = _mm256_sin_pd(v.avx);
#else
		v.c[0] = sin(v.c[0]);
		v.c[1] = sin(v.c[1]);
		v.c[2] = sin(v.c[2]);
		v.c[3] = sin(v.c[3]);
#endif
		return this;
	}

	inline void operator*=(double scalar) {
		v.c[0] *= scalar;
		v.c[1] *= scalar;
		v.c[2] *= scalar;
		v.c[3] *= scalar;
	}

	// for each multiply
	inline void operator*=(VectorDouble4D vec2) {
#if defined(USE_AVX)
		v.avx = _mm256_mul_pd(v.avx, vec2.v.avx);
#else
		v.c[0] *= vec2.v.c[0];
		v.c[1] *= vec2.v.c[1];
		v.c[2] *= vec2.v.c[2];
		v.c[3] *= vec2.v.c[3];
#endif
	}

	inline VectorDouble4D *forEachSqrt() {
#if defined(USE_AVX)
		v.avx = _mm256_sqrt_pd(v.avx);
#elif defined(USE_SSE)
		v.sse[0] = _mm_sqrt_pd(v.sse[0]);
		v.sse[1] = _mm_sqrt_pd(v.sse[1]);
#else
		v.c[0] = sqrt(v.c[0]);
		v.c[1] = sqrt(v.c[1]);
		v.c[2] = sqrt(v.c[2]);
		v.c[3] = sqrt(v.c[3]);
#endif
		return this;

	}

	inline VectorDouble4D operator*(double a) {
		VectorDouble4D ret;
#if defined(USE_AVX)
		doublevec4 b = {a, a, a, a};
		ret.v.avx = _mm256_mul_pd(v.avx, b.avx);
#else
		ret.c[0] *= v.c[0] *a;
		ret.c[1] *= v.c[1] * a ;
		ret.c[2] *= v.c[2] * a;
		ret.c[3] *= v.c[3] * a;
#endif
		return ret;
	}

	inline VectorDouble4D operator/(double a) {
		VectorDouble4D ret;
#if defined(USE_AVX)
		doublevec4 b = {a, a, a, a};
		ret.v.avx = _mm256_div_pd(v.avx, b.avx);
#else
		ret.c[0] /= v.c[0] *a;
		ret.c[1] /= v.c[1] * a ;
		ret.c[2] /= v.c[2] * a;
		ret.c[3] /= v.c[3] * a;
#endif
		return ret;
	}


	inline VectorDouble4D *capBetween1_0() {
#if defined(USE_AVX)
		doublevec4 one = {1.0f, 1.0f, 1.0f, 1.0f};
		v.avx = _mm256_max_pd(v.avx, one.avx);
		v.avx = _mm256_min_pd(v.avx, _mm256_setzero_pd());
#else
		if(v.c[0] > 1){
			v.c[0] = 1;
		} else if(v.c[0] < 0){
			v.c[0] = 0;
		}
		if(v.c[1] > 1){
			v.c[1] = 1;
		} else if(v.c[1] < 0){
			v.c[1] = 0;
		}
		if(v.c[2] > 1){
			v.c[2] = 1;
		} else if(v.c[2] < 0){
			v.c[2] = 0;
		}
		if(v.c[3] > 1){
			v.c[3] = 1;
		} else if(v.c[3] < 0){
			v.c[3] = 0;
		}
#endif
		return this;
	}

	inline VectorDouble4D *capBetweenX_Y(double upperBoundary, double lowerBoundary) {
#if defined(USE_AVX)
		doublevec4 upper = {upperBoundary, upperBoundary, upperBoundary, upperBoundary};
		doublevec4 lower = {lowerBoundary, lowerBoundary, lowerBoundary, lowerBoundary};
		v.avx = _mm256_max_pd(v.avx, upper.avx);
		v.avx = _mm256_min_pd(v.avx, lower.avx);
#else
		if(v.c[0] > upperBoundary){
			v.c[0] = upperBoundary;
		} else if(v.c[0] < lowerBoundary){
			v.c[0] = lowerBoundary;
		}
		if(v.c[1] > upperBoundary){
			v.c[1] = upperBoundary;
		} else if(v.c[1] < lowerBoundary){
			v.c[1] = lowerBoundary;
		}
		if(v.c[2] > upperBoundary){
			v.c[2] = upperBoundary;
		} else if(v.c[2] < lowerBoundary){
			v.c[2] = lowerBoundary;
		}
		if(v.c[3] > upperBoundary){
			v.c[3] = upperBoundary;
		} else if(v.c[3] < lowerBoundary){
			v.c[3] = lowerBoundary;
		}
#endif
		return this;
	}

	inline VectorDouble4D *map(double upperInput, double lowerInput, double upperOutput, double lowerOutput) {
#if defined(USE_AVX)
		doublevec4 a = {lowerInput, lowerInput, lowerInput, lowerInput};
		a.avx = _mm256_sub_pd(v.avx, a.avx);
		double factor = (upperOutput - lowerOutput) / (upperInput - lowerInput);
		doublevec4 b = {factor, factor, factor, factor};
		a.avx = _mm256_mul_pd(a.avx, b.avx);
		doublevec4 c = {lowerOutput, lowerOutput, lowerOutput, lowerOutput};
		v.avx = _mm256_add_pd(a.avx, c.avx);
#else
		v.c[0] =((v.c[0] - lowerInput) * ((upperOutput - lowerOutput) / (upperInput - lowerInput))) + lowerOutput;
		v.c[1] =((v.c[1] - lowerInput) * ((upperOutput - lowerOutput) / (upperInput - lowerInput))) + lowerOutput;
		v.c[2] =((v.c[2] - lowerInput) * ((upperOutput - lowerOutput) / (upperInput - lowerInput))) + lowerOutput;
		v.c[3] =((v.c[3] - lowerInput) * ((upperOutput - lowerOutput) / (upperInput - lowerInput))) + lowerOutput;
#endif
		return this;
	}


	inline VectorDouble4D(double a, double b, double c, double d) {
		v.c[0] = a;
		v.c[1] = b;
		v.c[2] = c;
		v.c[3] = d;
	}

	inline VectorDouble4D() {
		v.c[0] = 0;
		v.c[1] = 0;
		v.c[2] = 0;
		v.c[3] = 0;
	}

	inline explicit VectorDouble4D(double a) {
		v.c[0] = a;
		v.c[1] = a;
		v.c[2] = a;
		v.c[3] = a;
	}

#ifdef USE_AVX

	inline explicit VectorDouble4D(__m256d a) {
		v.avx = a;
	}

#endif


};

class MatrixDouble4X4 {
	doublemat4x4 m;

	inline VectorDouble4D operator[](uint32_t column) {
		VectorDouble4D ret;
#if defined(USE_AVX)
		ret.v.avx = m.avx[column];
#else
		ret.v.c[0] = m.c[column*4];
		ret.v.c[1] = m.c[column*4+1];
		ret.v.c[2] = m.c[column*4+2];
		ret.v.c[3] = m.c[column*4+3];
#endif
		return ret;//TODO maybe error handling
	}

	inline MatrixDouble4X4 *identity() {
		m = {{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0}};
	}
};


#endif //MATH_LIB_A_MATH_LIB_H
