//
// Created by af on 16.01.21.
//

#ifndef MATH_LIB_A_MATH_LIB_H
#define MATH_LIB_A_MATH_LIB_H

#define DEBUG_TO_INDEX(row, column) ((column - 1) * 4 + (row-1))

//#define USE_AVX512
//#define USE_AVX
//#define USE_SSE

//#define USE_SSE
//#define USE_SSE2
//#define USE_SSE3
//#define USE_SSE41
//#define USE_SSE42

//#define USE_AVX
//#define USE_AVX2
//#define USE_FMA


//#define USE_AVX512
//#define USE_AVX512F
//#define USE_KNC

#ifndef DEBUG

#if defined(__x86_64__)
#ifdef __SSE__
#define USE_SSE
#endif


#ifdef __SSE2__
#define USE_SSE2
#endif


#ifdef __SSE2__
#define USE_SSE2
#endif


#ifdef __SSE4_1__
#define USE_SSE41
#endif


#ifdef __SSE4_2__
#define USE_SSE42
#endif

#ifdef __AVX__
#define USE_AVX
#endif


#ifdef __AVX2__
#define USE_AVX2
#endif


#ifdef __FMA__
#define USE_FMA
#endif
#endif // INTEL

#endif //NDEBUG

#if defined(DEBUG)
#if defined(__x86_64__)
#if defined(USE_AVX512F)
#define USE_AVX512
#endif
#if defined(USE_AVX512)
#define USE_AVX512F
#define USE_FMA
#endif
#if defined(USE_FMA)
#define USE_AVX2
#endif
#if defined(USE_AVX2)
#define USE_AVX
#endif
#if defined(USE_AVX)
#define USE_SSE42
#endif
#if defined(USE_SSE42)
#define USE_SSE41
#endif
#if defined(USE_SSE41)
#define USE_SSE3
#endif
#if defined(USE_SSE3)
#define USE_SSE2
#endif
#if defined(USE_SSE2)
#define USE_SSE1
#endif
#if defined(USE_SSE1)
#define USE_SSE
#endif

#endif // DEBUG

#if defined(USE_AVX512)
#include <immintrin.h>
#endif


#if defined(USE_AVX)
#include <immintrin.h>
#endif
#ifdef USE_SSE

#include <emmintrin.h>

#endif //INTEL
#endif //DEBUG

#if defined(__EMSCRIPTEN__)
#include <wasm_simd128.h>
#endif

#if defined(__EMSCRIPTEN__)
#include <wasm_simd128.h>
#define USE_WASM_SIMD
#endif

#if defined(__ARM_NEON)

#include <arm_neon.h>

#define USE_NEON

#include <math.h>
#include <stdint.h>

#else

#include <cstdint>

#include <cmath>

#endif
union doublevec4 {
#ifdef USE_AVX
	__m256d avx;
#endif
#ifdef USE_SSE
	__m128d sse[2];
#endif
#ifdef USE_NEON
	float64x2_t neon[2];
#endif
	double c[4];
};

union doublevec8 {
#ifdef USE_AVX512
	__m512d avx512;
#endif
#ifdef USE_AVX
	__m256d avx[2];
#endif
#ifdef USE_SSE
	__m128d sse[4];
#endif
#ifdef USE_NEON
	float64x2_t neon[4];
#endif
	double c[8];
};

union doublemat4x4 {
#ifdef USE_AVX512
	__m512d avx512[2];
#endif
#ifdef USE_AVX
	__m256d avx[4];
#endif
#ifdef USE_SSE
	__m128d sse[8];
#endif
#ifdef USE_NEON
	float64x2_t neon[8];
#endif
	double c[16];
};

union doublevec2 {
#ifdef USE_SSE
	__m128d sse;
#endif
#ifdef USE_NEON
	float64x2_t neon;
#endif
	double c[2];
};

union u8vec2 {
	uint8_t c[2];
};


union u8vec3 {
	uint8_t c[3];
};


union u8vec4 {
	uint8_t c[4];
};


union u8vec8 {
	uint8_t c[8];
};


union u8vec16 {
#if defined(USE_SSE)
	__m128i sse;
#endif
	uint8_t c[16];
};


union u8vec32 {
#if defined(USE_AVX)
	__m256i avx;
#endif
#if defined(USE_SSE)
	__m128i sse[4];
#endif
	uint8_t c[32];
};


union u8vec64 {

#if defined(USE_AVX512)
	__m512i avx512;
#endif
#if defined(USE_AVX)
	__m256i avx[2];
#endif
#if defined(USE_SSE)
	__m128i sse[8];
#endif
	uint8_t c[64];
};

union u16vec2 {
	uint16_t c[2];
};


union u16vec3 {
	uint16_t c[3];
};


union u16vec4 {
#if defined(USE_SSE)
	__m128i sse;
#endif
	uint16_t c[4];
};


union u16vec8 {
#if defined(USE_SSE)
	__m128i sse[2];
#endif
	uint16_t c[8];
};


union u16vec16 {
#if defined(USE_SSE)
	__m128i sse[4];
#endif
	uint16_t c[16];
};


union u16vec32 {
#if defined(USE_AVX512)
	__m512i avx512;
#endif
#if defined(USE_AVX)
	__m256i avx[2];
#endif
#if defined(USE_SSE)
	__m128i sse[8];
#endif
	uint16_t c[32];
};


union u16vec64 {

#if defined(USE_AVX512)
	__m512i avx512[2];
#endif
#if defined(USE_AVX)
	__m256i avx[4];
#endif
#if defined(USE_SSE)
	__m128i sse[16];
#endif
	uint16_t c[64];
};

class VectorU16_2D {
public:
	u16vec2 v;

	VectorU16_2D() {
		v.c[0] = 0;
		v.c[1] = 0;
	}


	VectorU16_2D(uint16_t a, uint16_t b) {
		v.c[0] = a;
		v.c[1] = b;
	}
};

class VectorU8_4D {
public:
	u8vec4 v;

	VectorU8_4D() {
		v.c[0] = 0;
		v.c[1] = 0;
	}


	VectorU8_4D(uint8_t a, uint8_t b, uint8_t c, uint8_t d) {
		v.c[0] = a;
		v.c[1] = b;
		v.c[2] = c;
		v.c[3] = d;
	}
};


class VectorDouble4D {
private:

public:
	doublevec4 v{};

	inline double operator[](uint32_t position) {
		return v.c[position];
	}

	inline void operator+=(VectorDouble4D vec2) {
#if defined(USE_AVX)
		v.avx = _mm256_add_pd(v.avx, vec2.v.avx);
#elif defined(USE_SSE) // SSE2
		v.sse[0] = _mm_add_pd(v.sse[0], vec2.v.sse[0]);
		v.sse[1] = _mm_add_pd(v.sse[1], vec2.v.sse[1]);
#elif defined(USE_NEON)
		v.neon[0] = vaddq_f64(v.neon[0], vec2.v.neon[0]);
		v.neon[1] = vaddq_f64(v.neon[1], vec2.v.neon[1]);
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
#elif defined(USE_SSE2)
		VectorDouble4D ret;
		ret.v.sse[0] = _mm_add_pd(v.sse[0], vec2.v.sse[0]);
		ret.v.sse[1] = _mm_add_pd(v.sse[1], vec2.v.sse[1]);
		return ret;
#elif defined(USE_NEON)
		VectorDouble4D ret;
		ret.v.neon[0] = vaddq_f64(v.neon[0], vec2.v.neon[0]);
		ret.v.neon[1] = vaddq_f64(v.neon[1], vec2.v.neon[1]);
		return ret;
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
#elif defined(USE_SSE2)
		VectorDouble4D ret(a);
		ret.v.sse[0] = _mm_add_pd(v.sse[0], ret.v.sse[0]);
		ret.v.sse[1] = _mm_add_pd(v.sse[1], ret.v.sse[1]);
		return ret;
#elif defined(USE_NEON)
		VectorDouble4D ret(a);
		ret.v.neon[0] = vaddq_f64(v.neon[0], ret.v.neon[0]);
		ret.v.neon[1] = vaddq_f64(v.neon[1], ret.v.neon[1]);
		return ret;
#else
		VectorDouble4D ret(v.c[0] + a, v.c[1] + a, v.c[2] + a, v.c[3] + a);
		return ret;
#endif


	}

	inline VectorDouble4D *add(VectorDouble4D a) {
#if defined(USE_AVX)
		v.avx = _mm256_add_pd(v.avx, a.v.avx);
#elif defined(USE_SSE) // SSE2
		v.sse[0] = _mm_add_pd(v.sse[0], a.v.sse[0]);
		v.sse[1] = _mm_add_pd(v.sse[1], a.v.sse[1]);
#elif defined(USE_NEON)
		VectorDouble4D ret(a);
		v.neon[0] = vaddq_f64(v.neon[0], ret.v.neon[0]);
		v.neon[1] = vaddq_f64(v.neon[1], ret.v.neon[1]);
#else
		v.c[0] += a.v.c[0];
		v.c[1] += a.v.c[1];
		v.c[2] += a.v.c[2];
		v.c[3] += a.v.c[3];
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
#elif defined(USE_NEON)
		v.neon[0] = vmulq_f64(v.neon[0], vec2.v.neon[0]);
		v.neon[1] = vmulq_f64(v.neon[1], vec2.v.neon[1]);
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
		ret.v.c[0] *= v.c[0] * a;
		ret.v.c[1] *= v.c[1] * a;
		ret.v.c[2] *= v.c[2] * a;
		ret.v.c[3] *= v.c[3] * a;
#endif
		return ret;
	}

	inline VectorDouble4D operator/(double a) {
		VectorDouble4D ret;
#if defined(USE_AVX)
		doublevec4 b = {a, a, a, a};
		ret.v.avx = _mm256_div_pd(v.avx, b.avx);
#else
		ret.v.c[0] /= v.c[0] * a;
		ret.v.c[1] /= v.c[1] * a;
		ret.v.c[2] /= v.c[2] * a;
		ret.v.c[3] /= v.c[3] * a;
#endif
		return ret;
	}


	inline VectorDouble4D *capBetween1_0() {
#if defined(USE_AVX)
		doublevec4 one = {1.0f, 1.0f, 1.0f, 1.0f};
		v.avx = _mm256_max_pd(v.avx, one.avx);
		v.avx = _mm256_min_pd(v.avx, _mm256_setzero_pd());
#else
		if (v.c[0] > 1) {
			v.c[0] = 1;
		} else if (v.c[0] < 0) {
			v.c[0] = 0;
		}
		if (v.c[1] > 1) {
			v.c[1] = 1;
		} else if (v.c[1] < 0) {
			v.c[1] = 0;
		}
		if (v.c[2] > 1) {
			v.c[2] = 1;
		} else if (v.c[2] < 0) {
			v.c[2] = 0;
		}
		if (v.c[3] > 1) {
			v.c[3] = 1;
		} else if (v.c[3] < 0) {
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
		if (v.c[0] > upperBoundary) {
			v.c[0] = upperBoundary;
		} else if (v.c[0] < lowerBoundary) {
			v.c[0] = lowerBoundary;
		}
		if (v.c[1] > upperBoundary) {
			v.c[1] = upperBoundary;
		} else if (v.c[1] < lowerBoundary) {
			v.c[1] = lowerBoundary;
		}
		if (v.c[2] > upperBoundary) {
			v.c[2] = upperBoundary;
		} else if (v.c[2] < lowerBoundary) {
			v.c[2] = lowerBoundary;
		}
		if (v.c[3] > upperBoundary) {
			v.c[3] = upperBoundary;
		} else if (v.c[3] < lowerBoundary) {
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
		v.c[0] = ((v.c[0] - lowerInput) * ((upperOutput - lowerOutput) / (upperInput - lowerInput))) + lowerOutput;
		v.c[1] = ((v.c[1] - lowerInput) * ((upperOutput - lowerOutput) / (upperInput - lowerInput))) + lowerOutput;
		v.c[2] = ((v.c[2] - lowerInput) * ((upperOutput - lowerOutput) / (upperInput - lowerInput))) + lowerOutput;
		v.c[3] = ((v.c[3] - lowerInput) * ((upperOutput - lowerOutput) / (upperInput - lowerInput))) + lowerOutput;
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
public:
	doublemat4x4 m;

	inline VectorDouble4D operator[](uint32_t column) {
		VectorDouble4D ret;
#if defined(USE_AVX)
		ret.v.avx = m.avx[column];
#else
		ret.v.c[0] = m.c[column * 4];
		ret.v.c[1] = m.c[column * 4 + 1];
		ret.v.c[2] = m.c[column * 4 + 2];
		ret.v.c[3] = m.c[column * 4 + 3];
#endif
		return ret;
	}

	inline MatrixDouble4X4 *identity() {
		m = (doublemat4x4) {1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
		return this;
	}

	inline MatrixDouble4X4 operator*(MatrixDouble4X4 b) {
		MatrixDouble4X4 ret;
#if defined(USE_AVX512F)
		__m512d O0 = (__m512d) {b.m.c[0],b.m.c[0],b.m.c[0],b.m.c[0],b.m.c[4],b.m.c[4],b.m.c[4],b.m.c[4]};
		__m512d O1 = (__m512d) {b.m.c[1],b.m.c[1],b.m.c[1],b.m.c[1],b.m.c[5],b.m.c[5],b.m.c[5],b.m.c[5]};
		__m512d O2 = (__m512d) {b.m.c[2],b.m.c[2],b.m.c[2],b.m.c[2],b.m.c[6],b.m.c[6],b.m.c[6],b.m.c[6]};
		__m512d O3 = (__m512d) {b.m.c[3],b.m.c[3],b.m.c[3],b.m.c[3],b.m.c[7],b.m.c[7],b.m.c[7],b.m.c[7]};

		__m512d T0 =_mm512_insertf64x4(_mm512_castpd256_pd512(m.avx[0]),m.avx[0],1);
		__m512d T1 =_mm512_insertf64x4(_mm512_castpd256_pd512(m.avx[1]),m.avx[1],1);
		__m512d T2 =_mm512_insertf64x4(_mm512_castpd256_pd512(m.avx[2]),m.avx[2],1);
		__m512d T3 =_mm512_insertf64x4(_mm512_castpd256_pd512(m.avx[3]),m.avx[3],1);
		ret.m.avx512[0] = _mm512_mul_pd(T0, O0);
		ret.m.avx512[0] = _mm512_fmadd_pd(T1,O1,ret.m.avx512[0]);
		ret.m.avx512[0] = _mm512_fmadd_pd(T2,O2,ret.m.avx512[0]);
		ret.m.avx512[0] = _mm512_fmadd_pd(T3,O3,ret.m.avx512[0]);

		__m512d O4 = (__m512d) {b.m.c[8],b.m.c[8],b.m.c[8],b.m.c[8],b.m.c[12],b.m.c[12],b.m.c[12],b.m.c[12]};
		__m512d O5 = (__m512d) {b.m.c[9],b.m.c[9],b.m.c[9],b.m.c[9],b.m.c[13],b.m.c[13],b.m.c[13],b.m.c[13]};
		__m512d O6 = (__m512d) {b.m.c[10],b.m.c[10],b.m.c[10],b.m.c[10],b.m.c[14],b.m.c[14],b.m.c[14],b.m.c[14]};
		__m512d O7 = (__m512d) {b.m.c[11],b.m.c[11],b.m.c[11],b.m.c[11],b.m.c[15],b.m.c[15],b.m.c[15],b.m.c[15]};
		ret.m.avx512[1] = _mm512_mul_pd(T0, O4);
		ret.m.avx512[1] = _mm512_fmadd_pd(T1,O5,ret.m.avx512[1]);
		ret.m.avx512[1] = _mm512_fmadd_pd(T2,O6,ret.m.avx512[1]);
		ret.m.avx512[1] = _mm512_fmadd_pd(T3,O7,ret.m.avx512[1]);
#elif defined(USE_FMA)
		/*
		 * m0 * bcst 0
		 * m0 * bcst 4
		 * m0 * bcst 8
		 * m0 * bcst 12
		 */
		__m256d O0 = _mm256_broadcastsd_pd((__m128d) {b.m.c[0], 0.0f});
		__m256d O1 = _mm256_broadcastsd_pd((__m128d) {b.m.c[1], 0.0f});
		__m256d O2 = _mm256_broadcastsd_pd((__m128d) {b.m.c[2], 0.0f});
		__m256d O3 = _mm256_broadcastsd_pd((__m128d) {b.m.c[3], 0.0f});

		ret.m.avx[0] = _mm256_mul_pd(m.avx[0], O0);
		ret.m.avx[0] = _mm256_fmadd_pd(m.avx[1], O1, ret.m.avx[0]);
		ret.m.avx[0] = _mm256_fmadd_pd(m.avx[2], O2, ret.m.avx[0]);
		ret.m.avx[0] = _mm256_fmadd_pd(m.avx[3], O3, ret.m.avx[0]);

		__m256d O4 = _mm256_broadcastsd_pd((__m128d) {b.m.c[4], 0.0f});
		__m256d O5 = _mm256_broadcastsd_pd((__m128d) {b.m.c[5], 0.0f});
		__m256d O6 = _mm256_broadcastsd_pd((__m128d) {b.m.c[6], 0.0f});
		__m256d O7 = _mm256_broadcastsd_pd((__m128d) {b.m.c[7], 0.0f});

		ret.m.avx[1] = _mm256_mul_pd(m.avx[0], O4);
		ret.m.avx[1] = _mm256_fmadd_pd(m.avx[1], O5, ret.m.avx[1]);
		ret.m.avx[1] = _mm256_fmadd_pd(m.avx[2], O6, ret.m.avx[1]);
		ret.m.avx[1] = _mm256_fmadd_pd(m.avx[3], O7, ret.m.avx[1]);

		__m256d O8 = _mm256_broadcastsd_pd((__m128d) {b.m.c[8], 0.0f});
		__m256d O9 = _mm256_broadcastsd_pd((__m128d) {b.m.c[9], 0.0f});
		__m256d O10 = _mm256_broadcastsd_pd((__m128d) {b.m.c[10], 0.0f});
		__m256d O11 = _mm256_broadcastsd_pd((__m128d) {b.m.c[11], 0.0f});

		ret.m.avx[2] = _mm256_mul_pd(m.avx[0], O8);
		ret.m.avx[2] = _mm256_fmadd_pd(m.avx[1], O9, ret.m.avx[2]);
		ret.m.avx[2] = _mm256_fmadd_pd(m.avx[2], O10, ret.m.avx[2]);
		ret.m.avx[2] = _mm256_fmadd_pd(m.avx[3], O11, ret.m.avx[2]);

		__m256d O12 = _mm256_broadcastsd_pd((__m128d) {b.m.c[12], 0.0f});
		__m256d O13 = _mm256_broadcastsd_pd((__m128d) {b.m.c[13], 0.0f});
		__m256d O14 = _mm256_broadcastsd_pd((__m128d) {b.m.c[14], 0.0f});
		__m256d O15 = _mm256_broadcastsd_pd((__m128d) {b.m.c[15], 0.0f});

		ret.m.avx[3] = _mm256_mul_pd(m.avx[0], O12);
		ret.m.avx[3] = _mm256_fmadd_pd(m.avx[1], O13, ret.m.avx[3]);
		ret.m.avx[3] = _mm256_fmadd_pd(m.avx[2], O14, ret.m.avx[3]);
		ret.m.avx[3] = _mm256_fmadd_pd(m.avx[3], O15, ret.m.avx[3]);

#elif defined(USE_SSE2)

		ret.m.sse[0] = _mm_mul_pd(m.sse[0], (__m128d) {b.m.c[0], b.m.c[0]});
		__m128d cache = _mm_mul_pd(m.sse[2], (__m128d) {b.m.c[1], b.m.c[1]});
		ret.m.sse[0] = _mm_add_pd(cache, ret.m.sse[0]);
		cache = _mm_mul_pd(m.sse[4], (__m128d) {b.m.c[2], b.m.c[2]});
		ret.m.sse[0] = _mm_add_pd(cache, ret.m.sse[0]);
		cache = _mm_mul_pd(m.sse[6], (__m128d) {b.m.c[3], b.m.c[3]});
		ret.m.sse[0] = _mm_add_pd(cache, ret.m.sse[0]);
		//
		ret.m.sse[1] = _mm_mul_pd(m.sse[1], (__m128d) {b.m.c[0], b.m.c[0]});
		cache = _mm_mul_pd(m.sse[3], (__m128d) {b.m.c[1], b.m.c[1]});
		ret.m.sse[1] = _mm_add_pd(cache, ret.m.sse[1]);
		cache = _mm_mul_pd(m.sse[5], (__m128d) {b.m.c[2], b.m.c[2]});
		ret.m.sse[1] = _mm_add_pd(cache, ret.m.sse[1]);
		cache = _mm_mul_pd(m.sse[7], (__m128d) {b.m.c[3], b.m.c[3]});
		ret.m.sse[1] = _mm_add_pd(cache, ret.m.sse[1]);
		//

		ret.m.sse[2] = _mm_mul_pd(m.sse[0], (__m128d) {b.m.c[4], b.m.c[4]});
		cache = _mm_mul_pd(m.sse[2], (__m128d) {b.m.c[5], b.m.c[5]});
		ret.m.sse[2] = _mm_add_pd(cache, ret.m.sse[2]);
		cache = _mm_mul_pd(m.sse[4], (__m128d) {b.m.c[6], b.m.c[6]});
		ret.m.sse[2] = _mm_add_pd(cache, ret.m.sse[2]);
		cache = _mm_mul_pd(m.sse[6], (__m128d) {b.m.c[7], b.m.c[7]});
		ret.m.sse[2] = _mm_add_pd(cache, ret.m.sse[2]);
		//
		ret.m.sse[3] = _mm_mul_pd(m.sse[1], (__m128d) {b.m.c[4], b.m.c[4]});
		cache = _mm_mul_pd(m.sse[3], (__m128d) {b.m.c[5], b.m.c[5]});
		ret.m.sse[3] = _mm_add_pd(cache, ret.m.sse[3]);
		cache = _mm_mul_pd(m.sse[5], (__m128d) {b.m.c[6], b.m.c[6]});
		ret.m.sse[3] = _mm_add_pd(cache, ret.m.sse[3]);
		cache = _mm_mul_pd(m.sse[7], (__m128d) {b.m.c[7], b.m.c[7]});
		ret.m.sse[3] = _mm_add_pd(cache, ret.m.sse[3]);
		//

		ret.m.sse[4] = _mm_mul_pd(m.sse[0], (__m128d) {b.m.c[8], b.m.c[8]});
		cache = _mm_mul_pd(m.sse[2], (__m128d) {b.m.c[9], b.m.c[9]});
		ret.m.sse[4] = _mm_add_pd(cache, ret.m.sse[4]);
		cache = _mm_mul_pd(m.sse[4], (__m128d) {b.m.c[10], b.m.c[10]});
		ret.m.sse[4] = _mm_add_pd(cache, ret.m.sse[4]);
		cache = _mm_mul_pd(m.sse[6], (__m128d) {b.m.c[11], b.m.c[11]});
		ret.m.sse[4] = _mm_add_pd(cache, ret.m.sse[4]);
		//
		ret.m.sse[5] = _mm_mul_pd(m.sse[1], (__m128d) {b.m.c[8], b.m.c[8]});
		cache = _mm_mul_pd(m.sse[3], (__m128d) {b.m.c[9], b.m.c[9]});
		ret.m.sse[5] = _mm_add_pd(cache, ret.m.sse[5]);
		cache = _mm_mul_pd(m.sse[5], (__m128d) {b.m.c[10], b.m.c[10]});
		ret.m.sse[5] = _mm_add_pd(cache, ret.m.sse[5]);
		cache = _mm_mul_pd(m.sse[7], (__m128d) {b.m.c[11], b.m.c[11]});
		ret.m.sse[5] = _mm_add_pd(cache, ret.m.sse[5]);
		//

		ret.m.sse[6] = _mm_mul_pd(m.sse[0], (__m128d) {b.m.c[12], b.m.c[12]});
		cache = _mm_mul_pd(m.sse[2], (__m128d) {b.m.c[13], b.m.c[13]});
		ret.m.sse[6] = _mm_add_pd(cache, ret.m.sse[6]);
		cache = _mm_mul_pd(m.sse[4], (__m128d) {b.m.c[14], b.m.c[14]});
		ret.m.sse[6] = _mm_add_pd(cache, ret.m.sse[6]);
		cache = _mm_mul_pd(m.sse[6], (__m128d) {b.m.c[15], b.m.c[15]});
		ret.m.sse[6] = _mm_add_pd(cache, ret.m.sse[6]);
		//
		ret.m.sse[7] = _mm_mul_pd(m.sse[1], (__m128d) {b.m.c[12], b.m.c[12]});
		cache = _mm_mul_pd(m.sse[3], (__m128d) {b.m.c[13], b.m.c[13]});
		ret.m.sse[7] = _mm_add_pd(cache, ret.m.sse[7]);
		cache = _mm_mul_pd(m.sse[5], (__m128d) {b.m.c[14], b.m.c[14]});
		ret.m.sse[7] = _mm_add_pd(cache, ret.m.sse[7]);
		cache = _mm_mul_pd(m.sse[7], (__m128d) {b.m.c[15], b.m.c[15]});
		ret.m.sse[7] = _mm_add_pd(cache, ret.m.sse[7]);
#elif defined(USE_NEON)
		ret.m.neon[0] = vmulq_f64(m.neon[0], (float64x2_t) {b.m.c[0], b.m.c[0]});
		ret.m.neon[0] = vfmaq_f64(ret.m.neon[0], m.neon[2], (float64x2_t) {b.m.c[1], b.m.c[1]});
		ret.m.neon[0] = vfmaq_f64(ret.m.neon[0], m.neon[4], (float64x2_t) {b.m.c[2], b.m.c[2]});
		ret.m.neon[0] = vfmaq_f64(ret.m.neon[0], m.neon[6], (float64x2_t) {b.m.c[3], b.m.c[3]});
		//
		ret.m.neon[1] = vmulq_f64(m.neon[1], (float64x2_t) {b.m.c[0], b.m.c[0]});
		ret.m.neon[1] = vfmaq_f64(ret.m.neon[1], m.neon[3], (float64x2_t) {b.m.c[1], b.m.c[1]});
		ret.m.neon[1] = vfmaq_f64(ret.m.neon[1], m.neon[5], (float64x2_t) {b.m.c[2], b.m.c[2]});
		ret.m.neon[1] = vfmaq_f64(ret.m.neon[1], m.neon[7], (float64x2_t) {b.m.c[3], b.m.c[3]});
		//
		ret.m.neon[2] = vmulq_f64(m.neon[0], (float64x2_t) {b.m.c[4], b.m.c[4]});
		ret.m.neon[2] = vfmaq_f64(ret.m.neon[2], m.neon[2], (float64x2_t) {b.m.c[5], b.m.c[5]});
		ret.m.neon[2] = vfmaq_f64(ret.m.neon[2], m.neon[4], (float64x2_t) {b.m.c[6], b.m.c[6]});
		ret.m.neon[2] = vfmaq_f64(ret.m.neon[2], m.neon[6], (float64x2_t) {b.m.c[7], b.m.c[7]});

		//
		ret.m.neon[3] = vmulq_f64(m.neon[1], (float64x2_t) {b.m.c[4], b.m.c[4]});
		ret.m.neon[3] = vfmaq_f64(ret.m.neon[3], m.neon[3], (float64x2_t) {b.m.c[5], b.m.c[5]});
		ret.m.neon[3] = vfmaq_f64(ret.m.neon[3], m.neon[5], (float64x2_t) {b.m.c[6], b.m.c[6]});
		ret.m.neon[3] = vfmaq_f64(ret.m.neon[3], m.neon[7], (float64x2_t) {b.m.c[7], b.m.c[7]});

		//
		ret.m.neon[4] = vmulq_f64(m.neon[0], (float64x2_t) {b.m.c[8], b.m.c[8]});
		ret.m.neon[4] = vfmaq_f64(ret.m.neon[4], m.neon[2], (float64x2_t) {b.m.c[9], b.m.c[9]});
		ret.m.neon[4] = vfmaq_f64(ret.m.neon[4], m.neon[4], (float64x2_t) {b.m.c[10], b.m.c[10]});
		ret.m.neon[4] = vfmaq_f64(ret.m.neon[4], m.neon[6], (float64x2_t) {b.m.c[11], b.m.c[11]});

		//
		ret.m.neon[5] = vmulq_f64(m.neon[1], (float64x2_t) {b.m.c[8], b.m.c[8]});
		ret.m.neon[5] = vfmaq_f64(ret.m.neon[5], m.neon[3], (float64x2_t) {b.m.c[9], b.m.c[9]});
		ret.m.neon[5] = vfmaq_f64(ret.m.neon[5], m.neon[5], (float64x2_t) {b.m.c[10], b.m.c[10]});
		ret.m.neon[5] = vfmaq_f64(ret.m.neon[5], m.neon[7], (float64x2_t) {b.m.c[11], b.m.c[11]});

		//

		ret.m.neon[6] = vmulq_f64(m.neon[0], (float64x2_t) {b.m.c[12], b.m.c[12]});
		ret.m.neon[6] = vfmaq_f64(ret.m.neon[6], m.neon[2], (float64x2_t) {b.m.c[13], b.m.c[13]});
		ret.m.neon[6] = vfmaq_f64(ret.m.neon[6], m.neon[4], (float64x2_t) {b.m.c[14], b.m.c[14]});
		ret.m.neon[6] = vfmaq_f64(ret.m.neon[6], m.neon[6], (float64x2_t) {b.m.c[15], b.m.c[15]});

		//
		ret.m.neon[7] = vmulq_f64(m.neon[1], (float64x2_t) {b.m.c[12], b.m.c[12]});
		ret.m.neon[7] = vfmaq_f64(ret.m.neon[7], m.neon[3], (float64x2_t) {b.m.c[13], b.m.c[13]});
		ret.m.neon[7] = vfmaq_f64(ret.m.neon[7], m.neon[5], (float64x2_t) {b.m.c[14], b.m.c[14]});
		ret.m.neon[7] = vfmaq_f64(ret.m.neon[7], m.neon[7], (float64x2_t) {b.m.c[15], b.m.c[15]});


#else
		ret.m.c[0] = m.c[0] * b.m.c[0] + m.c[4] * b.m.c[1] + m.c[8] * b.m.c[2] +
					 m.c[12] * b.m.c[3];// c11 = a11 * b11 + a12 * b21 + ...
		ret.m.c[1] = m.c[1] * b.m.c[0] + m.c[5] * b.m.c[1] + m.c[9] * b.m.c[2] +
					 m.c[13] * b.m.c[3];// c12 = a21 + b11 + a22 + b21
		ret.m.c[2] = m.c[2] * b.m.c[0] + m.c[6] * b.m.c[1] + m.c[10] * b.m.c[2] + m.c[14] * b.m.c[3];
		ret.m.c[3] = m.c[3] * b.m.c[0] + m.c[7] * b.m.c[1] + m.c[11] * b.m.c[2] + m.c[15] * b.m.c[3];
		ret.m.c[4] = m.c[0] * b.m.c[4] + m.c[4] * b.m.c[5] + m.c[8] * b.m.c[6] +
					 m.c[12] * b.m.c[7];// c21 = a11 * b12 + b12 * b22 + ...
		ret.m.c[5] = m.c[1] * b.m.c[4] + m.c[5] * b.m.c[5] + m.c[9] * b.m.c[6] + m.c[13] * b.m.c[7];
		ret.m.c[6] = m.c[2] * b.m.c[4] + m.c[6] * b.m.c[5] + m.c[10] * b.m.c[6] + m.c[14] * b.m.c[7];
		ret.m.c[7] = m.c[3] * b.m.c[4] + m.c[7] * b.m.c[5] + m.c[11] * b.m.c[6] + m.c[15] * b.m.c[7];
		ret.m.c[8]  = m.c[0] * b.m.c[8]  + m.c[4] * b.m.c[9]  + m.c[8]  * b.m.c[10] + m.c[12] * b.m.c[11];
		ret.m.c[9]  = m.c[1] * b.m.c[8]  + m.c[5] * b.m.c[9]  + m.c[9]  * b.m.c[10] + m.c[13] * b.m.c[11];
		ret.m.c[10] = m.c[2] * b.m.c[8]  + m.c[6] * b.m.c[9]  + m.c[10] * b.m.c[10] + m.c[14] * b.m.c[11];
		ret.m.c[11] = m.c[3] * b.m.c[8]  + m.c[7] * b.m.c[9]  + m.c[11] * b.m.c[10] + m.c[15] * b.m.c[11];
		ret.m.c[12] = m.c[0] * b.m.c[12] + m.c[4] * b.m.c[13] + m.c[8]  * b.m.c[14] + m.c[12] * b.m.c[15];
		ret.m.c[13] = m.c[1] * b.m.c[12] + m.c[5] * b.m.c[13] + m.c[9]  * b.m.c[14] + m.c[13] * b.m.c[15];
		ret.m.c[14] = m.c[2] * b.m.c[12] + m.c[6] * b.m.c[13] + m.c[10] * b.m.c[14] + m.c[14] * b.m.c[15];
		ret.m.c[15] = m.c[3] * b.m.c[12] + m.c[7] * b.m.c[13] + m.c[11] * b.m.c[14] + m.c[15] * b.m.c[15];
#endif
		return ret;
	}

	inline VectorDouble4D operator*(VectorDouble4D b) {
		VectorDouble4D ret;
		ret.v.c[0] = m.c[0] * b.v.c[0] + m.c[4] * b.v.c[1] + m.c[8] * b.v.c[2] + m.c[12] * b.v.c[3];
		ret.v.c[1] = m.c[1] * b.v.c[0] + m.c[5] * b.v.c[1] + m.c[9] * b.v.c[2] + m.c[13] * b.v.c[3];
		ret.v.c[2] = m.c[2] * b.v.c[0] + m.c[6] * b.v.c[1] + m.c[10] * b.v.c[2] + m.c[14] * b.v.c[3];
		ret.v.c[3] = m.c[3] * b.v.c[0] + m.c[7] * b.v.c[1] + m.c[11] * b.v.c[2] + m.c[15] * b.v.c[3];
		return ret;
	}

	inline MatrixDouble4X4() {
		m = (doublemat4x4) {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	}

	inline MatrixDouble4X4(VectorDouble4D a, VectorDouble4D b, VectorDouble4D c, VectorDouble4D d) {
		m.c[0] = a.v.c[0];
		m.c[1] = a.v.c[1];
		m.c[2] = a.v.c[2];
		m.c[3] = a.v.c[3];
		m.c[4] = b.v.c[0];
		m.c[5] = b.v.c[1];
		m.c[6] = b.v.c[2];
		m.c[7] = b.v.c[3];
		m.c[8] = c.v.c[0];
		m.c[9] = c.v.c[1];
		m.c[10] = c.v.c[2];
		m.c[11] = c.v.c[3];
		m.c[12] = d.v.c[0];
		m.c[13] = d.v.c[1];
		m.c[14] = d.v.c[2];
		m.c[15] = d.v.c[3];
	}
};

class VectorDouble8D {
public:
	doublevec8 v;

	inline double operator[](uint32_t position) {
		return v.c[position];
	}

	inline void operator+=(VectorDouble8D vec2) {
#if defined(USE_AVX512F) || defined(KNCNI)
		v.avx512 = _mm512_add_pd(v.avx512, vec2.v.avx512);
#elif defined(USE_AVX)
		v.avx[0] = _mm256_add_pd(v.avx[0], vec2.v.avx[0]);
		v.avx[1] = _mm256_add_pd(v.avx[1], vec2.v.avx[1]);
#elif defined(USE_SSE) // SSE2
		v.sse[0] = _mm_add_pd(v.sse[0], vec2.v.sse[0]);
		v.sse[1] = _mm_add_pd(v.sse[1], vec2.v.sse[1]);
		v.sse[2] = _mm_add_pd(v.sse[2], vec2.v.sse[2]);
		v.sse[3] = _mm_add_pd(v.sse[3], vec2.v.sse[3]);
#else
		v.c[0] += vec2[0];
		v.c[1] += vec2[1];
		v.c[2] += vec2[2];
		v.c[3] += vec2[3];
		v.c[4] += vec2[4];
		v.c[5] += vec2[5];
		v.c[6] += vec2[6];
		v.c[7] += vec2[7];
#endif


	}

	inline VectorDouble8D(double a, double b, double c, double d, double e, double f, double g, double h) {
		v.c[0] = a;
		v.c[1] = b;
		v.c[2] = c;
		v.c[3] = d;
		v.c[4] = e;
		v.c[5] = f;
		v.c[6] = g;
		v.c[7] = h;
	}

	inline VectorDouble8D(VectorDouble4D a, VectorDouble4D b) {
#if defined(USE_AVX)
		v.avx[0] = a.v.avx;
		v.avx[1] = b.v.avx;
#else
		v.c[0] = a.v.c[0];
		v.c[1] = a.v.c[1];
		v.c[2] = a.v.c[2];
		v.c[3] = a.v.c[3];
		v.c[4] = b.v.c[0];
		v.c[5] = b.v.c[1];
		v.c[6] = b.v.c[2];
		v.c[7] = b.v.c[3];
#endif
	}

	inline VectorDouble8D() {
		v.c[0] = 0.0f;
		v.c[1] = 0.0f;
		v.c[2] = 0.0f;
		v.c[3] = 0.0f;
		v.c[4] = 0.0f;
		v.c[5] = 0.0f;
		v.c[6] = 0.0f;
		v.c[7] = 0.0f;
	}
};


class Complex64 {
public:
	doublevec2 c;


	inline Complex64(double real, double img) {
		c.c[0] = real;
		c.c[1] = img;
	}

	inline Complex64() {
		c.c[0] = 0;
		c.c[1] = 0;
	}

	inline Complex64 *add(Complex64 a) {
		c.c[0] += a.c.c[0];
		c.c[1] += a.c.c[1];
		return this;
	}

	inline Complex64 *subtract(Complex64 a) {
		c.c[0] -= a.c.c[0];
		c.c[1] -= a.c.c[1];
		return this;
	}

	inline Complex64 *conjugate() {
		c.c[1] = -c.c[1];
		return this;
	}


	inline Complex64 *multiply(Complex64 a) {
		double d1 = c.c[0] * a.c.c[0] - c.c[1] * a.c.c[1];
		double d2 = c.c[0] * a.c.c[1] + c.c[1] * a.c.c[0];
		c.c[0] = d1;
		c.c[1] = d2;
		return this;
	}

	inline Complex64 *divide(Complex64 a) {
		double d1 = (c.c[0] * a.c.c[0] + c.c[1] * a.c.c[1]) / (a.c.c[0] * a.c.c[0] + a.c.c[1] * a.c.c[1]);
		double d2 = (c.c[1] * a.c.c[0] - c.c[0] * a.c.c[1]) / (a.c.c[0] * a.c.c[0] + a.c.c[1] * a.c.c[1]);

		c.c[0] = d1;
		c.c[1] = d2;

		return this;
	}

	inline Complex64 *sqrt() {
		double d2 = ::sqrt((-c.c[0] + ::sqrt(c.c[0] * c.c[0] + c.c[1] * c.c[1])) / (2));
		double d1;
		if (d2 == 0) {
			d1 = ::sqrt(c.c[0]);
		} else {
			d1 = c.c[1] / (2 * d2);
		}
		c.c[0] = d1;
		c.c[1] = d2;
		return this;
	}

	inline Complex64 *sin() {
		double d1 = ::sin(c.c[0]) * ::cosh(c.c[1]);
		double d2 = ::cos(c.c[1]) * ::sinh(c.c[0]);

		c.c[0] = d1;
		c.c[1] = d2;
		return this;
	}


	inline Complex64 *exp() {
		double d1 = ::exp(c.c[0]) * ::cos(c.c[1]);
		double d2 = ::exp(c.c[0]) * ::sin(c.c[1]);


		c.c[0] = d1;
		c.c[1] = d2;
		return this;
	}

	inline Complex64 *exp(double n) {
		double d1 = ::pow(n, c.c[0]) * ::cos(c.c[1] * ::log(n));
		double d2 = ::pow(n, c.c[0]) * ::sin(c.c[1] * ::log(n));


		c.c[0] = d1;
		c.c[1] = d2;
		return this;
	}

	inline double abs() {
		return ::sqrt(c.c[0] * c.c[0] + c.c[1] * c.c[1]);
	}

	inline bool abs_gt(double a) {
		return a * a < c.c[0] * c.c[0] + c.c[1] * c.c[1];
	}

	inline bool abs_lt(double a) {
		return a * a > c.c[0] * c.c[0] + c.c[1] * c.c[1];
	}

	inline bool abs_eq(double a) {
		return a * a == c.c[0] * c.c[0] + c.c[1] * c.c[1];
	}

	inline double imaginary() {
		return c.c[1];
	}

	inline double real() {
		return c.c[0];
	}


};

#endif //MATH_LIB_A_MATH_LIB_H
