/*
 * simd.hpp
 *
 *  Created on: Jul 3, 2020
 *      Author: dmarce1
 */

#ifndef TIGERGRAV_SIMD_HPP_
#define TIGERGRAV_SIMD_HPP_

#include <tigersph/defs.hpp>

#include <immintrin.h>

#include <cmath>

#ifdef USE_AVX512
#define SIMDALIGN                  __attribute__((aligned(64)))
#define SIMD_VLEN 16
#define _simd_float                 __m512
#define _simd_int                   __m512i
#define _mmx_add_ps(a,b)            _mm512_add_ps((a),(b))
#define _mmx_sub_ps(a,b)            _mm512_sub_ps((a),(b))
#define _mmx_mul_ps(a,b)            _mm512_mul_ps((a),(b))
#define _mmx_div_ps(a,b)            _mm512_div_ps((a),(b))
#define _mmx_sqrt_ps(a)             _mm512_sqrt_ps(a)
#define _mmx_min_ps(a, b)           _mm512_min_ps((a),(b))
#define _mmx_max_ps(a, b)           _mm512_max_ps((a),(b))
#define _mmx_or_ps(a, b)            _mm512_or_ps((a),(b))
#define _mmx_and_ps(a, b)           _mm512_and_ps((a),(b))
#define _mmx_andnot_ps(a, b)        _mm512_andnot_ps((a),(b))
#define _mmx_rsqrt_ps(a)            _mm512_rsqrt14_ps(a)
#define _mmx_add_epi32(a,b)         _mm512_add_epi32((a),(b))
#define _mmx_sub_epi32(a,b)         _mm512_sub_epi32((a),(b))
#define _mmx_mul_epi32(a,b)         _mm512_mullo_epi32((a),(b))
#define _mmx_cvtps_epi32(a)         _mm512_cvtps_epi32((a))
#define _mmx_fmadd_ps(a,b,c)        _mm512_fmadd_ps ((a),(b),(c))
#define _mmx_cmp_ps(a,b,c)       _mm512_cmp_ps_mask(a,b,c)
#endif

#ifdef USE_AVX2
#define SIMDALIGN                  __attribute__((aligned(32)))
#define SIMD_VLEN 8
#define _simd_float                 __m256
#define _simd_int                   __m256i
#define _mmx_add_ps(a,b)            _mm256_add_ps((a),(b))
#define _mmx_sub_ps(a,b)            _mm256_sub_ps((a),(b))
#define _mmx_mul_ps(a,b)            _mm256_mul_ps((a),(b))
#define _mmx_div_ps(a,b)            _mm256_div_ps((a),(b))
#define _mmx_sqrt_ps(a)             _mm256_sqrt_ps(a)
#define _mmx_min_ps(a, b)           _mm256_min_ps((a),(b))
#define _mmx_max_ps(a, b)           _mm256_max_ps((a),(b))
#define _mmx_or_ps(a, b)            _mm256_or_ps((a),(b))
#define _mmx_and_ps(a, b)           _mm256_and_ps((a),(b))
#define _mmx_andnot_ps(a, b)        _mm256_andnot_ps((a),(b))
#define _mmx_rsqrt_ps(a)            _mm256_rsqrt_ps(a)
#define _mmx_add_epi32(a,b)         _mm256_add_epi32((a),(b))
#define _mmx_sub_epi32(a,b)         _mm256_sub_epi32((a),(b))
#define _mmx_mul_epi32(a,b)         _mm256_mullo_epi32((a),(b))
#define _mmx_cvtps_epi32(a)         _mm256_cvtps_epi32((a))
#define _mmx_fmadd_ps(a,b,c)        _mm256_fmadd_ps ((a),(b),(c))
#define _mmx_cmp_ps(a,b,c)        	_mm256_cmp_ps(a,b,c)
#endif

#ifdef USE_AVX
#define SIMDALIGN                  __attribute__((aligned(16)))
#define SIMD_VLEN 4
#define _simd_float                 __m128
#define _simd_int                   __m128i
#define _mmx_add_ps(a,b)            _mm_add_ps((a),(b))
#define _mmx_sub_ps(a,b)            _mm_sub_ps((a),(b))
#define _mmx_mul_ps(a,b)            _mm_mul_ps((a),(b))
#define _mmx_div_ps(a,b)            _mm_div_ps((a),(b))
#define _mmx_sqrt_ps(a)             _mm_sqrt_ps(a)
#define _mmx_min_ps(a, b)           _mm_min_ps((a),(b))
#define _mmx_max_ps(a, b)           _mm_max_ps((a),(b))
#define _mmx_or_ps(a, b)            _mm_or_ps((a),(b))
#define _mmx_and_ps(a, b)           _mm_and_ps((a),(b))
#define _mmx_andnot_ps(a, b)        _mm_andnot_ps((a),(b))
#define _mmx_rsqrt_ps(a)            _mm_rsqrt_ps(a)
#define _mmx_add_epi32(a,b)         _mm_add_epi32((a),(b))
#define _mmx_sub_epi32(a,b)         _mm_sub_epi32((a),(b))
#define _mmx_mul_epi32(a,b)         _mm_mullo_epi32((a),(b))
#define _mmx_cvtps_epi32(a)         _mm_cvtps_epi32((a))
#define _mmx_cmp_ps(a,b,c)        	_mm_cmp_ps(a,b,c)
#endif

class simd_int;
class simd_float;

class simd_float {
private:
	_simd_float v;
public:
	static constexpr std::size_t size() {
		return SIMD_VLEN;
	}
	simd_float() = default;
	inline ~simd_float() = default;
	simd_float(const simd_float&) = default;
	inline simd_float(float d) {
#ifdef USE_AVX512
        v = _mm512_set_ps(d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d);
#endif
#ifdef USE_AVX2
		v = _mm256_set_ps(d, d, d, d, d, d, d, d);
#endif
#ifdef USE_AVX
		v = _mm_set_ps(d, d, d, d);
#endif
	}
	union union_mm {
#ifdef USE_AVX512
		__m512 m16[1];
		__m256 m8[2];
		__m128 m4[4];
		float m1[16];
#endif
#ifdef USE_AVX2
		__m256 m8[1];
		__m128 m4[2];
		float m1[8];
#endif
#ifdef USE_AVX
		__m128 m4[2];
		float m1[4];
#endif
	};
	inline float sum() const {
		union_mm s;
#ifdef USE_AVX512
		s.m16[0] = v;
		s.m8[0] = _mm256_add_ps(s.m8[0],s.m8[1]);
		s.m4[0] = _mm_add_ps(s.m4[0], s.m4[1]);
#endif
#ifdef USE_AVX2
		s.m8[0] = v;
		s.m4[0] = _mm_add_ps(s.m4[0], s.m4[1]);
#endif
#ifdef USE_AVX
		s.m4[0] = v;
#endif
		s.m1[4] = s.m1[2];
		s.m1[5] = s.m1[3];
		s.m4[0] = _mm_add_ps(s.m4[0], s.m4[1]);
		return s.m1[0] + s.m1[1];
	}
	inline simd_float(const simd_int &other);
	inline simd_int to_int() const;
	inline simd_float& operator=(const simd_float &other) = default;
	simd_float& operator=(simd_float &&other) = default;
	inline simd_float operator+(const simd_float &other) const {
		simd_float r;
		r.v = _mmx_add_ps(v, other.v);
		return r;
	}
	inline simd_float operator-(const simd_float &other) const {
		simd_float r;
		r.v = _mmx_sub_ps(v, other.v);
		return r;
	}
	inline simd_float operator*(const simd_float &other) const {
		simd_float r;
		r.v = _mmx_mul_ps(v, other.v);
		return r;
	}
	inline simd_float operator/(const simd_float &other) const {
		simd_float r;
		r.v = _mmx_div_ps(v, other.v);
		return r;
	}
	inline simd_float operator+() const {
		return *this;
	}
	inline simd_float operator-() const {
		return simd_float(0.0) - *this;
	}
	inline simd_float& operator+=(const simd_float &other) {
		*this = *this + other;
		return *this;
	}
	inline simd_float& operator-=(const simd_float &other) {
		*this = *this - other;
		return *this;
	}
	inline simd_float& operator*=(const simd_float &other) {
		*this = *this * other;
		return *this;
	}
	inline simd_float& operator/=(const simd_float &other) {
		*this = *this / other;
		return *this;
	}

	inline simd_float operator*(float d) const {
		const simd_float other = d;
		return other * *this;
	}
	inline simd_float operator/(float d) const {
		const simd_float other = 1.0 / d;
		return *this * other;
	}

	inline simd_float operator*=(float d) {
		*this = *this * d;
		return *this;
	}
	inline simd_float operator/=(float d) {
		*this = *this * (1.0 / d);
		return *this;
	}
	inline float& operator[](std::size_t i) {
		float *a = reinterpret_cast<float*>(&v[i / SIMD_VLEN]);
		return a[i];
	}
	inline float operator[](std::size_t i) const {
		const float *a = reinterpret_cast<const float*>(&v[i / SIMD_VLEN]);
		return a[i];
	}
	friend simd_float sqrt(const simd_float&);
	friend simd_float rsqrt(const simd_float&);
	friend simd_float operator*(float, const simd_float &other);
	friend simd_float operator/(float, const simd_float &other);
	friend simd_float max(const simd_float &a, const simd_float &b);
	friend simd_float min(const simd_float &a, const simd_float &b);
	friend simd_float fmadd(const simd_float &a, const simd_float &b, const simd_float &c);
	friend simd_float abs(const simd_float &a);
	simd_float operator<(simd_float other) const {
		simd_float rc;
		static const simd_float one(1);
		static const simd_float zero(0);
		auto mask = _mmx_cmp_ps(v, other.v, _CMP_LT_OQ);
#ifdef USE_AVX512
			rc.v = _mm512_mask_mov_ps(zero.v,mask,one.v);
#else
		rc.v = _mmx_and_ps(mask, one.v);
#endif
		return rc;
	}
	simd_float operator>(simd_float other) const {
		simd_float rc;
		static const simd_float one(1);
		static const simd_float zero(0);
		auto mask = _mmx_cmp_ps(v, other.v, _CMP_GT_OQ);
#ifdef USE_AVX512
			rc.v = _mm512_mask_mov_ps(zero.v,mask,one.v);
#else
		rc.v = _mmx_and_ps(mask, one.v);
#endif
		return rc;
	}
}
SIMDALIGN;

inline simd_float fmadd(const simd_float &a, const simd_float &b, const simd_float &c) {

	simd_float v;
#ifdef USE_AVX
	v = a*b + c;
#else
	v.v = _mmx_fmadd_ps(a.v, b.v, c.v);
#endif
	return v;
}

inline simd_float sqrt(const simd_float &vec) {
	simd_float r;
	r.v = _mmx_sqrt_ps(vec.v);
	return r;

}

inline simd_float rsqrt(const simd_float &vec) {
	simd_float r;
	r.v = _mmx_rsqrt_ps(vec.v);
	return r;
}

inline simd_float operator*(float d, const simd_float &other) {
	const simd_float a = d;
	return a * other;
}

inline simd_float operator/(float d, const simd_float &other) {
	const simd_float a = d;
	return a / other;
}

inline simd_float max(const simd_float &a, const simd_float &b) {
	simd_float r;
	r.v = _mmx_max_ps(a.v, b.v);
	return r;
}

inline simd_float min(const simd_float &a, const simd_float &b) {
	simd_float r;
	r.v = _mmx_min_ps(a.v, b.v);
	return r;
}

inline simd_float abs(const simd_float &a) {
	return max(a, -a);
}

#endif /* TIGERGRAV_SIMD_HPP_ */
