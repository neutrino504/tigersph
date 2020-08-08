/*
 * general_vect.hpp
 *
 *  Created on: Nov 30, 2019
 *      Author: dmarce1
 */

#ifndef VECT_HPP_
#define VECT_HPP_

#include <tigersph/defs.hpp>
#include <tigersph/cuda_export.hpp>

#include <array>
#include <atomic>
#include <cmath>
#include <limits>

template<class T, int N>
class general_vect {
#ifdef __CUDA_ARCH__
	T v[N];
#else
	std::array<T, N> v;
#endif
public:
	CUDA_EXPORT general_vect() {
		for( int d = 0; d < NDIM; d++) {
			v[d] = std::numeric_limits<T>::signaling_NaN();
		}
	}
	CUDA_EXPORT general_vect(std::array<T, N> a);
	CUDA_EXPORT general_vect(T a);
	CUDA_EXPORT T& operator[](int i);
	CUDA_EXPORT T operator[](int i) const;
	CUDA_EXPORT general_vect operator-() const;
	CUDA_EXPORT general_vect operator-(const general_vect &other) const;
	CUDA_EXPORT general_vect operator+(const general_vect &other) const;
	CUDA_EXPORT general_vect operator*(T r) const;
	CUDA_EXPORT general_vect operator/(T r) const;
	CUDA_EXPORT general_vect operator-=(const general_vect &other) ;
	CUDA_EXPORT general_vect operator+=(const general_vect &other) ;
	CUDA_EXPORT general_vect operator*=(T r);
	CUDA_EXPORT general_vect operator/=(T r);
	CUDA_EXPORT bool operator<(const general_vect &other) const;
	CUDA_EXPORT bool operator<=(const general_vect &other) const;
	CUDA_EXPORT bool operator>(const general_vect &other) const;
	CUDA_EXPORT bool operator>=(const general_vect &other) const;
	CUDA_EXPORT bool operator==(const general_vect &other) const;
	CUDA_EXPORT bool operator!=(const general_vect &other) const;
	CUDA_EXPORT T dot(const general_vect &other) const;
	template<class U>
	CUDA_EXPORT operator general_vect<U,N>() const;
	template<class Arc>
	void serialize(Arc &&a, unsigned) {
		a & v;
	}

};

template<class T, int N>
template<class U>
CUDA_EXPORT inline general_vect<T,N>::operator general_vect<U,N>() const {
	general_vect<U,N> a;
	for( int i = 0; i < N; i++) {
		a[i] = (U) v[i];
	}
	return a;
}

template<class T, int N>
CUDA_EXPORT bool inline general_vect<T, N>::operator<(const general_vect &other) const {
	for (int n = 0; n < N; n++) {
		if ((*this)[n] < other[n]) {
			return true;
		} else if ((*this)[n] > other[n]) {
			return false;
		}
	}
	return false;
}

template<class T, int N>
CUDA_EXPORT bool inline general_vect<T, N>::operator<=(const general_vect &other) const {
	return *this < other || *this == other;
}

template<class T, int N>
CUDA_EXPORT bool inline general_vect<T, N>::operator>(const general_vect &other) const {
	return !(*this <= other);
}

template<class T, int N>
CUDA_EXPORT bool inline general_vect<T, N>::operator>=(const general_vect &other) const {
	return !(*this < other);
}

template<class T, int N>
CUDA_EXPORT inline general_vect<T, N>::general_vect(std::array<T, N> a) :
		v(a) {
}

template<class T, int N>
CUDA_EXPORT inline general_vect<T, N>::general_vect(T a) {

	for (int i = 0; i < N; i++) {
		v[i] = a;
	}
}

template<class T, int N>
CUDA_EXPORT inline bool general_vect<T, N>::operator==(const general_vect<T, N> &other) const {

	for (int dim = 0; dim < NDIM; dim++) {
		if ((*this)[dim] != other[dim]) {
			return false;
		}
	}
	return true;
}

template<class T, int N>
CUDA_EXPORT inline bool general_vect<T, N>::operator!=(const general_vect<T, N> &other) const {
	return !((*this) == other);
}

template<class T, int N>
CUDA_EXPORT inline T& general_vect<T, N>::operator[](int i) {
	return v[i];
}

template<class T, int N>
CUDA_EXPORT inline T general_vect<T, N>::operator[](int i) const {
	return v[i];
}

template<class T, int N>
CUDA_EXPORT inline general_vect<T, N> general_vect<T, N>::operator-() const {
	general_vect<T, N> result;

	for (int dim = 0; dim < N; dim++) {
		result[dim] = -v[dim];
	}
	return result;
}

template<class T, int N>
CUDA_EXPORT inline general_vect<T, N> general_vect<T, N>::operator-(const general_vect<T, N> &other) const {
	general_vect<T, N> result;

	for (int dim = 0; dim < N; dim++) {
		result[dim] = v[dim] - other[dim];
	}
	return result;
}

template<class T, int N>
CUDA_EXPORT inline general_vect<T, N> general_vect<T, N>::operator-=(const general_vect<T, N> &other) {

	for (int dim = 0; dim < N; dim++) {
		v[dim] -= other[dim];
	}
	return *this;
}

template<class T, int N>
CUDA_EXPORT inline general_vect<T, N> general_vect<T, N>::operator+=(const general_vect<T, N> &other) {

	for (int dim = 0; dim < N; dim++) {
		v[dim] += other[dim];
	}
	return *this;
}

template<class T, int N>
CUDA_EXPORT inline general_vect<T, N> general_vect<T, N>::operator*=(T r)  {

	for (int dim = 0; dim < N; dim++) {
		v[dim] *= r;
	}
	return *this;
}

template<class T, int N>
CUDA_EXPORT inline general_vect<T, N> general_vect<T, N>::operator/=(T r)  {

	for (int dim = 0; dim < N; dim++) {
		v[dim] /= r;
	}
	return *this;
}

template<class T, int N>
CUDA_EXPORT inline general_vect<T, N> general_vect<T, N>::operator+(const general_vect<T, N> &other) const {
	general_vect<T, N> result;

	for (int dim = 0; dim < N; dim++) {
		result[dim] = v[dim] + other[dim];
	}
	return result;
}

template<class T, int N>
CUDA_EXPORT inline general_vect<T, N> general_vect<T, N>::operator*(T r) const {
	general_vect<T, N> result;

	for (int dim = 0; dim < N; dim++) {
		result[dim] = v[dim] * r;
	}
	return result;
}

template<class T, int N>
CUDA_EXPORT inline general_vect<T, N> general_vect<T, N>::operator/(T r) const {
	general_vect<T, N> result;

	for (int dim = 0; dim < N; dim++) {
		result[dim] = v[dim] / r;
	}
	return result;
}

template<class T, int N>
CUDA_EXPORT inline T general_vect<T, N>::dot(const general_vect<T, N> &other) const {
	T result = v[0] * other[0];

	for (int dim = 1; dim < N; dim++) {
		result += v[dim] * other[dim];
	}
	return result;
}

template<class T, int N>
CUDA_EXPORT inline T abs(const general_vect<T, N> &v) {
	return sqrt(v.dot(v));
}

template<class T, int N>
CUDA_EXPORT inline general_vect<T, N> abs(const general_vect<T, N> &a, const general_vect<T, N> &b) {
	general_vect<T, N> c;

	for (int i = 0; i < N; i++) {
		c[i] = abs(a[i] - b[i]);
	}
	return c;
}

template<class T, int N>
CUDA_EXPORT inline general_vect<T, N> max(const general_vect<T, N> &a, const general_vect<T, N> &b) {
	general_vect<T, N> c;

	for (int i = 0; i < N; i++) {
		c[i] = max(a[i], b[i]);
	}
	return c;
}

template<class T, int N>
CUDA_EXPORT inline general_vect<T, N> min(const general_vect<T, N> &a, const general_vect<T, N> &b) {
	general_vect<T, N> c;

	for (int i = 0; i < N; i++) {
		c[i] = min(a[i], b[i]);
	}
	return c;
}

template<class T>
using vect = general_vect<T, NDIM>;

#endif /* VECT_HPP_ */
