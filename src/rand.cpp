/*
 * rand.cpp
 *
 *  Created on: Dec 6, 2019
 *      Author: dmarce1
 */

#include <cstdlib>

#include <tigersph/rand.hpp>

double rand_unit_box() {
	return (rand() + 0.5) / (RAND_MAX + 1.0) - 0.5;
}

double rand1() {
	return (rand() + 0.5) / (RAND_MAX + 1.0);
}

double rand_sign() {
	if( rand() % 2 ) {
		return +1.0;
	} else {
		return -1.0;
	}
}

vect<double> rand_unit_vect() {
	vect<double> n;
	double sum = 0.0;
	for (int dim = 0; dim < NDIM; dim++) {
		n[dim] = rand_unit_box();
		sum += n[dim] * n[dim];
	}
	sum = 1.0 / sqrt(sum);
	for (int dim = 0; dim < NDIM; dim++) {
		n[dim] *= sum;
	}
	return n;
}

double rand_normal() {
	const auto u1 = rand1();
	const auto u2 = rand1();
	return std::sqrt(-2.0*log(u1)) * cos(2.0 * M_PI * u2);
}
