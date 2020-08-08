/*
 * position.hpp
 *
 *  Created on: Jul 11, 2020
 *      Author: dmarce1
 */

#ifndef TIGERGRAV_POSITION_HPP_
#define TIGERGRAV_POSITION_HPP_

#include <tigersph/vect.hpp>

using pos_type = std::uint32_t;



inline double pos_to_double(pos_type x) {
	return ((double) x + (double) 0.5) / ((double) std::numeric_limits<pos_type>::max() + (double) 1.0);
}


inline vect<double> pos_to_double(vect<pos_type> x) {
	vect<double> f;
	for (int dim = 0; dim < NDIM; dim++) {
		f[dim] = pos_to_double(x[dim]);
	}
	return f;
}

inline pos_type double_to_pos(double x) {
	return x * ((double) std::numeric_limits<pos_type>::max() + (double) 1.0) + 0.5;
}

inline vect<pos_type> double_to_pos(vect<double> d) {
	vect<pos_type> x;
	for (int dim = 0; dim < NDIM; dim++) {
		x[dim] = double_to_pos(d[dim]);
	}
	return x;
}



#endif /* TIGERGRAV_POSITION_HPP_ */
