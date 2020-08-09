/*
 * particle.hpp
 *
 *  Created on: Jun 2, 2020
 *      Author: dmarce1
 */

#ifndef TIGERGRAV_PARTICLE_HPP_
#define TIGERGRAV_PARTICLE_HPP_

#include <tigersph/range.hpp>
#include <tigersph/time.hpp>
#include <tigersph/vect.hpp>

#include <vector>



struct particle {
	vect<pos_type> x;
	vect<float> v;
	rung_type rung;
	struct {
		std::uint8_t out :1;
	} flags;
	template<class A>
	void serialize(A &&arc, unsigned) {
		int tmp;
		arc & x;
		arc & v;
		arc & rung;
		tmp = flags.out;
		arc & tmp;
		flags.out = tmp;
	}

};

using part_vect = std::vector<particle>;

template<class I, class F>
I bisect(I begin, I end, F &&below) {
	auto lo = begin;
	auto hi = end - 1;
	while (lo < hi) {
		if (!below(*lo)) {
			while (lo != hi) {
				if (below(*hi)) {
					auto tmp = *lo;
					*lo = *hi;
					*hi = tmp;
					break;
				}
				hi--;
			}

		}
		lo++;
	}
	return hi;
}

#endif /* TIGERGRAV_PARTICLE_HPP_ */
