/*
 * math.hpp
 *
 *  Created on: Nov 29, 2019
 *      Author: dmarce1
 */

#ifndef RANGE_HPP_
#define RANGE_HPP_

#include <tigersph/vect.hpp>

#include <limits>


using pos_type = float;


struct range {
	vect<pos_type> min;
	vect<pos_type> max;
	template<class Arc>
	void serialize(Arc& arc,unsigned) {
		arc & min;
		arc & max;
	}
};

range reflect_range(const range&, int dim, pos_type axis);
vect<pos_type> range_center(const range &r);
range shift_range(const range& r, const vect<pos_type>&);
range scale_range(const range& , pos_type);
vect<pos_type> range_span(const range&);
bool in_range(const vect<pos_type>&, const range&);
bool in_range(const range&, const range&);
bool ranges_intersect(const range&, const range&);
range range_around(const vect<pos_type>&, pos_type);
range expand_range(const range&, pos_type);
pos_type range_volume(const range&);
range null_range();
bool operator==(const range&, const range&);
bool operator!=(const range&, const range&);

#endif /* MATH_HPP_ */

