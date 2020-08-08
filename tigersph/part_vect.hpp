#pragma once

#include <tigersph/particle.hpp>
#include <tigersph/output.hpp>
#include <tigersph/options.hpp>
#include <tigersph/range.hpp>
#include <tigersph/gravity.hpp>

#ifdef HPX_LITE
#include <hpx/hpx_lite.hpp>
#else
#include <hpx/include/async.hpp>
#endif

using part_iter = std::uint64_t;
using const_part_iter = std::uint64_t;


struct statistics {
	vect<float> g;
	vect<float> p;
	float pot;
	float kin;

	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & g;
		arc & p;
		arc & pot;
		arc & kin;
	}

	void zero() {
		pot = kin = 0.0;
		g = p = vect<float>(0);
	}
	statistics operator+(const statistics &other) const {
		statistics C;
		C.g = g + other.g;
		C.p = p + other.p;
		C.pot = pot + other.pot;
		C.kin = kin + other.kin;
		return C;
	}
};



struct kick_return {
	statistics stats;
	rung_type rung;
	std::vector<output> out;
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & stats;
		arc & rung;
		arc & out;
	}
};


void part_vect_write(part_iter b, part_iter e, std::vector<particle> these_parts);
hpx::future<std::vector<particle>> part_vect_read(part_iter b, part_iter e);
void part_vect_init();
hpx::future<std::vector<vect<pos_type>>> part_vect_read_position(part_iter b, part_iter e);
part_iter part_vect_sort(part_iter b, part_iter e, double mid, int dim);
range part_vect_range(part_iter b, part_iter e);
int part_vect_locality_id(part_iter);
void part_vect_cache_reset();
std::pair<float, vect<float>> part_vect_center_of_mass(part_iter b, part_iter e);
multipole_info part_vect_multipole_info(vect<float> com, rung_type mrung, part_iter b, part_iter e);
void part_vect_drift(float dt);
std::vector<vect<float>> part_vect_read_active_positions(part_iter b, part_iter e, rung_type rung);
kick_return part_vect_kick(part_iter b, part_iter e, rung_type rung, bool do_out, std::vector<force>&& f);

