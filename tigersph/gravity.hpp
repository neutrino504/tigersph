#pragma once

#include <tigersph/defs.hpp>
#include <tigersph/simd.hpp>
#include <tigersph/expansion.hpp>
#include <tigersph/multipole.hpp>
#include <tigersph/particle.hpp>

#include <memory>
#include <map>

struct multi_src {
	multipole<ireal> m;
	vect<ireal> x;
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & m;
		arc & x;
	}

};


struct multipole_info {
	multipole<ireal> m;
	vect<ireal> x;
	ireal r;
	bool has_active;
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & m;
		arc & x;
		arc & r;
		arc & has_active;
	}
};


std::uint64_t gravity_PP_direct(std::vector<force> &g, const std::vector<vect<float>> &x, std::vector<vect<float>> &y);
std::uint64_t gravity_PC_direct(std::vector<force> &g, const std::vector<vect<float>> &x, std::vector<multi_src> &y);
std::uint64_t gravity_CC_direct(expansion<float>&, const vect<ireal> &x, std::vector<multi_src> &y);
std::uint64_t gravity_CP_direct(expansion<float> &L, const vect<ireal> &x, std::vector<vect<float>> &y);
std::uint64_t gravity_PP_ewald(std::vector<force> &g, const std::vector<vect<float>> &x, std::vector<vect<float>> &y);
std::uint64_t gravity_PC_ewald(std::vector<force> &g, const std::vector<vect<float>> &x, std::vector<multi_src> &y);
std::uint64_t gravity_CC_ewald(expansion<float>&, const vect<ireal> &x, std::vector<multi_src> &y);
std::uint64_t gravity_CP_ewald(expansion<float> &L, const vect<ireal> &x, std::vector<vect<float>> &y);


float ewald_near_separation(const vect<float> x);
float ewald_far_separation(const vect<float> x, float r);
void init_ewald();
