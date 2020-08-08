#pragma once

#include <tigersph/vect.hpp>

#include <vector>

struct output {
	vect<double> x;
	vect<float> v;
	vect<float> g;
	float phi;
	int rung;
	bool operator<(const output &other) {
		for (int dim = 0; dim < NDIM; dim++) {
			if (x[dim] < other.x[dim]) {
				return true;
			} else if (x[dim] > other.x[dim]) {
				return false;
			}
		}
		return false;
	}
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & x;
		arc & v;
		arc & g;
		arc & phi;
		arc & rung;
	}
};

struct error {
	double err;
	double err99;
	vect<double> g;
};

void output_particles(const std::vector<output> &parts, const std::string filename);

error compute_error(std::vector<output> test, std::vector<output> direct);
