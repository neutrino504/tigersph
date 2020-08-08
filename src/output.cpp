/*
 * bin2silo.cpp
 *
 *  Created on: Jul 5, 2020
 *      Author: dmarce1
 */

#include <tigersph/options.hpp>
#include <tigersph/output.hpp>
#include <silo.h>

#include <set>
#include <thread>
#include <fenv.h>

error compute_error(std::vector<output> test, std::vector<output> direct) {
	error err;
	std::multiset<double> top_errs;
	err.err = 0.0;
	err.g = vect<float>(0);
	if (test.size() != direct.size()) {
		ERROR();
	}
	for (int i = 0; i < test.size(); i++) {
		for (int dim = 0; dim < NDIM; dim++) {
			if (test[i].x[dim] != direct[i].x[dim]) {
				printf("Test set does not match reference\n");
				abort();
			}
		}
		const auto dg = abs(test[i].g - direct[i].g);
		const auto gd = abs(test[i].g);
		const auto this_err = dg / gd;
		err.g = err.g + test[i].g;
		err.err += this_err;
		if (top_errs.size() < test.size() / 100) {
			top_errs.insert(this_err);
		} else {
			if (this_err > *(top_errs.begin())) {
				top_errs.erase(top_errs.begin());
				top_errs.insert(this_err);
			}
		}
	}
	err.g = err.g / test.size();
	err.err /= test.size();
	err.err99 = *(top_errs.begin());
	return err;
}

void output_particles(const std::vector<output> &parts, const std::string filename) {
	std::thread([&] {
		static const auto opts = options::get();
		static const float m = 1.0 / opts.problem_size;
		std::array<std::vector<double>, NDIM> x;
		std::array<std::vector<float>, NDIM> g;
		std::array<std::vector<float>, NDIM> v;
		std::vector<float> phi;
		std::vector<int> rung;

		printf( "SILO out\n");
		for (auto i = parts.begin(); i != parts.end(); i++) {
			for (int dim = 0; dim < NDIM; dim++) {
				x[dim].push_back(i->x[dim]);
				g[dim].push_back(i->g[dim]);
				v[dim].push_back(i->v[dim]);
			}
			rung.push_back(i->rung);
			phi.push_back(i->phi);
		}
		DBfile *db = DBCreateReal(filename.c_str(), DB_CLOBBER, DB_LOCAL, "Meshless", DB_PDB);
		const int nparts = phi.size();
		double *coords[NDIM] = {x[0].data(), x[1].data(), x[2].data()};
		DBPutPointmesh(db, "points", NDIM, coords, nparts, DB_DOUBLE, NULL);
		for (int dim = 0; dim < NDIM; dim++) {
			std::string nm = std::string() + "v_" + char('x' + char(dim));
			DBPutPointvar1(db, nm.c_str(), "points", v[dim].data(), nparts, DB_FLOAT, NULL);
		}
		for (int dim = 0; dim < NDIM; dim++) {
			std::string nm = std::string() + "g_" + char('x' + char(dim));
			DBPutPointvar1(db, nm.c_str(), "points", g[dim].data(), nparts, DB_FLOAT, NULL);
		}
		DBPutPointvar1(db, "phi", "points", phi.data(), nparts, DB_FLOAT, NULL);
		DBPutPointvar1(db, "rung", "points", rung.data(), nparts, DB_INT, NULL);
		DBClose(db);}
	).join();
}
