#include <tigersph/initialize.hpp>
#include <tigersph/options.hpp>
#include <tigersph/rand.hpp>

#ifdef HPX_LITE
#include <hpx/hpx_lite.hpp>
#else
#include <hpx/include/components.hpp>
#endif

part_vect initial_particle_set(std::string pn, int N, int Nout) {
	part_vect parts;
	parts.reserve(N);
	srand(hpx::get_locality_id() * 0xA3CF98A7);
	if (pn == "collapse") {
		for (int i = 0; i < N; i++) {
			particle p;
			do {
				for (int dim = 0; dim < NDIM; dim++) {
					p.x[dim] = 2.0 * rand1() - 1.0;
					p.v[dim] = 0.0;
				}
			} while (abs(p.x) > 1.0);
			parts.push_back(std::move(p));
		}
	} else {
		printf("Problem %s unknown\n", pn.c_str());
		abort();
	}
	int j = 0;
	for (auto i = parts.begin(); i != parts.end(); i++, j++) {
		i->rung = 0;
		i->flags.out = j < Nout;
	}
	return parts;
}
