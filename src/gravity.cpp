#include <tigersph/expansion.hpp>
#include <tigersph/gravity.hpp>
#include <tigersph/options.hpp>
#include <tigersph/simd.hpp>

static const auto one = simd_float(1.0);
static const auto half = simd_float(0.5);
static const simd_float eps = simd_float(std::numeric_limits<float>::min());

std::uint64_t gravity_PP_direct(std::vector<force> &f, const std::vector<vect<float>> &x, std::vector<vect<float>> &y) {
	if (x.size() == 0) {
		return 0;
	}
	std::uint64_t flop = 0;
	static const auto opts = options::get();
	static const simd_float M(opts.m_tot / opts.problem_size);
	static const auto h = opts.soft_len;
	static const auto h2 = h * h;
	static const simd_float huge = std::numeric_limits<float>::max() / 10.0;
	static const simd_float tiny = std::numeric_limits<float>::min() * 10.0;
	static const simd_float H(h);
	static const simd_float H2(h * h);
	static const simd_float Hinv(1.0 / h);
	static const simd_float H3inv(1.0 / h / h / h);
	static const simd_float H5inv(1.0 / h / h / h / h / h);
	static const simd_float H4(h2 * h2);
	static const auto _15o8 = simd_float(15.0 / 8.0);
	static const auto _10o8 = simd_float(10.0 / 8.0);
	static const auto _3o8 = simd_float(3.0 / 8.0);
	static const auto _2p5 = simd_float(2.5);
	static const auto _1p5 = simd_float(1.5);
	static const auto zero = simd_float(0);
	vect<simd_float> X, Y;
	std::vector<vect<simd_float>> G(x.size());
	std::vector<simd_float> Phi(x.size());
	for (int i = 0; i < x.size(); i++) {
		G[i] = vect<simd_float>(0.0);
		Phi[i] = 0.0;
	}
	const auto cnt1 = y.size();
	const auto cnt2 = ((cnt1 - 1 + simd_float::size()) / simd_float::size()) * simd_float::size();
	y.resize(cnt2);
	for (int j = cnt1; j < cnt2; j++) {
		y[j] = vect<float>(1.0e+10);
	}
	for (int j = 0; j < cnt1; j += simd_float::size()) {
		for (int dim = 0; dim < NDIM; dim++) {
			for (int k = 0; k < simd_float::size(); k++) {
				Y[dim][k] = y[j + k][dim];
			}
		}
		for (int i = 0; i < x.size(); i++) {
			for (int dim = 0; dim < NDIM; dim++) {
				X[dim] = simd_float(x[i][dim]);
			}

			vect<simd_float> dX = X - Y;             																		// 3 OP
			const simd_float r2 = dX.dot(dX);																				// 1 OP
			const simd_float zero_mask = r2 > simd_float(0);
			const simd_float rinv = zero_mask * rsqrt(r2 + tiny);       													// 2 OP
			const simd_float r = sqrt(r2);																					// 1 OP
			const simd_float rinv3 = rinv * rinv * rinv;   																	// 2 OP
			const simd_float sw_far = H < r;   																				// 1 OP
			const simd_float sw_near = simd_float(1) - sw_far;																// 1 OP
			const simd_float roverh = min(r * Hinv, 1);																		// 2 OP
			const simd_float roverh2 = roverh * roverh;																		// 1 OP
			const simd_float roverh4 = roverh2 * roverh2;																	// 1 OP
			const simd_float fnear = (_2p5 - _1p5 * roverh2) * H3inv;														// 3 OP
			for (int dim = 0; dim < NDIM; dim++) {
				G[i][dim] -= dX[dim] * M * (sw_far * rinv3 + sw_near * fnear);  											// 18 OP
			}
			const auto tmp = M * zero_mask; 	            																// 1 OP
			const auto near = (_15o8 - _10o8 * roverh2 + _3o8 * roverh4) * Hinv;											// 5 OP
			Phi[i] -= (sw_far * rinv + sw_near * near) * tmp;		        												// 5 OP
		}
	}
	for (int i = 0; i < x.size(); i++) {
		for (int dim = 0; dim < NDIM; dim++) {
			f[i].g[dim] += G[i][dim].sum();
		}
		f[i].phi += Phi[i].sum();
	}
	y.resize(cnt1);
	return (65 * cnt1 + simd_float::size() * 4) * x.size();
}

std::uint64_t gravity_PC_direct(std::vector<force> &f, const std::vector<vect<float>> &x, std::vector<multi_src> &y) {
	if (x.size() == 0) {
		return 0;
	}
	std::uint64_t flop = 0;
	static const auto opts = options::get();
	static const auto h = opts.soft_len;
	static const auto h2 = h * h;
	static const simd_float H(h);
	static const simd_float H2(h2);
	vect<simd_float> X, Y;
	multipole<simd_float> M;
	std::vector<vect<simd_float>> G(x.size(), vect<simd_float>(simd_float(0)));
	std::vector<simd_float> Phi(x.size(), simd_float(0.0));
	const auto cnt1 = y.size();
	const auto cnt2 = ((cnt1 - 1 + simd_float::size()) / simd_float::size()) * simd_float::size();
	y.resize(cnt2);
	for (int j = cnt1; j < cnt2; j++) {
		y[j].m = 0.0;
		y[j].x = vect<float>(1.0);
	}
	for (int j = 0; j < cnt1; j += simd_float::size()) {
		for (int n = 0; n < MP; n++) {
			for (int k = 0; k < simd_float::size(); k++) {
				M[n][k] = y[j + k].m[n];
			}
		}
		for (int dim = 0; dim < NDIM; dim++) {
			for (int k = 0; k < simd_float::size(); k++) {
				Y[dim][k] = y[j + k].x[dim];
			}
		}
		for (int i = 0; i < x.size(); i++) {
			for (int dim = 0; dim < NDIM; dim++) {
				X[dim] = simd_float(x[i][dim]);
			}

			vect<simd_float> dX = X - Y;             		// 3 OP
			auto this_f = multipole_interaction(M, dX); // 517 OP

			for (int dim = 0; dim < NDIM; dim++) {
				G[i][dim] += this_f.second[dim];  // 6 OP
			}
			Phi[i] += this_f.first;		        // 2 OP
		}
	}
	for (int i = 0; i < x.size(); i++) {
		for (int dim = 0; dim < NDIM; dim++) {
			f[i].g[dim] += G[i][dim].sum();
		}
		f[i].phi += Phi[i].sum();
	}
	y.resize(cnt1);
	return (546 * cnt1 + simd_float::size() * 4) * x.size();
}

std::uint64_t gravity_CC_direct(expansion<float> &L, const vect<ireal> &x, std::vector<multi_src> &y) {
	if (y.size() == 0) {
		return 0;
	}
	static const auto one = simd_real(1.0);
	static const auto half = simd_real(0.5);
	std::uint64_t flop = 0;
	static const auto opts = options::get();
	static const auto h = opts.soft_len;
	static const auto h2 = h * h;
	static const simd_real H(h);
	static const simd_real H2(h2);
	vect<simd_real> X, Y;
	multipole<simd_real> M;
	expansion<simd_real> Lacc;
	Lacc = simd_real(0);
	const auto cnt1 = y.size();
	const auto cnt2 = ((cnt1 - 1 + simd_real::size()) / simd_real::size()) * simd_real::size();
	y.resize(cnt2);
	for (int j = cnt1; j < cnt2; j++) {
		y[j].m = 0.0;
		y[j].x = vect<ireal>(1.0);
	}

	for (int dim = 0; dim < NDIM; dim++) {
		X[dim] = simd_real(x[dim]);
	}
	for (int j = 0; j < cnt1; j += simd_real::size()) {
		for (int n = 0; n < MP; n++) {
			for (int k = 0; k < simd_real::size(); k++) {
				M[n][k] = y[j + k].m[n];
			}
		}
		for (int dim = 0; dim < NDIM; dim++) {
			for (int k = 0; k < simd_real::size(); k++) {
				Y[dim][k] = y[j + k].x[dim];
			}
		}
		vect<simd_real> dX = X - Y;             										// 3 OP
		multipole_interaction(Lacc, M, dX);												// 670 OP
	}

	for (int i = 0; i < LP; i++) {
		L[i] += Lacc[i].sum();
	}
	y.resize(cnt1);
	return 691 * cnt1 + LP * simd_float::size();
}

std::uint64_t gravity_CP_direct(expansion<float> &L, const vect<ireal> &x, std::vector<vect<float>> &y) {
	if (y.size() == 0) {
		return 0;
	}
	static const auto one = simd_real(1.0);
	static const auto half = simd_real(0.5);
	std::uint64_t flop = 0;
	static const auto opts = options::get();
	static const auto m = opts.m_tot / opts.problem_size;
	static const auto h = opts.soft_len;
	static const auto h2 = h * h;
	static const simd_real H(h);
	static const simd_real H2(h2);
	vect<simd_real> X, Y;
	simd_real M;
	expansion<simd_real> Lacc;
	Lacc = simd_real(0);
	const auto cnt1 = y.size();
	const auto cnt2 = ((cnt1 - 1 + simd_real::size()) / simd_real::size()) * simd_real::size();
	y.resize(cnt2);
	for (int j = cnt1; j < cnt2; j++) {
		y[j] = vect<ireal>(1.0);
	}

	for (int dim = 0; dim < NDIM; dim++) {
		X[dim] = simd_real(x[dim]);
	}
	M = m;
	for (int j = 0; j < cnt1; j += simd_real::size()) {
		for (int dim = 0; dim < NDIM; dim++) {
			for (int k = 0; k < simd_real::size(); k++) {
				Y[dim][k] = y[j + k][dim];
			}
		}
		if (j + simd_real::size() > cnt1) {
			for (int k = cnt1; k < cnt2; k++) {
				M[k - j] = 0.0;
			}
		}
		vect<simd_real> dX = X - Y;             										// 3 OP
		multipole_interaction(Lacc, M, dX);												// 390 OP
	}

	for (int i = 0; i < LP; i++) {
		L[i] += Lacc[i].sum();
	}
	y.resize(cnt1);
	return 408 * cnt1 + simd_float::size() * LP;
}
