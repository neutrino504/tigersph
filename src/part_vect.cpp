#include <tigersph/part_vect.hpp>
#include <tigersph/initialize.hpp>
#include <tigersph/simd.hpp>

#ifdef HPX_LITE
#include <hpx/hpx_lite.hpp>
#else
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/include/async.hpp>
#include <hpx/include/serialization.hpp>
#include <hpx/synchronization/spinlock.hpp>
#include <hpx/local_lcos/promise.hpp>
#endif

#include <unistd.h>
#include <thread>
#include <unordered_map>

using mutex_type = hpx::lcos::local::spinlock;

static std::vector<particle> particles;
static part_iter part_begin;
static part_iter part_end;
static std::vector<hpx::id_type> localities;
static int myid;

#define POS_CACHE_SIZE 1024
static std::unordered_map<part_iter, hpx::shared_future<std::vector<vect<pos_type>>>> pos_cache[POS_CACHE_SIZE];
static mutex_type pos_cache_mtx[POS_CACHE_SIZE];

HPX_PLAIN_ACTION (part_vect_init);
HPX_PLAIN_ACTION (part_vect_read_position);
HPX_PLAIN_ACTION (part_vect_range);
HPX_PLAIN_ACTION (part_vect_cache_reset);
HPX_PLAIN_ACTION (part_vect_center_of_mass);
HPX_PLAIN_ACTION (part_vect_multipole_info);
HPX_PLAIN_ACTION (part_vect_drift);
HPX_PLAIN_ACTION (part_vect_read_active_positions);
HPX_PLAIN_ACTION (part_vect_kick);

inline particle& parts(part_iter i) {
	const int j = i - part_begin;
#ifndef NDEBUG
	if (j < 0 || j >= particles.size()) {
		printf("Index out of bounds! %i should be between 0 and %i\n", j, particles.size());
		abort();
	}
#endif
	return particles[j];
}

int round_robin(int me, int round, int N) {
	int nrounds = N % 2 == 0 ? N : N + 1;
	if (N % 2 == 0) {
		return (me + round % nrounds) % N;
	} else {
		int tmp = (me + round % nrounds) % (N + 1);
		if (tmp == N) {
			return -1;
		} else {
			return tmp;
		}
	}
}

void part_vect_sort_begin(part_iter b, part_iter e, part_iter mid, float xmid, int dim);
std::vector<particle> part_vect_sort_end(part_iter b, part_iter e, part_iter mid, float xmid, int dim, std::vector<particle>);
pos_type part_vect_find_median_helper(part_iter b, part_iter e, part_iter median, pos_type xmid, int dim);

HPX_PLAIN_ACTION (part_vect_sort_begin);
HPX_PLAIN_ACTION (part_vect_sort_end);
HPX_PLAIN_ACTION (part_vect_avg_pos);
HPX_PLAIN_ACTION (part_vect_find_median_helper);

mutex_type sort_mutex;

avg_pos_return part_vect_avg_pos(part_iter b, part_iter e, int dim) {
	static const auto opts = options::get();
	std::vector < hpx::future < avg_pos_return >> futs;
	const part_iter N = localities.size();
	const part_iter M = opts.problem_size;
	pos_type min = +std::numeric_limits < pos_type > ::max();
	pos_type max = -std::numeric_limits < pos_type > ::max();
	if (part_end < e) {
		for (int n = part_vect_locality_id(part_end); n * M / N < e; n++) {
			auto fut = hpx::async < part_vect_avg_pos_action > (localities[n], n * M / N, std::min(e, (n + 1) * M / N), dim);
			futs.push_back(std::move(fut));
		}
	}
	const auto this_e = std::min(part_end, e);
	pos_type avg = 0.0;
	std::uint64_t count = this_e - b;
	for (auto i = b; i < this_e; i++) {
		const auto x = parts(i).x[dim];
		avg += x;
		min = std::min(min, x);
		max = std::max(max, x);
	}
	int j = 0;
	if (part_end < e) {
		for (int n = part_vect_locality_id(part_end); n * M / N < e; n++) {
			const auto rc = futs[j++].get();
			const auto this_avg = rc.avg;
			const auto this_count = std::min(e, (n + 1) * M / N) - n * M / N;
			count += this_count;
			avg += this_avg * this_count;
			min = std::min(min, rc.min);
			max = std::max(max, rc.max);
		}
	}
	avg_pos_return rc;
	avg /= count;
	rc.avg = avg;
	rc.max = max;
	rc.min = min;
	return rc;
}

pos_type part_vect_find_median_helper(part_iter b, part_iter e, part_iter median, pos_type xmid, int dim) {
	pos_type r;
//	printf("%i %i %i %.8e\n", b, e, median, xmid);
//	if (e - b == 2) {
//		printf( "!\n");
//		for (auto i = b; i < std::min(part_end, e); i++) {
//			printf("%.8e\n", parts(i).x[dim]);
//		}
//
//	}
	if (e - b <= 1) {
		r = xmid;
	} else {
		auto m = part_vect_sort(b, e, xmid, dim);
		if (median > m && median < e) {
			const auto new_b = m;
			const auto new_e = e;
			const auto new_loc = localities[part_vect_locality_id(new_b)];
			const auto rc = hpx::async < part_vect_avg_pos_action > (new_loc, new_b, new_e, dim).get();
//			printf("%i %i %i %i %.8e\n", b, e, m, median, xmid);
//			printf( "- %.8e %.8e %.8e\n", rc.min, rc.max, rc.avg);
			if (almost_equal(rc.max, rc.min)) {
				r = rc.max;
			} else {
				const auto new_xmid = rc.avg;
				r = hpx::async < part_vect_find_median_helper_action > (new_loc, new_b, new_e, median, new_xmid, dim).get();
			}
		} else if (median >= b && median < m) {
			const auto new_b = b;
			const auto new_e = m;
			const auto rc = part_vect_avg_pos(new_b, new_e, dim);
//			printf("%i %i %i %i %e\n", b, e, m, median, xmid);
//			printf( "- %.8e %.8e %.8e\n", rc.min, rc.max, rc.avg);
			if (almost_equal(rc.max, rc.min)) {
				r = rc.max;
			} else {
				const auto new_xmid = rc.avg;
				r = part_vect_find_median_helper(new_b, new_e, median, new_xmid, dim);
			}
		} else if (median == m) {
			return xmid;
		} else {
			assert(false);
		}
	}
	return r;
}

pos_type part_vect_find_median(part_iter b, part_iter e, int dim) {
	const auto xmid = part_vect_avg_pos(b, e, dim).avg;
	return part_vect_find_median_helper(b, e, (b + e) / 2, xmid, dim);
}

void part_vect_sort_begin(part_iter b, part_iter e, part_iter mid, float xmid, int dim) {
	static const auto opts = options::get();
	if (b == mid || e == mid) {
//		printf("%i %i %i\n", b, mid, e);
		return;
	}
	std::vector<hpx::future<void>> futs;
	if (part_vect_locality_id(b) == myid) {
		for (int n = myid + 1; n <= part_vect_locality_id(mid); n++) {
			futs.push_back(hpx::async < part_vect_sort_begin_action > (localities[n], b, e, mid, xmid, dim));
		}
	}
	int nproc = part_vect_locality_id(e - 1) - part_vect_locality_id(b) + 1;
	int low_proc = part_vect_locality_id(b);
	int chunk_size = opts.problem_size / localities.size() / 10 + 1;
	bool done = false;
	part_iter begin = std::max(b, part_begin);
	for (int round = 0; !done; round++) {
		int other = round_robin(myid - low_proc, round, nproc) + low_proc;
//		printf("%i %i %i %i %i %i\n", myid, b, mid, e, round % nproc, other);
		if (other >= 0 && other >= part_vect_locality_id(mid)) {
			std::vector<particle> send;
			send.reserve(chunk_size);
			done = true;
			std::unique_lock<mutex_type> lock(sort_mutex);
//			printf( "--%i %i %i\n",myid, begin,std::min(part_end, mid) );
			for (part_iter i = begin; i < std::min(part_end, mid); i++) {
				if (parts(i).x[dim] >= xmid) {
					send.push_back(parts(i));
					done = false;
					if (send.size() >= chunk_size) {
						break;
					}
				}
			}
			if (!done) {
				//			printf("%i Sending %i to %i\n", myid, send.size(), other);
				lock.unlock();
				auto recv = part_vect_sort_end_action()(localities[other], b, e, mid, xmid, dim, std::move(send));
				lock.lock();
				int j = 0;
				if (recv.size()) {
					for (part_iter i = begin; i < std::min(part_end, mid); i++) {
						begin = i;
						if (parts(i).x[dim] >= xmid) {
							parts(i) = recv[j++];
							if (recv.size() == j) {
								break;
							}
						}
					}
				}
			}
		}
	}
	hpx::wait_all(futs.begin(), futs.end());
}

std::vector<particle> part_vect_sort_end(part_iter b, part_iter e, part_iter mid, float xmid, int dim, std::vector<particle> low) {
	std::lock_guard<mutex_type> lock(sort_mutex);
	int j = 0;
	for (part_iter i = std::max(mid, part_begin); i < std::min(part_end, e); i++) {
		if (parts(i).x[dim] < xmid) {
			auto tmp = parts(i);
			parts(i) = low[j];
			low[j++] = tmp;
			if (j == low.size()) {
				break;
			}
		}
	}
	low.resize(j);
//	if(myid==1)
//	printf("%i used %i\n", myid, j);
	return low;
}

HPX_PLAIN_ACTION (part_vect_read);
HPX_PLAIN_ACTION (part_vect_write);

void part_vect_write(part_iter b, part_iter e, std::vector<particle> these_parts) {
	const auto id = part_vect_locality_id(b);
	if (id == myid) {
		int i = b;
		for (auto this_part : these_parts) {
			parts(i) = this_part;
			i++;
			if (i == part_end) {
				break;
			}
		}
		std::vector<particle> next_parts;
		while (i < e) {
			next_parts.push_back(these_parts[i - b]);
			i++;
		}
		if (next_parts.size()) {
			part_vect_write_action()(localities[myid + 1], e - next_parts.size(), e, std::move(next_parts));
		}
	} else {
		part_vect_write_action()(localities[id], b, e, std::move(these_parts));
	}
}

hpx::future<std::vector<particle>> part_vect_read(part_iter b, part_iter e) {
	const auto id = part_vect_locality_id(b);
	std::vector<particle> these_parts;
	if (id == myid) {
		these_parts.reserve(e - b);
		const auto this_e = std::min(e, part_end);
		for (int i = b; i < this_e; i++) {
			these_parts.push_back(parts(i));
		}
		if (these_parts.size() != e - b) {
			auto fut = part_vect_read_action()(localities[myid + 1], b + these_parts.size(), e);
			auto next_parts = fut.get();
			for (const auto &p : next_parts) {
				these_parts.push_back(p);
			}
			if (these_parts.size() != e - b) {
				printf("Error in part_vect_read\n");
				abort();
			}
		}
	} else {
		auto fut = part_vect_read_action()(localities[id], b, e);
		these_parts = fut.get();
	}
	return hpx::make_ready_future(these_parts);
}
kick_return part_vect_kick(part_iter b, part_iter e, rung_type min_rung, bool do_out, std::vector<force> &&f) {
	kick_return rc;
	const auto opts = options::get();
	const float eps = 10.0 * std::numeric_limits<float>::min();
	const float m = opts.m_tot / opts.problem_size;
	rc.rung = 0;
	part_iter j = 0;
	if (do_out) {
		rc.stats.zero();
	}
	for (auto i = b; i != std::min(e, part_end); i++) {
		if (parts(i).rung >= min_rung || do_out) {
			if (parts(i).rung >= min_rung) {
				if (parts(i).rung != 0) {
					const float dt = rung_to_dt(parts(i).rung);
					auto &v = parts(i).v;
					v = v + f[j].g * 0.5 * dt;
				}
				const float a = abs(f[j].g);
				float dt = std::min(opts.dt_max, opts.eta * std::sqrt(opts.soft_len / (a + eps)));
				rung_type rung = dt_to_rung(dt);
				rung = std::max(rung, min_rung);
				rc.rung = std::max(rc.rung, rung);
				dt = rung_to_dt(rung);
				parts(i).rung = std::max(std::max(rung, rung_type(parts(i).rung - 1)), (rung_type) 1);
				auto &v = parts(i).v;
				v = v + f[j].g * 0.5 * dt;
			}
			if (do_out) {
				rc.stats.g = rc.stats.g + f[j].g * m;
				rc.stats.p = rc.stats.p + parts(i).v * m;
				rc.stats.pot += 0.5 * m * f[j].phi;
				rc.stats.kin += 0.5 * m * parts(i).v.dot(parts(i).v);
				if (parts(i).flags.out) {
					output out;
					out.x = parts(i).x;
					out.v = parts(i).v;
					out.g = f[j].g;
					out.phi = f[j].phi;
					out.rung = parts(i).rung;
					rc.out.push_back(out);
				}
			}
			j++;
		}
	}
	const auto k = j;
	for (; j < f.size(); j++) {
		f[j - k] = f[j];
	}
	f.resize(f.size() - k);
	if (f.size()) {
		const auto other = part_vect_kick_action()(localities[myid + 1], part_end, e, min_rung, do_out, std::move(f));
		if (do_out) {
			rc.stats = rc.stats + other.stats;
			rc.rung = std::max(rc.rung, other.rung);
			rc.out.insert(rc.out.end(), other.out.begin(), other.out.end());
		}
	}
	return rc;
}

std::vector<vect<float>> part_vect_read_active_positions(part_iter b, part_iter e, rung_type rung) {
	std::vector<vect<float>> x;
	hpx::future < std::vector<vect<float>> > fut;
	x.reserve(e - b);
	if (e > part_end) {
		fut = hpx::async < part_vect_read_active_positions_action > (localities[myid + 1], part_end, e, rung);
	}
	for (part_iter i = b; i < std::min(e, part_end); i++) {
		if (parts(i).rung >= rung) {
			x.push_back(parts(i).x);
		}
	}
	if (e > part_end) {
		auto other = fut.get();
		for (const auto &o : other) {
			x.push_back(o);
		}
	}
	return x;
}

void part_vect_drift(float dt) {
	static const auto opts = options::get();
	std::vector<hpx::future<void>> futs;
	if (myid == 0) {
		for (int i = 1; i < localities.size(); i++) {
			futs.push_back(hpx::async < part_vect_drift_action > (localities[i], dt));
		}
	}
	const part_iter chunk_size = std::max(part_iter(1), (part_end - part_begin) / std::thread::hardware_concurrency());
	for (part_iter i = part_begin; i < part_end; i += chunk_size) {
		auto func = [i, chunk_size, dt]() {
			const auto end = std::min(part_end, i + chunk_size);
			for (int j = i; j < end; j++) {
				const vect<double> dx = parts(j).v * dt;
				vect<double> x = parts(j).x;
				x += dx;
				parts(j).x = x;
			}
		};
		futs.push_back(hpx::async(std::move(func)));
	}
	hpx::wait_all(futs.begin(), futs.end());
}

std::pair<float, vect<float>> part_vect_center_of_mass(part_iter b, part_iter e) {
	static const auto m = options::get().m_tot / options::get().problem_size;
	std::pair<float, vect<float>> rc;
	hpx::future < std::pair<float, vect<float>> > fut;
	if (e > part_end) {
		fut = hpx::async < part_vect_center_of_mass_action > (localities[myid + 1], part_end, e);
	}
	const auto this_end = std::min(part_end, e);
	rc.first = 0.0;
	rc.second = vect<float>(0.0);
	for (part_iter i = b; i < this_end; i++) {
		rc.first += m;
		rc.second += parts(i).x * m;
	}
	if (e > part_end) {
		auto tmp = fut.get();
		rc.first += tmp.first;
		rc.second = rc.second + tmp.second * tmp.first;
	}
	if (rc.first > 0.0) {
		rc.second = rc.second / rc.first;
	}
	return rc;
}

multipole_info part_vect_multipole_info(vect<float> com, rung_type mrung, part_iter b, part_iter e) {
	static const auto opts = options::get();
	static const auto m = opts.m_tot / opts.problem_size;
	multipole_info rc;
	hpx::future<multipole_info> fut;
	if (e > part_end) {
		fut = hpx::async < part_vect_multipole_info_action > (localities[myid + 1], com, mrung, part_end, e);
	}
	const auto this_end = std::min(part_end, e);
	rc.m = 0.0;
	rc.x = com;
	multipole<simd_float> M;
	vect<simd_float> Xcom;
	M = simd_float(0.0);
	for (int dim = 0; dim < NDIM; dim++) {
		Xcom[dim] = simd_float(com[dim]);
	}
	for (part_iter i = b; i < this_end; i += simd_float::size()) {
		vect<simd_float> X;
		simd_float mass;
		for (int k = 0; k < simd_float::size(); k++) {
			if (i + k < this_end) {
				for (int dim = 0; dim < NDIM; dim++) {
					X[dim][k] = parts(i + k).x[dim];
				}
				mass[k] = m;
			} else {
				for (int dim = 0; dim < NDIM; dim++) {
					X[dim][k] = 0.0;
				}
				mass[k] = 0.0;
			}
		}
		const auto dx = X - Xcom;
		M() += mass;
		for (int j = 0; j < NDIM; j++) {
			for (int k = 0; k <= j; k++) {
				const auto Xjk = dx[j] * dx[k];
				M(j, k) += mass * Xjk;
				for (int l = 0; l <= k; l++) {
					M(j, k, l) += mass * Xjk * dx[l];
				}
			}
		}
	}
	for (int i = 0; i < MP; i++) {
		rc.m[i] = M[i].sum();
	}
	rc.r = 0.0;
	rc.has_active = false;
	for (part_iter i = b; i < this_end; i++) {
		rc.r = std::max(rc.r, (ireal) abs(parts(i).x - rc.x));
		if (parts(i).rung >= mrung) {
			rc.has_active = true;
		}
	}
	if (e > part_end) {
		auto tmp = fut.get();
		for (int i = 0; i < MP; i++) {
			rc.m[i] += tmp.m[i];
		}
		rc.r = std::max(rc.r, tmp.r);
		rc.has_active = rc.has_active || tmp.has_active;
	}
	return rc;
}

void part_vect_cache_reset() {
	std::vector<hpx::future<void>> futs;
	if (myid == 0) {
		for (int i = 1; i < localities.size(); i++) {
			futs.push_back(hpx::async < part_vect_cache_reset_action > (localities[i]));
		}
	}
	for (int i = 0; i < POS_CACHE_SIZE; i++) {
		pos_cache[i].clear();
	}
	hpx::wait_all(futs.begin(), futs.end());
}

void part_vect_init() {
	localities = hpx::find_all_localities();
	myid = hpx::get_locality_id();

	const auto opts = options::get();

	const part_iter N = localities.size();
	const part_iter n = myid;
	const part_iter M = opts.problem_size;
	std::vector<hpx::future<void>> futs;
	if (n == 0) {
		for (int i = 1; i < N; i++) {
			futs.push_back(hpx::async < part_vect_init_action > (localities[i]));
		}
	}
	part_begin = n * M / N;
	part_end = (n + 1) * M / N;
	particles = initial_particle_set(opts.problem, part_end - part_begin, (n + 1) * opts.out_parts / N - n * opts.out_parts / N);
	hpx::wait_all(futs.begin(), futs.end());
}

inline hpx::future<std::vector<vect<pos_type>>> part_vect_read_pos_cache(part_iter b, part_iter e) {
	const int index = (b / sizeof(particle)) % POS_CACHE_SIZE;
	std::unique_lock<mutex_type> lock(pos_cache_mtx[index]);
	auto iter = pos_cache[index].find(b);
	if (iter == pos_cache[index].end()) {
		hpx::lcos::local::promise < hpx::future < std::vector<vect<pos_type>> >> promise;
		auto fut = promise.get_future();
		pos_cache[index][b] = fut.then([b](decltype(fut) f) {
			return f.get().get();
		});
		lock.unlock();
		promise.set_value(hpx::async < part_vect_read_position_action > (localities[part_vect_locality_id(b)], b, e));
	}
	return hpx::async(hpx::launch::deferred, [b, index]() {
		std::unique_lock < mutex_type > lock(pos_cache_mtx[index]);
		auto future = pos_cache[index][b];
		lock.unlock();
		return future.get();
	});
}

hpx::future<std::vector<vect<pos_type>>> part_vect_read_position(part_iter b, part_iter e) {
	if (b >= part_begin && b < part_end) {
		if (e <= part_end) {
			return hpx::async(hpx::launch::deferred, [=]() {
				std::vector<vect<pos_type>> these_parts;
				these_parts.reserve(e - b);
				for (int i = b; i < e; i++) {
					these_parts.push_back(parts(i).x);
				}
				return these_parts;
			});
		} else {
			auto fut = part_vect_read_pos_cache(part_end, e);
			return hpx::async(hpx::launch::deferred, [=](decltype(fut) f) {
				std::vector<vect<pos_type>> these_parts;
				these_parts.reserve(e - b);
				for (int i = b; i < part_end; i++) {
					these_parts.push_back(parts(i).x);
				}
				auto next_parts = f.get();
				for (const auto &p : next_parts) {
					these_parts.push_back(p);
				}
				return these_parts;
			}, std::move(fut));

		}
	} else {
		return part_vect_read_pos_cache(b, e);
	}
}

part_iter part_vect_count_lo(part_iter, part_iter, double xmid, int dim);

HPX_PLAIN_ACTION (part_vect_count_lo);

part_iter part_vect_count_lo(part_iter b, part_iter e, double xmid, int dim) {
//	printf( "part_vect_count_lo\n");
	const part_iter N = localities.size();
	const part_iter M = options::get().problem_size;
	part_iter count = 0;
	std::vector < hpx::future < part_iter >> futs;
	if (part_end < e) {
		for (int n = part_vect_locality_id(part_end); n * M / N < e; n++) {
			auto fut = hpx::async < part_vect_count_lo_action > (localities[n], n * M / N, std::min(e, (n + 1) * M / N), xmid, dim);
			futs.push_back(std::move(fut));

		}
	}
	auto this_begin = std::max(b, part_begin);
	auto this_end = std::min(e, part_end);
	const part_iter chunk_size = std::max(part_iter(1), (part_end - part_begin) / std::thread::hardware_concurrency());
	for (part_iter i = this_begin; i < this_end; i += chunk_size) {
		auto func = [i, this_end, dim, xmid, chunk_size]() {
			const auto end = std::min(i + chunk_size, this_end);
			part_iter count = 0;
			for (part_iter j = i; j < end; j++) {
				if (parts(j).x[dim] < xmid) {
					count++;
				}
			}
			return count;
		};
		if (chunk_size >= this_end - this_begin) {
			count += func();
		} else {
			futs.push_back(hpx::async(func));
		}
	}
	for (auto &fut : futs) {
		count += fut.get();
	}
	return count;
}

part_iter part_vect_sort(part_iter b, part_iter e, double xmid, int dim) {
	if (e == b) {
		return e;
	} else {
		if (e <= part_end) {
			auto lo = b;
			auto hi = e;
			while (lo < hi) {
				if (parts(lo).x[dim] >= xmid) {
					while (lo != hi) {
						hi--;
						if (parts(hi).x[dim] < xmid) {
							auto tmp = parts(lo);
							parts(lo) = parts(hi);
							parts(hi) = tmp;
							break;
						}
					}
				}
				lo++;
			}
			return hi;
		} else {
			auto count = part_vect_count_lo(b, e, xmid, dim);
			part_vect_sort_begin_action()(localities[part_vect_locality_id(b)], b, e, b + count, xmid, dim);
			return b + count;
		}
	}
}

range part_vect_range(part_iter b, part_iter e) {
	const auto id = part_vect_locality_id(b);
	range r;
	if (id == myid) {
		for (int dim = 0; dim < NDIM; dim++) {
			r.max[dim] = -std::numeric_limits<pos_type>::max();
			r.min[dim] = +std::numeric_limits<pos_type>::max();
		}
		const auto this_e = std::min(e, part_end);
		hpx::future<range> fut;
		if (this_e != e) {
			fut = hpx::async < part_vect_range_action > (localities[myid + 1], this_e, e);
		}
		for (part_iter i = b; i < this_e; i++) {
			const auto &p = parts(i);
			for (int dim = 0; dim < NDIM; dim++) {
				const auto x = p.x[dim];
				r.max[dim] = std::max(r.max[dim], x);
				r.min[dim] = std::min(r.min[dim], x);
			}
		}
		if (this_e != e) {
			const auto tmp = fut.get();
			for (int dim = 0; dim < NDIM; dim++) {
				r.max[dim] = std::max(r.max[dim], tmp.max[dim]);
				r.min[dim] = std::min(r.min[dim], tmp.max[dim]);
			}
		}
	} else {
		auto fut = hpx::async < part_vect_range_action > (localities[id], b, e);
		r = fut.get();
	}
	return r;
}

int part_vect_locality_id(part_iter i) {
	const part_iter N = localities.size();
	const part_iter M = options::get().problem_size;
	part_iter n = i * N / M;
	while (i < n * M / N) {
		n--;
	}
	while (i >= (n + 1) * M / N) {
		n++;
	}
	return std::min(n, N - 1);
}
