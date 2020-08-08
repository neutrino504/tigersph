/*
 * mutex.hpp
 *
 *  Created on: Dec 11, 2015
 *      Author: dmarce1
 */

#ifndef HPX_MUTEX_HPP_
#define HPX_MUTEX_HPP_

#include <atomic>
#include <mutex>

namespace hpx {


class mutex {
private:
	std::atomic<int> locked;
public:
	constexpr mutex() :
			locked(0) {
	}
	mutex(const mutex&) = delete;
	mutex(mutex&&) = delete;
	~mutex() = default;
	mutex& operator=(const mutex&) = delete;
	mutex& operator=(mutex&&) = delete;
	void lock();
	void unlock();
};

class spinlock {
private:
	std::atomic<int> locked;
public:
	constexpr spinlock() :
			locked(0) {
	}
	spinlock(const spinlock&) = delete;
	spinlock(spinlock&&) = delete;
	~spinlock() = default;
	spinlock& operator=(const spinlock&) = delete;
	spinlock& operator=(spinlock&&) = delete;
	void lock();
	void unlock();
};

namespace lcos {
namespace local {
using spinlock = hpx::spinlock;
using mutex =hpx::mutex;
}
}
}

#include "hpx/thread.hpp"

namespace hpx {

inline void mutex::lock() {
	while (locked++ != 0) {
		locked--;
		hpx::this_thread::yield();
	}
}

inline void mutex::unlock() {
	locked--;
}

inline void spinlock::lock() {
	while (locked++ != 0) {
		locked--;
		hpx::this_thread::yield();
	}
}

inline void spinlock::unlock() {
	locked--;
}

}

#endif /* MUTEX_HPP_ */
