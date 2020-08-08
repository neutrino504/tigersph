#pragma once

#include <tigersph/gravity.hpp>
#include <tigersph/multipole.hpp>
#include <tigersph/output.hpp>
#include <tigersph/options.hpp>
#include <tigersph/part_vect.hpp>
#include <tigersph/range.hpp>

#include <array>
#include <atomic>
#include <memory>

#ifndef HPX_LITE
#include <hpx/synchronization/spinlock.hpp>
#include <hpx/include/async.hpp>
#include <hpx/include/threads.hpp>
#include <hpx/include/components.hpp>
#else
#include <hpx/hpx_lite.hpp>
#endif

class tree;

using mutex_type = hpx::lcos::local::spinlock;
template<class T>
using future_type = hpx::future<T>;

class node_attr;
class multi_src;
class check_item;
class multipole_return;

struct raw_id_type {
	int loc_id;
	std::uint64_t ptr;
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & loc_id;
		arc & ptr;
	}
	bool operator==(const raw_id_type &other) const {
		return ptr == other.ptr && loc_id == other.loc_id;
	}
};
using id_type = hpx::id_type;

class tree_client {
	id_type ptr;
public:
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & ptr;
	}
	tree_client() = default;
	tree_client(id_type ptr_);
	check_item get_check_item() const;
	multipole_return compute_multipoles(rung_type min_rung, bool do_out, int stack_cnt) const;
	void drift(float dt) const;
	kick_return kick_fmm(std::vector<check_item> checklist, const vect<ireal> &Lcom, expansion<float> L, rung_type min_rung, bool do_output,
			int stack_cnt) const;
	bool refine(int) const;
};

class raw_tree_client {
	raw_id_type ptr;
public:
	raw_tree_client() = default;
	raw_tree_client(raw_id_type ptr_);
	hpx::future<node_attr> get_node_attributes() const;
	hpx::future<multi_src> get_multi_srcs() const;
	int get_locality() const {
		return ptr.loc_id;
	}
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & ptr;
	}
};

struct check_item {
	bool opened;
	raw_tree_client node;
	float r;
	vect<float> x;
	part_iter pbegin;
	part_iter pend;
	bool is_leaf;
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & opened;
		arc & node;
		arc & r;
		arc & x;
		arc & pbegin;
		arc & pend;
		arc & is_leaf;
	}

};

struct node_attr {
	std::array<check_item, NCHILD> children;
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & children;
	}
};

struct multipole_return {
	multipole_info m;
	range r;
	check_item c;
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & m;
		arc & r;
		arc & c;
	}
};

class tree: public hpx::components::managed_component_base<tree> {
	multipole_info multi;
	part_iter part_begin;
	part_iter part_end;
	int level;
	range box;
	bool leaf;
	std::array<tree_client, NCHILD> children;
	std::array<check_item, NCHILD> child_check;

	static float theta_inv;
	static std::atomic<std::uint64_t> flop;

public:
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & multi;
		arc & part_begin;
		arc & part_end;
		arc & children;
		arc & child_check;
		arc & level;
		arc & box;
		arc & leaf;
	}
	static std::uint64_t get_flop();
	tree() = default;
	static void set_theta(float);
	static void reset_flop();
	tree(range, part_iter, part_iter, int level);
	bool refine(int);
	bool is_leaf() const;
	multipole_return compute_multipoles(rung_type min_rung, bool do_out, int stack_cnt);
	node_attr get_node_attributes() const;
	multi_src get_multi_srcs() const;
	void drift(float);
	kick_return kick_fmm(std::vector<check_item> checklist, const vect<ireal> &Lcom, expansion<float> L, rung_type min_rung, bool do_output, int stack_ccnt);
	check_item get_check_item() const; //
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,refine);//
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,compute_multipoles);//
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,kick_fmm);//
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,drift);//
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_check_item);
	//
};

inline tree_client::tree_client(id_type ptr_) {
	ptr = ptr_;
}

inline raw_tree_client::raw_tree_client(raw_id_type ptr_) {
	ptr = ptr_;
}

inline check_item tree_client::get_check_item() const {
	return tree::get_check_item_action()(ptr);
}

inline multipole_return tree_client::compute_multipoles(rung_type min_rung, bool do_out, int stack_cnt) const {
	return tree::compute_multipoles_action()(ptr, min_rung, do_out, stack_cnt);
}

inline void tree_client::drift(float dt) const {
	return tree::drift_action()(ptr, dt);
}

inline bool tree_client::refine(int stack_cnt) const {
	return tree::refine_action()(ptr, stack_cnt);
}

inline kick_return tree_client::kick_fmm(std::vector<check_item> checklist, const vect<ireal> &Lcom, expansion<float> L, rung_type min_rung, bool do_output,
		int stack_cnt) const {
	return tree::kick_fmm_action()(ptr, std::move(checklist), Lcom, L, min_rung, do_output, stack_cnt);
}
