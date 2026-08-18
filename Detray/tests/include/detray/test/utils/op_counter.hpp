// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace detray::test {

// ---------------------------------------------------------------------------
//  op-count instrumentation
// ---------------------------------------------------------------------------
//
// A scalar that carries, alongside its value, the number of arithmetic
// operators that produced it. Feeding it through the library as the scalar
// type lets us tally the emitted op count. The count is purely structural
// (driven by the compile-time zero/one elisions), so the fill values are
// irrelevant: we measure default-constructed inputs.
struct op_counts {
  long mul = 0, add = 0, sub = 0, neg = 0, div = 0;
  long flops() const { return mul + add + sub + neg + div; }
  friend op_counts operator+(op_counts a, op_counts b) {
    return {a.mul + b.mul, a.add + b.add, a.sub + b.sub, a.neg + b.neg,
            a.div + b.div};
  }
  op_counts &operator+=(op_counts o) { return *this = *this + o; }
};

// Every intermediate value is materialized through exactly one
// counted::operator= (the `rv.at<I,J>() = ...` in fill_with, or a local `auto t
// = ...`). We harvest the producing expression's ops there -- once -- and zero
// the stored copy, so a value read by many consumers is counted once at its
// store, not once per reader. That makes the tally DAG-correct (each executed
// op counted once) instead of a per-output-cell expression-tree sum.
inline op_counts g_executed{};

template <typename T>
struct counted {
  T v{};
  op_counts ops{};
  counted() = default;
  counted(const counted &) = default;
  counted(T x) : v(x) {}  // constant/literal -> ops stays {}
  counted(T x, op_counts o) : v(x), ops(o) {}
  counted &operator=(const counted &o) {
    g_executed = g_executed + o.ops;  // realize this store's cost, once
    v = o.v;
    ops = {};  // already counted -- don't let consumers re-count it
    return *this;
  }
  friend counted operator*(counted a, counted b) {
    return {a.v * b.v, a.ops + b.ops + op_counts{.mul = 1}};
  }
  friend counted operator+(counted a, counted b) {
    return {a.v + b.v, a.ops + b.ops + op_counts{.add = 1}};
  }
  friend counted operator-(counted a, counted b) {
    return {a.v - b.v, a.ops + b.ops + op_counts{.sub = 1}};
  }
  friend counted operator/(counted a, counted b) {
    return {a.v / b.v, a.ops + b.ops + op_counts{.div = 1}};
  }
  counted operator-() const { return {-v, ops + op_counts{.neg = 1}}; }
};

// Run a kernel with the op tally zeroed and return the ops it executed. Only
// stored cells (symmetric results store just their upper triangle) and shared
// intermediates are counted, each exactly once -- see counted::operator=.
template <typename F>
op_counts measure(F &&run) {
  g_executed = op_counts{};
  run();
  return g_executed;
}

}  // namespace detray::test
