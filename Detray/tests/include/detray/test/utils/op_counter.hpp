// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace detray::test {

/// Struct for recording floating point operation counts
struct op_counts {
  long mul = 0, add = 0, sub = 0, neg = 0, div = 0;
  long flops() const { return mul + add + sub + neg + div; }
  friend op_counts operator+(op_counts a, op_counts b) {
    return {a.mul + b.mul, a.add + b.add, a.sub + b.sub, a.neg + b.neg,
            a.div + b.div};
  }
  op_counts &operator+=(op_counts o) { return *this = *this + o; }
};

/// Global object that is accessed by @c measure to provide accurate op
/// counting for a code block. Do not use in a multi-threaded context.
inline op_counts g_executed{};

/// A scalar type annotated with an operation count.
template <typename T>
struct counted {
  T v{};
  op_counts ops{};
  counted() = default;
  counted(const counted &) = default;
  counted(T x) : v(x) {}
  counted(T x, op_counts o) : v(x), ops(o) {}
  counted &operator=(const counted &o) {
    // Realize this store's cost once.
    g_executed = g_executed + o.ops;
    v = o.v;
    // Reset the count because the global state has been updated.
    ops = {};
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

/// Run a kernel with the op tally zeroed and return the ops it executed.
template <typename F>
op_counts measure(F &&run) {
  g_executed = op_counts{};
  run();
  return g_executed;
}

}  // namespace detray::test
