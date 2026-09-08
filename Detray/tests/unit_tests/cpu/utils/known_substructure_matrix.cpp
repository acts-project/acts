// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s).
#include "detray/algebra/common/known_substructure_matrix.hpp"

#include "detray/test/utils/op_counter.hpp"

// Google Test include(s).
#include <gtest/gtest.h>

// System include(s)
#include <array>
#include <cstddef>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

namespace {

using namespace detray::ksm;
using detray::test::counted;
using detray::test::op_counts;

template <typename A, typename B>
concept addable = (A::rows == B::rows) && (A::columns == B::columns);

template <typename A, typename B>
concept multipliable = (A::columns == B::rows);

template <typename A, typename B>
concept congruable = concepts::is_symmetric<typename B::canonical_type> &&
                     (A::columns == B::rows);

// Count the structural zeros (cells whose value type is `zero`) of a
// substructure. Inversion permutes the cofactor pattern, so the inverse must
// carry exactly as many structural zeros as its input -- see the inverse block
// in test_single.
template <typename Sub, std::size_t... Ks>
constexpr std::size_t count_structural_zeros(std::index_sequence<Ks...>) {
  constexpr std::size_t C = Sub::columns;
  return (std::size_t{0} + ... +
          (std::is_same_v<typename Sub::template value_at<Ks / C, Ks % C>, zero>
               ? std::size_t{1}
               : std::size_t{0}));
}

// A matrix whose every cell is the same integer constant. Only used to test
// that a uniform-integer add/subtract preserves the variable count (see
// test_single); not a generally useful shape, so it lives here rather than in
// the header's builders.
template <int N>
struct const_cell {
  template <std::size_t I, std::size_t J>
  struct cell {
    using type = integral_value<N>;
  };
};
template <std::size_t Rows, std::size_t Cols, int N>
using make_const_substructure =
    make_substructure<const_cell<N>::template cell, Rows, Cols>;

// One operation's structure-aware cost paired with its structure-blind (dense)
// baseline, which the dense-flop invariant is asserted against.
struct op_stat {
  op_counts ops{};    // structure-aware
  op_counts naive{};  // dense baseline of the same shape
};

// The structure-blind substructure of the same shape as S.
template <typename S>
using dense_of = make_dense_substructure<S::rows, S::columns>;

// Build a counted matrix for each substructure in Subs, apply op to them, and
// tally the arithmetic ops the op executed (harvested at the result's stores by
// counted::operator=, so shared intermediates are billed once -- see
// op_counter.hpp).
template <typename... Subs, typename Op>
op_counts measure(Op op) {
  return detray::test::measure(
      [&] { (void)op(matrix<Subs, counted<double>>{}...); });
}

// The same op run on the structured Subs and on their dense equivalents: the
// structure-aware count paired with the structure-blind baseline.
template <typename... Subs, typename Op>
op_stat measure_vs_dense(Op op) {
  return {measure<Subs...>(op), measure<dense_of<Subs>...>(op)};
}

// --- runtime checking ------------------------------------------------------

void check(bool cond, const std::string &name) {
  EXPECT_TRUE(cond) << name;
}

// Hard invariant on the structure-blind baseline: a fully dense computation
// must cost exactly the textbook flop count -- M*N for add/sub, M*P*(2N-1) for
// an M*N by N*P multiply. If the dense path disagrees, the evaluator is doing
// redundant work even with no structure to exploit, so raise an error.
void require_dense_flops(long got, long expected, const std::string &what) {
  EXPECT_EQ(got, expected) << what
                           << ": dense flop count must be the textbook "
                              "count with no structure to exploit";
}

template <typename Sub, std::size_t I, std::size_t J, typename S>
void set_cell(matrix<Sub, S> &m, S v) {
  using value_t = typename Sub::template value_at<I, J>;
  if constexpr (std::is_same_v<value_t, variable>) {
    m.template at<I, J>() = S{v + static_cast<S>((I + 1) * 412 + J)};
  } else if constexpr (concepts::is_index_variable<value_t>) {
    m.template at<I, J>() = S{v + static_cast<S>(value_t::index)};
  }
}

template <typename Sub, std::size_t I, std::size_t J, typename S>
void set_truth_cell(std::array<std::array<S, Sub::columns>, Sub::rows> &m,
                    S v) {
  using value_t = typename Sub::template value_at<I, J>;
  if constexpr (std::is_same_v<value_t, variable>) {
    m.at(I).at(J) = S{v + static_cast<S>((I + 1) * 412 + J)};
  } else if constexpr (concepts::is_index_variable<value_t>) {
    m.at(I).at(J) = S{v + static_cast<S>(value_t::index)};
  } else {
    m.at(I).at(J) = S{static_cast<S>(value_t::value)};
  }
}

template <typename Sub, typename S, std::size_t... Ks>
void fill_seq(matrix<Sub, S> &m, S base, std::index_sequence<Ks...>) {
  (set_cell<Sub, Ks / Sub::columns, Ks % Sub::columns>(m, base), ...);
}

template <typename Sub, typename S, std::size_t... Ks>
void fill_truth_seq(std::array<std::array<S, Sub::columns>, Sub::rows> &m,
                    S base, std::index_sequence<Ks...>) {
  (set_truth_cell<Sub, Ks / Sub::columns, Ks % Sub::columns>(m, base), ...);
}

template <typename S>
bool check_eq(std::size_t i, std::size_t j, S a, S b) {
  if (a == b) {
    return true;
  } else {
    ADD_FAILURE() << "Different values at (" << i << ", " << j << "): " << a
                  << " vs " << b;
    return false;
  }
}

template <typename MA, std::size_t... Ks>
bool verify_truth_values(
    const MA &a,
    const std::array<std::array<typename MA::scalar_type, MA::columns>,
                     MA::rows> &a_t,
    std::index_sequence<Ks...>) {
  constexpr std::size_t C = MA::substructure_type::columns;

  return (... && check_eq(Ks / C, Ks % C, a_t.at(Ks / C).at(Ks % C),
                          a.template at<Ks / C, Ks % C>()));
}

template <typename MA, typename MB, typename MR, std::size_t... Ks>
bool verify_addition_values(
    [[maybe_unused]] const MA &a, [[maybe_unused]] const MB &b, const MR &r,
    const std::array<std::array<typename MA::scalar_type, MA::columns>,
                     MA::rows> &a_t,
    const std::array<std::array<typename MB::scalar_type, MB::columns>,
                     MB::rows> &b_t,
    std::index_sequence<Ks...>) {
  static_assert(MR::substructure_type::columns ==
                MA::substructure_type::columns);
  static_assert(MR::substructure_type::columns ==
                MB::substructure_type::columns);
  static_assert(MR::substructure_type::rows == MA::substructure_type::rows);
  static_assert(MR::substructure_type::rows == MB::substructure_type::rows);
  constexpr std::size_t C = MR::substructure_type::columns;
  std::array<std::array<typename MA::scalar_type, MA::columns>, MA::rows> r_t{};

  ((r_t.at(Ks / C).at(Ks % C) =
        a_t.at(Ks / C).at(Ks % C) + b_t.at(Ks / C).at(Ks % C)),
   ...);

  return (... && check_eq(Ks / C, Ks % C, r_t.at(Ks / C).at(Ks % C),
                          r.template at<Ks / C, Ks % C>()));
}

template <typename MA, typename MB, typename MR, std::size_t... Ks>
bool verify_subtraction_values(const MA &a, const MB &b, const MR &r,
                               std::index_sequence<Ks...>) {
  constexpr std::size_t C = MR::substructure_type::columns;
  return (... &&
          (r.template at<Ks / C, Ks % C>() ==
           a.template at<Ks / C, Ks % C>() - b.template at<Ks / C, Ks % C>()));
}

template <typename MA, typename MR, std::size_t... Ks>
bool verify_negation_values(const MA &a, const MR &r,
                            std::index_sequence<Ks...>) {
  constexpr std::size_t C = MR::substructure_type::columns;
  return (... && (r.template at<Ks / C, Ks % C>() ==
                  -a.template at<Ks / C, Ks % C>()));
}

// transpose(a)(j, i) must read back the same scalar as a(i, j). The index walk
// is over a's M*N cells; the transposed matrix t is N*M, so cell (i, j) of a is
// compared against cell (j, i) of t.
template <typename MA, typename MT, std::size_t... Ks>
bool verify_transpose_values(const MA &a, const MT &t,
                             std::index_sequence<Ks...>) {
  constexpr std::size_t C = MA::substructure_type::columns;
  return (... && check_eq(Ks / C, Ks % C, a.template at<Ks / C, Ks % C>(),
                          t.template at<Ks % C, Ks / C>()));
}

// transpose is an involution: transposing twice returns the original cell by
// cell. tt has the same shape as a, so the comparison is position for position.
template <typename MA, typename MTT, std::size_t... Ks>
bool verify_involution_values(const MA &a, const MTT &tt,
                              std::index_sequence<Ks...>) {
  constexpr std::size_t C = MA::substructure_type::columns;
  return (... && check_eq(Ks / C, Ks % C, a.template at<Ks / C, Ks % C>(),
                          tt.template at<Ks / C, Ks % C>()));
}

// Two same-shape matrices agree cell by cell. Used to check a structure-aware
// result against an independent reference matrix of the same dimensions.
template <typename MA, typename MB, std::size_t... Ks>
bool verify_equal_matrices(const MA &a, const MB &b,
                           std::index_sequence<Ks...>) {
  static_assert(MA::rows == MB::rows && MA::columns == MB::columns);
  constexpr std::size_t C = MA::substructure_type::columns;
  return (... && check_eq(Ks / C, Ks % C, a.template at<Ks / C, Ks % C>(),
                          b.template at<Ks / C, Ks % C>()));
}

// Approximate cell comparison. Inverse introduces division, so its results are
// no longer exact even for integer-valued fills; an absolute epsilon suffices
// here because the only inverse value check is against the identity (0s and
// 1s).
bool check_close(std::size_t i, std::size_t j, double a, double b, double eps) {
  double d = a - b;
  if (d < 0)
    d = -d;
  if (d <= eps)
    return true;
  ADD_FAILURE() << "Not close at (" << i << ", " << j << "): " << a << " vs "
                << b;
  return false;
}

// A square product matrix equals the identity: 1 on the diagonal, 0 off it.
template <typename MP, std::size_t... Ks>
bool verify_identity(const MP &p, std::index_sequence<Ks...>) {
  constexpr std::size_t C = MP::substructure_type::columns;
  constexpr double eps = 1e-9;
  return (... && check_close(Ks / C, Ks % C, p.template at<Ks / C, Ks % C>(),
                             (Ks / C == Ks % C) ? 1.0 : 0.0, eps));
}

// Reference (fold-based) dot product of row I of `a` with column J of `b`.
template <std::size_t I, std::size_t J, typename S, typename MA, typename MB,
          std::size_t... Ks>
S ref_dot(const MA &a, const MB &b, std::index_sequence<Ks...>) {
  return (S{0} + ... + (a.template at<I, Ks>() * b.template at<Ks, J>()));
}

template <typename S, typename MA, typename MB, typename MR, std::size_t... Ks>
bool verify_multiplication_values(const MA &a, const MB &b, const MR &r,
                                  std::index_sequence<Ks...>) {
  constexpr std::size_t C = MR::substructure_type::columns;
  constexpr std::size_t K = MA::substructure_type::columns;
  return (... &&
          (r.template at<Ks / C, Ks % C>() ==
           ref_dot<Ks / C, Ks % C, S>(a, b, std::make_index_sequence<K>{})));
}

// --- per-pair driver -------------------------------------------------------

template <typename A, typename B>
void test_pair(const std::string &an, const std::string &bn) {
  using S = double;  // integer-valued fills keep == comparisons exact
  const std::string pair = an + ", " + bn;
  SCOPED_TRACE(pair);

  [[maybe_unused]] constexpr auto a_indices =
      std::make_index_sequence<A::rows * A::columns>{};
  [[maybe_unused]] constexpr auto b_indices =
      std::make_index_sequence<B::rows * B::columns>{};

  // Each matrix must store exactly one scalar per free value of A and of B.
  static_assert(verify_stored_count<typename A::canonical_type>());
  static_assert(verify_stored_count<typename B::canonical_type>());

  matrix<A, S> ma;
  matrix<B, S> mb;
  fill_seq(ma, S{1}, a_indices);
  fill_seq(mb, S{1000}, b_indices);
  std::array<std::array<S, A::columns>, A::rows> A_truth{};
  std::array<std::array<S, B::columns>, B::rows> B_truth{};
  fill_truth_seq<A>(A_truth, S{1}, a_indices);
  fill_truth_seq<B>(B_truth, S{1000}, b_indices);

  bool vta = verify_truth_values(ma, A_truth, a_indices);
  check(vta, an + ": A");
  bool vtb = verify_truth_values(mb, B_truth, b_indices);
  check(vtb, bn + ": B");

  if constexpr (addable<A, B>) {
    using ADD = typename A::canonical_type::template addition_type<
        typename B::canonical_type>;
    using SUB = typename A::canonical_type::template subtraction_type<
        typename B::canonical_type>;
    static_assert(ADD::rows == A::rows && ADD::columns == A::columns,
                  "addition result has wrong shape");
    static_assert(SUB::rows == A::rows && SUB::columns == A::columns,
                  "subtraction result has wrong shape");

    auto add = [](auto a, auto b) { return a + b; };
    auto sub = [](auto a, auto b) { return a - b; };
    op_stat as = measure_vs_dense<A, B>(add);
    op_stat ss = measure_vs_dense<A, B>(sub);

    const long elems =
        static_cast<long>(A::rows) * static_cast<long>(A::columns);
    require_dense_flops(as.naive.flops(), elems, pair + ": A + B dense");
    require_dense_flops(ss.naive.flops(), elems, pair + ": A - B dense");

    bool va =
        verify_addition_values(ma, mb, ma + mb, A_truth, B_truth, a_indices);
    check(va, pair + ": A + B (ops=" + std::to_string(as.ops.flops()) +
                  ", naive=" + std::to_string(as.naive.flops()) + ")");
    bool vb = verify_subtraction_values(ma, mb, ma - mb, a_indices);
    check(vb, pair + ": A - B (ops=" + std::to_string(ss.ops.flops()) +
                  ", naive=" + std::to_string(ss.naive.flops()) + ")");
  }

  if constexpr (multipliable<A, B>) {
    using MUL = typename A::canonical_type::template multiplication_type<
        typename B::canonical_type>;
    static_assert(MUL::rows == A::rows && MUL::columns == B::columns,
                  "multiplication result has wrong shape");
    constexpr auto r_indices =
        std::make_index_sequence<MUL::rows * MUL::columns>{};

    auto mul = [](auto a, auto b) { return a * b; };
    op_stat ms = measure_vs_dense<A, B>(mul);

    const long M = static_cast<long>(A::rows),
               N = static_cast<long>(A::columns),
               P = static_cast<long>(B::columns);
    require_dense_flops(ms.naive.flops(), M * P * (2 * N - 1),
                        pair + ": A * B dense");

    bool vm = verify_multiplication_values<S>(ma, mb, ma * mb, r_indices);
    check(vm, pair + ": A * B (add=" + std::to_string(ms.ops.add) +
                  ", mul=" + std::to_string(ms.ops.mul) +
                  ", ops=" + std::to_string(ms.ops.flops()) +
                  ", naive=" + std::to_string(ms.naive.flops()) + ")");
  }

  // Congruence A*B*A^T: only when B is symmetric (so the result is symmetric).
  if constexpr (congruable<A, B>) {
    using CONG = typename A::canonical_type::template congruence_type<
        typename B::canonical_type>;
    static_assert(CONG::rows == A::rows && CONG::columns == A::rows,
                  "congruence result must be A.rows x A.rows");
    constexpr auto c_indices =
        std::make_index_sequence<CONG::rows * CONG::columns>{};

    // Oracle: the plain triple product from already-tested ops. mb is
    // numerically symmetric (its mirror cells share storage), so this is the
    // genuine A*B*A^T.
    auto cong = ma.congruence(mb);
    auto naive = (ma * mb) * ma.transpose();

    check(verify_equal_matrices(cong, naive, c_indices),
          pair + ": congruence == A*B*transpose(A)");

    // Flop pin against the naive implementation of the same product. Both paths
    // form A*B and transpose(A); congruence() then fills only the stored cells
    // of a result *declared* symmetric, where the plain triple product computes
    // both mirror halves as independent dot products. Computing a subset of the
    // same dot products can never cost more, so this is derived rather than a
    // snapshot -- no constant to bump when the algebra improves. It is <=, not
    // <: when the result's symmetry buys nothing (a 1x1 congruence, or a shape
    // whose mirror cells were already shared) the two coincide exactly.
    auto congruence_op = [](auto a, auto b) { return a.congruence(b); };
    auto triple = [](auto a, auto b) { return (a * b) * a.transpose(); };
    const long cong_flops = measure<A, B>(congruence_op).flops();
    const long triple_flops = measure<A, B>(triple).flops();
    check(cong_flops <= triple_flops,
          pair + ": congruence costs <= naive A*B*transpose(A) (" +
              std::to_string(cong_flops) +
              " <= " + std::to_string(triple_flops) + ")");

    // Declared structure: A*B*A^T is symmetric for symmetric B, the property a
    // plain triple product drops (its mirror cells are distinct dot products).
    check(concepts::is_symmetric<typename CONG::canonical_type>,
          pair + ": congruence(A, symmetric B) is symmetric");

    // Storage payoff: a symmetric m*m result stores at most the upper triangle,
    // m(m+1)/2, versus up to m*m for the structure-blind product.
    constexpr std::size_t m = CONG::rows;
    constexpr std::size_t tri = m * (m + 1) / 2;
    check(CONG::num_variables <= tri, pair +
                                          ": congruence stores <= m(m+1)/2 (" +
                                          std::to_string(CONG::num_variables) +
                                          " <= " + std::to_string(tri) + ")");
  }
}

// --- per-matrix driver -----------------------------------------------------
//
// Operations that act on a single matrix (currently negation) don't fit the
// pair driver, so they get their own pass over the substructure list. Mirrors
// test_pair: compile-time shape + stored-count checks, an op-count vs
// dense-baseline measurement, and a runtime value check.

template <typename A>
void test_single(const std::string &an) {
  using S = double;
  SCOPED_TRACE(an);

  constexpr auto a_indices = std::make_index_sequence<A::rows * A::columns>{};

  matrix<A, S> ma;
  fill_seq(ma, S{1}, a_indices);
  std::array<std::array<S, A::columns>, A::rows> A_truth{};
  fill_truth_seq<A>(A_truth, S{1}, a_indices);
  check(verify_truth_values(ma, A_truth, a_indices), an + ": A");

  using NEG =
      typename detail::negated_substructure<typename A::canonical_type>::type;
  static_assert(NEG::rows == A::rows && NEG::columns == A::columns,
                "negation result has wrong shape");

  auto neg = [](auto a) { return -a; };
  op_stat ns = measure_vs_dense<A>(neg);

  // Negating an M*N matrix is one negation per element.
  const long elems = static_cast<long>(A::rows) * static_cast<long>(A::columns);
  require_dense_flops(ns.naive.flops(), elems, an + ": -A dense");

  bool vn = verify_negation_values(ma, -ma, a_indices);
  check(vn, an + ": -A (ops=" + std::to_string(ns.ops.flops()) +
                ", naive=" + std::to_string(ns.naive.flops()) + ")");

  // Structure invariant: negation is a per-cell sign flip, so it preserves
  // every equality between cells. -A must therefore store exactly as many
  // values as A. (Runtime check, not static_assert: the current algebra
  // collapses negated index_variables to anonymous variables and violates
  // this for any matrix with shared variables.)
  constexpr std::size_t a_vars = A::canonical_type::num_variables;
  constexpr std::size_t neg_vars = NEG::num_variables;
  check(neg_vars == a_vars, an + ": -A preserves stored count (" +
                                std::to_string(a_vars) + " -> " +
                                std::to_string(neg_vars) + ")");

  // Structure invariant: if A is symmetric then A + A is symmetric with the
  // same sharing, so A + A must store exactly as many values as A. Only checked
  // for symmetric A; same caveat as above about the additive algebra.
  if constexpr (concepts::is_symmetric<typename A::canonical_type>) {
    using APA = typename A::canonical_type::template addition_type<
        typename A::canonical_type>;
    constexpr std::size_t apa_vars = APA::num_variables;
    check(apa_vars == a_vars,
          an + ": A+A preserves stored count when A symmetric (" +
              std::to_string(a_vars) + " -> " + std::to_string(apa_vars) + ")");
  }

  // Structure invariant: adding or subtracting a matrix whose every cell is the
  // *same* integer constant shifts values but touches no structure -- a
  // variable cell stays that same variable (+/- c), a structural zero/constant
  // stays constant, and equal cells stay equal -- so the stored-variable count
  // is unchanged. Both conditions are load-bearing: a *non-uniform* constant
  // would split a shared variable into distinct value_(add|sub)<var, c> slots,
  // and an all-*variable* matrix would turn structural zeros into stored
  // variables.
  {
    using ConstM = typename make_const_substructure<A::rows, A::columns,
                                                    7>::canonical_type;
    using APC = typename A::canonical_type::template addition_type<ConstM>;
    using AMC = typename A::canonical_type::template subtraction_type<ConstM>;
    constexpr std::size_t apc_vars = APC::num_variables;
    constexpr std::size_t amc_vars = AMC::num_variables;
    check(apc_vars == a_vars, an + ": A + const preserves stored count (" +
                                  std::to_string(a_vars) + " -> " +
                                  std::to_string(apc_vars) + ")");
    check(amc_vars == a_vars, an + ": A - const preserves stored count (" +
                                  std::to_string(a_vars) + " -> " +
                                  std::to_string(amc_vars) + ")");
  }

  // -- transpose: the other single-matrix operation -------------------------
  // Shape swaps M*N -> N*M and each cell mirrors the input. Since transpose
  // only permutes cells it preserves the stored-value count and symmetry.
  using TR = typename detail::transposed_substructure<
      typename A::canonical_type>::type;
  static_assert(TR::rows == A::columns && TR::columns == A::rows,
                "transpose result has wrong shape");

  auto mt = ma.transpose();
  using MT = std::decay_t<decltype(mt)>;
  static_assert(MT::rows == A::columns && MT::columns == A::rows,
                "transpose() matrix has wrong shape");

  check(verify_transpose_values(ma, mt, a_indices),
        an + ": transpose(A)(j,i) == A(i,j)");
  check(verify_involution_values(ma, mt.transpose(), a_indices),
        an + ": transpose(transpose(A)) == A");

  // Transpose permutes cells, so it stores exactly as many values as A.
  constexpr std::size_t t_vars = TR::num_variables;
  check(t_vars == a_vars, an + ": transpose preserves stored count (" +
                              std::to_string(a_vars) + " -> " +
                              std::to_string(t_vars) + ")");

  // The transpose of a symmetric matrix is symmetric (in fact equal to the
  // original, which the value check above already pins).
  if constexpr (concepts::is_symmetric<typename A::canonical_type>) {
    check(concepts::is_symmetric<typename TR::canonical_type>,
          an + ": transpose(symmetric) stays symmetric");
  }

  // -- inverse: closed-form for small (1x1, 2x2) invertible matrices --------
  // Gated on the library's is_invertible: singularity is semantic, not a shape
  // predicate, so structurally-singular matrices (a det that folds to a
  // compile-time zero) are excluded here rather than divided by. The defining
  // property is S*inverse(S) == I; the inverse of a symmetric matrix is
  // symmetric, reusing the upper-triangle storage; and 1/det introduces a
  // division.
  if constexpr ((A::rows == A::columns) && (A::rows <= 2) &&
                concepts::is_invertible<typename A::canonical_type>) {
    using CA = typename A::canonical_type;
    using INV = typename detail::inverted_substructure<
        typename CA::canonical_type>::type;
    static_assert(INV::rows == CA::rows && INV::columns == CA::columns,
                  "inverse has wrong shape");

    // Defining property: S * inverse(S) == I.
    check(verify_identity(ma * ma.inverse(), a_indices),
          an + ": S * inverse(S) == I");

    // The structure-aware inverse, paired with its naive baseline: the same
    // inverse of a fully-dense matrix (no zeros, ones, or shared variables to
    // exploit), which it saves against.
    auto inv = [](auto a) { return a.inverse(); };
    op_stat is = measure_vs_dense<A>(inv);

    // 1/det normally introduces a division -- except when det folds to 1, where
    // the inverse is just the adjugate and no division is performed. Pin both.
    if constexpr (std::same_as<
                      typename detail::get_determinant_helper<CA>::type, one>) {
      check(is.ops.div == 0, an + ": inverse(det==1) uses no division (div=" +
                                 std::to_string(is.ops.div) + ")");
    } else {
      check(is.ops.div >= 1, an + ": inverse uses division (div=" +
                                 std::to_string(is.ops.div) + ")");
    }

    // For symmetric input the inverse is symmetric, so it stores only the upper
    // triangle (m(m+1)/2), not m*m.
    if constexpr (concepts::is_symmetric<CA>) {
      check(concepts::is_symmetric<typename INV::canonical_type>,
            an + ": inverse(symmetric) is symmetric");
      constexpr std::size_t m = CA::rows;
      constexpr std::size_t tri = m * (m + 1) / 2;
      check(INV::num_variables <= tri, an + ": inverse stores <= m(m+1)/2 (" +
                                           std::to_string(INV::num_variables) +
                                           " <= " + std::to_string(tri) + ")");
    }

    // Structure preservation: inversion permutes the cofactor pattern, so the
    // inverse keeps exactly as many structural zeros as the input -- a diagonal
    // matrix inverts to a diagonal one. (When value_division was fed the
    // cofactor *metafunctions* instead of their ::type, the operands never
    // matched zero/one, every cell fell back to a dense variable, and the
    // structural zeros were silently lost.)
    constexpr std::size_t a_zeros = count_structural_zeros<CA>(
        std::make_index_sequence<CA::rows * CA::columns>{});
    constexpr std::size_t inv_zeros =
        count_structural_zeros<typename INV::canonical_type>(
            std::make_index_sequence<INV::rows * INV::columns>{});
    check(inv_zeros == a_zeros,
          an + ": inverse preserves structural-zero count (" +
              std::to_string(a_zeros) + " -> " + std::to_string(inv_zeros) +
              ")");
  }
}

// --- direct dot-product flop-count checks ----------------------------------
//
// The cartesian harness measures savings against the dense ceiling but has no
// "minimal" oracle, so these pin the *exact* op count for dot products that
// mix constant and variable terms. The constant terms must fold to a single
// compile-time offset, so however many constants a dot product contains they
// cost at most one runtime add. Each matrix below is all-zero except row 0 of
// A and column 0 of B, so the only non-trivial result cell is (0,0) and the
// whole result's op count is that single dot product.

// Each pair's (0,0) dot product realises one constant/variable evaluation-tree
// shape; everything else is zero so (0,0) is the only non-trivial result cell.
// These are reused below as S32-S37 of the cartesian list.

// const + const + var :  row [2, 3, var] . col [4, 5, var]
using cfd_ccv_a =
    substructure<row<integral_value<2>, integral_value<3>, variable>,
                 row<zero, zero, zero>, row<zero, zero, zero>>;
using cfd_ccv_b =
    substructure<row<integral_value<4>, zero, zero>,
                 row<integral_value<5>, zero, zero>, row<variable, zero, zero>>;
// const + var + const :  row [2, var, 3] . col [4, var, 5]
using cfd_cvc_a =
    substructure<row<integral_value<2>, variable, integral_value<3>>,
                 row<zero, zero, zero>, row<zero, zero, zero>>;
using cfd_cvc_b =
    substructure<row<integral_value<4>, zero, zero>, row<variable, zero, zero>,
                 row<integral_value<5>, zero, zero>>;
// var + const + const + var :  row [var, 2, 3, var] . col [var, 4, 5, var]
using cfd_vccv_a =
    substructure<row<variable, integral_value<2>, integral_value<3>, variable>,
                 row<zero, zero, zero, zero>, row<zero, zero, zero, zero>,
                 row<zero, zero, zero, zero>>;
using cfd_vccv_b = substructure<
    row<variable, zero, zero, zero>, row<integral_value<4>, zero, zero, zero>,
    row<integral_value<5>, zero, zero, zero>, row<variable, zero, zero, zero>>;

template <typename A, typename B>
void check_dot_ops(const std::string &name, long exp_mul, long exp_add) {
  using S = double;
  SCOPED_TRACE(name);
  constexpr auto a_indices = std::make_index_sequence<A::rows * A::columns>{};
  constexpr auto b_indices = std::make_index_sequence<B::rows * B::columns>{};
  constexpr auto r_indices = std::make_index_sequence<A::rows * B::columns>{};

  // Values must stay correct regardless of how the dot product is grouped.
  matrix<A, S> ma;
  matrix<B, S> mb;
  fill_seq(ma, S{1}, a_indices);
  fill_seq(mb, S{1000}, b_indices);
  check(verify_multiplication_values<S>(ma, mb, ma * mb, r_indices),
        name + " values");

  // Op count must hit the structural minimum (constants folded to one offset).
  op_counts mo = detray::test::measure([] {
    (void)(matrix<A, counted<double>>{} * matrix<B, counted<double>>{});
  });
  check(mo.mul == exp_mul && mo.add == exp_add,
        name + " ops (mul=" + std::to_string(mo.mul) + " exp " +
            std::to_string(exp_mul) + ", add=" + std::to_string(mo.add) +
            " exp " + std::to_string(exp_add) + ")");
}

/// A structural cell of a result has no storage, so nothing may be executed to
/// produce it.
///
/// Only a class-type scalar can show this. For a built-in scalar, assigning to
/// the prvalue the const accessor returns is ill formed, so a structural cell
/// is skipped whatever the fill gate asks. For a class type the same
/// assignment is well formed, because it is a call to @c operator=, so a gate
/// that asks whether the assignment compiles -- rather than whether the cell
/// is stored -- evaluates one cell per distinct structural constant and throws
/// the result away. That is invisible in the values and shows up only here.
///
/// @c I - B is the smallest case that exposes it. Where @p B leaves a column
/// structurally zero, that column of the result stays the identity's own
/// constants, and every cell the result does store costs exactly one operation
/// -- a subtraction on the diagonal, a negation off it. So the executed count
/// must equal the stored count, for every @p K.
template <std::size_t N, std::size_t K>
void check_no_structural_evaluation(const std::string &name) {
  SCOPED_TRACE(name);
  using lhs_t = typename make_identity_substructure<N>::canonical_type;
  using rhs_t =
      typename make_left_columns_substructure<N, N, K>::canonical_type;
  using result_t = typename lhs_t::template subtraction_type<rhs_t>;

  auto sub = [](auto a, auto b) { return a - b; };
  const long executed = measure<lhs_t, rhs_t>(sub).flops();
  const long stored = static_cast<long>(result_t::num_variables);

  check(executed == stored, name + ": executed " + std::to_string(executed) +
                                " ops for " + std::to_string(stored) +
                                " stored cells");
}

// --- the substructure list -------------------------------------------------

using substructure_tuple = std::tuple<

    substructure<row<one, variable>, row<zero, variable>>,
    substructure<row<variable, variable>, row<variable, one>>,
    substructure<row<integral_value<2>, variable>, row<one, integral_value<3>>>,
    substructure<row<one, zero, variable>, row<variable, variable, zero>>,
    substructure<row<variable, one>, row<zero, variable>, row<one, variable>>,
    substructure<row<one, zero, zero>, row<zero, one, zero>,
                 row<variable, variable, one>>,
    substructure<row<variable, integral_value<-1>, zero>,
                 row<zero, variable, variable>, row<one, zero, variable>>,
    substructure<row<one, zero, zero, zero, zero, zero, zero, variable>,
                 row<zero, one, zero, zero, zero, zero, variable, zero>,
                 row<zero, zero, one, zero, zero, variable, zero, zero>,
                 row<zero, zero, zero, one, variable, zero, zero, zero>,
                 row<zero, zero, zero, variable, one, zero, zero, zero>,
                 row<zero, zero, variable, zero, zero, one, zero, zero>,
                 row<zero, variable, zero, zero, zero, zero, one, zero>,
                 row<variable, zero, zero, zero, zero, zero, zero, one>>,
    substructure<row<variable, variable, zero, zero, zero, zero, zero, one>,
                 row<zero, variable, variable, zero, zero, zero, one, zero>,
                 row<zero, zero, variable, variable, zero, one, zero, zero>,
                 row<zero, zero, zero, variable, variable, zero, zero, zero>,
                 row<one, zero, zero, zero, variable, variable, zero, zero>>,
    substructure<row<variable, zero, zero, zero, one>,
                 row<zero, variable, zero, one, zero>,
                 row<zero, zero, variable, zero, zero>,
                 row<zero, one, zero, variable, zero>,
                 row<one, zero, zero, zero, variable>,
                 row<variable, zero, one, zero, zero>,
                 row<zero, variable, zero, zero, one>,
                 row<zero, zero, variable, variable, zero>>,
    substructure<row<one, variable, zero, zero, zero, zero>,
                 row<variable, one, variable, zero, zero, zero>,
                 row<zero, variable, one, variable, zero, zero>,
                 row<zero, zero, variable, one, variable, zero>,
                 row<zero, zero, zero, variable, one, variable>,
                 row<zero, zero, zero, zero, variable, one>>,
    substructure<row<variable, zero, zero, zero, zero, zero, one>,
                 row<zero, variable, zero, zero, zero, one, zero>,
                 row<zero, zero, variable, zero, one, zero, zero>,
                 row<zero, zero, zero, variable, zero, zero, zero>,
                 row<zero, zero, one, zero, variable, zero, zero>,
                 row<zero, one, zero, zero, zero, variable, zero>,
                 row<one, zero, zero, zero, zero, zero, variable>>,
    substructure<row<variable, variable, variable, variable>>,
    substructure<row<variable>, row<one>, row<variable>, row<zero>>,
    substructure<row<variable, zero, variable, one, variable, variable>>,
    substructure<row<variable>, row<variable>, row<one>, row<variable>,
                 row<zero>, row<variable>>,
    substructure<row<one, variable, zero, variable>>,
    substructure<row<variable>, row<variable>, row<variable>, row<variable>>,
    substructure<row<index_variable<0>, index_variable<1>, index_variable<2>>,
                 row<index_variable<1>, index_variable<3>, index_variable<4>>,
                 row<index_variable<2>, index_variable<4>, index_variable<5>>>,
    substructure<row<index_variable<0>, index_variable<0>>,
                 row<index_variable<0>, index_variable<0>>>,
    substructure<row<index_variable<0>, zero, zero, zero>,
                 row<zero, index_variable<0>, zero, zero>,
                 row<zero, zero, index_variable<0>, zero>,
                 row<zero, zero, zero, index_variable<0>>>,
    substructure<row<index_variable<0>, variable, one>,
                 row<zero, index_variable<1>, variable>,
                 row<index_variable<0>, zero, index_variable<1>>>,
    substructure<
        row<index_variable<0>, variable, index_variable<1>, index_variable<0>>>,
    substructure<row<index_variable<0>>, row<index_variable<1>>,
                 row<index_variable<2>>, row<index_variable<3>>>,
    substructure<row<index_variable<0>, index_variable<1>, index_variable<2>,
                     index_variable<3>>,
                 row<index_variable<0>, index_variable<1>, index_variable<2>,
                     index_variable<3>>>,
    substructure<row<index_variable<7>, variable, one>,
                 row<zero, index_variable<0>, variable>,
                 row<index_variable<0>, zero, index_variable<7>>>,
    substructure<row<variable>>,                 // S26: 1x1, free variable
    substructure<row<zero>>,                     // S27: 1x1, fixed 0 (singular)
    substructure<row<integral_value<2>>>,        // S28: 1x1, fixed non-0
    substructure<row<variable, variable>>,       // S29: 1x2
    substructure<row<variable>, row<variable>>,  // S30: 2x1
    substructure<row<integral_value<1>, integral_value<2>>,
                 row<integral_value<2>, integral_value<4>>>,
    cfd_ccv_a, cfd_ccv_b,    // S32, S33: 3x3
    cfd_cvc_a, cfd_cvc_b,    // S34, S35: 3x3
    cfd_vccv_a, cfd_vccv_b,  // S36, S37: 4x4
    substructure<row<index_variable<0>, index_variable<1>>,
                 row<index_variable<1>, index_variable<2>>>,
    substructure<row<index_variable<0>, index_variable<1>, zero>,
                 row<index_variable<1>, integral_value<2>, index_variable<2>>,
                 row<zero, index_variable<2>, index_variable<3>>>,
    substructure<row<variable, variable>, row<variable, variable>>,
    substructure<row<index_variable<0>, zero>, row<zero, index_variable<0>>>,
    substructure<row<index_variable<0>, index_variable<0>>,
                 row<index_variable<0>, index_variable<3>>>>;

inline constexpr std::size_t num_substructures =
    std::tuple_size_v<substructure_tuple>;

std::string sub_name(std::size_t i) {
  return "S" + std::to_string(i);
}

/// A substructure paired with its position in @c substructure_tuple. The typed
/// tests below run one case per substructure, and the wrapper is what lets a
/// case name itself the way the original standalone harness did (S0, S1, ...)
/// and recover its own index to drive the pair sweep. Wrapping also keeps the
/// list duplicate-safe: two structurally identical entries stay distinct types.
template <std::size_t I, typename Sub>
struct numbered {
  static constexpr std::size_t index = I;
  using type = Sub;
};

template <typename Tuple, typename Seq>
struct numbered_types;

template <typename... Subs, std::size_t... Is>
struct numbered_types<std::tuple<Subs...>, std::index_sequence<Is...>> {
  using type = ::testing::Types<numbered<Is, Subs>...>;
};

using substructure_list =
    typename numbered_types<substructure_tuple,
                            std::make_index_sequence<num_substructures>>::type;

/// Names the typed-test cases S0..S42 rather than 0..42.
struct substructure_name {
  template <typename T>
  static std::string GetName(int i) {
    return sub_name(static_cast<std::size_t>(i));
  }
};

/// One row of the cartesian sweep: A against every B in the list.
template <typename A, std::size_t Ai, std::size_t... Bi>
void test_pair_row(std::index_sequence<Bi...>) {
  (test_pair<A, std::tuple_element_t<Bi, substructure_tuple>>(sub_name(Ai),
                                                              sub_name(Bi)),
   ...);
}

}  // namespace

template <typename T>
class detray_ksm_substructure : public ::testing::Test {};

TYPED_TEST_SUITE(detray_ksm_substructure, substructure_list, substructure_name);

// Single-matrix operations: negation, transpose, inverse, and the structural
// invariants those must preserve.
TYPED_TEST(detray_ksm_substructure, unary_ops) {
  test_single<typename TypeParam::type>(sub_name(TypeParam::index));
}

// This substructure paired against every substructure in the list: addition,
// subtraction, multiplication and congruence, wherever the shapes allow.
TYPED_TEST(detray_ksm_substructure, pairwise_ops) {
  test_pair_row<typename TypeParam::type, TypeParam::index>(
      std::make_index_sequence<num_substructures>{});
}

// Constant folding in a dot product: however many constant terms a dot product
// contains, they must collapse to a single compile-time offset.
GTEST_TEST(detray_ksm, constant_fold_dot) {
  // 23 + v0*v1                 -> 1 mul, 1 add
  check_dot_ops<cfd_ccv_a, cfd_ccv_b>("const+const+var", 1, 1);
  // 23 + v0*v1                 -> 1 mul, 1 add
  check_dot_ops<cfd_cvc_a, cfd_cvc_b>("const+var+const", 1, 1);
  // v0*v1 + v2*v3 + 23         -> 2 mul, 2 add
  check_dot_ops<cfd_vccv_a, cfd_vccv_b>("var+const+const+var", 2, 2);
}

// A cell the substructure fixes has no storage, so producing it must cost
// nothing. K is how many columns of the subtrahend are stored: the rest are
// structurally zero, so those columns of I - B remain the identity's constants.
GTEST_TEST(detray_ksm, structural_cells_are_not_evaluated) {
  // Every column structural: the whole result is constant and costs nothing.
  check_no_structural_evaluation<6, 0>("6x6, no stored column");
  check_no_structural_evaluation<6, 1>("6x6, one stored column");
  check_no_structural_evaluation<6, 2>("6x6, two stored columns");
  check_no_structural_evaluation<6, 3>("6x6, three stored columns");
  // No column structural: the baseline, which was always counted correctly.
  check_no_structural_evaluation<6, 6>("6x6, fully dense");
}
