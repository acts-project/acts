// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <optional>

namespace Acts::detail {

/// A scalar objective in @p N parameters, as consumed by the routines below
///
/// The value may be infinite to mark a point the objective does not accept.
template <typename Callable, int N>
concept ScalarObjective =
    requires(const Callable& objective, const Vector<N>& point) {
      { objective(point) } -> std::convertible_to<double>;
    };

/// Stopping criteria and iteration cap for the simplex search
struct NelderMeadOptions {
  /// Iteration cap; a well behaved low-dimensional problem converges in far
  /// fewer steps, so hitting this means something is wrong
  std::size_t maxIterations = 2000;
  /// Convergence threshold on the objective, relative to its value at the
  /// best vertex
  double valueTolerance = 1e-12;
  /// Convergence threshold on the simplex extent in parameter space
  double parameterTolerance = 1e-10;
};

/// Nelder-Mead simplex minimisation of a scalar objective in @p N parameters
///
/// @param objective The function to minimise
/// @param start Initial vertex
/// @param steps Per-coordinate displacement spanning the initial simplex
/// @param options Convergence and iteration settings
/// @return The minimising point, or `std::nullopt` if the search did not
///         converge within the iteration cap or no vertex ever evaluated
///         to a finite value
template <int N, ScalarObjective<N> Callable>
std::optional<Vector<N>> nelderMead(const Callable& objective,
                                    const Vector<N>& start,
                                    const Vector<N>& steps,
                                    const NelderMeadOptions& options = {}) {
  // Standard Nelder-Mead coefficients
  constexpr double reflection = 1.0;
  constexpr double expansion = 2.0;
  constexpr double contraction = 0.5;
  constexpr double shrink = 0.5;

  constexpr int nVertices = N + 1;

  // Initial simplex: displace each coordinate on its own natural scale
  std::array<Vector<N>, nVertices> simplex;
  simplex[0] = start;
  for (int i = 0; i < N; ++i) {
    simplex[i + 1] = start;
    simplex[i + 1](i) += steps(i);
  }

  std::array<double, nVertices> values;
  for (int v = 0; v < nVertices; ++v) {
    values[v] = objective(simplex[v]);
  }
  if (std::ranges::none_of(values, [](double v) { return std::isfinite(v); })) {
    return std::nullopt;
  }

  const auto sortSimplex = [&simplex, &values]() {
    std::array<int, nVertices> order{};
    for (int v = 0; v < nVertices; ++v) {
      order[v] = v;
    }
    std::ranges::sort(
        order, [&values](int a, int b) { return values[a] < values[b]; });
    const std::array<Vector<N>, nVertices> sortedSimplex = [&]() {
      std::array<Vector<N>, nVertices> result;
      for (int v = 0; v < nVertices; ++v) {
        result[v] = simplex[order[v]];
      }
      return result;
    }();
    const std::array<double, nVertices> sortedValues = [&]() {
      std::array<double, nVertices> result;
      for (int v = 0; v < nVertices; ++v) {
        result[v] = values[order[v]];
      }
      return result;
    }();
    simplex = sortedSimplex;
    values = sortedValues;
  };

  for (std::size_t iteration = 0; iteration < options.maxIterations;
       ++iteration) {
    sortSimplex();

    // Converged once the simplex is tiny both in objective and in parameters.
    // Requiring both avoids stopping early on a flat stretch.
    const double valueSpread = std::abs(values[N] - values[0]);
    double parameterSpread = 0;
    for (int v = 1; v < nVertices; ++v) {
      parameterSpread = std::max(
          parameterSpread,
          (simplex[v] - simplex[0]).template lpNorm<Eigen::Infinity>());
    }
    const double scale = std::max(1.0, std::abs(values[0]));
    if (valueSpread < options.valueTolerance * scale &&
        parameterSpread < options.parameterTolerance) {
      return simplex[0];
    }

    // Centroid of all but the worst vertex
    Vector<N> centroid = Vector<N>::Zero();
    for (int v = 0; v + 1 < nVertices; ++v) {
      centroid += simplex[v] / static_cast<double>(nVertices - 1);
    }

    const auto along = [&centroid, &simplex](double factor) {
      return Vector<N>(centroid + factor * (centroid - simplex[N]));
    };

    const Vector<N> reflected = along(reflection);
    const double reflectedValue = objective(reflected);

    if (reflectedValue < values[0]) {
      // Better than the best so far: try to push further in that direction
      const Vector<N> expanded = along(expansion);
      const double expandedValue = objective(expanded);
      if (expandedValue < reflectedValue) {
        simplex[N] = expanded;
        values[N] = expandedValue;
      } else {
        simplex[N] = reflected;
        values[N] = reflectedValue;
      }
      continue;
    }

    if (reflectedValue < values[N - 1]) {
      simplex[N] = reflected;
      values[N] = reflectedValue;
      continue;
    }

    const Vector<N> contracted = along(-contraction);
    const double contractedValue = objective(contracted);
    if (contractedValue < values[N]) {
      simplex[N] = contracted;
      values[N] = contractedValue;
      continue;
    }

    // Nothing helped: pull every vertex towards the best one
    for (int v = 1; v < nVertices; ++v) {
      simplex[v] = simplex[0] + shrink * (simplex[v] - simplex[0]);
      values[v] = objective(simplex[v]);
    }
  }

  return std::nullopt;
}

/// Covariance from the inverse of a finite-difference Hessian at @p point
///
/// @param objective The function whose curvature is probed
/// @param point The point to evaluate the Hessian at, typically the minimum
///              of @p objective
/// @param steps Per-coordinate finite-difference step
/// @return The inverse Hessian, or `std::nullopt` if it is not positive
///         definite, i.e. @p point is not a genuine minimum
template <int N, ScalarObjective<N> Callable>
std::optional<SquareMatrix<N>> numericalCovariance(const Callable& objective,
                                                   const Vector<N>& point,
                                                   const Vector<N>& steps) {
  for (int i = 0; i < N; ++i) {
    if (!(steps(i) > 0)) {
      return std::nullopt;
    }
  }

  const double centre = objective(point);

  SquareMatrix<N> hessian = SquareMatrix<N>::Zero();
  for (int i = 0; i < N; ++i) {
    Vector<N> plusI = point;
    plusI(i) += steps(i);
    Vector<N> minusI = point;
    minusI(i) -= steps(i);

    hessian(i, i) = (objective(plusI) - 2 * centre + objective(minusI)) /
                    (steps(i) * steps(i));

    for (int j = i + 1; j < N; ++j) {
      Vector<N> plusIplusJ = point;
      plusIplusJ(i) += steps(i);
      plusIplusJ(j) += steps(j);
      Vector<N> plusIminusJ = point;
      plusIminusJ(i) += steps(i);
      plusIminusJ(j) -= steps(j);
      Vector<N> minusIplusJ = point;
      minusIplusJ(i) -= steps(i);
      minusIplusJ(j) += steps(j);
      Vector<N> minusIminusJ = point;
      minusIminusJ(i) -= steps(i);
      minusIminusJ(j) -= steps(j);

      const double mixed = (objective(plusIplusJ) - objective(plusIminusJ) -
                            objective(minusIplusJ) + objective(minusIminusJ)) /
                           (4 * steps(i) * steps(j));
      hessian(i, j) = mixed;
      hessian(j, i) = mixed;
    }
  }

  if (!hessian.allFinite()) {
    return std::nullopt;
  }

  const Eigen::LLT<SquareMatrix<N>> llt(hessian);
  if (llt.info() != Eigen::Success) {
    // Not positive definite, so not a minimum we can put errors on
    return std::nullopt;
  }

  const SquareMatrix<N> covariance = llt.solve(SquareMatrix<N>::Identity());
  if (!covariance.allFinite()) {
    return std::nullopt;
  }

  return covariance;
}

}  // namespace Acts::detail
