// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace Acts {

/// @brief Coefficients of an explicit Runge-Kutta method
///
/// The tableau defines the stages
///
///   Y_i = y + h * sum_{j < i} a(i, j) * K_j,  K_i = f(Y_i),
///
/// the solution y_1 = y + h * sum_i b(i) * K_i and optionally an embedded
/// solution with the weights @ref bEmbedded for an error estimate.
///
/// Applied to the first-order system of the free parameters, a tableau is the
/// same as the Runge-Kutta-Nyström method with the coefficients A^2 and b A,
/// so any tableau for a first-order system works for the equation of motion.
class ButcherTableau {
 public:
  /// Construct and validate a tableau.
  ///
  /// @param name Name of the method
  /// @param order Order of the solution
  /// @param embeddedOrder Order of the embedded solution, ignored without
  ///        embedded weights
  /// @param c Stage nodes
  /// @param a Stage coefficients, one row per stage, each row with the
  ///        coefficients of the previous stages
  /// @param b Weights of the solution
  /// @param bEmbedded Weights of the embedded solution, or empty
  ///
  /// @throws std::invalid_argument if the sizes do not match, if the first
  ///         node is not zero, or if a row of @p a does not sum to its node
  ButcherTableau(std::string name, unsigned order, unsigned embeddedOrder,
                 std::vector<double> c, std::vector<std::vector<double>> a,
                 std::vector<double> b, std::vector<double> bEmbedded);

  /// @return the name of the method
  const std::string& name() const { return m_name; }

  /// @return the number of stages
  std::size_t stages() const { return m_c.size(); }

  /// @return the order of the solution
  unsigned order() const { return m_order; }

  /// @return the order of the embedded solution
  unsigned embeddedOrder() const { return m_embeddedOrder; }

  /// @return true if the tableau has embedded weights
  bool hasEmbedded() const { return !m_bEmbedded.empty(); }

  /// @param i Stage index
  /// @return the node of stage @p i
  double c(std::size_t i) const { return m_c[i]; }

  /// @param i Stage index
  /// @param j Index of a previous stage, j < i
  /// @return the coefficient of stage @p j in stage @p i
  double a(std::size_t i, std::size_t j) const { return m_a[i * stages() + j]; }

  /// @param i Stage index
  /// @return the weight of stage @p i in the solution
  double b(std::size_t i) const { return m_b[i]; }

  /// @param i Stage index
  /// @return the weight of stage @p i in the embedded solution
  double bEmbedded(std::size_t i) const { return m_bEmbedded[i]; }

  /// The classical fourth-order method without an embedded solution
  /// @return the shared tableau
  static std::shared_ptr<const ButcherTableau> classicalRk4();

  /// The Dormand-Prince 5(4) pair (Dormand, Prince 1980)
  /// @return the shared tableau
  static std::shared_ptr<const ButcherTableau> dormandPrince54();

  /// The "most efficient" Verner 9(8) pair (J. H. Verner, Numerical
  /// Algorithms 53, 2010), with the coefficients of OrdinaryDiffEq.jl `Vern9`
  /// @return the shared tableau
  static std::shared_ptr<const ButcherTableau> verner98();

 private:
  std::string m_name;
  unsigned m_order;
  unsigned m_embeddedOrder;
  std::vector<double> m_c;
  std::vector<double> m_a;
  std::vector<double> m_b;
  std::vector<double> m_bEmbedded;
};

}  // namespace Acts
