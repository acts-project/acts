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
/// The equation of motion is a second-order equation for the position. A
/// tableau with the coefficients A and b, applied to the first-order system of
/// the free parameters, is the same as the general Runge-Kutta-Nyström method
/// with the position coefficients A^2 and b A (Hairer, Nørsett, Wanner,
/// Solving Ordinary Differential Equations I, 2nd ed., Springer 1993, Section
/// II.14).
///
/// @note The Runge-Kutta-Nyström step of @ref EigenStepper and of the ATLAS
///       stepper is a different fourth-order method than @ref classicalRk4
///       on the first-order system, so the two do not give the same result.
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
  /// (W. Kutta, Z. Math. Phys. 46 (1901) 435-453; Hairer, Nørsett, Wanner,
  /// Solving Ordinary Differential Equations I, Section II.1)
  /// @return the shared tableau
  static std::shared_ptr<const ButcherTableau> classicalRk4();

  /// The Dormand-Prince 5(4) pair RK5(4)7M (J. R. Dormand, P. J. Prince,
  /// A family of embedded Runge-Kutta formulae, J. Comput. Appl. Math. 6
  /// (1980) 19-26, Table 2; Hairer, Nørsett, Wanner, Solving Ordinary
  /// Differential Equations I, Section II.5, DOPRI5). The stepper propagates
  /// the fifth-order solution and uses the fourth-order solution for the
  /// error estimate.
  /// @return the shared tableau
  static std::shared_ptr<const ButcherTableau> dormandPrince54();

  /// The "most efficient" Verner 9(8) pair. The method is described in
  /// J. H. Verner, Numerically optimal Runge-Kutta pairs with interpolants,
  /// Numerical Algorithms 53 (2010) 383-396, doi:10.1007/s11075-009-9290-3.
  /// The coefficients are from J. H. Verner, https://www.sfu.ca/~jverner/,
  /// file RKV98.IIa.Efficient.000000349.081210.CoeffsOnlyRADandFLOATS. The
  /// stepper propagates the ninth-order solution and uses the eighth-order
  /// solution for the error estimate.
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
