// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include <array>
#include <tuple>
#include <vector>

namespace Acts {
namespace Experimental {

/// Axis generators are used to allow defining different grid types
/// for indexed geometry objects.
namespace GridAxisGenerators {

/// @brief  Templated base generator for equidistant axis as a tuple - 1D
///
/// @tparam aType the type of the axis (Bound, Closed, Open)
template <detail::AxisBoundaryType aType>
struct Eq {
  /// Broadcast the return_type
  using return_type =
      std::tuple<detail::Axis<detail::AxisType::Equidistant, aType>>;

  /// Broadcast the grid type
  template <typename T>
  using grid_type =
      detail::Grid<T, detail::Axis<detail::AxisType::Equidistant, aType>>;

  std::array<ActsScalar, 2u> range = {};
  std::size_t nBins = 0u;

  /// Call operator that generates the Axis
  return_type operator()() const {
    detail::Axis<detail::AxisType::Equidistant, aType> eAxis(range[0u],
                                                             range[1u], nBins);
    return {eAxis};
  }
};

// All 1D equidistant options
using EqBound = Eq<detail::AxisBoundaryType::Bound>;
using EqOpen = Eq<detail::AxisBoundaryType::Open>;
using EqClosed = Eq<detail::AxisBoundaryType::Closed>;

/// @brief  Templated base generator for vairable axis as a tuple - 1D
///
/// @tparam aType the type of the axis (Bound, Closed, Open)
template <detail::AxisBoundaryType aType>
struct Var {
  /// Broadcast the return_type
  using return_type =
      std::tuple<detail::Axis<detail::AxisType::Variable, aType>>;

  /// Broadcast the grid type
  template <typename T>
  using grid_type =
      detail::Grid<T, detail::Axis<detail::AxisType::Variable, aType>>;

  std::vector<ActsScalar> edges = {};

  /// Call operator that generates the Axis
  return_type operator()() const {
    detail::Axis<detail::AxisType::Variable, aType> vAxis(edges);
    return {vAxis};
  }
};

// All 1D variable options
using VarBound = Var<detail::AxisBoundaryType::Bound>;
using VarOpen = Var<detail::AxisBoundaryType::Open>;
using VarClosed = Var<detail::AxisBoundaryType::Closed>;

/// @brief  Templated base generator for two equidistant axes as a tuple - 2D
///
/// @tparam aType the type of the first axis (Bound, Closed, Open)
/// @tparam bType the type of the second axis (Bound, Closed, Open)
///
template <detail::AxisBoundaryType aType, detail::AxisBoundaryType bType>
struct EqEq {
  /// Broadcast the return_type
  using return_type =
      std::tuple<detail::Axis<detail::AxisType::Equidistant, aType>,
                 detail::Axis<detail::AxisType::Equidistant, bType>>;

  /// Broadcast the grid type
  template <typename T>
  using grid_type =
      detail::Grid<T, detail::Axis<detail::AxisType::Equidistant, aType>,
                   detail::Axis<detail::AxisType::Equidistant, bType>>;

  std::array<ActsScalar, 2u> range0 = {};
  std::size_t nBins0 = 0u;
  std::array<ActsScalar, 2u> range1 = {};
  std::size_t nBins1 = 1u;

  /// Call operator that generates the Axis
  return_type operator()() const {
    // Create the two axis
    detail::Axis<detail::AxisType::Equidistant, aType> aEq(range0[0u],
                                                           range0[1u], nBins0);
    detail::Axis<detail::AxisType::Equidistant, bType> bEq(range1[0u],
                                                           range1[1u], nBins1);
    return std::tie(aEq, bEq);
  }
};

// All 2D EqEq options
using EqBoundEqBound =
    EqEq<detail::AxisBoundaryType::Bound, detail::AxisBoundaryType::Bound>;
using EqBoundEqOpen =
    EqEq<detail::AxisBoundaryType::Bound, detail::AxisBoundaryType::Open>;
using EqBoundEqClosed =
    EqEq<detail::AxisBoundaryType::Bound, detail::AxisBoundaryType::Closed>;
using EqOpenEqBound =
    EqEq<detail::AxisBoundaryType::Open, detail::AxisBoundaryType::Bound>;
using EqOpenEqOpen =
    EqEq<detail::AxisBoundaryType::Open, detail::AxisBoundaryType::Open>;
using EqOpenEqClosed =
    EqEq<detail::AxisBoundaryType::Open, detail::AxisBoundaryType::Closed>;
using EqClosedEqBound =
    EqEq<detail::AxisBoundaryType::Closed, detail::AxisBoundaryType::Bound>;
using EqClosedEqOpen =
    EqEq<detail::AxisBoundaryType::Closed, detail::AxisBoundaryType::Open>;
using EqClosedEqClosed =
    EqEq<detail::AxisBoundaryType::Closed, detail::AxisBoundaryType::Closed>;

/// @brief  Templated base generator for equidistant / variable axes as a tuple - 2D
///
/// @tparam aType the type of the first axis (Bound, Closed, Open)
/// @tparam bType the type of the second axis (Bound, Closed, Open)
///
template <detail::AxisBoundaryType aType, detail::AxisBoundaryType bType>
struct EqVar {
  /// Broadcast the return_type
  using return_type =
      std::tuple<detail::Axis<detail::AxisType::Equidistant, aType>,
                 detail::Axis<detail::AxisType::Variable, bType>>;

  /// Broadcast the grid type
  template <typename T>
  using grid_type =
      detail::Grid<T, detail::Axis<detail::AxisType::Equidistant, aType>,
                   detail::Axis<detail::AxisType::Variable, bType>>;

  std::array<ActsScalar, 2u> range = {};
  std::size_t nBins = 0u;
  std::vector<ActsScalar> edges = {};

  /// Call operator that generates the Axis
  return_type operator()() const {
    detail::Axis<detail::AxisType::Equidistant, aType> eqA(range[0u], range[1u],
                                                           nBins);
    detail::Axis<detail::AxisType::Variable, bType> varB(edges);
    return std::tie(eqA, varB);
  }
};

// All 2D EqVar options
using EqBoundVarBound =
    EqVar<detail::AxisBoundaryType::Bound, detail::AxisBoundaryType::Bound>;
using EqBoundVarOpen =
    EqVar<detail::AxisBoundaryType::Bound, detail::AxisBoundaryType::Open>;
using EqBoundVarClosed =
    EqVar<detail::AxisBoundaryType::Bound, detail::AxisBoundaryType::Closed>;
using EqOpenVarBound =
    EqVar<detail::AxisBoundaryType::Open, detail::AxisBoundaryType::Bound>;
using EqOpenVarOpen =
    EqVar<detail::AxisBoundaryType::Open, detail::AxisBoundaryType::Open>;
using EqOpenVarClosed =
    EqVar<detail::AxisBoundaryType::Open, detail::AxisBoundaryType::Closed>;
using EqClosedVarBound =
    EqVar<detail::AxisBoundaryType::Closed, detail::AxisBoundaryType::Bound>;
using EqClosedVarOpen =
    EqVar<detail::AxisBoundaryType::Closed, detail::AxisBoundaryType::Open>;
using EqClosedVarClosed =
    EqVar<detail::AxisBoundaryType::Closed, detail::AxisBoundaryType::Closed>;

/// @brief  Templated base generator for a variable, equidistant axes tuple - 2D
///
/// @tparam aType the type of the first axis (Bound, Closed, Open)
/// @tparam bType the type of the second axis (Bound, Closed, Open)
///
template <detail::AxisBoundaryType aType, detail::AxisBoundaryType bType>
struct VarEq {
  /// Broadcast the return_type
  using return_type =
      std::tuple<detail::Axis<detail::AxisType::Variable, aType>,
                 detail::Axis<detail::AxisType::Equidistant, bType>>;

  /// Broadcast the grid type
  template <typename T>
  using grid_type =
      detail::Grid<T, detail::Axis<detail::AxisType::Variable, aType>,
                   detail::Axis<detail::AxisType::Equidistant, bType>>;

  std::vector<ActsScalar> edges = {};
  std::array<ActsScalar, 2u> range = {};
  std::size_t nBins = 0u;

  /// Call operator that generates the Axis
  return_type operator()() const {
    detail::Axis<detail::AxisType::Variable, aType> varA(edges);
    detail::Axis<detail::AxisType::Equidistant, bType> eqB(range[0u], range[1u],
                                                           nBins);
    return std::tie(varA, eqB);
  }
};

// All 2D VarEq options
using VarBoundEqBound =
    VarEq<detail::AxisBoundaryType::Bound, detail::AxisBoundaryType::Bound>;
using VarBoundEqOpen =
    VarEq<detail::AxisBoundaryType::Bound, detail::AxisBoundaryType::Open>;
using VarBoundEqClosed =
    VarEq<detail::AxisBoundaryType::Bound, detail::AxisBoundaryType::Closed>;
using VarOpenEqBound =
    VarEq<detail::AxisBoundaryType::Open, detail::AxisBoundaryType::Bound>;
using VarOpenEqOpen =
    VarEq<detail::AxisBoundaryType::Open, detail::AxisBoundaryType::Open>;
using VarOpenEqClosed =
    VarEq<detail::AxisBoundaryType::Open, detail::AxisBoundaryType::Closed>;
using VarClosedEqBound =
    VarEq<detail::AxisBoundaryType::Closed, detail::AxisBoundaryType::Bound>;
using VarClosedEqOpen =
    VarEq<detail::AxisBoundaryType::Closed, detail::AxisBoundaryType::Open>;
using VarClosedEqClosed =
    VarEq<detail::AxisBoundaryType::Closed, detail::AxisBoundaryType::Closed>;

/// @brief  Templated base generator for a two variable axes tuple - 2D
///
/// @tparam aType the type of the first axis (Bound, Closed, Open)
/// @tparam bType the type of the second axis (Bound, Closed, Open)
///
template <detail::AxisBoundaryType aType, detail::AxisBoundaryType bType>
struct VarVar {
  /// Broadcast the return_type
  using return_type =
      std::tuple<detail::Axis<detail::AxisType::Variable, aType>,
                 detail::Axis<detail::AxisType::Variable, bType>>;

  /// Broadcast the grid type
  template <typename T>
  using grid_type =
      detail::Grid<T, detail::Axis<detail::AxisType::Variable, aType>,
                   detail::Axis<detail::AxisType::Variable, bType>>;

  std::vector<ActsScalar> edges0 = {};
  std::vector<ActsScalar> edges1 = {};

  /// Call operator that generates the Axis
  return_type operator()() const {
    detail::Axis<detail::AxisType::Variable, aType> varA(edges0);
    detail::Axis<detail::AxisType::Variable, bType> varB(edges1);
    return std::tie(varA, varB);
  }
};

// All 2D VarVar options
using VarBoundVarBound =
    VarVar<detail::AxisBoundaryType::Bound, detail::AxisBoundaryType::Bound>;
using VarBoundVarOpen =
    VarVar<detail::AxisBoundaryType::Bound, detail::AxisBoundaryType::Open>;
using VarBoundVarClosed =
    VarVar<detail::AxisBoundaryType::Bound, detail::AxisBoundaryType::Closed>;
using VarOpenVarBound =
    VarVar<detail::AxisBoundaryType::Open, detail::AxisBoundaryType::Bound>;
using VarOpenVarOpen =
    VarVar<detail::AxisBoundaryType::Open, detail::AxisBoundaryType::Open>;
using VarOpenVarClosed =
    VarVar<detail::AxisBoundaryType::Open, detail::AxisBoundaryType::Closed>;
using VarClosedVarBound =
    VarVar<detail::AxisBoundaryType::Closed, detail::AxisBoundaryType::Bound>;
using VarClosedVarOpen =
    VarVar<detail::AxisBoundaryType::Closed, detail::AxisBoundaryType::Open>;
using VarClosedVarClosed =
    VarVar<detail::AxisBoundaryType::Closed, detail::AxisBoundaryType::Closed>;

}  // namespace GridAxisGenerators

}  // namespace Experimental
}  // namespace Acts
