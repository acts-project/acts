// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/TypeList.hpp"

#include <array>
#include <tuple>
#include <vector>

/// Axis generators are used to allow defining different grid types
/// for indexed geometry objects.
///
/// The call operator() API allows to plug axis generators into
/// dedicated code snippets and create fitting axis types on the fly
/// which then turn into concrete Grid types.
namespace Acts::GridAxisGenerators {

/// @brief  Templated base generator for equidistant axis as a tuple - 1D
///
/// @tparam aType the type of the axis (Bound, Closed, Open)
template <AxisBoundaryType aType>
struct Eq {
  /// Broadcast the return_type
  using return_type = std::tuple<Axis<AxisType::Equidistant, aType>>;

  /// Broadcast the grid type
  template <typename T>
  using grid_type = Grid<T, Axis<AxisType::Equidistant, aType>>;

  std::array<double, 2u> range = {};
  std::size_t nBins = 0u;

  /// Call operator that generates the Axis
  return_type operator()() const {
    Axis<AxisType::Equidistant, aType> eAxis(range[0u], range[1u], nBins);
    return {eAxis};
  }
};

// All 1D equidistant options
using EqBound = Eq<AxisBoundaryType::Bound>;
using EqOpen = Eq<AxisBoundaryType::Open>;
using EqClosed = Eq<AxisBoundaryType::Closed>;

/// @brief  Templated base generator for variable axis as a tuple - 1D
///
/// @tparam aType the type of the axis (Bound, Closed, Open)
template <AxisBoundaryType aType>
struct Var {
  /// Broadcast the return_type
  using return_type = std::tuple<Axis<AxisType::Variable, aType>>;

  /// Broadcast the grid type
  template <typename T>
  using grid_type = Grid<T, Axis<AxisType::Variable, aType>>;

  std::vector<double> edges = {};

  /// Call operator that generates the Axis
  return_type operator()() const {
    Axis<AxisType::Variable, aType> vAxis(edges);
    return {vAxis};
  }
};

// All 1D variable options
using VarBound = Var<AxisBoundaryType::Bound>;
using VarOpen = Var<AxisBoundaryType::Open>;
using VarClosed = Var<AxisBoundaryType::Closed>;

/// @brief  Templated base generator for two equidistant axes as a tuple - 2D
///
/// @tparam aType the type of the first axis (Bound, Closed, Open)
/// @tparam bType the type of the second axis (Bound, Closed, Open)
///
template <AxisBoundaryType aType, AxisBoundaryType bType>
struct EqEq {
  /// Broadcast the return_type
  using return_type = std::tuple<Axis<AxisType::Equidistant, aType>,
                                 Axis<AxisType::Equidistant, bType>>;

  /// Broadcast the grid type
  template <typename T>
  using grid_type = Grid<T, Axis<AxisType::Equidistant, aType>,
                         Axis<AxisType::Equidistant, bType>>;

  std::array<double, 2u> range0 = {};
  std::size_t nBins0 = 0u;
  std::array<double, 2u> range1 = {};
  std::size_t nBins1 = 1u;

  /// Call operator that generates the Axis
  return_type operator()() const {
    // Create the two axis
    Axis<AxisType::Equidistant, aType> aEq(range0[0u], range0[1u], nBins0);
    Axis<AxisType::Equidistant, bType> bEq(range1[0u], range1[1u], nBins1);
    return {aEq, bEq};
  }
};

// All 2D EqEq options
using EqBoundEqBound = EqEq<AxisBoundaryType::Bound, AxisBoundaryType::Bound>;
using EqBoundEqOpen = EqEq<AxisBoundaryType::Bound, AxisBoundaryType::Open>;
using EqBoundEqClosed = EqEq<AxisBoundaryType::Bound, AxisBoundaryType::Closed>;
using EqOpenEqBound = EqEq<AxisBoundaryType::Open, AxisBoundaryType::Bound>;
using EqOpenEqOpen = EqEq<AxisBoundaryType::Open, AxisBoundaryType::Open>;
using EqOpenEqClosed = EqEq<AxisBoundaryType::Open, AxisBoundaryType::Closed>;
using EqClosedEqBound = EqEq<AxisBoundaryType::Closed, AxisBoundaryType::Bound>;
using EqClosedEqOpen = EqEq<AxisBoundaryType::Closed, AxisBoundaryType::Open>;
using EqClosedEqClosed =
    EqEq<AxisBoundaryType::Closed, AxisBoundaryType::Closed>;

/// @brief  Templated base generator for equidistant / variable axes as a tuple - 2D
///
/// @tparam aType the type of the first axis (Bound, Closed, Open)
/// @tparam bType the type of the second axis (Bound, Closed, Open)
///
template <AxisBoundaryType aType, AxisBoundaryType bType>
struct EqVar {
  /// Broadcast the return_type
  using return_type = std::tuple<Axis<AxisType::Equidistant, aType>,
                                 Axis<AxisType::Variable, bType>>;

  /// Broadcast the grid type
  template <typename T>
  using grid_type = Grid<T, Axis<AxisType::Equidistant, aType>,
                         Axis<AxisType::Variable, bType>>;

  std::array<double, 2u> range = {};
  std::size_t nBins = 0u;
  std::vector<double> edges = {};

  /// Call operator that generates the Axis
  return_type operator()() const {
    Axis<AxisType::Equidistant, aType> eqA(range[0u], range[1u], nBins);
    Axis<AxisType::Variable, bType> varB(edges);
    return {eqA, varB};
  }
};

// All 2D EqVar options
using EqBoundVarBound = EqVar<AxisBoundaryType::Bound, AxisBoundaryType::Bound>;
using EqBoundVarOpen = EqVar<AxisBoundaryType::Bound, AxisBoundaryType::Open>;
using EqBoundVarClosed =
    EqVar<AxisBoundaryType::Bound, AxisBoundaryType::Closed>;
using EqOpenVarBound = EqVar<AxisBoundaryType::Open, AxisBoundaryType::Bound>;
using EqOpenVarOpen = EqVar<AxisBoundaryType::Open, AxisBoundaryType::Open>;
using EqOpenVarClosed = EqVar<AxisBoundaryType::Open, AxisBoundaryType::Closed>;
using EqClosedVarBound =
    EqVar<AxisBoundaryType::Closed, AxisBoundaryType::Bound>;
using EqClosedVarOpen = EqVar<AxisBoundaryType::Closed, AxisBoundaryType::Open>;
using EqClosedVarClosed =
    EqVar<AxisBoundaryType::Closed, AxisBoundaryType::Closed>;

/// @brief  Templated base generator for a variable, equidistant axes tuple - 2D
///
/// @tparam aType the type of the first axis (Bound, Closed, Open)
/// @tparam bType the type of the second axis (Bound, Closed, Open)
///
template <AxisBoundaryType aType, AxisBoundaryType bType>
struct VarEq {
  /// Broadcast the return_type
  using return_type = std::tuple<Axis<AxisType::Variable, aType>,
                                 Axis<AxisType::Equidistant, bType>>;

  /// Broadcast the grid type
  template <typename T>
  using grid_type = Grid<T, Axis<AxisType::Variable, aType>,
                         Axis<AxisType::Equidistant, bType>>;

  std::vector<double> edges = {};
  std::array<double, 2u> range = {};
  std::size_t nBins = 0u;

  /// Call operator that generates the Axis
  return_type operator()() const {
    Axis<AxisType::Variable, aType> varA(edges);
    Axis<AxisType::Equidistant, bType> eqB(range[0u], range[1u], nBins);
    return {varA, eqB};
  }
};

// All 2D VarEq options
using VarBoundEqBound = VarEq<AxisBoundaryType::Bound, AxisBoundaryType::Bound>;
using VarBoundEqOpen = VarEq<AxisBoundaryType::Bound, AxisBoundaryType::Open>;
using VarBoundEqClosed =
    VarEq<AxisBoundaryType::Bound, AxisBoundaryType::Closed>;
using VarOpenEqBound = VarEq<AxisBoundaryType::Open, AxisBoundaryType::Bound>;
using VarOpenEqOpen = VarEq<AxisBoundaryType::Open, AxisBoundaryType::Open>;
using VarOpenEqClosed = VarEq<AxisBoundaryType::Open, AxisBoundaryType::Closed>;
using VarClosedEqBound =
    VarEq<AxisBoundaryType::Closed, AxisBoundaryType::Bound>;
using VarClosedEqOpen = VarEq<AxisBoundaryType::Closed, AxisBoundaryType::Open>;
using VarClosedEqClosed =
    VarEq<AxisBoundaryType::Closed, AxisBoundaryType::Closed>;

/// @brief  Templated base generator for a two variable axes tuple - 2D
///
/// @tparam aType the type of the first axis (Bound, Closed, Open)
/// @tparam bType the type of the second axis (Bound, Closed, Open)
///
template <AxisBoundaryType aType, AxisBoundaryType bType>
struct VarVar {
  /// Broadcast the return_type
  using return_type = std::tuple<Axis<AxisType::Variable, aType>,
                                 Axis<AxisType::Variable, bType>>;

  /// Broadcast the grid type
  template <typename T>
  using grid_type =
      Grid<T, Axis<AxisType::Variable, aType>, Axis<AxisType::Variable, bType>>;

  std::vector<double> edges0 = {};
  std::vector<double> edges1 = {};

  /// Call operator that generates the Axis
  return_type operator()() const {
    Axis<AxisType::Variable, aType> varA(edges0);
    Axis<AxisType::Variable, bType> varB(edges1);
    return {varA, varB};
  }
};

// All 2D VarVar options
using VarBoundVarBound =
    VarVar<AxisBoundaryType::Bound, AxisBoundaryType::Bound>;
using VarBoundVarOpen = VarVar<AxisBoundaryType::Bound, AxisBoundaryType::Open>;
using VarBoundVarClosed =
    VarVar<AxisBoundaryType::Bound, AxisBoundaryType::Closed>;
using VarOpenVarBound = VarVar<AxisBoundaryType::Open, AxisBoundaryType::Bound>;
using VarOpenVarOpen = VarVar<AxisBoundaryType::Open, AxisBoundaryType::Open>;
using VarOpenVarClosed =
    VarVar<AxisBoundaryType::Open, AxisBoundaryType::Closed>;
using VarClosedVarBound =
    VarVar<AxisBoundaryType::Closed, AxisBoundaryType::Bound>;
using VarClosedVarOpen =
    VarVar<AxisBoundaryType::Closed, AxisBoundaryType::Open>;
using VarClosedVarClosed =
    VarVar<AxisBoundaryType::Closed, AxisBoundaryType::Closed>;

// Generate the possible axes in this case
using PossibleAxes =
    TypeList<EqBound, EqOpen, EqClosed,
             // All 1D Var  options
             VarBound, VarOpen, VarClosed,
             // All 2D EqEq options
             EqBoundEqBound, EqBoundEqOpen, EqBoundEqClosed, EqOpenEqBound,
             EqOpenEqOpen, EqOpenEqClosed, EqClosedEqBound, EqClosedEqOpen,
             EqClosedEqClosed,
             // All 2D EqVar options
             EqBoundVarBound, EqBoundVarOpen, EqBoundVarClosed, EqOpenVarBound,
             EqOpenVarOpen, EqOpenVarClosed, EqClosedVarBound, EqClosedVarOpen,
             EqClosedVarClosed,
             // All 2D VarEq options
             VarBoundEqBound, VarBoundEqOpen, VarBoundEqClosed, VarOpenEqBound,
             VarOpenEqOpen, VarOpenEqClosed, VarClosedEqBound, VarClosedEqOpen,
             VarClosedEqClosed,
             // All 2D VarEq options
             VarBoundVarBound, VarBoundVarOpen, VarBoundVarClosed,
             VarOpenVarBound, VarOpenVarOpen, VarOpenVarClosed,
             VarClosedVarBound, VarClosedVarOpen, VarClosedVarClosed>;

}  // namespace Acts::GridAxisGenerators
