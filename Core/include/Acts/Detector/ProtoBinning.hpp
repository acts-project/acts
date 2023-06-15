// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/detail/AxisFwd.hpp"

#include <stdexcept>
#include <vector>

namespace Acts {

namespace Experimental {

/// @brief  Simple helper class to define a binning structure
///
/// @note no checks are performed on the consistency, this is
/// only for convenience that the binning can be defined and then
/// translated into concrete axis types
struct ProtoBinning {
  /// The binning value of this
  BinningValue binValue = BinningValue::binValues;
  /// The axis type: equidistant or variable
  Acts::detail::AxisType axisType = Acts::detail::AxisType::Equidistant;
  /// The axis boundary type: Open, Bound or Closed
  Acts::detail::AxisBoundaryType boundaryType =
      Acts::detail::AxisBoundaryType::Bound;
  /// The binning edges
  std::vector<ActsScalar> edges = {};
  /// An expansion for the filling (in bins)
  std::size_t expansion = 0u;

  /// Convenience constructors - recommended for variable binning
  ///
  /// @param bValue the value/cast in which this is binned
  /// @param bType the axis boundary type
  /// @param e the bin edges (variable binning)
  /// @param exp the expansion (in bins)
  ProtoBinning(BinningValue bValue, Acts::detail::AxisBoundaryType bType,
               const std::vector<ActsScalar>& e, std::size_t exp = 0u)
      : binValue(bValue),
        axisType(Acts::detail::AxisType::Variable),
        boundaryType(bType),
        edges(e),
        expansion(exp) {
    if (edges.size() < 2u) {
      throw std::invalid_argument(
          "ProtoBinning: Invalid binning, at least two edges are needed.");
    }
  }

  /// Convenience constructors - recommended for variable binning
  ///
  /// @param bValue the value/cast in which this is binned
  /// @param bType the axis boundary type
  /// @param minE the lowest edge of the binning
  /// @param maxE the highest edge of the binning
  /// @param nbins the number of bins
  /// @param exp the expansion (in bins)
  ProtoBinning(BinningValue bValue, Acts::detail::AxisBoundaryType bType,
               ActsScalar minE, ActsScalar maxE, std::size_t nbins,
               std::size_t exp = 0u)
      : binValue(bValue), boundaryType(bType), expansion(exp) {
    if (minE >= maxE) {
      throw std::invalid_argument(
          "ProtoBinning: Invalid binning, min edge needs to be smaller than "
          "max edge.");
    }
    if (nbins < 1u) {
      throw std::invalid_argument(
          "ProtoBinning: Invalid binning, at least one bin is needed.");
    }

    ActsScalar stepE = (maxE - minE) / nbins;
    edges.reserve(nbins + 1);
    for (std::size_t i = 0; i <= nbins; i++) {
      edges.push_back(minE + i * stepE);
    }
  }

  std::size_t bins() const { return edges.size() - 1u; }
};

}  // namespace Experimental
}  // namespace Acts
