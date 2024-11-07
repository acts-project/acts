// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/BinUtility.hpp"

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Acts::Experimental {

/// @brief  Simple helper class to define a binning structure
///
/// @note no checks are performed on the consistency, this is
/// only for convenience that the binning can be defined and then
/// translated into concrete axis types
struct ProtoBinning {
  /// The binning value of this
  AxisDirection axisDirection;
  /// The axis type: equidistant or variable
  AxisType axisType = AxisType::Equidistant;
  /// The axis boundary type: Open, Bound or Closed
  AxisBoundaryType axisBoundaryType = AxisBoundaryType::Bound;
  /// The binning edges
  std::vector<double> edges = {};
  /// An expansion for the filling (in bins)
  std::size_t expansion = 0u;
  /// Indication if this is an auto-range binning
  bool autorange = false;

  /// Convenience constructors - for variable binning
  ///
  /// @param aDescr The axis description/direction for the binning
  /// @param bType the axis boundary type
  /// @param e the bin edges (variable binning)
  /// @param exp the expansion (in bins)
  ProtoBinning(AxisDirection aDescr, AxisBoundaryType bType,
               const std::vector<double>& e, std::size_t exp = 0u)
      : axisDirection(aDescr),
        axisType(AxisType::Variable),
        axisBoundaryType(bType),
        edges(e),
        expansion(exp) {
    if (edges.size() < 2u) {
      throw std::invalid_argument(
          "ProtoBinning: Invalid binning, at least two edges are needed.");
    }
  }

  /// Convenience constructors - for equidistant binning
  ///
  /// @param aDescr The axis description/direction for the binning
  /// @param bType the axis boundary type
  /// @param minE the lowest edge of the binning
  /// @param maxE the highest edge of the binning
  /// @param nbins the number of bins
  /// @param exp the expansion (in bins)
  ProtoBinning(AxisDirection aDescr, AxisBoundaryType bType, double minE,
               double maxE, std::size_t nbins, std::size_t exp = 0u)
      : axisDirection(aDescr), axisBoundaryType(bType), expansion(exp) {
    if (minE >= maxE) {
      std::string msg = "ProtoBinning: Invalid binning for value '";
      msg += axisDirectionToString(aDescr);
      msg += "', min edge (" + std::to_string(minE) + ") ";
      msg += " needs to be smaller than max edge (";
      msg += std::to_string(maxE) + ").";
      throw std::invalid_argument(msg);
    }
    if (nbins < 1u) {
      throw std::invalid_argument(
          "ProtoBinning: Invalid binning, at least one bin is needed.");
    }

    double stepE = (maxE - minE) / nbins;
    edges.reserve(nbins + 1);
    for (std::size_t i = 0; i <= nbins; i++) {
      edges.push_back(minE + i * stepE);
    }
  }

  /// Placeholder constructors - for equidistant binning
  ///
  /// @note this is designed to give a binning prescription
  /// when the actual extent is not yet evaluated, only works
  /// for equidistant binning obviously
  ///
  /// @param aDescr The axis description/direction for the binning
  /// @param bType the axis boundary type
  /// @param nbins the number of bins
  /// @param exp the expansion (in bins)
  ProtoBinning(AxisDirection aDescr, AxisBoundaryType bType, std::size_t nbins,
               std::size_t exp = 0u)
      : axisDirection(aDescr),
        axisBoundaryType(bType),
        edges(nbins + 1, 0.),
        expansion(exp),
        autorange(true) {}

  // Return the number of bins
  std::size_t bins() const { return edges.size() - 1u; }

  // Screen output
  std::string toString() const {
    std::stringstream ss;
    ss << "ProtoBinning: " << bins() << " bins in " << axisDirection;
    ss << axisType;
    if (!autorange) {
      ss << "within [" << edges.front() << ", " << edges.back() << "] ";
    } else {
      ss << "within automatic range";
    }
    return ss.str();
  }
};

/// @brief A binning description, it helps for screen output
struct BinningDescription {
  /// Convert the binning description into a bin utility
  ///
  /// @param binUtility the bin utility to be converted into a BinningDescription
  static BinningDescription fromBinUtility(const BinUtility& binUtility) {
    BinningDescription bDesc;
    for (const auto& bData : binUtility.binningData()) {
      std::vector<double> edges;
      if (bData.axisType == AxisType::Equidistant) {
        bDesc.binning.push_back(ProtoBinning(bData.axisDirection,
                                             bData.axisBoundaryType, bData.min,
                                             bData.max, bData.bins(), 0u));

      } else {
        std::for_each(bData.boundaries().begin(), bData.boundaries().end(),
                      [&](double edge) { edges.push_back(edge); });
        bDesc.binning.push_back(ProtoBinning(
            bData.axisDirection, bData.axisBoundaryType, edges, 0u));
      }
    }
    return bDesc;
  }

  /// Convert to a BinUtility - only basic types are supported
  ///
  BinUtility toBinUtility() const {
    BinUtility binUtility;
    for (const auto& b : binning) {
      if (b.axisType == AxisType::Equidistant) {
        binUtility += BinUtility(b.bins(), b.edges.front(), b.edges.back(),
                                 b.axisBoundaryType, b.axisDirection);
      } else {
        std::vector<float> edges;
        std::for_each(b.edges.begin(), b.edges.end(),
                      [&](double edge) { edges.push_back(edge); });
        binUtility += BinUtility(edges, b.axisBoundaryType, b.axisDirection);
      }
    }
    return binUtility;
  }

  /// The contained binnings
  std::vector<ProtoBinning> binning;

  // Screen output
  std::string toString() const {
    std::stringstream ss;
    ss << "BinningDescription: " << binning.size() << "D" << std::endl;
    for (const auto& b : binning) {
      ss << "  " << b.toString() << std::endl;
    }
    return ss.str();
  }
};

}  // namespace Acts::Experimental
