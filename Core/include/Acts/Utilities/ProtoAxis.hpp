// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/IAxis.hpp"

namespace Acts {
/// @brief Description of a ProtoAxis which extends the IAxis
/// and lets the user to define a certain axis type, boundary type
/// and associated binning.
class ProtoAxis : public IAxis {
 public:
  /// Convenience constructors - for variable binning
  ///
  /// @param aDir the value/cast in which this is binned
  /// @param abType the axis boundary type
  /// @param edges the bin edges (variable binning)
  ProtoAxis(AxisDirection aDir, Acts::AxisBoundaryType abType,
            const std::vector<double>& edges);

  /// Convenience constructors - for equidistant binning
  ///
  /// @param aDir the value/cast in which this is binned
  /// @param abType the axis boundary type
  /// @param minE the lowest edge of the binning
  /// @param maxE the highest edge of the binning
  /// @param nbins the number of bins
  ProtoAxis(AxisDirection aDir, AxisBoundaryType abType, double minE,
            double maxE, std::size_t nbins);

  /// Placeholder constructor for auto-range binning
  ///
  /// @param aDir the value/cast in which this is binned
  /// @param abType the axis boundary type
  /// @param nbins the number of bins
  ///
  /// @note that auto-range is only supported for equidistant binning
  ProtoAxis(AxisDirection aDir, AxisBoundaryType abType,
                           std::size_t nbins);

  /// @brief returns whether the axis is equidistant
  ///
  /// @return bool is equidistant
  bool isEquidistant() const final;

  /// @brief returns whether the axis is variable
  ///
  /// @return bool is variable
  bool isVariable() const final;

  /// @brief returns the axis direction
  ///
  /// @return @c AxisDirection of this axis
  AxisDirection getAxisDirection() const;

  /// @brief returns the type of the axis
  ///
  /// @return @c AxisType of this axis
  AxisType getType() const final;

  /// @brief returns the boundary type set in the template param
  ///
  /// @return @c AxisBoundaryType of this axis
  AxisBoundaryType getBoundaryType() const final;

  /// @brief Return a vector of bin edges
  /// @return Vector which contains the bin edges
  std::vector<double> getBinEdges() const final;

  /// @brief get minimum of binning range
  ///
  /// @return minimum of binning range
  double getMin() const final;

  /// @brief get maximum of binning range
  ///
  /// @return maximum of binning range
  double getMax() const final;

  /// @brief get total number of bins
  ///
  /// @return total number of bins (excluding under-/overflow bins)
  std::size_t getNBins() const final;

  /// @brief check if this is an auto-range binning
  bool isAutorange() const;

 protected:
  /// Dispatch to the correct stream operator
  /// @param os output stream
  void toStream(std::ostream& os) const final;

  /// Dump into a string
  /// @return the string representation
  std::string toString() const;

  /// The axis direction
  AxisDirection m_axisDir = AxisDirection::AxisX;

  /// The axis type: equidistant or variable
  AxisType m_axisType = AxisType::Equidistant;

  /// The axis boundary type: Open, Bound or Closed
  AxisBoundaryType m_axisBoundaryType = AxisBoundaryType::Bound;

  /// The bin edges
  std::vector<double> m_edges = {};

  /// Indicate if this is a place holder auto-range binning
  bool m_autorange = false;
};

}  // namespace Acts
