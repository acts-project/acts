// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/IAxis.hpp"

namespace Acts {
/// @brief Description of a ProtoAxis which holds an IAxis
/// and lets the user to define a certain axis type, boundary type
/// and associated binning.
///
/// In addition to a simple axis definitions, it holds also a description
/// of the axis direction.
class ProtoAxis {
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
  ProtoAxis(AxisDirection aDir, AxisBoundaryType abType, std::size_t nbins);

  /// @brief returns the axis direction
  ///
  /// @return @c AxisDirection of this axis
  AxisDirection getAxisDirection() const;

  /// @brief Return the IAxis representation
  ///
  /// @return @c AxisType of this axis
  const IAxis& getAxis() const;

  /// @brief check if this is an auto-range binning
  bool isAutorange() const;

  /// Dump into a string
  /// @return the string representation
  std::string toString() const;

 protected:
  /// Dispatch to the correct stream operator
  /// @param os output stream
  void toStream(std::ostream& os) const;

  /// The axis direction
  AxisDirection m_axisDir = AxisDirection::AxisX;

  /// The axis representation
  std::unique_ptr<IAxis> m_axis = nullptr;

  /// Indicate if this is a place holder auto-range binning
  bool m_autorange = false;
};

}  // namespace Acts
