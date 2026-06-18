// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/IAxis.hpp"

#include "Acts/Utilities/Axis.hpp"

#include <algorithm>
#include <stdexcept>

namespace Acts {

std::unique_ptr<IAxis> IAxis::createEquidistant(
    AxisBoundaryType aBoundaryType, double min, double max, std::size_t nbins,
    std::optional<AxisDirection> direction) {
  using enum AxisType;
  using enum AxisBoundaryType;

  switch (aBoundaryType) {
    case Open:
      return std::make_unique<Axis<Equidistant, Open>>(min, max, nbins,
                                                       direction);
    case Bound:
      return std::make_unique<Axis<Equidistant, Bound>>(min, max, nbins,
                                                        direction);
    case Closed:
      return std::make_unique<Axis<Equidistant, Closed>>(min, max, nbins,
                                                         direction);
    default:  // should never happen
      throw std::logic_error("Unknown axis boundary type");
  }
}

std::unique_ptr<IAxis> IAxis::createVariable(
    AxisBoundaryType aBoundaryType, const std::vector<double>& edges,
    std::optional<AxisDirection> direction) {
  using enum AxisType;
  using enum AxisBoundaryType;

  switch (aBoundaryType) {
    case Open:
      return std::make_unique<Axis<Variable, Open>>(edges, direction);
    case Bound:
      return std::make_unique<Axis<Variable, Bound>>(edges, direction);
    case Closed:
      return std::make_unique<Axis<Variable, Closed>>(edges, direction);
    default:  // should never happen
      throw std::logic_error("Unknown axis boundary type");
  }
}

}  // namespace Acts
