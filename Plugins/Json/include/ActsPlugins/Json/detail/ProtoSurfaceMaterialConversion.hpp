// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Acts::detail {

/// Convert legacy proto binning at the I/O boundary. Proto ranges and
/// transforms were replaced by the surface geometry during mapping, so keep
/// ranges deferred. Flatten subdivisions into normalized variable edges and pad
/// an omitted axis with one bin. Cylinder azimuth is expressed in the local
/// rphi coordinate.
/// @param binning the legacy binning description
/// @return the two-dimensional binning spec
inline MultiAxisSpec2D protoSurfaceMaterialBinning(const BinUtility& binning) {
  const auto& data = binning.binningData();
  if (data.size() > 2) {
    throw std::invalid_argument(
        "Proto surface material needs at most two axes");
  }
  std::array specs{AxisSpec::DeferredEquidistant(1),
                   AxisSpec::DeferredEquidistant(1)};
  const bool cylinder = std::ranges::any_of(data, [](const auto& axis) {
    return axis.binvalue == AxisDirection::AxisZ ||
           axis.binvalue == AxisDirection::AxisRPhi;
  });
  for (std::size_t i = 0; i < data.size(); ++i) {
    const auto& axis = data[i];
    auto direction = axis.binvalue;
    if (cylinder && direction == AxisDirection::AxisPhi) {
      direction = AxisDirection::AxisRPhi;
    }
    const auto boundary = axis.option == closed ? AxisBoundaryType::Closed
                                                : AxisBoundaryType::Bound;
    if (axis.type == equidistant && !axis.subBinningData) {
      specs[i] = AxisSpec::Equidistant(axis.bins(), std::nullopt, std::nullopt,
                                       boundary, direction);
    } else {
      const auto& edges = axis.boundaries();
      if (edges.size() < 2 || !(axis.min < axis.max)) {
        throw std::invalid_argument("Invalid legacy proto variable axis");
      }
      std::vector<double> normalized;
      normalized.reserve(edges.size());
      for (float edge : edges) {
        normalized.push_back((static_cast<double>(edge) - axis.min) /
                             (static_cast<double>(axis.max) - axis.min));
      }
      normalized.front() = 0.;
      normalized.back() = 1.;
      specs[i] = AxisSpec::DeferredVariable(std::move(normalized), boundary,
                                            direction);
    }
  }
  if (data.size() == 1) {
    AxisDirection other;
    switch (*specs[0].direction()) {
      case AxisDirection::AxisX:
        other = AxisDirection::AxisY;
        break;
      case AxisDirection::AxisY:
        other = AxisDirection::AxisX;
        break;
      case AxisDirection::AxisR:
        other = AxisDirection::AxisPhi;
        break;
      case AxisDirection::AxisPhi:
        other = AxisDirection::AxisR;
        break;
      case AxisDirection::AxisRPhi:
        other = AxisDirection::AxisZ;
        break;
      case AxisDirection::AxisZ:
        other = AxisDirection::AxisRPhi;
        break;
      default:
        throw std::invalid_argument("Unsupported legacy proto axis direction");
    }
    specs[1] = AxisSpec::DeferredEquidistant(1, other);
  }
  return MultiAxisSpec2D(std::move(specs));
}

}  // namespace Acts::detail
