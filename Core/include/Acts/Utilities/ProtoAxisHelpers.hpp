// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"

// #include <array>
// #include <cstddef>
// #include <iostream>
// #include <iterator>
// #include <memory>
// #include <stdexcept>
#include <string>
#include <vector>

namespace Acts {

/// @brief Get the number of bins from a ProtoAxis
/// @param axis DirectedProtoAxis object
/// @return Number of bins in the axis
inline std::size_t binsOfProtoAxis(DirectedProtoAxis& axis) {
  return axis.getAxis().getNBins();
}

/// @brief Get the total number of bins from multiple ProtoAxes
/// @param axes Vector of DirectedProtoAxis objects
/// @return Total number of bins across all axes
inline std::size_t totalBinsFromProtoAxes(
    const std::vector<DirectedProtoAxis>& axes) {
  return axes[0].getAxis().getNBins() *
         (axes.size() > 1 ? axes[1].getAxis().getNBins() : 1) *
         (axes.size() > 2 ? axes[2].getAxis().getNBins() : 1);
}

/// @brief Get the number of bins from a specific ProtoAxis in a collection
/// @param axes DirectedProtoAxis vector
/// @param ba Bin axis index
/// @return Number of bins in the specified axis
inline std::size_t binsFromProtoAxes(std::vector<DirectedProtoAxis> axes,
                                     std::size_t ba) {
  BinningData bd(axes[ba]);
  return bd.bins();
}

/// @brief Get the bin index from a ProtoAxis using local coordinates
/// @param axis DirectedProtoAxis object
/// @param lp Local position vector
/// @return Bin index corresponding to the local position
inline std::size_t binFromProtoAxis(const DirectedProtoAxis& axis,
                                    const Vector2& lp) {
  BinningData bd(axis);
  return bd.searchLocal(lp);
}

/// @brief Get the bin index from a ProtoAxis using global coordinates
/// @param axis DirectedProtoAxis object
/// @param gp Global position vector
/// @return Bin index corresponding to the global position
inline std::size_t binFromProtoAxis(const DirectedProtoAxis& axis,
                                    const Vector3& gp) {
  BinningData bd(axis);
  const Transform3& invTransform3 = Transform3::Identity().inverse();
  return bd.searchGlobal(invTransform3 * gp);
}

/// @brief Get the bin triple from multiple ProtoAxes using global coordinates
/// @param axes Vector of DirectedProtoAxis objects
/// @param gp Global position vector
/// @return Array of bin indices corresponding to the global position for each axis
inline std::array<std::size_t, 3> binTripleFromProtoAxes(
    const std::vector<DirectedProtoAxis>& axes, const Vector3& gp) {
  const Transform3& invTransform3 = Transform3::Identity().inverse();
  const Vector3& bPosition = invTransform3 * gp;
  std::array<std::size_t, 3> bTriple = {0, 0, 0};
  if (axes.size() > 0) {
    BinningData bd0(axes[0]);
    bTriple[0] = bd0.searchGlobal(bPosition);
  }
  if (axes.size() > 1) {
    BinningData bd1(axes[1]);
    bTriple[1] = bd1.searchGlobal(bPosition);
  }
  if (axes.size() > 2) {
    BinningData bd2(axes[2]);
    bTriple[2] = bd2.searchGlobal(bPosition);
  }
  return bTriple;
}

/// @brief Get the maximum bin index from a specific ProtoAxis in a collection
/// @param axes DirectedProtoAxis vector
/// @param ba Bin axis index
/// @return Maximum bin index in the specified axis
inline std::size_t maxBin(std::vector<DirectedProtoAxis>& axes,
                          std::size_t ba = 0) {
  std::vector<BinningData> binningDataVec;
  binningDataVec.reserve(axes.size());
  for (const auto& axis : axes) {
    binningDataVec.emplace_back(BinningData(axis));
  }
  if (ba >= binningDataVec.size()) {
    return 0;
  }
  return (binningDataVec[ba].bins() - 1);
}

}  // namespace Acts