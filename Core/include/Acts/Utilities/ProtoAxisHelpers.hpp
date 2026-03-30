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

#include <span>
#include <string>
#include <vector>

namespace Acts::ProtoAxisHelpers {

/// @brief Get the number of bins from a ProtoAxis
/// @param axis DirectedProtoAxis object
/// @return Number of bins in the axis
inline std::size_t binsOfProtoAxis(const Acts::DirectedProtoAxis& axis) {
  return axis.getAxis().getNBins();
}

/// @brief Get the total number of bins from multiple ProtoAxes
/// @param axes Span of DirectedProtoAxis objects
/// @return Total number of bins across all axes
inline std::size_t totalBinsFromProtoAxes(
    std::span<const Acts::DirectedProtoAxis> axes) {
  if (axes.empty() || axes.size() > 3) {
    throw std::runtime_error(
        "Unsupported number of axes, must be 1-3, instead got " +
        std::to_string(axes.size()) + ")");
  }
  return axes[0].getAxis().getNBins() *
         (axes.size() > 1 ? axes[1].getAxis().getNBins() : 1) *
         (axes.size() > 2 ? axes[2].getAxis().getNBins() : 1);
}

/// @brief Get the number of bins from a specific ProtoAxis in a collection
/// @param axes DirectedProtoAxis span
/// @param ba Bin axis index
/// @return Number of bins in the specified axis
inline std::size_t binsFromProtoAxes(
    std::span<const Acts::DirectedProtoAxis> axes, std::size_t ba) {
  if (axes.empty() || axes.size() > 3) {
    throw std::runtime_error(
        "Unsupported number of axes, must be 1-3, instead got " +
        std::to_string(axes.size()) + ")");
  }
  Acts::BinningData bd(axes[ba]);
  return bd.bins();
}

/// @brief Get the bin index from a ProtoAxis using local coordinates
/// @param axis DirectedProtoAxis object
/// @param lp Local position vector
/// @return Bin index corresponding to the local position
inline std::size_t binFromProtoAxis(const Acts::DirectedProtoAxis& axis,
                                    const Acts::Vector2& lp) {
  Acts::BinningData bd(axis);
  return bd.searchLocal(lp);
}

/// @brief Get the bin index from a ProtoAxis using global coordinates
/// @param axis DirectedProtoAxis object
/// @param gp Global position vector
/// @return Bin index corresponding to the global position
inline std::size_t binFromProtoAxis(const Acts::DirectedProtoAxis& axis,
                                    const Acts::Vector3& gp) {
  Acts::BinningData bd(axis);
  return bd.searchGlobal(gp);
}

/// @brief Get the bin triple from multiple ProtoAxes using global coordinates
/// @param axes Span of DirectedProtoAxis objects
/// @param gp Global position vector
/// @return Array of bin indices corresponding to the global position for each axis
inline std::array<std::size_t, 3> binTripleFromProtoAxes(
    std::span<const Acts::DirectedProtoAxis> axes, const Acts::Vector3& gp) {
  const Acts::Vector3& bPosition = gp;
  std::array<std::size_t, 3> bTriple = {0, 0, 0};
  if (axes.empty() || axes.size() > 3) {
    throw std::runtime_error(
        "Unsupported number of axes, must be 1-3, instead got " +
        std::to_string(axes.size()) + ")");
  }
  if (axes.size() == 1) {
    Acts::BinningData bd0(axes[0]);
    bTriple[0] = bd0.searchGlobal(bPosition);
  }
  if (axes.size() == 2) {
    Acts::BinningData bd1(axes[1]);
    bTriple[1] = bd1.searchGlobal(bPosition);
  }
  if (axes.size() == 3) {
    Acts::BinningData bd2(axes[2]);
    bTriple[2] = bd2.searchGlobal(bPosition);
  }
  return bTriple;
}

/// @brief Get the maximum bin index from a specific ProtoAxis in a collection
/// @param axes DirectedProtoAxis span
/// @param ba Bin axis index
/// @return Maximum bin index in the specified axis
inline std::size_t maxBin(std::span<const Acts::DirectedProtoAxis> axes,
                          std::size_t ba = 0) {
  if (axes.empty() || axes.size() > 3) {
    throw std::runtime_error(
        "Unsupported number of axes, must be 1-3, instead got " +
        std::to_string(axes.size()) + ")");
  }
  std::vector<Acts::BinningData> binningDataVec;
  binningDataVec.reserve(axes.size());
  for (const auto& axis : axes) {
    binningDataVec.emplace_back(axis);
  }
  if (ba >= binningDataVec.size()) {
    return 0;
  }
  return (binningDataVec.at(ba).bins() - 1);
}

}  // namespace Acts::ProtoAxisHelpers
