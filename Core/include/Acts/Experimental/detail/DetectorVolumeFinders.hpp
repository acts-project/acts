// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Experimental/Detector.hpp"
#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

namespace Acts {
namespace Experimental {
namespace detail {

/// Trial and error volume finder
///
/// @param gctx is the geometry context of this call
/// @param detector is the detector
/// @param position is the position
///
/// @return a detector volume pointer or null
const DetectorVolume* trialAndError(const Acts::GeometryContext& gctx,
                                    const Detector& detector,
                                    const Acts::Vector3& position) {
  for (const auto v : detector.volumes()) {
    if (v->inside(gctx, position)) {
      return v;
    }
  }
  return nullptr;
}

/// Generate a volume finder with trial and error
///
/// @return a managed detector volume finder
inline static ManagedDetectorVolumeFinder tryAllVolumes() {
  ManagedDetectorVolumeFinder managedFinder;
  DetectorVolumeFinder volumeFinder;
  volumeFinder.connect<&trialAndError>();
  managedFinder.delegate = std::move(volumeFinder);
  managedFinder.implementation = nullptr;
  return managedFinder;
}

/// @brief  This is a templated grid surface attacher, it relies on
/// an indexed based grid lookup and a range interation
///
/// @tparam grid_type the grid type
///
template <typename grid_type>
class GridDetectorVolumeFinder : public IDelegateImpl {
 public:
  /// @brief The grid
  grid_type grid;

  std::array<BinningValue, grid_type::DIM> casts;

  /// @brief  Constructor for a grid based surface attacher
  /// @param igrid the grid that is moved into this attacher
  /// @param icasts is the cast values array
  GridDetectorVolumeFinder(
      grid_type&& igrid, const std::array<BinningValue, grid_type::DIM>& icasts)
      : grid(std::move(igrid)), casts(icasts) {}

  /// Trial and error volume finder
  ///
  /// @param gctx is the geometry context of this call
  /// @param detector is the detector
  /// @param position is the position
  ///
  /// @return a detector volume pointer or null
  const DetectorVolume* find(const Acts::GeometryContext& gctx,
                             const Detector& detector,
                             const Acts::Vector3& position) const {
    // Get the surfaces containerd in this volume
    const auto& volumes = detector.volumes();
    // Retrieve the grid indices and fill the surface candidates
    auto vIndex = grid.atPosition(castPosition(position));
    if (vIndex >= 0 and vIndex < volumes.size()){
      return volumes[vIndex];
    }
    return nullptr;
  }

  /// Cast into a lookup position
  ///
  /// @param position is the position of the update call
  std::array<ActsScalar, grid_type::DIM> castPosition(
      const Vector3& position) const {
    std::array<ActsScalar, grid_type::DIM> casted;
    fillCasts(position, casted,
              std::make_integer_sequence<std::size_t, grid_type::DIM>{});
    return casted;
  }

 private:
  /// Unroll cast loop
  /// @param position is the position of the update call
  /// @param a is the array to be filled
  template <typename Array, std::size_t... idx>
  void fillCasts(const Vector3& position, Array& a,
                 std::index_sequence<idx...>) const {
    ((a[idx] = VectorHelpers::cast(position, casts[idx])), ...);
  }
};

}  // namespace detail
}  // namespace Experimental
}  // namespace Acts