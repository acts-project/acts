// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointColumnProxy.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/EventData/StripSpacePointCalibrationDetails.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Seeding/GbtsLayerDescription.hpp"

#include <cstdint>
#include <memory>
#include <optional>

namespace Acts {

class GraphBasedTrackSeeder;
namespace detail {
class GbtsNodeStorageImpl;
}

/// Storage for the GBTS graph nodes.
///
/// Nodes go in one at a time through `insert`, which takes plain scalars so
/// that a caller can fill the storage from its own space point EDM. `finalize`
/// then orders the nodes by (eta bin, phi) into a space point container, with
/// the derived per-node data in dynamic columns on that container.
class GbtsNodeStorage final {
 public:
  /// Move the handle; node and column addresses remain stable.
  GbtsNodeStorage(GbtsNodeStorage&&) noexcept;
  /// Move-assign the handle.
  GbtsNodeStorage& operator=(GbtsNodeStorage&&) noexcept;
  GbtsNodeStorage(const GbtsNodeStorage&) = delete;
  GbtsNodeStorage& operator=(const GbtsNodeStorage&) = delete;
  /// Destroy the event storage.
  ~GbtsNodeStorage();

  //! [gbts insert]
  /// Insert a space point, deriving r and phi from the global position.
  /// @param index Index of the space point in the caller's own collection
  /// @param x Global x coordinate
  /// @param y Global y coordinate
  /// @param z Global z coordinate
  /// @param layerIndex Dense GBTS layer index
  /// @param clusterWidth Pixel cluster width
  /// @param localPositionY Local y cluster position
  /// @return The eta bin the node was placed in, or nullopt if it was rejected
  std::optional<std::uint32_t> insert(SpacePointIndex index, float x, float y,
                                      float z, GbtsLayerIndex layerIndex,
                                      float clusterWidth = 0.f,
                                      float localPositionY = 0.f);
  //! [gbts insert]

  /// Insert a space point for callers that already have r and phi.
  /// @param index Index of the space point in the caller's own collection
  /// @param x Global x coordinate
  /// @param y Global y coordinate
  /// @param z Global z coordinate
  /// @param r Transverse distance from the beamline
  /// @param phi Azimuthal angle in the xy plane
  /// @param layerIndex Dense GBTS layer index
  /// @param clusterWidth Pixel cluster width
  /// @param localPositionY Local y cluster position
  /// @param strip Stereo pair the point was formed from, null for a pixel one
  /// @return The eta bin the node was placed in, or nullopt if it was rejected
  std::optional<std::uint32_t> insert(
      SpacePointIndex index, float x, float y, float z, float r, float phi,
      GbtsLayerIndex layerIndex, float clusterWidth = 0.f,
      float localPositionY = 0.f,
      const OuterStripSpacePointCalibrationDetails* strip = nullptr);

  /// Insert a space point from an ACTS space point container.
  /// @param sp The space point to insert
  /// @param layerColumn Column holding the dense GBTS layer index
  /// @param clusterWidthColumn Column holding the pixel cluster width
  /// @param localPositionYColumn Column holding the local y cluster position
  /// @param strips Whether the container carries the stereo pairs
  /// @return The eta bin the node was placed in, or nullopt if it was rejected
  std::optional<std::uint32_t> insert(
      const ConstSpacePointProxy& sp,
      const ConstSpacePointColumnProxy<GbtsLayerIndex>& layerColumn,
      const ConstSpacePointColumnProxy<float>& clusterWidthColumn,
      const ConstSpacePointColumnProxy<float>& localPositionYColumn,
      bool strips = false) {
    return insert(sp.index(), sp.x(), sp.y(), sp.z(), sp.r(), sp.phi(),
                  sp.extra(layerColumn), sp.extra(clusterWidthColumn),
                  sp.extra(localPositionYColumn),
                  strips ? &sp.outerStripCalibrationDetails() : nullptr);
  }

  /// Insert every space point of a container. A container carrying
  /// `SpacePointColumns::StripCalibrationDetails` has its stereo pairs taken
  /// with it, for the layers the configuration marks as strip layers.
  ///
  /// @param spacePoints The space points to insert
  /// @param layerColumn Column holding the dense GBTS layer index
  /// @param clusterWidthColumn Column holding the pixel cluster width
  /// @param localPositionYColumn Column holding the local y cluster position
  void extend(const SpacePointContainer& spacePoints,
              const ConstSpacePointColumnProxy<GbtsLayerIndex>& layerColumn,
              const ConstSpacePointColumnProxy<float>& clusterWidthColumn,
              const ConstSpacePointColumnProxy<float>& localPositionYColumn);

  /// Sort the nodes by (eta bin, phi) and build the derived per-node data.
  /// Must be called once after all inserts and before the storage is read.
  void finalize();

  /// Get the total number of nodes
  /// @return Total number of nodes
  std::uint32_t numberOfNodes() const;

  /// Map a node index back to the index the caller used when inserting it.
  /// @param node Node index
  /// @return The caller's space point index
  SpacePointIndex spacePointIndex(SpacePointIndex node) const;

  /// Whether any node carries a stereo pair, i.e. whether the graph has a
  /// strip path to take. A caller that fed strip space points in can check
  /// here that their pairs arrived.
  /// @return Whether there are any
  bool hasStrips() const;

 private:
  friend class GraphBasedTrackSeeder;
  explicit GbtsNodeStorage(std::unique_ptr<detail::GbtsNodeStorageImpl> impl);
  std::unique_ptr<detail::GbtsNodeStorageImpl> m_impl;
};

}  // namespace Acts
