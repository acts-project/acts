// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>

namespace Acts {

class TrackingVolume;
class PortalShellBase;
class Volume;

namespace Experimental {

struct GeometryIdentifierBlueprintNodeImpl;

/// @brief Blueprint node for configuring and applying geometry identifiers to volumes
///
/// This node must have exactly one child and applies geometry identifier
/// configurations to the volumes in its subtree during finalization. Multiple
/// configurations can be chained using the fluent interface.
class GeometryIdentifierBlueprintNode : public BlueprintNode {
 public:
  /// Virtual destructor pushed to cpp file to avoid defining implementation
  /// details
  ~GeometryIdentifierBlueprintNode() override;

  /// Default constructor
  GeometryIdentifierBlueprintNode();

  /// @brief Build the volume hierarchy under this node
  /// @param options Blueprint build options
  /// @param gctx The geometry context
  /// @param logger The logger instance
  /// @return Reference to the built volume
  /// @note Requires exactly one child node
  Volume& build(const BlueprintOptions& options, const GeometryContext& gctx,
                const Logger& logger = Acts::getDummyLogger()) override;

  /// @brief Connect portals in the volume hierarchy
  /// @param options Blueprint build options
  /// @param gctx The geometry context
  /// @param logger The logger instance
  /// @return Reference to the connected portal shell
  PortalShellBase& connect(
      const BlueprintOptions& options, const GeometryContext& gctx,
      const Logger& logger = Acts::getDummyLogger()) override;

  /// @brief Finalize the volume hierarchy and apply geometry identifier configurations
  /// @param options Blueprint build options
  /// @param gctx The geometry context
  /// @param parent The parent tracking volume
  /// @param logger The logger instance
  /// @note Applies all configured geometry ID assignments to new volumes in the subtree
  void finalize(const BlueprintOptions& options, const GeometryContext& gctx,
                TrackingVolume& parent,
                const Logger& logger = Acts::getDummyLogger()) override;

  /// @brief Set a fixed layer ID for volumes in this subtree
  /// @param layer The layer ID value to set
  /// @return Reference to this node for method chaining
  /// @note Will throw if volumes already have layer IDs assigned
  GeometryIdentifierBlueprintNode& setLayerIdTo(
      GeometryIdentifier::Value layer);

  /// @brief Incrementally assign layer IDs to volumes in this subtree
  /// @param start The starting layer ID value (default: 0)
  /// @return Reference to this node for method chaining
  /// @note Will throw if volumes already have layer IDs assigned
  /// @note Layer IDs are assigned sequentially starting from the given value
  GeometryIdentifierBlueprintNode& incrementLayerIds(
      GeometryIdentifier::Value start = 0);

  /// @brief Set the same volume ID for all volumes in this subtree
  /// @param volumeId The volume ID to set
  /// @return Reference to this node for method chaining
  /// @note Will throw if volumes already have volume IDs assigned
  /// @note Applies recursively to all descendant volumes
  GeometryIdentifierBlueprintNode& setAllVolumeIdsTo(
      GeometryIdentifier::Value volumeId);

  /// Predicate function to compare two @ref Acts::TrackingVolume with each other to determine their *closure order*.
  using CompareVolumes =
      std::function<bool(const TrackingVolume&, const TrackingVolume&)>;

  /// @brief Configure this node to order eligible tracking volumes using the provided
  /// function
  /// @param compare Function to use for sorting volumes
  /// @return Reference to this node for method chaining
  GeometryIdentifierBlueprintNode& sortBy(const CompareVolumes& compare);

  /// @brief Get the name of this node
  /// @return String containing concatenated configuration names
  const std::string& name() const override;

 private:
  std::unique_ptr<GeometryIdentifierBlueprintNodeImpl> m_impl;
};

}  // namespace Experimental
}  // namespace Acts
