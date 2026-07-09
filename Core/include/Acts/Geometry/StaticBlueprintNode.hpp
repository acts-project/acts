// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/PortalShell.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"

namespace Acts::Experimental {

/// The static blueprint node wraps a single already-constructred @c TrackingVolume.
/// The node will present this volume to its hierarchy. The volume is given as
/// mutable, and will be potentially enlarged in order to connect to neighboring
/// volumes.
/// - In case the volume already has child volumes, they will be retained.
/// - In case the volume already has a registered navigation policy, it will be
///   overwritten with the one configured on this node, regardless of content.
class StaticBlueprintNode : public BlueprintNode {
 public:
  /// Construct the static node from an existing volume
  /// @param volume The volume to wrap
  explicit StaticBlueprintNode(std::unique_ptr<TrackingVolume> volume);

  /// Get the name of this node. It is automatically taken from the wrapped
  /// volume
  /// @return The name of the volume
  const std::string& name() const override;

  /// @copydoc BlueprintNode::build
  /// Build-phase of the blueprint construction. Returns the wrapped volume for
  /// sizing.
  Volume& build(const BlueprintOptions& options, const GeometryContext& gctx,
                const Logger& logger = Acts::getDummyLogger()) override;

  /// @copydoc BlueprintNode::connect
  PortalShellBase& connect(
      const BlueprintOptions& options, const GeometryContext& gctx,
      const Logger& logger = Acts::getDummyLogger()) override;

  /// @copydoc BlueprintNode::finalize
  void finalize(const BlueprintOptions& options, const GeometryContext& gctx,
                TrackingVolume& parent,
                const Logger& logger = Acts::getDummyLogger()) override;

  /// Set the navigation policy factory for this node
  /// @param navigationPolicyFactory Shared pointer to navigation policy factory
  /// @return Reference to this node for chaining
  virtual StaticBlueprintNode& setNavigationPolicyFactory(
      std::shared_ptr<NavigationPolicyFactory> navigationPolicyFactory);

  /// Release ownership of the wrapped tracking volume
  /// @return The tracking volume as a unique_ptr (m_volume becomes nullptr)
  std::unique_ptr<TrackingVolume> releaseVolume() {
    return std::move(m_volume);
  }

  /// Get the navigation policy factory for this node
  /// @return Pointer to the navigation policy factory (may be nullptr)
  const NavigationPolicyFactory* navigationPolicyFactory() const;

  /// Take over the ownership over a Surface or VolumePlacement and pass it on
  /// to the tracking geometry
  /// @tparam Obj_t Either the PlacementOwnPtr variant or any other pointer
  ///               type where the object inhherits from the Surface or
  ///               VolumePlacement base class
  /// @param placement Pointer to the placement to be managed by the
  ///                   tracking volume
  template <typename Obj_t>
  void retainPlacement(Obj_t placement)
    requires(std::is_constructible_v<TrackingVolume::PlacementOwnPtr, Obj_t>)
  {
    m_placements.emplace_back(std::move(placement));
  }

 protected:
  void addToGraphviz(std::ostream& os) const override;

  /// The wrapped tracking volume managed by this blueprint node
  std::unique_ptr<TrackingVolume> m_volume;
  /// Vector of volume or surface plaements to be owned by the tracking geometry
  std::vector<TrackingVolume::PlacementOwnPtr> m_placements{};

  /// Portal shell representation for geometry connection
  std::unique_ptr<PortalShellBase> m_shell;

  /// Factory for creating navigation policies for this volume
  std::shared_ptr<NavigationPolicyFactory> m_navigationPolicyFactory;
};

}  // namespace Acts::Experimental
