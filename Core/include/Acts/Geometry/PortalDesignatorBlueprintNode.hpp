// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"

#include <string_view>

namespace Acts::Experimental {

namespace detail {
class PortalDesignatorBlueprintNodeImpl;
}

/// This node type assigns string tags to specific portal faces of its child
/// volume during the blueprint construction. The tags can be used to look up
/// the corresponding portals from the final @ref Acts::TrackingGeometry, e.g.
/// "the portal connecting the tracker and the calorimeter".
///
/// Tagging happens in the *finalize* phase, after all portal merging and
/// fusing has completed, so the tagged portal is the final, shared portal that
/// ends up in the geometry. The position of this node in the blueprint tree
/// determines which face is meant, which makes the lookup robust against volume
/// subdivision and portal ordering.
/// @note This node can only have a single child. This is not an error during
///       tree building, but during geometry construction.
class PortalDesignatorBlueprintNode final : public BlueprintNode {
 public:
  /// Main constructor for the portal designator node.
  /// @param name The name of the node (for debug only)
  explicit PortalDesignatorBlueprintNode(std::string_view name);

  ~PortalDesignatorBlueprintNode() override;

  /// @copydoc BlueprintNode::name
  const std::string& name() const override;

  /// @copydoc BlueprintNode::toStream
  void toStream(std::ostream& os) const override;

  /// This method participates in the geometry construction.
  /// It checks that this node only has a single child and forwards the call.
  /// @param options The global blueprint options
  /// @param gctx The geometry context (nominal usually)
  /// @param logger The logger to use
  /// @return The child volume
  Volume& build(const BlueprintOptions& options, const GeometryContext& gctx,
                const Logger& logger = Acts::getDummyLogger()) override;

  /// This method participates in the geometry construction.
  /// It captures the populated portal shell from its only child so the tags can
  /// be applied in the finalize phase (after all merging/fusing), and forwards
  /// the call.
  /// @param options The global blueprint options
  /// @param gctx The geometry context (nominal usually)
  /// @param logger The logger to use
  /// @return The portal shell of the child
  PortalShellBase& connect(
      const BlueprintOptions& options, const GeometryContext& gctx,
      const Logger& logger = Acts::getDummyLogger()) override;

  /// This method participates in the geometry construction.
  /// It applies the configured tags to the (now final) portals of the captured
  /// shell, then passes the call through to its only child.
  /// @param options The global blueprint options
  /// @param gctx The geometry context (nominal usually)
  /// @param parent The parent volume
  /// @param logger The logger to use during construction
  void finalize(const BlueprintOptions& options, const GeometryContext& gctx,
                TrackingVolume& parent, const Logger& logger) override;

  /// Tag a cylinder face with a string label.
  /// @note This method can be called multiple times to tag different faces.
  /// @param face The face of the cylinder to tag
  /// @param label The tag to assign to the portal at that face
  /// @return The portal designator node
  /// @note If this node has previously been configured with a different volume
  ///       shape, this will throw an exception during geometry construction.
  PortalDesignatorBlueprintNode& tagFace(CylinderVolumeBounds::Face face,
                                         const std::string& label);

  /// Tag a cuboid face with a string label.
  /// @note This method can be called multiple times to tag different faces.
  /// @param face The face of the cuboid to tag
  /// @param label The tag to assign to the portal at that face
  /// @return The portal designator node
  /// @note If this node has previously been configured with a different volume
  ///       shape, this will throw an exception during geometry construction.
  PortalDesignatorBlueprintNode& tagFace(CuboidVolumeBounds::Face face,
                                         const std::string& label);

 private:
  /// @copydoc BlueprintNode::addToGraphviz
  void addToGraphviz(std::ostream& os) const override;

  detail::PortalDesignatorBlueprintNodeImpl& impl();
  const detail::PortalDesignatorBlueprintNodeImpl& impl() const;

  std::unique_ptr<detail::PortalDesignatorBlueprintNodeImpl> m_impl;
};

}  // namespace Acts::Experimental
