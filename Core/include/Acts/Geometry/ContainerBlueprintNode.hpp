// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/BlueprintOptions.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Geometry/VolumeResizeStrategy.hpp"
#include "Acts/Geometry/VolumeStack.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/GraphViz.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

#include <map>

namespace Acts::Experimental {

/// This class handles the case of wrapping a set of children
/// and stacking them in a configured direction.
/// The container does not result in an extra volume in the hierarchy, as all
/// input volumes and any gap volumes produced are directly registered in the
/// volume of the parent of this node.
class ContainerBlueprintNode : public BlueprintNode {
 public:
  /// Main constructor for the container node.
  /// @param name The name of the node (for debug only)
  /// @param axis The stacking axis direction in local reference frame
  /// @param attachmentStrategy The attachment strategy for the stack
  /// @param resizeStrategy The resize strategy
  ContainerBlueprintNode(
      const std::string& name, AxisDirection axis,
      VolumeAttachmentStrategy attachmentStrategy =
          VolumeAttachmentStrategy::Midpoint,
      VolumeResizeStrategy resizeStrategy = VolumeResizeStrategy::Expand);

  /// @copydoc BlueprintNode::name
  const std::string& name() const override;

  /// This participates in the construction of the geometry via the blueprint
  /// tree. The steps are approximately as follows:
  /// -# Collect all child volumes
  /// -# Package them into a VolumeStack (cuboid or cylinder), which performs
  ///    sizing and/or gap creation
  /// -# Return the VolumeStack as a volume up the tree
  ///
  /// @param options The global blueprint options
  /// @param gctx The geometry context (nominal usually)
  /// @param logger The logger to use
  /// @return The combined VolumeStack
  Volume& build(const Experimental::BlueprintOptions& options,
                const GeometryContext& gctx,
                const Logger& logger = Acts::getDummyLogger()) override;

  /// This participates in the construction of the geometry via the blueprint
  /// tree. The steps are approximately as follows:
  /// -# Register portals created for gap volumes, as they're not handled by
  ///    dedicated nodes
  /// -# Register gap volumes in the @p parent volume
  /// -# Create a configured @ref Acts::INavigationPolicy for the gap
  /// -# Call `finalize` on all children while passing through @p parent.
  ///
  /// @param options The global blueprint options
  /// @param gctx The geometry context (nominal usually)
  /// @param parent The parent volume
  /// @param logger The logger to use
  void finalize(const Experimental::BlueprintOptions& options,
                const GeometryContext& gctx, TrackingVolume& parent,
                const Logger& logger) override;

  /// Setter for the stacking direction
  /// @param direction The stacking direction
  /// @return This node for chaining
  ContainerBlueprintNode& setDirection(AxisDirection direction);

  /// Setter for the attachment strategy
  /// @param attachmentStrategy The attachment strategy
  /// @return This node for chaining
  ContainerBlueprintNode& setAttachmentStrategy(
      VolumeAttachmentStrategy attachmentStrategy);

  /// Setter for the resize strategy
  /// @param resizeStrategy The resize strategy
  /// @return This node for chaining
  ContainerBlueprintNode& setResizeStrategy(
      VolumeResizeStrategy resizeStrategy);

  /// Accessor to the stacking direction
  /// @return The stacking direction
  AxisDirection direction() const;

  /// Accessor to the attachment strategy
  /// @return The attachment strategy
  VolumeAttachmentStrategy attachmentStrategy() const;

  /// Accessor to the resize strategy
  /// @return The resize strategy
  VolumeResizeStrategy resizeStrategy() const;

 private:
  /// @copydoc BlueprintNode::addToGraphviz
  void addToGraphviz(std::ostream& os) const override;

 protected:
  virtual std::unique_ptr<VolumeStack> makeStack(std::vector<Volume*>& volumes,
                                                 const Logger& logger) = 0;

  virtual const std::string& typeName() const = 0;

  template <typename BaseShell, typename SingleShell>
  std::vector<BaseShell*> collectChildShells(
      const Experimental::BlueprintOptions& options,
      const GeometryContext& gctx, VolumeStack& stack,
      const std::string& prefix, const Logger& logger);

  template <typename BaseShell, typename SingleShell, typename ShellStack>
  PortalShellBase& connectImpl(const Experimental::BlueprintOptions& options,
                               const GeometryContext& gctx, VolumeStack* stack,
                               const std::string& prefix, const Logger& logger);

  std::string m_name;
  AxisDirection m_direction = AxisDirection::AxisZ;
  VolumeAttachmentStrategy m_attachmentStrategy{
      VolumeAttachmentStrategy::Midpoint};
  VolumeResizeStrategy m_resizeStrategy{VolumeResizeStrategy::Expand};

  std::vector<Volume*> m_childVolumes;
  // This is going to be an instance of a *stack* of volumes, which is created
  // by the derived classes
  std::unique_ptr<VolumeStack> m_stack{nullptr};
  std::map<const Volume*, BlueprintNode*> m_volumeToNode;

  std::unique_ptr<PortalShellBase> m_shell{nullptr};
  std::vector<std::pair<std::unique_ptr<PortalShellBase>,
                        std::unique_ptr<TrackingVolume>>>
      m_gaps;
};

class CylinderContainerBlueprintNode final : public ContainerBlueprintNode {
 public:
  using ContainerBlueprintNode::ContainerBlueprintNode;

  /// This participates in the construction of the geometry via the blueprint
  /// tree. The steps are approximately as follows:
  /// -# Walk through all volumes that were created by the build phase
  /// -# Check if they are: *real* child volumes or gap volumes
  ///   - If gap volume: produce a @ref Acts::TrackingVolume, and wrap it in a single use shell
  ///   - If child volume: locate the right child node it came from, call
  ///   ` connect` and collect the returned shell
  /// -# Produce a combined StackPortalShell (cuboid or cylinder) from all the
  /// shells
  /// -# Return that shell representation
  ///
  /// @param options The global blueprint options
  /// @param gctx The geometry context (nominal usually)
  /// @param logger The logger to use
  /// @return The combined StackPortalShell (cuboid or cylinder)
  PortalShellBase& connect(
      const Experimental::BlueprintOptions& options,
      const GeometryContext& gctx,
      const Logger& logger = Acts::getDummyLogger()) override;

  std::unique_ptr<VolumeStack> makeStack(std::vector<Volume*>& volumes,
                                         const Logger& logger) override;

 protected:
  const std::string& typeName() const override;
};

class CuboidContainerBlueprintNode final : public ContainerBlueprintNode {
 public:
  using ContainerBlueprintNode::ContainerBlueprintNode;

  /// This participates in the construction of the geometry via the blueprint
  /// tree. The steps are approximately as follows:
  /// -# Walk through all volumes that were created by the build phase
  /// -# Check if they are: *real* child volumes or gap volumes
  ///   - If gap volume: produce a @ref Acts::TrackingVolume, and wrap it in a single use shell
  ///   - If child volume: locate the right child node it came from, call
  ///   ` connect` and collect the returned shell
  /// -# Produce a combined StackPortalShell (cuboid or cylinder) from all the
  /// shells
  /// -# Return that shell representation
  ///
  /// @param options The global blueprint options
  /// @param gctx The geometry context (nominal usually)
  /// @param logger The logger to use
  /// @return The combined StackPortalShell (cuboid or cylinder)
  PortalShellBase& connect(
      const Experimental::BlueprintOptions& options,
      const GeometryContext& gctx,
      const Logger& logger = Acts::getDummyLogger()) override;

  std::unique_ptr<VolumeStack> makeStack(std::vector<Volume*>& volumes,
                                         const Logger& logger) override;

 protected:
  const std::string& typeName() const override;
};

}  // namespace Acts::Experimental
