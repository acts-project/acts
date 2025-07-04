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

/// @class ContainerBlueprintNode
///
/// A blueprint node that can contain multiple child volumes. It is responsible
/// for managing the child volumes and their shells. The child volumes can be
/// either gap volumes or volumes from child nodes.
///
/// The container node is responsible for:
/// 1. Managing the child volumes and their shells
/// 2. Creating gap volumes between child volumes
/// 3. Collecting shells from child nodes and gap volumes
/// 4. Building the volume stack
///
/// The container node is an abstract base class. Derived classes must
/// implement:
/// 1. makeStack - to create the appropriate volume stack
/// 2. typeName - to provide the type name for debug output
///
class ContainerBlueprintNode : public BlueprintNode {
 public:
  /// Main constructor for the container node.
  /// @param name The name of the node (for debug only)
  /// @param axis The stacking axis direction in local reference frame
  /// @param attachmentStrategy The attachment strategy for the stack
  /// @param resizeStrategy The resize strategy for the stack
  ContainerBlueprintNode(
      const std::string& name, AxisDirection axis,
      VolumeAttachmentStrategy attachmentStrategy =
          VolumeAttachmentStrategy::Midpoint,
      VolumeResizeStrategy resizeStrategy = VolumeResizeStrategy::Expand);

  /// Main constructor for the container node.
  /// @param name The name of the node (for debug only)
  /// @param axis The stacking axis direction in local reference frame
  /// @param attachmentStrategy The attachment strategy for the stack
  /// @param resizeStrategies The resize strategies for the stack
  ContainerBlueprintNode(
      const std::string& name, AxisDirection axis,
      VolumeAttachmentStrategy attachmentStrategy,
      std::pair<VolumeResizeStrategy, VolumeResizeStrategy> resizeStrategies);

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
  /// @note @p resizeStrategy is used for both sides of the container
  /// @return This node for chaining
  ContainerBlueprintNode& setResizeStrategy(
      VolumeResizeStrategy resizeStrategy);

  /// Setter for the resize strategies
  /// @param inner The inner resize strategy
  /// @param outer The outer resize strategy
  /// @return This node for chaining
  ContainerBlueprintNode& setResizeStrategies(VolumeResizeStrategy inner,
                                              VolumeResizeStrategy outer);

  /// Accessor to the stacking direction
  /// @return The stacking direction
  AxisDirection direction() const;

  /// Accessor to the attachment strategy
  /// @return The attachment strategy
  VolumeAttachmentStrategy attachmentStrategy() const;

  /// Accessor to the resize strategy
  /// @return The resize strategy
  [[deprecated("Use resizeStrategies() instead")]]
  VolumeResizeStrategy resizeStrategy() const;

  /// Accessor to the resize strategies
  /// @return The resize strategies
  std::pair<VolumeResizeStrategy, VolumeResizeStrategy> resizeStrategies()
      const;

  /// @copydoc BlueprintNode::addToGraphviz
  void addToGraphviz(std::ostream& os) const override;

 protected:
  /// Make the volume stack for the container. This is called by the build
  /// method and is implemented by the derived classes.
  /// @param volumes The volumes to stack
  /// @param logger The logger to use
  /// @return The volume stack
  virtual std::unique_ptr<VolumeStack> makeStack(std::vector<Volume*>& volumes,
                                                 const Logger& logger) = 0;

  /// Get the type name of the container. This is used for the debug output
  /// of the container and encoding the volume shape in the dot graph.
  /// @return The type name
  virtual const std::string& typeName() const = 0;

  /// Collect shells from child nodes and gap volumes
  ///
  /// This function is responsible for collecting shells from child nodes and
  /// creating shells for gap volumes. It is used by the connect method to
  /// prepare the shells for the volume stack.
  ///
  /// The function processes each volume in m_childVolumes in two ways:
  /// 1. For gap volumes:
  ///    - Creates a TrackingVolume from the gap volume
  ///    - Assigns a unique name (ContainerName::GapN)
  ///    - Creates a single shell for the gap volume
  ///    - Stores both the shell and gap volume in m_gaps for later use
  ///
  /// 2. For child volumes:
  ///    - Looks up the corresponding child node in m_volumeToNode
  ///    - Calls connect() on the child node to get its shell
  ///    - Validates that the shell type matches the expected type
  ///    - Ensures the shell is valid
  ///
  /// The function maintains the order of volumes as they appear in
  /// m_childVolumes, which is important for the final stack shell construction.
  ///
  /// @tparam BaseShell The base shell type (e.g. CylinderPortalShell)
  /// @tparam SingleShell The single shell type (e.g. SingleCylinderPortalShell)
  /// @param options The blueprint options
  /// @param gctx The geometry context
  /// @param stack The volume stack
  /// @param prefix The prefix for debug output
  /// @param logger The logger to use
  /// @return A vector of shells in the same order as m_childVolumes
  template <typename BaseShell, typename SingleShell>
  std::vector<BaseShell*> collectChildShells(
      const Experimental::BlueprintOptions& options,
      const GeometryContext& gctx, VolumeStack& stack,
      const std::string& prefix, const Logger& logger);

  /// Implementation of the connect method for container nodes
  ///
  /// This method is responsible for:
  /// 1. Collecting shells from child nodes and gap volumes
  /// 2. Validating that the number of shells matches the number of child
  /// volumes
  /// 3. Ensuring all shells are valid
  /// 4. Creating a merged stack shell from all collected shells
  ///
  /// @tparam BaseShell The base shell type (e.g. CylinderPortalShell)
  /// @tparam SingleShell The single shell type (e.g. SingleCylinderPortalShell)
  /// @tparam ShellStack The stack shell type (e.g. StackCylinderPortalShell)
  /// @param options The blueprint options
  /// @param gctx The geometry context
  /// @param stack The volume stack
  /// @param prefix The prefix for debug output
  /// @param logger The logger to use
  /// @return The merged stack shell
  template <typename BaseShell, typename SingleShell, typename ShellStack>
  PortalShellBase& connectImpl(const Experimental::BlueprintOptions& options,
                               const GeometryContext& gctx, VolumeStack* stack,
                               const std::string& prefix, const Logger& logger);

  std::string m_name;
  AxisDirection m_direction = AxisDirection::AxisZ;
  VolumeAttachmentStrategy m_attachmentStrategy{
      VolumeAttachmentStrategy::Midpoint};

  std::pair<VolumeResizeStrategy, VolumeResizeStrategy> m_resizeStrategies{
      VolumeResizeStrategy::Expand, VolumeResizeStrategy::Expand};

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
  inline static const std::string s_typeName = "Cylinder";
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
  inline static const std::string s_typeName = "Cuboid";
  const std::string& typeName() const override;
};

}  // namespace Acts::Experimental
