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
#include "Acts/Geometry/detail/ContainerBlueprintNodeTraits.hpp"
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
template <class Bounds>
class ContainerBlueprintNode final : public Experimental::BlueprintNode {
 public:
  using VolumeStack =
      typename ContainerBlueprintNodeTraits<Bounds>::VolumeStack;

  using BaseShell = typename ContainerBlueprintNodeTraits<Bounds>::BaseShell;
  using SingleShell =
      typename ContainerBlueprintNodeTraits<Bounds>::SingleShell;
  using ShellStack = typename ContainerBlueprintNodeTraits<Bounds>::ShellStack;

  /// Main constructor for the container node.
  /// @param name The name of the node (for debug only)
  /// @param axis The stacking axis direction in local reference frame
  /// @param attachmentStrategy The attachment strategy for the stack
  /// @param resizeStrategy The resize strategy
  ContainerBlueprintNode(
      const std::string& name, AxisDirection axis,
      VolumeAttachmentStrategy attachmentStrategy =
          VolumeAttachmentStrategy::Midpoint,
      VolumeResizeStrategy resizeStrategy = VolumeResizeStrategy::Expand)
      : m_name(name),
        m_direction(axis),
        m_attachmentStrategy(attachmentStrategy),
        m_resizeStrategy(resizeStrategy) {}

  /// @copydoc BlueprintNode::name
  const std::string& name() const override { return m_name; }

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
                const Logger& logger = Acts::getDummyLogger()) override {
    ACTS_DEBUG(prefix() << "container build (dir=" << m_direction << ")");

    if (m_stack != nullptr) {
      ACTS_ERROR(prefix() << "Volume is already built");
      throw std::runtime_error("Volume is already built");
    }

    for (auto& child : children()) {
      Volume& volume = child.build(options, gctx, logger);
      m_childVolumes.push_back(&volume);
      // We need to remember which volume we got from which child, so we can
      // assemble a crrect portal shell later
      m_volumeToNode[&volume] = &child;
    }
    ACTS_VERBOSE(prefix() << "-> Collected " << m_childVolumes.size()
                          << " child volumes");
    ACTS_VERBOSE(prefix() << "-> Building the stack");
    m_stack = std::make_unique<VolumeStack>(m_childVolumes, m_direction,
                                            m_attachmentStrategy,
                                            m_resizeStrategy, logger);
    ACTS_DEBUG(prefix() << "-> Stack bounds are: " << m_stack->volumeBounds());

    ACTS_DEBUG(prefix() << " *** build complete ***");

    return *m_stack;
  }

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
  ShellStack& connect(const Experimental::BlueprintOptions& options,
                      const GeometryContext& gctx,
                      const Logger& logger = Acts::getDummyLogger()) override {
    ACTS_DEBUG(prefix() << "Container connect");
    if (m_stack == nullptr) {
      ACTS_ERROR(prefix() << "Volume is not built");
      throw std::runtime_error("Volume is not built");
    }

    ACTS_DEBUG(prefix() << "Collecting child shells from " << children().size()
                        << " children");

    // We have child volumes and gaps as bare Volumes in `m_childVolumes` after
    // `build()` has completed. For the stack shell, we need TrackingVolumes in
    // the right order.

    std::vector<BaseShell*> shells = collectChildShells(options, gctx, logger);

    // Sanity checks
    throw_assert(shells.size() == m_childVolumes.size(),
                 "Number of shells does not match number of child volumes");

    throw_assert(
        std::ranges::none_of(
            shells, [](const auto* shell) { return shell == nullptr; }),
        "Invalid shell pointer");

    throw_assert(
        std::ranges::all_of(shells,
                            [](const auto* shell) { return shell->isValid(); }),
        "Invalid shell");

    ACTS_DEBUG(prefix() << "Producing merged stack shell in " << m_direction
                        << " direction from " << shells.size() << " shells");
    m_shell = std::make_unique<ShellStack>(gctx, std::move(shells), m_direction,
                                           logger);

    assert(m_shell != nullptr && "No shell was built at the end of connect");
    assert(m_shell->isValid() && "Shell is not valid at the end of connect");
    return *m_shell;
  }

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
                const Logger& logger) override {
    ACTS_DEBUG(prefix() << "Finalizing container");

    if (m_stack == nullptr) {
      ACTS_ERROR(prefix() << "Volume is not built");
      throw std::runtime_error("Volume is not built");
    }

    if (m_shell == nullptr) {
      ACTS_ERROR(prefix() << "Volume is not connected");
      throw std::runtime_error("Volume is not connected");
    }

    const auto* policyFactory = options.defaultNavigationPolicyFactory.get();

    ACTS_DEBUG(prefix() << "Registering " << m_gaps.size()
                        << " gap volumes with parent");
    for (auto& [shell, gap] : m_gaps) {
      auto* gapPtr = gap.get();
      parent.addVolume(std::move(gap));
      shell->applyToVolume();
      auto policy = policyFactory->build(gctx, *gapPtr, logger);
      gapPtr->setNavigationPolicy(std::move(policy));
    }

    ACTS_DEBUG(prefix() << "Finalizing " << children().size() << " children");

    for (auto& child : children()) {
      child.finalize(options, gctx, parent, logger);
    }
  }

  /// Setter for the stacking direction
  /// @param direction The stacking direction
  /// @return This node for chaining
  ContainerBlueprintNode& setDirection(AxisDirection direction) {
    if (m_stack != nullptr) {
      throw std::runtime_error("Cannot change direction after build");
    }
    m_direction = direction;
    return *this;
  }

  /// Setter for the attachment strategy
  /// @param attachmentStrategy The attachment strategy
  /// @return This node for chaining
  ContainerBlueprintNode& setAttachmentStrategy(
      VolumeAttachmentStrategy attachmentStrategy) {
    if (m_stack != nullptr) {
      throw std::runtime_error("Cannot change direction after build");
    }
    m_attachmentStrategy = attachmentStrategy;
    return *this;
  }

  /// Setter for the resize strategy
  /// @param resizeStrategy The resize strategy
  /// @return This node for chaining
  ContainerBlueprintNode& setResizeStrategy(
      VolumeResizeStrategy resizeStrategy) {
    if (m_stack != nullptr) {
      throw std::runtime_error("Cannot change direction after build");
    }
    m_resizeStrategy = resizeStrategy;
    return *this;
  }

  /// Accessor to the stacking direction
  /// @return The stacking direction
  AxisDirection direction() const { return m_direction; }

  /// Accessor to the attachment strategy
  /// @return The attachment strategy
  VolumeAttachmentStrategy attachmentStrategy() const {
    return m_attachmentStrategy;
  }

  /// Accessor to the resize strategy
  /// @return The resize strategy
  VolumeResizeStrategy resizeStrategy() const { return m_resizeStrategy; }

 private:
  /// @copydoc BlueprintNode::addToGraphviz
  void addToGraphviz(std::ostream& os) const override {
    std::stringstream ss;
    ss << "<b>" + name() + "</b>";
    ss << "<br/>Container";
    ss << "<br/>dir: " << m_direction;
    GraphViz::Node node{.id = name(),
                        .label = ss.str(),
                        .shape = GraphViz::Shape::DoubleOctagon};
    os << node << std::endl;
    for (const auto& child : children()) {
      os << indent() << GraphViz::Edge{{.id = name()}, {.id = child.name()}}
         << std::endl;
      child.addToGraphviz(os);
    }
  }

  /// Helper function to check if a volume was created as a gap volume.
  /// @param volume The volume to check
  /// @return True if the volume is a gap volume, false otherwise
  bool isGapVolume(const Volume& volume) const {
    assert(m_stack != nullptr);
    return std::ranges::any_of(
        m_stack->gaps(), [&](const auto& gap) { return gap.get() == &volume; });
  }

  std::vector<BaseShell*> collectChildShells(
      const Experimental::BlueprintOptions& options,
      const GeometryContext& gctx, const Logger& logger) {
    std::vector<BaseShell*> shells;
    ACTS_DEBUG(prefix() << "Have " << m_childVolumes.size()
                        << " child volumes");
    for (Volume* volume : m_childVolumes) {
      if (isGapVolume(*volume)) {
        // We need to create a TrackingVolume from the gap and put it in the
        // shell
        auto gap = std::make_unique<TrackingVolume>(*volume);
        gap->setVolumeName(name() + "::Gap" +
                           std::to_string(m_gaps.size() + 1));
        ACTS_DEBUG(prefix() << " ~> Gap volume (" << gap->volumeName()
                            << "): " << gap->volumeBounds());
        auto shell = std::make_unique<SingleShell>(*gap);
        assert(shell->isValid());
        shells.push_back(shell.get());

        m_gaps.emplace_back(std::move(shell), std::move(gap));

      } else {
        // Figure out which child we got this volume from
        auto it = m_volumeToNode.find(volume);
        if (it == m_volumeToNode.end()) {
          throw std::runtime_error("Volume not found in child volumes");
        }
        BlueprintNode& child = *it->second;

        ACTS_DEBUG(prefix() << " ~> Child (" << child.name()
                            << ") volume: " << volume->volumeBounds());

        auto* shell =
            dynamic_cast<BaseShell*>(&child.connect(options, gctx, logger));
        if (shell == nullptr) {
          ACTS_ERROR(prefix() << "Child volume stack type mismatch");
          throw std::runtime_error("Child volume stack type mismatch");
        }
        assert(shell->isValid());

        shells.push_back(shell);
      }
    }
    return shells;
  }

  std::string m_name;

  AxisDirection m_direction = AxisDirection::AxisZ;

  VolumeAttachmentStrategy m_attachmentStrategy{
      VolumeAttachmentStrategy::Midpoint};

  VolumeResizeStrategy m_resizeStrategy{VolumeResizeStrategy::Expand};

  // Is only initialized during `build`
  std::vector<Volume*> m_childVolumes;
  std::unique_ptr<VolumeStack> m_stack{nullptr};
  std::map<const Volume*, BlueprintNode*> m_volumeToNode;
  std::vector<
      std::pair<std::unique_ptr<SingleShell>, std::unique_ptr<TrackingVolume>>>
      m_gaps;
  std::unique_ptr<ShellStack> m_shell{nullptr};
};

}  // namespace Acts::Experimental
