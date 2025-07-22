// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/ContainerBlueprintNode.hpp"

#include "Acts/Geometry/CuboidPortalShell.hpp"
#include "Acts/Geometry/CuboidVolumeStack.hpp"
#include "Acts/Geometry/CylinderPortalShell.hpp"
#include "Acts/Geometry/CylinderVolumeStack.hpp"

namespace Acts::Experimental {

ContainerBlueprintNode::ContainerBlueprintNode(
    const std::string& name, AxisDirection axis,
    VolumeAttachmentStrategy attachmentStrategy,
    VolumeResizeStrategy resizeStrategy)
    : ContainerBlueprintNode(name, axis, attachmentStrategy,
                             {resizeStrategy, resizeStrategy}) {}

ContainerBlueprintNode::ContainerBlueprintNode(
    const std::string& name, AxisDirection axis,
    VolumeAttachmentStrategy attachmentStrategy,
    std::pair<VolumeResizeStrategy, VolumeResizeStrategy> resizeStrategies)
    : m_name(name),
      m_direction(axis),
      m_attachmentStrategy(attachmentStrategy),
      m_resizeStrategies(resizeStrategies) {}

const std::string& ContainerBlueprintNode::name() const {
  return m_name;
}

Volume& ContainerBlueprintNode::build(
    const Experimental::BlueprintOptions& options, const GeometryContext& gctx,
    const Logger& logger) {
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
  m_stack = makeStack(m_childVolumes, logger);
  ACTS_DEBUG(prefix() << "-> Stack bounds are: " << m_stack->volumeBounds());

  ACTS_DEBUG(prefix() << " *** build complete ***");

  return *m_stack;
}

void ContainerBlueprintNode::finalize(
    const Experimental::BlueprintOptions& options, const GeometryContext& gctx,
    TrackingVolume& parent, const Logger& logger) {
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

ContainerBlueprintNode& ContainerBlueprintNode::setDirection(
    AxisDirection direction) {
  if (m_stack != nullptr) {
    throw std::runtime_error("Cannot change direction after build");
  }
  m_direction = direction;
  return *this;
}

ContainerBlueprintNode& ContainerBlueprintNode::setAttachmentStrategy(
    VolumeAttachmentStrategy attachmentStrategy) {
  if (m_stack != nullptr) {
    throw std::runtime_error("Cannot change direction after build");
  }
  m_attachmentStrategy = attachmentStrategy;
  return *this;
}

ContainerBlueprintNode& ContainerBlueprintNode::setResizeStrategy(
    VolumeResizeStrategy resizeStrategy) {
  if (m_stack != nullptr) {
    throw std::runtime_error("Cannot change direction after build");
  }
  m_resizeStrategies = {resizeStrategy, resizeStrategy};
  return *this;
}

ContainerBlueprintNode& ContainerBlueprintNode::setResizeStrategies(
    VolumeResizeStrategy inner, VolumeResizeStrategy outer) {
  if (m_stack != nullptr) {
    throw std::runtime_error("Cannot change direction after build");
  }
  m_resizeStrategies = {inner, outer};
  return *this;
}

AxisDirection ContainerBlueprintNode::direction() const {
  return m_direction;
}

VolumeAttachmentStrategy ContainerBlueprintNode::attachmentStrategy() const {
  return m_attachmentStrategy;
}

VolumeResizeStrategy ContainerBlueprintNode::resizeStrategy() const {
  if (m_resizeStrategies.first != m_resizeStrategies.second) {
    throw std::runtime_error(
        "Resize strategy is not the same for inner and outer. Use "
        "resizeStrategies() instead.");
  }
  return m_resizeStrategies.first;
}

std::pair<VolumeResizeStrategy, VolumeResizeStrategy>
ContainerBlueprintNode::resizeStrategies() const {
  return m_resizeStrategies;
}

void ContainerBlueprintNode::addToGraphviz(std::ostream& os) const {
  std::stringstream ss;
  ss << "<b>" + name() + "</b>";
  ss << "<br/>" << typeName() << "Container";
  ss << "<br/>dir: " << m_direction;
  GraphViz::Node node{
      .id = name(), .label = ss.str(), .shape = GraphViz::Shape::DoubleOctagon};
  os << node << std::endl;
  for (const auto& child : children()) {
    os << indent() << GraphViz::Edge{{.id = name()}, {.id = child.name()}}
       << std::endl;
    child.addToGraphviz(os);
  }
}

template <typename BaseShell, typename SingleShell>
std::vector<BaseShell*> ContainerBlueprintNode::collectChildShells(
    const Experimental::BlueprintOptions& options, const GeometryContext& gctx,
    VolumeStack& stack, const std::string& prefix, const Logger& logger) {
  std::vector<BaseShell*> shells;
  ACTS_DEBUG(prefix << "Have " << m_childVolumes.size() << " child volumes");
  std::size_t nGaps = 0;
  for (Volume* volume : m_childVolumes) {
    if (stack.isGapVolume(*volume)) {
      // We need to create a TrackingVolume from the gap and put it in the
      // shell
      auto gap = std::make_unique<TrackingVolume>(*volume);
      gap->setVolumeName(name() + "::Gap" + std::to_string(nGaps + 1));
      nGaps++;
      ACTS_DEBUG(prefix << " ~> Gap volume (" << gap->volumeName()
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

      ACTS_DEBUG(prefix << " ~> Child (" << child.name()
                        << ") volume: " << volume->volumeBounds());

      auto* shell =
          dynamic_cast<BaseShell*>(&child.connect(options, gctx, logger));
      if (shell == nullptr) {
        ACTS_ERROR(prefix << "Child volume stack type mismatch");
        throw std::runtime_error("Child volume stack type mismatch");
      }
      assert(shell->isValid());

      shells.push_back(shell);
    }
  }
  return shells;
}

template <typename BaseShell, typename SingleShell, typename ShellStack>
PortalShellBase& ContainerBlueprintNode::connectImpl(
    const Experimental::BlueprintOptions& options, const GeometryContext& gctx,
    VolumeStack* stack, const std::string& prefix, const Logger& logger) {
  ACTS_DEBUG(prefix << "Container connect");
  if (stack == nullptr) {
    ACTS_ERROR(prefix << "Volume is not built");
    throw std::runtime_error("Volume is not built");
  }
  ACTS_DEBUG(prefix << "Collecting child shells from " << children().size()
                    << " children");

  // We have child volumes and gaps as bare Volumes in `m_childVolumes` after
  // `build()` has completed. For the stack shell, we need TrackingVolumes in
  // the right order.

  std::vector<BaseShell*> shells = collectChildShells<BaseShell, SingleShell>(
      options, gctx, *stack, prefix, logger);

  // Sanity checks
  throw_assert(shells.size() == m_childVolumes.size(),
               "Number of shells does not match number of child volumes");

  throw_assert(std::ranges::none_of(
                   shells, [](const auto* shell) { return shell == nullptr; }),
               "Invalid shell pointer");

  throw_assert(std::ranges::all_of(
                   shells, [](const auto* shell) { return shell->isValid(); }),
               "Invalid shell");

  ACTS_DEBUG(prefix << "Producing merged stack shell in " << direction()
                    << " direction from " << shells.size() << " shells");
  m_shell = std::make_unique<ShellStack>(gctx, std::move(shells), direction(),
                                         logger);

  assert(m_shell != nullptr && "No shell was built at the end of connect");
  assert(m_shell->isValid() && "Shell is not valid at the end of connect");
  return *m_shell;
}

PortalShellBase& CylinderContainerBlueprintNode::connect(
    const Experimental::BlueprintOptions& options, const GeometryContext& gctx,
    const Logger& logger) {
  return connectImpl<CylinderPortalShell, SingleCylinderPortalShell,
                     CylinderStackPortalShell>(options, gctx, m_stack.get(),
                                               prefix(), logger);
}

const std::string& CylinderContainerBlueprintNode::typeName() const {
  return s_typeName;
}

std::unique_ptr<VolumeStack> CylinderContainerBlueprintNode::makeStack(
    std::vector<Volume*>& volumes, const Logger& logger) {
  return std::make_unique<CylinderVolumeStack>(
      volumes, m_direction, m_attachmentStrategy, m_resizeStrategies, logger);
}

PortalShellBase& CuboidContainerBlueprintNode::connect(
    const Experimental::BlueprintOptions& options, const GeometryContext& gctx,
    const Logger& logger) {
  return connectImpl<CuboidPortalShell, SingleCuboidPortalShell,
                     CuboidStackPortalShell>(options, gctx, m_stack.get(),
                                             prefix(), logger);
}

const std::string& CuboidContainerBlueprintNode::typeName() const {
  return s_typeName;
}

std::unique_ptr<VolumeStack> CuboidContainerBlueprintNode::makeStack(
    std::vector<Volume*>& volumes, const Logger& logger) {
  return std::make_unique<CuboidVolumeStack>(volumes, m_direction,
                                             m_attachmentStrategy,
                                             m_resizeStrategies.first, logger);
}

}  // namespace Acts::Experimental
