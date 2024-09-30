// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/RootBlueprintNode.hpp"

#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/PortalShell.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Utilities/GraphViz.hpp"

namespace Acts {

RootBlueprintNode::RootBlueprintNode(const Config &cfg) : m_cfg(cfg) {}

const std::string &RootBlueprintNode::name() const {
  static const std::string root = "root";
  return root;
}

Volume &RootBlueprintNode::build(const Options & /*options*/,
                                 const Logger & /*logger*/) {
  throw std::logic_error("Root node cannot be built");
}

PortalShellBase &RootBlueprintNode::connect(const Options & /*options*/,
                                            const GeometryContext & /*gctx*/,
                                            const Logger & /*logger*/) {
  throw std::logic_error("Root node cannot be connected");
}

void RootBlueprintNode::finalize(const Options & /*options*/,
                                 TrackingVolume & /*parent*/,
                                 const Logger & /*logger*/) {
  throw std::logic_error("Root node cannot be finalized");
}

void RootBlueprintNode::addToGraphviz(std::ostream &os) const {
  GraphViz::Node node{
      .id = name(), .label = "Root", .shape = GraphViz::Shape::House};

  os << node;
  BlueprintNode::addToGraphviz(os);
}

std::unique_ptr<TrackingGeometry> RootBlueprintNode::construct(
    const Options &options, const GeometryContext &gctx, const Logger &logger) {
  using enum BinningValue;

  ACTS_INFO(prefix() << "Building tracking geometry from blueprint tree");

  options.validate();

  if (m_cfg.envelope == ExtentEnvelope::Zero()) {
    ACTS_WARNING(prefix() << "Root node is configured with zero envelope. This "
                             "might lead to navigation issues");
  }

  if (children().size() != 1) {
    ACTS_ERROR(prefix() << "Root node must have exactly one child");
    throw std::logic_error("Root node must have exactly one child");
  }

  auto &child = children().at(0);

  ACTS_DEBUG(prefix() << "Executing building on tree");
  Volume &topVolume = child.build(options, logger);
  const auto &bounds = topVolume.volumeBounds();

  std::stringstream ss;
  bounds.toStream(ss);
  ACTS_DEBUG(prefix() << "have top volume: " << ss.str() << "\n"
                      << topVolume.transform().matrix());

  std::shared_ptr<VolumeBounds> worldBounds;

  if (const auto *cyl = dynamic_cast<const CylinderVolumeBounds *>(&bounds)) {
    using enum CylinderVolumeBounds::BoundValues;

    // Make a copy that we'll modify
    auto newBounds = std::make_shared<CylinderVolumeBounds>(*cyl);

    const auto &zEnv = m_cfg.envelope[binZ];
    if (zEnv[0] != zEnv[1]) {
      ACTS_ERROR(
          prefix() << "Root node cylinder envelope for z must be symmetric");
    }

    const auto &rEnv = m_cfg.envelope[binR];

    newBounds->set({
        {eHalfLengthZ, newBounds->get(eHalfLengthZ) + zEnv[0]},
        {eMinR, std::max(0.0, newBounds->get(eMinR) - rEnv[0])},
        {eMaxR, newBounds->get(eMaxR) + rEnv[1]},
    });

    worldBounds = std::move(newBounds);

  } else if (const auto *box =
                 dynamic_cast<const CuboidVolumeBounds *>(&bounds)) {
    throw std::logic_error{"Not implemented"};
  } else {
    throw std::logic_error{"Unsupported volume bounds type"};
  }

  ACTS_DEBUG(prefix() << "New root volume bounds are: " << *worldBounds);

  auto world = std::make_unique<TrackingVolume>(
      topVolume.transform(), std::move(worldBounds), "World");

  // @TODO: This needs to become configurable
  world->setNavigationPolicy(
      options.defaultNavigationPolicyFactory->build(*world));

  // Need one-sided portal shell that connects outwards to nullptr
  SingleCylinderPortalShell worldShell{*world};
  worldShell.applyToVolume();

  auto &shell = child.connect(options, gctx, logger);

  shell.connectOuter(*world);

  child.finalize(options, *world, logger);

  // @TODO: Handle material decorator, geo id hook

  return std::make_unique<TrackingGeometry>(
      std::move(world), nullptr, m_cfg.geometryIdentifierHook, logger);
}

}  // namespace Acts
