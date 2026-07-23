// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Blueprint.hpp"

#include "Acts/Geometry/CuboidPortalShell.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderPortalShell.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/PadBlueprintNode.hpp"
#include "Acts/Geometry/PortalShell.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Geometry/detail/AlignablePortalVisitor.hpp"
#include "Acts/Geometry/detail/BoundDeduplicator.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "Acts/Utilities/GraphViz.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>

namespace {
const std::string s_rootName = "Root";
}

namespace Acts::Experimental {

///@class BlueprintVisitor
/// A class for visiting blueprint hierarchy and apply the geometry identifiers
class BlueprintVisitor : public TrackingGeometryMutableVisitor {
 public:
  explicit BlueprintVisitor(
      const Logger &logger,
      std::array<const TrackingVolume *, GeometryIdentifier::getMaxVolume()>
          &volumesById)
      : TrackingGeometryMutableVisitor(true),
        m_volumesById(volumesById),
        m_logger(logger) {}

  void visitVolume(TrackingVolume &volume) override {
    GeometryIdentifier::Value iportal = 0;
    GeometryIdentifier::Value isensitive = 0;

    auto id = volume.geometryId();

    if (id == GeometryIdentifier{}) {
      auto it = std::ranges::find(m_volumesById, nullptr);
      if (it == m_volumesById.end()) {
        ACTS_ERROR("No free volume IDs left, all " << m_volumesById.size()
                                                   << " are used");
        // @TODO: Maybe link to documentation about this
        throw std::logic_error("No free volume IDs left");
      }

      id = GeometryIdentifier().withVolume(
          std::distance(m_volumesById.begin(), it) + 1);

      ACTS_VERBOSE("Assigning volume ID " << id << " for "
                                          << volume.volumeName());
      volume.assignGeometryId(id);
      *it = &volume;
    }

    for (auto &portal : volume.portals()) {
      if (portal.surface().geometryId() != GeometryIdentifier{}) {
        continue;
      }
      iportal += 1;
      auto portalId = id.withBoundary(iportal);
      ACTS_VERBOSE("Assigning portal ID: " << portalId);
      portal.surface().assignGeometryId(portalId);
    }
    for (auto &surface : volume.surfaces()) {
      if (surface.geometryId() != GeometryIdentifier{}) {
        continue;
      }
      isensitive += 1;
      auto surfaceId = id.withSensitive(isensitive);
      ACTS_VERBOSE("Assigning surface ID: " << surfaceId);
      surface.assignGeometryId(surfaceId);
    }
  }

 private:
  std::array<const TrackingVolume *, GeometryIdentifier::getMaxVolume()>
      &m_volumesById;
  const Logger &m_logger;
  const Acts::Logger &logger() const { return m_logger; }
};

Blueprint::Blueprint(const Config &config) : m_cfg(config) {}

const std::string &Blueprint::name() const {
  return s_rootName;
}

Volume &Blueprint::build(const BlueprintOptions & /*options*/,
                         const GeometryContext & /*gctx*/,
                         const Logger & /*logger*/) {
  throw std::logic_error("Root node cannot be built");
}

PortalShellBase &Blueprint::connect(const BlueprintOptions & /*options*/,
                                    const GeometryContext & /*gctx*/,
                                    const Logger & /*logger*/) {
  throw std::logic_error("Root node cannot be connected");
}

void Blueprint::finalize(const BlueprintOptions & /*options*/,
                         const GeometryContext & /*gctx*/,
                         TrackingVolume & /*parent*/,
                         const Logger & /*logger*/) {
  throw std::logic_error("Root node cannot be finalized");
}

void Blueprint::addToGraphviz(std::ostream &os) const {
  GraphViz::Node node{
      .id = name(), .label = "World", .shape = GraphViz::Shape::House};

  os << node;
  BlueprintNode::addToGraphviz(os);
}

std::unique_ptr<TrackingGeometry> Blueprint::construct(
    const BlueprintOptions &options, const GeometryContext &gctx,
    const Logger &logger) {
  using enum AxisDirection;

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

  auto child = childPtr()[0];
  clearChildren();

  ACTS_DEBUG(prefix() << "Executing building on tree");

  auto autoSizingNode = std::make_shared<Acts::Experimental::PadBlueprintNode>(
      "World", m_cfg.envelope);

  autoSizingNode->addChild(std::move(child));

  autoSizingNode->build(options, gctx, logger);

  auto world = autoSizingNode->trackingVolume();

  ACTS_DEBUG(prefix() << "New root volume bounds are: "
                      << world->volumeBounds());

  world->setNavigationPolicy(std::make_unique<Acts::TryAllNavigationPolicy>(
      gctx, *world, logger, Acts::TryAllNavigationPolicy::Config{}));

  autoSizingNode->connect(options, gctx, logger);

  if (m_cfg.boundDeduplication) {
    ACTS_DEBUG("Deduplicate equivalent bounds");
    detail::BoundDeduplicator deduplicator{};
    world->apply(deduplicator);
  }
  autoSizingNode->finalize(options, gctx, *world, logger);

  std::set<std::string, std::less<>> volumeNames;
  std::array<const TrackingVolume *, GeometryIdentifier::getMaxVolume()>
      volumesById{};
  volumesById.fill(nullptr);

  // @TODO: Take this from GeometryIdentifier instead of hard-coding

  world->apply([&, this](TrackingVolume &volume) {
    if (volumeNames.contains(volume.volumeName())) {
      ACTS_ERROR(prefix() << "Duplicate volume name: " << volume.volumeName());
      throw std::logic_error("Duplicate volume name");
    }
    volumeNames.insert(volume.volumeName());

    if (volume.geometryId() != GeometryIdentifier{}) {
      // We can have multiple volumes with the same volume ID component, but
      // they should differ in other components like "layer"
      if (volumesById.at(volume.geometryId().volume() - 1) == nullptr) {
        volumesById.at(volume.geometryId().volume() - 1) = &volume;
      }
    }

    // Clear boundary surfaces!
    volume.clearBoundarySurfaces();
  });

  std::size_t unusedVolumeIds = std::ranges::count(volumesById, nullptr);
  ACTS_DEBUG(prefix() << "Number of unused volume IDs: " << unusedVolumeIds);

  ACTS_DEBUG(prefix() << "Assigning volume IDs for remaining volumes");

  BlueprintVisitor visitor{logger, volumesById};
  world->apply(visitor);

  Acts::detail::AlignablePortalVisitor alignPortals{gctx, logger};
  world->apply(alignPortals);

  // All navigation policies are attached at this point: initialize each one's
  // statelessness cache (probing whether it pushes only default states) so the
  // navigator can skip the per-volume-entry state creation for volumes with
  // stateless policies.
  world->apply([&](TrackingVolume &volume) {
    if (INavigationPolicy *policy = volume.navigationPolicy();
        policy != nullptr) {
      policy->initializeStatelessCache(gctx, logger);
    }
  });

  return std::make_unique<TrackingGeometry>(
      std::shared_ptr<TrackingVolume>(autoSizingNode->releaseVolume()), nullptr,
      GeometryIdentifierHook{}, logger, false);
}

}  // namespace Acts::Experimental
