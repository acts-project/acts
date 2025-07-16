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
#include "Acts/Geometry/PortalShell.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
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

      ACTS_DEBUG("Assigning volume ID " << id << " for "
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
      ACTS_DEBUG("Assigning portal ID: " << portalId);
      portal.surface().assignGeometryId(portalId);
    }
    for (auto &surface : volume.surfaces()) {
      if (surface.geometryId() != GeometryIdentifier{}) {
        continue;
      }
      isensitive += 1;
      auto surfaceId = id.withSensitive(isensitive);
      ACTS_DEBUG("Assigning surface ID: " << surfaceId);
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

  auto &child = children().at(0);

  ACTS_DEBUG(prefix() << "Executing building on tree");
  Volume &topVolume = child.build(options, gctx, logger);
  const auto &bounds = topVolume.volumeBounds();

  std::stringstream ss;
  bounds.toStream(ss);
  ACTS_DEBUG(prefix() << "have top volume: " << ss.str() << "\n"
                      << topVolume.transform().matrix());

  std::unique_ptr<TrackingVolume> world;
  static const std::string worldName = "World";

  if (const auto *cyl = dynamic_cast<const CylinderVolumeBounds *>(&bounds);
      cyl != nullptr) {
    ACTS_VERBOSE(prefix() << "Expanding cylinder bounds");
    using enum CylinderVolumeBounds::BoundValues;

    // Make a copy that we'll modify
    auto newBounds = std::make_shared<CylinderVolumeBounds>(*cyl);

    const auto &zEnv = m_cfg.envelope[AxisZ];
    if (zEnv[0] != zEnv[1]) {
      ACTS_ERROR(
          prefix() << "Root node cylinder envelope for z must be symmetric");
      throw std::logic_error(
          "Root node cylinder envelope for z must be "
          "symmetric");
    }

    const auto &rEnv = m_cfg.envelope[AxisR];

    newBounds->set({
        {eHalfLengthZ, newBounds->get(eHalfLengthZ) + zEnv[0]},
        {eMinR, std::max(0.0, newBounds->get(eMinR) - rEnv[0])},
        {eMaxR, newBounds->get(eMaxR) + rEnv[1]},
    });

    ACTS_DEBUG(prefix() << "Applied envelope to cylinder: Z=" << zEnv[0]
                        << ", Rmin=" << rEnv[0] << ", Rmax=" << rEnv[1]);

    world = std::make_unique<TrackingVolume>(topVolume.transform(),
                                             std::move(newBounds), worldName);

    // Need one-sided portal shell that connects outwards to nullptr
    SingleCylinderPortalShell worldShell{*world};
    worldShell.applyToVolume();

  } else if (const auto *box =
                 dynamic_cast<const CuboidVolumeBounds *>(&bounds);
             box != nullptr) {
    ACTS_VERBOSE(prefix() << "Expanding cuboid bounds");
    // Make a copy that we'll modify
    auto newBounds = std::make_shared<CuboidVolumeBounds>(*box);

    // Get the current half lengths
    double halfX = newBounds->get(CuboidVolumeBounds::eHalfLengthX);
    double halfY = newBounds->get(CuboidVolumeBounds::eHalfLengthY);
    double halfZ = newBounds->get(CuboidVolumeBounds::eHalfLengthZ);

    // Apply envelope to each dimension
    const auto &xEnv = m_cfg.envelope[AxisX];
    const auto &yEnv = m_cfg.envelope[AxisY];
    const auto &zEnv = m_cfg.envelope[AxisZ];

    // Check if envelopes are symmetric for all dimensions
    if (xEnv[0] != xEnv[1]) {
      ACTS_ERROR(
          prefix() << "Root node cuboid envelope for X must be symmetric");
      throw std::logic_error(
          "Root node cuboid envelope for X must be symmetric");
    }

    if (yEnv[0] != yEnv[1]) {
      ACTS_ERROR(
          prefix() << "Root node cuboid envelope for Y must be symmetric");
      throw std::logic_error(
          "Root node cuboid envelope for Y must be symmetric");
    }

    if (zEnv[0] != zEnv[1]) {
      ACTS_ERROR(
          prefix() << "Root node cuboid envelope for Z must be symmetric");
      throw std::logic_error(
          "Root node cuboid envelope for Z must be symmetric");
    }

    newBounds->set({
        {CuboidVolumeBounds::eHalfLengthX, halfX + xEnv[0]},
        {CuboidVolumeBounds::eHalfLengthY, halfY + yEnv[0]},
        {CuboidVolumeBounds::eHalfLengthZ, halfZ + zEnv[0]},
    });

    ACTS_DEBUG(prefix() << "Applied envelope to cuboid: X=" << xEnv[0]
                        << ", Y=" << yEnv[0] << ", Z=" << zEnv[0]);

    world = std::make_unique<TrackingVolume>(topVolume.transform(),
                                             std::move(newBounds), worldName);

    // Need one-sided portal shell that connects outwards to nullptr
    SingleCuboidPortalShell worldShell{*world};
    worldShell.applyToVolume();

  } else {
    throw std::logic_error{"Unsupported volume bounds type"};
  }

  ACTS_DEBUG(prefix() << "New root volume bounds are: "
                      << world->volumeBounds());

  world->setNavigationPolicy(std::make_unique<Acts::TryAllNavigationPolicy>(
      gctx, *world, logger, Acts::TryAllNavigationPolicy::Config{}));

  auto &shell = child.connect(options, gctx, logger);

  shell.fill(*world);

  child.finalize(options, gctx, *world, logger);

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

  return std::make_unique<TrackingGeometry>(
      std::move(world), nullptr, GeometryIdentifierHook{}, logger, false);
}

}  // namespace Acts::Experimental
