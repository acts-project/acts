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
#include "Acts/Utilities/GraphViz.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>

namespace {
const std::string s_rootName = "Root";
}

namespace Acts::Experimental {

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

  world->setNavigationPolicy(
      options.defaultNavigationPolicyFactory->build(gctx, *world, logger));

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

  TrackingVolume *currentVolume = nullptr;
  GeometryIdentifier::Value iportal = 0;
  GeometryIdentifier::Value isensitive = 0;
  world->apply(overloaded{
      [&](TrackingVolume &volume) {
        iportal = 0;
        isensitive = 0;
        currentVolume = &volume;

        if (volume.geometryId() != GeometryIdentifier{}) {
          return;
        }

        auto it = std::ranges::find(volumesById, nullptr);
        if (it == volumesById.end()) {
          ACTS_ERROR(prefix() << "No free volume IDs left, all "
                              << volumesById.size() << " are used");
          // @TODO: Maybe link to documentation about this
          throw std::logic_error("No free volume IDs left");
        }

        auto id = GeometryIdentifier().withVolume(
            std::distance(volumesById.begin(), it) + 1);

        ACTS_DEBUG(prefix() << "Assigning volume ID " << id << " for "
                            << volume.volumeName());
        volume.assignGeometryId(id);

        *it = &volume;
      },
      [&](::Acts::Portal &portal) {
        if (currentVolume == nullptr) {
          // This should not really happen
          ACTS_ERROR(prefix() << "No current volume found");
          throw std::logic_error("No current volume found");
        }
        if (portal.surface().geometryId() != GeometryIdentifier{}) {
          return;
        }

        iportal += 1;
        auto id = currentVolume->geometryId().withBoundary(iportal);
        ACTS_VERBOSE(prefix() << "Assigning portal ID: " << id);
        portal.surface().assignGeometryId(id);
      },
      [&](Surface &surface) {
        if (currentVolume == nullptr) {
          ACTS_ERROR(prefix() << "No current volume found");
          throw std::logic_error("No current volume found");
        }

        if (surface.geometryId() != GeometryIdentifier{}) {
          return;
        }

        isensitive += 1;
        auto id = currentVolume->geometryId().withSensitive(isensitive);
        ACTS_VERBOSE(prefix() << "Assigning surface ID: " << id);
        surface.assignGeometryId(id);
      }});

  return std::make_unique<TrackingGeometry>(
      std::move(world), nullptr, GeometryIdentifierHook{}, logger, false);
}

}  // namespace Acts::Experimental
