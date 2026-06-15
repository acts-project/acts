// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/PadBlueprintNode.hpp"

#include "Acts/Geometry/CuboidPortalShell.hpp"
#include "Acts/Geometry/CylinderPortalShell.hpp"
#include "Acts/Geometry/Extent.hpp"

namespace Acts::Experimental {

PadBlueprintNode::PadBlueprintNode(const std::string &name,
                                   const ExtentEnvelope &envelope)
    : StaticBlueprintNode(nullptr), m_envelope(envelope), m_name(name) {}

Volume &PadBlueprintNode::build(const BlueprintOptions &options,
                                const GeometryContext &gctx,
                                const Logger &logger) {
  using enum AxisDirection;

  if (children().size() != 1) {
    ACTS_ERROR(prefix() << "PadBlueprintNode must have exactly one child, "
                           "but has "
                        << children().size());
    throw std::invalid_argument("PadBlueprintNode must have exactly one child");
  }

  Volume &topVolume = children().at(0).build(options, gctx, logger);

  const auto &bounds = topVolume.volumeBounds();

  std::stringstream ss;
  bounds.toStream(ss);
  ACTS_DEBUG(prefix() << "have top volume: " << ss.str() << "\n"
                      << topVolume.localToGlobalTransform(gctx).matrix());

  if (const auto *cyl = dynamic_cast<const CylinderVolumeBounds *>(&bounds);
      cyl != nullptr) {
    ACTS_VERBOSE(prefix() << "Expanding cylinder bounds");
    using enum CylinderVolumeBounds::BoundValues;

    // Make a copy that we'll modify
    auto newBounds = std::make_shared<CylinderVolumeBounds>(*cyl);

    const auto &zEnv = m_envelope[AxisZ];
    if (zEnv[0] != zEnv[1]) {
      ACTS_ERROR(
          prefix() << "Root node cylinder envelope for z must be symmetric");
      throw std::logic_error(
          "Root node cylinder envelope for z must be "
          "symmetric");
    }

    const auto &rEnv = m_envelope[AxisR];

    newBounds->set({
        {eHalfLengthZ, newBounds->get(eHalfLengthZ) + zEnv[0]},
        {eMinR, std::max(0.0, newBounds->get(eMinR) - rEnv[0])},
        {eMaxR, newBounds->get(eMaxR) + rEnv[1]},
    });

    ACTS_DEBUG(prefix() << "Applied envelope to cylinder: Z=" << zEnv[0]
                        << ", Rmin=" << rEnv[0] << ", Rmax=" << rEnv[1]);

    m_volume = std::make_unique<TrackingVolume>(
        topVolume.localToGlobalTransform(gctx), std::move(newBounds), m_name);

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
    const auto &xEnv = m_envelope[AxisX];
    const auto &yEnv = m_envelope[AxisY];
    const auto &zEnv = m_envelope[AxisZ];

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

    m_volume = std::make_unique<TrackingVolume>(
        topVolume.localToGlobalTransform(gctx), std::move(newBounds), m_name);

  } else {
    throw std::logic_error{"Unsupported volume bounds type"};
  }

  return topVolume;
}

TrackingVolume *PadBlueprintNode::trackingVolume() const {
  return m_volume.get();
}

}  // namespace Acts::Experimental
