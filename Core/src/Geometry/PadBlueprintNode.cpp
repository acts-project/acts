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

namespace {

/// Throw if @p env is asymmetric, naming @p direction in the message.
void requireSymmetric(const std::array<double, 2> &env, AxisDirection direction,
                      const Logger &logger) {
  if (env[0] == env[1]) {
    return;
  }
  std::stringstream ss;
  ss << "Pad envelope for " << direction << " must be symmetric";
  ACTS_ERROR(ss.str());
  throw std::logic_error(ss.str());
}

}  // namespace

PadBlueprintNode::PadBlueprintNode(const std::string &name,
                                   const ExtentEnvelope &envelope)
    : StaticBlueprintNode(nullptr), m_envelope(envelope), m_name(name) {}

std::unique_ptr<TrackingVolume> PadBlueprintNode::padded(
    const GeometryContext &gctx, const Volume &inner,
    const ExtentEnvelope &envelope, const std::string &name,
    const Logger &logger) {
  using enum AxisDirection;

  const auto &bounds = inner.volumeBounds();

  std::stringstream ss;
  bounds.toStream(ss);
  ACTS_DEBUG("Padding volume: " << ss.str() << "\n"
                                << inner.localToGlobalTransform(gctx).matrix());

  std::shared_ptr<VolumeBounds> newBounds;

  if (const auto *cyl = dynamic_cast<const CylinderVolumeBounds *>(&bounds);
      cyl != nullptr) {
    ACTS_VERBOSE("Expanding cylinder bounds");
    using enum CylinderVolumeBounds::BoundValues;

    const auto &zEnv = envelope[AxisZ];
    const auto &rEnv = envelope[AxisR];
    requireSymmetric(zEnv, AxisZ, logger);

    // Make a copy that we'll modify
    auto cylBounds = std::make_shared<CylinderVolumeBounds>(*cyl);
    cylBounds->set({
        {eHalfLengthZ, cylBounds->get(eHalfLengthZ) + zEnv[0]},
        {eMinR, std::max(0.0, cylBounds->get(eMinR) - rEnv[0])},
        {eMaxR, cylBounds->get(eMaxR) + rEnv[1]},
    });

    ACTS_DEBUG("Applied envelope to cylinder: Z="
               << zEnv[0] << ", Rmin=" << rEnv[0] << ", Rmax=" << rEnv[1]);
    newBounds = std::move(cylBounds);

  } else if (const auto *box =
                 dynamic_cast<const CuboidVolumeBounds *>(&bounds);
             box != nullptr) {
    ACTS_VERBOSE("Expanding cuboid bounds");
    using enum CuboidVolumeBounds::BoundValues;

    // A cuboid is centered on its transform in every direction, so *all* of
    // the envelopes have to be symmetric.
    const auto &xEnv = envelope[AxisX];
    const auto &yEnv = envelope[AxisY];
    const auto &zEnv = envelope[AxisZ];
    requireSymmetric(xEnv, AxisX, logger);
    requireSymmetric(yEnv, AxisY, logger);
    requireSymmetric(zEnv, AxisZ, logger);

    // Make a copy that we'll modify
    auto boxBounds = std::make_shared<CuboidVolumeBounds>(*box);
    boxBounds->set({
        {eHalfLengthX, boxBounds->get(eHalfLengthX) + xEnv[0]},
        {eHalfLengthY, boxBounds->get(eHalfLengthY) + yEnv[0]},
        {eHalfLengthZ, boxBounds->get(eHalfLengthZ) + zEnv[0]},
    });

    ACTS_DEBUG("Applied envelope to cuboid: X=" << xEnv[0] << ", Y=" << yEnv[0]
                                                << ", Z=" << zEnv[0]);
    newBounds = std::move(boxBounds);

  } else {
    throw std::logic_error{"Unsupported volume bounds type"};
  }

  return std::make_unique<TrackingVolume>(inner.localToGlobalTransform(gctx),
                                          std::move(newBounds), name);
}

Volume &PadBlueprintNode::build(const BlueprintOptions &options,
                                const GeometryContext &gctx,
                                const Logger &logger) {
  if (children().size() != 1) {
    ACTS_ERROR(prefix() << "PadBlueprintNode must have exactly one child, "
                           "but has "
                        << children().size());
    throw std::invalid_argument("PadBlueprintNode must have exactly one child");
  }

  const Volume &inner = children().at(0).build(options, gctx, logger);

  m_volume = padded(gctx, inner, m_envelope, m_name, logger);

  return *m_volume;
}

}  // namespace Acts::Experimental
