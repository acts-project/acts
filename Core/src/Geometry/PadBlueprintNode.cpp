// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/PadBlueprintNode.hpp"

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Geometry/CuboidPortalShell.hpp"
#include "Acts/Geometry/CylinderPortalShell.hpp"
#include "Acts/Geometry/Extent.hpp"

namespace Acts {

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
                                   const ExtentEnvelope &envelope,
                                   std::optional<Transform3> axisTransform)
    : StaticBlueprintNode(nullptr),
      m_envelope(envelope),
      m_name(name),
      m_axisTransform(std::move(axisTransform)) {}

std::unique_ptr<TrackingVolume> PadBlueprintNode::padded(
    const GeometryContext &gctx, const Volume &inner,
    const ExtentEnvelope &envelope, const std::string &name,
    const std::optional<Transform3> &axisTransform, const Logger &logger) {
  using enum AxisDirection;

  const auto &bounds = inner.volumeBounds();
  const Transform3 childGlobal = inner.localToGlobalTransform(gctx);

  std::stringstream ss;
  bounds.toStream(ss);
  ACTS_DEBUG("Padding volume: " << ss.str() << "\n" << childGlobal.matrix());

  std::shared_ptr<VolumeBounds> newBounds;
  // By default the padded volume is centered on the child. An axis-aligned
  // enclosure overrides this to sit on the reference axis instead.
  Transform3 volumeTransform = childGlobal;

  if (const auto *cyl = dynamic_cast<const CylinderVolumeBounds *>(&bounds);
      cyl != nullptr) {
    ACTS_VERBOSE("Expanding cylinder bounds");
    using enum CylinderVolumeBounds::BoundValues;

    const auto &zEnv = envelope[AxisZ];
    const auto &rEnv = envelope[AxisR];
    requireSymmetric(zEnv, AxisZ, logger);

    // Make a copy that we'll modify
    auto cylBounds = std::make_shared<CylinderVolumeBounds>(*cyl);

    if (axisTransform.has_value()) {
      // Axis-aligned enclosure: express the child in the reference axis frame,
      // recenter the envelope onto the axis, and grow the radial bounds so the
      // displaced child stays fully enclosed.
      const Transform3 childInAxisFrame =
          axisTransform->inverse() * childGlobal;

      // The child cylinder axis must stay parallel to the reference axis,
      // otherwise it cannot be represented as an axis-aligned cylinder.
      constexpr auto tolerance = s_onSurfaceTolerance;
      if (std::abs(childInAxisFrame.rotation().col(eX)[eZ]) >= tolerance ||
          std::abs(childInAxisFrame.rotation().col(eY)[eZ]) >= tolerance) {
        ACTS_ERROR("Pad child cylinder tilts relative to the reference axis");
        throw std::logic_error(
            "PadBlueprintNode axis-aligned child tilts relative to the "
            "reference axis");
      }

      const double radialOffset =
          childInAxisFrame.translation().head<2>().norm();
      const double childMinR = cylBounds->get(eMinR);
      const double childMaxR = cylBounds->get(eMaxR);
      // Nearest child material to the axis: the offset either pushes the inner
      // hole outward (offset < minR), the outer edge toward the axis
      // (offset > maxR), or the annulus straddles the axis (result 0).
      const double minR = std::max(
          0.0, std::max(childMinR - radialOffset, radialOffset - childMaxR) -
                   rEnv[0]);
      const double maxR = childMaxR + radialOffset + rEnv[1];

      cylBounds->set({
          {eHalfLengthZ, cylBounds->get(eHalfLengthZ) + zEnv[0]},
          {eMinR, minR},
          {eMaxR, maxR},
      });

      // Drop the transverse offset so the envelope sits on the reference axis,
      // keeping the child's longitudinal position and orientation.
      Transform3 envelopeInAxisFrame = childInAxisFrame;
      envelopeInAxisFrame.translation()[eX] = 0.0;
      envelopeInAxisFrame.translation()[eY] = 0.0;
      volumeTransform = axisTransform.value() * envelopeInAxisFrame;

      ACTS_DEBUG("Applied axis-aligned envelope to cylinder: Z="
                 << zEnv[0] << ", Rmin=" << minR << ", Rmax=" << maxR
                 << " around child offset " << radialOffset);
    } else {
      cylBounds->set({
          {eHalfLengthZ, cylBounds->get(eHalfLengthZ) + zEnv[0]},
          {eMinR, std::max(0.0, cylBounds->get(eMinR) - rEnv[0])},
          {eMaxR, cylBounds->get(eMaxR) + rEnv[1]},
      });

      ACTS_DEBUG("Applied envelope to cylinder: Z="
                 << zEnv[0] << ", Rmin=" << rEnv[0] << ", Rmax=" << rEnv[1]);
    }
    newBounds = std::move(cylBounds);

  } else if (const auto *box =
                 dynamic_cast<const CuboidVolumeBounds *>(&bounds);
             box != nullptr) {
    if (axisTransform.has_value()) {
      ACTS_ERROR("Axis-aligned padding is only supported for cylinder volumes");
      throw std::logic_error(
          "PadBlueprintNode axis-aligned padding requires a cylinder child");
    }

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

  return std::make_unique<TrackingVolume>(volumeTransform, std::move(newBounds),
                                          name);
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

  m_volume = padded(gctx, inner, m_envelope, m_name, m_axisTransform, logger);

  return *m_volume;
}

PadBlueprintNode &PadBlueprintNode::setAxisTransform(
    const Transform3 &axisTransform) {
  m_axisTransform = axisTransform;
  return *this;
}

const std::optional<Transform3> &PadBlueprintNode::axisTransform() const {
  return m_axisTransform;
}

}  // namespace Acts
