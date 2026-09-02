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

/// Result of growing a symmetric half-length by a per-side envelope.
struct AxialGrowth {
  double halfLength;   ///< new half-length
  double centerShift;  ///< shift of the volume center relative to the child
};

/// Grow a symmetric half-length by a per-side envelope @p env = {low, high}.
/// In @ref PadBlueprintNode::Centering::Centered mode the envelope must be
/// symmetric, so the center does not move; in @ref
/// PadBlueprintNode::Centering::FitBounds mode an asymmetric envelope shifts
/// the center to the midpoint of the expanded extent.
AxialGrowth growAxial(double halfLength, const std::array<double, 2> &env,
                      PadBlueprintNode::Centering centering,
                      AxisDirection direction, const Logger &logger) {
  if (centering == PadBlueprintNode::Centering::Centered) {
    requireSymmetric(env, direction, logger);
  }
  return {halfLength + 0.5 * (env[0] + env[1]), 0.5 * (env[1] - env[0])};
}

}  // namespace

PadBlueprintNode::PadBlueprintNode(const std::string &name,
                                   const ExtentEnvelope &envelope)
    : StaticBlueprintNode(nullptr), m_envelope(envelope), m_name(name) {}

std::unique_ptr<TrackingVolume> PadBlueprintNode::padded(
    const GeometryContext &gctx, const Volume &inner,
    const ExtentEnvelope &envelope, const std::string &name,
    const Logger &logger) {
  return padded(gctx, inner, envelope, name, std::nullopt, Centering::Centered,
                logger);
}

std::unique_ptr<TrackingVolume> PadBlueprintNode::padded(
    const GeometryContext &gctx, const Volume &inner,
    const ExtentEnvelope &envelope, const std::string &name,
    const std::optional<Transform3> &referenceAxis, Centering centering,
    const Logger &logger) {
  using enum AxisDirection;

  const auto &bounds = inner.volumeBounds();
  const Transform3 childGlobal = inner.localToGlobalTransform(gctx);

  std::stringstream ss;
  bounds.toStream(ss);
  ACTS_DEBUG("Padding volume: " << ss.str() << "\n" << childGlobal.matrix());

  std::shared_ptr<VolumeBounds> newBounds;
  // By default the padded volume is anchored on the child. Asymmetric axial
  // growth (FitBounds) and reference-axis alignment shift it as needed below.
  Transform3 volumeTransform = childGlobal;

  if (const auto *cyl = dynamic_cast<const CylinderVolumeBounds *>(&bounds);
      cyl != nullptr) {
    ACTS_VERBOSE("Expanding cylinder bounds");
    using enum CylinderVolumeBounds::BoundValues;

    const auto &zEnv = envelope[AxisZ];
    const auto &rEnv = envelope[AxisR];

    // Make a copy that we'll modify
    auto cylBounds = std::make_shared<CylinderVolumeBounds>(*cyl);

    // Axial growth is shared between both frames and may shift the z center.
    const auto [halfLengthZ, zShift] =
        growAxial(cylBounds->get(eHalfLengthZ), zEnv, centering, AxisZ, logger);

    double minR = 0.;
    double maxR = 0.;

    if (referenceAxis.has_value()) {
      // Reference-axis enclosure: express the child in the axis frame, recenter
      // transversely onto the axis, and grow the radial bounds so the displaced
      // child stays fully enclosed.
      const Transform3 childInAxisFrame =
          referenceAxis->inverse() * childGlobal;

      // The child cylinder axis must stay parallel to the reference axis,
      // otherwise it cannot be represented as an axis-aligned cylinder.
      constexpr auto tolerance = s_onSurfaceTolerance;
      if (std::abs(childInAxisFrame.rotation().col(eX)[eZ]) >= tolerance ||
          std::abs(childInAxisFrame.rotation().col(eY)[eZ]) >= tolerance) {
        ACTS_ERROR("Pad child cylinder tilts relative to the reference axis");
        throw std::logic_error(
            "PadBlueprintNode reference-axis child tilts relative to the "
            "reference axis");
      }

      const double radialOffset =
          childInAxisFrame.translation().head<2>().norm();
      const double childMinR = cylBounds->get(eMinR);
      const double childMaxR = cylBounds->get(eMaxR);
      // Nearest child material to the axis: the offset either pushes the inner
      // hole outward (offset < minR), the outer edge toward the axis
      // (offset > maxR), or the annulus straddles the axis (result 0).
      minR = std::max(
          0.0, std::max(childMinR - radialOffset, radialOffset - childMaxR) -
                   rEnv[0]);
      maxR = childMaxR + radialOffset + rEnv[1];

      // Drop the transverse offset and apply the axial shift in the axis frame.
      Transform3 envelopeInAxisFrame = childInAxisFrame;
      envelopeInAxisFrame.translation()[eX] = 0.0;
      envelopeInAxisFrame.translation()[eY] = 0.0;
      envelopeInAxisFrame.translation()[eZ] += zShift;
      volumeTransform = referenceAxis.value() * envelopeInAxisFrame;

      ACTS_DEBUG("Applied reference-axis envelope to cylinder: Rmin="
                 << minR << ", Rmax=" << maxR << " around child offset "
                 << radialOffset);
    } else {
      minR = std::max(0.0, cylBounds->get(eMinR) - rEnv[0]);
      maxR = cylBounds->get(eMaxR) + rEnv[1];
      volumeTransform = childGlobal * Translation3{Vector3{0., 0., zShift}};

      ACTS_DEBUG("Applied envelope to cylinder: Rmin=" << rEnv[0]
                                                       << ", Rmax=" << rEnv[1]);
    }

    cylBounds->set({
        {eHalfLengthZ, halfLengthZ},
        {eMinR, minR},
        {eMaxR, maxR},
    });
    newBounds = std::move(cylBounds);

  } else if (const auto *box =
                 dynamic_cast<const CuboidVolumeBounds *>(&bounds);
             box != nullptr) {
    if (referenceAxis.has_value()) {
      ACTS_ERROR(
          "Reference-axis padding is only supported for cylinder volumes");
      throw std::logic_error(
          "PadBlueprintNode reference-axis padding requires a cylinder child");
    }

    ACTS_VERBOSE("Expanding cuboid bounds");
    using enum CuboidVolumeBounds::BoundValues;

    // Make a copy that we'll modify
    auto boxBounds = std::make_shared<CuboidVolumeBounds>(*box);

    const auto [halfX, xShift] =
        growAxial(boxBounds->get(eHalfLengthX), envelope[AxisX], centering,
                  AxisX, logger);
    const auto [halfY, yShift] =
        growAxial(boxBounds->get(eHalfLengthY), envelope[AxisY], centering,
                  AxisY, logger);
    const auto [halfZ, zShift] =
        growAxial(boxBounds->get(eHalfLengthZ), envelope[AxisZ], centering,
                  AxisZ, logger);

    boxBounds->set({
        {eHalfLengthX, halfX},
        {eHalfLengthY, halfY},
        {eHalfLengthZ, halfZ},
    });
    volumeTransform =
        childGlobal * Translation3{Vector3{xShift, yShift, zShift}};

    ACTS_DEBUG("Applied envelope to cuboid");
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

  m_volume = padded(gctx, inner, m_envelope, m_name, m_referenceAxis,
                    m_centering, logger);

  return *m_volume;
}

PadBlueprintNode &PadBlueprintNode::setEnvelope(
    const ExtentEnvelope &envelope) {
  m_envelope = envelope;
  return *this;
}

const ExtentEnvelope &PadBlueprintNode::envelope() const {
  return m_envelope;
}

PadBlueprintNode &PadBlueprintNode::setReferenceAxis(const Transform3 &axis) {
  m_referenceAxis = axis;
  return *this;
}

PadBlueprintNode &PadBlueprintNode::clearReferenceAxis() {
  m_referenceAxis.reset();
  return *this;
}

const std::optional<Transform3> &PadBlueprintNode::referenceAxis() const {
  return m_referenceAxis;
}

PadBlueprintNode &PadBlueprintNode::setCentering(Centering centering) {
  m_centering = centering;
  return *this;
}

PadBlueprintNode::Centering PadBlueprintNode::centering() const {
  return m_centering;
}

}  // namespace Acts
