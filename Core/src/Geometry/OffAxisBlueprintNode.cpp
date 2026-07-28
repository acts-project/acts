// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/OffAxisBlueprintNode.hpp"

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"

namespace Acts {

OffAxisBlueprintNode::OffAxisBlueprintNode(std::string_view name,
                                           Transform3 axisTransform)
    : StaticBlueprintNode(nullptr),
      m_name(name),
      m_axisTransform(std::move(axisTransform)) {}

const std::string& OffAxisBlueprintNode::name() const {
  return m_name;
}

Volume& OffAxisBlueprintNode::build(const BlueprintOptions& options,
                                    const GeometryContext& gctx,
                                    const Logger& logger) {
  ACTS_DEBUG(prefix() << "Off-axis build");

  if (children().size() != 1u) {
    ACTS_ERROR(prefix() << "OffAxisBlueprintNode requires exactly one child");
    throw std::invalid_argument(
        "OffAxisBlueprintNode requires exactly one child");
  }

  Volume& childVolume = children().at(0).build(options, gctx, logger);
  const auto* childBounds =
      dynamic_cast<const CylinderVolumeBounds*>(&childVolume.volumeBounds());
  if (childBounds == nullptr) {
    ACTS_ERROR(prefix()
               << "OffAxisBlueprintNode only supports cylinder child volumes");
    throw std::invalid_argument(
        "OffAxisBlueprintNode only supports cylinder child volumes");
  }

  const Transform3 childGlobal = childVolume.localToGlobalTransform(gctx);
  const Transform3 childInAxisFrame = m_axisTransform.inverse() * childGlobal;

  constexpr auto tolerance = s_onSurfaceTolerance;
  if (std::abs(childInAxisFrame.rotation().col(eX)[eZ]) >= tolerance ||
      std::abs(childInAxisFrame.rotation().col(eY)[eZ]) >= tolerance) {
    ACTS_ERROR(prefix() << "Child volume rotation tilts relative to axis frame");
    throw std::invalid_argument(
        "OffAxisBlueprintNode child rotation tilts relative to axis frame");
  }

  const Vector3 localTranslation = childInAxisFrame.translation();
  const double radialOffset = localTranslation.head<2>().norm();
  const double childMinR = childBounds->get(CylinderVolumeBounds::eMinR);
  const double childMaxR = childBounds->get(CylinderVolumeBounds::eMaxR);
  const double minR = std::max(
      0., std::max(childMinR - radialOffset, radialOffset - childMaxR));
  const double maxR =
      childMaxR + radialOffset;
  const double halfLengthZ =
      childBounds->get(CylinderVolumeBounds::eHalfLengthZ);

  Transform3 envelopeInAxisFrame = childInAxisFrame;
  envelopeInAxisFrame.translation()[eX] = 0.;
  envelopeInAxisFrame.translation()[eY] = 0.;
  const Transform3 envelopeGlobal = m_axisTransform * envelopeInAxisFrame;

  m_volume = std::make_unique<TrackingVolume>(
      envelopeGlobal,
      std::make_shared<CylinderVolumeBounds>(minR, maxR, halfLengthZ), m_name);

  ACTS_DEBUG(prefix() << "Built envelope volume " << m_name << " with r=["
                      << minR << ", " << maxR << "] around child offset "
                      << radialOffset);

  return *m_volume;
}

OffAxisBlueprintNode& OffAxisBlueprintNode::setAxisTransform(
    const Transform3& axisTransform) {
  m_axisTransform = axisTransform;
  return *this;
}

const Transform3& OffAxisBlueprintNode::axisTransform() const {
  return m_axisTransform;
}

}  // namespace Acts
