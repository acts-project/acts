// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Digitization/Channelizer.hpp"

namespace ActsFatras {

Acts::Result<std::vector<Segmentizer::ChannelSegment>> Channelizer::channelize(
    const Hit& hit, const Acts::Surface& surface,
    const Acts::GeometryContext& gctx, const Acts::Vector3& driftDir,
    const Acts::IMultiAxis& segmentation, double thickness,
    double minRelPerpDrift) const {
  // Drifted surface and scalor 2D to 3D segment
  // SurfaceDrift handles the surface-type-specific local frame internally
  // (plane/disc Cartesian, cylinder unrolled (rPhi, z))
  const auto atReadoutPlane = m_surfaceDrift.toReadout(
      gctx, surface, thickness, hit.position(), hit.direction(), driftDir);
  if (!atReadoutPlane.ok()) {
    return atReadoutPlane.error();
  }
  // The drifted and the full segment
  const auto& [driftedSegment, fullSegment] = *atReadoutPlane;

  // Applies the surface mask (also surface-type agnostic)
  const auto maskedSegmentRes = m_surfaceMask.apply(surface, driftedSegment);
  if (!maskedSegmentRes.ok()) {
    return maskedSegmentRes.error();
  }

  // Now Channelize, i.e. segments are mapped to the readout grid
  auto segments =
      m_segmentizer.segments(gctx, surface, segmentation, *maskedSegmentRes);

  const double driftedPathLength =
      (driftedSegment[1] - driftedSegment[0]).norm();
  // In case we have close-to-nominal incident, we could run into numerical
  // problems
  if (std::abs(driftedPathLength) < minRelPerpDrift * thickness &&
      segments.size() == 1) {
    segments[0].activation = thickness;
    return segments;
  }

  const double fullPathLength = (fullSegment[1] - fullSegment[0]).norm();
  const double scale2Dto3D = fullPathLength / driftedPathLength;
  // scale the activations
  for (auto& segment : segments) {
    segment.activation *= scale2Dto3D;
  }

  return segments;
}

}  // namespace ActsFatras
