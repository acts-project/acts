// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Digitization/Channelizer.hpp"

Acts::Result<std::vector<ActsFatras::Segmentizer::ChannelSegment>>
ActsFatras::Channelizer::channelize(const Hit& hit,
                                    const Acts::Surface& surface,
                                    const Acts::GeometryContext& gctx,
                                    const Acts::Vector3& driftDir,
                                    const Acts::BinUtility& segmentation,
                                    double thickness,
                                    double minRelPerpDrift) const {
  // Drifted surface and scalor 2D to 3D segment
  Acts::Result<
      std::tuple<PlanarSurfaceDrift::Segment2D, PlanarSurfaceDrift::Segment3D>>
      atReadoutPlane = m_surfaceDrift.toReadout(
          gctx, surface, thickness, hit.position(), hit.direction(), driftDir);
  if (!atReadoutPlane.ok()) {
    return atReadoutPlane.error();
  }
  // The drifted and the full segment
  const auto& [driftedSegment, fullSegement] = atReadoutPlane.value();

  // Applies the surface mask
  auto maskedSegmentRes = m_surfaceMask.apply(surface, driftedSegment);
  if (!maskedSegmentRes.ok()) {
    return maskedSegmentRes.error();
  }

  // Now Channelize, i.e. segments are mapped to the readout grid
  auto segments =
      m_segmentizer.segments(gctx, surface, segmentation, *maskedSegmentRes);

  double driftedPathLength = (driftedSegment[1] - driftedSegment[0]).norm();
  // In case we have close-to-nominal incident, we could run into numerical
  // problems, we'll take a 0.1 permille of the thickness as a threshold
  if (std::abs(driftedPathLength) < minRelPerpDrift * thickness &&
      segments.size() == 1) {
    segments[0].activation = thickness;
    return segments;
  }

  double fullPathLength = (fullSegement[1] - fullSegement[0]).norm();
  double sclale2Dto3D = fullPathLength / driftedPathLength;
  // scale the activations
  for (auto& segment : segments) {
    segment.activation *= sclale2Dto3D;
  }

  return segments;
}
