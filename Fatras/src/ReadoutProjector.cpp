// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/Digitization/ReadoutProjector.hpp"
#include "ActsFatras/Digitization/DigitizationError.hpp"
#include "ActsFatras/EventData/Hit.hpp"

Acts::Result<ActsFatras::ProjectedSegment>
ActsFatras::ReadoutProjector::project(
    const ActsFatras::DigitizationInput& dInput, const Acts::Vector3D& dDir,
    double thickness) const {
  using Plane = Eigen::Hyperplane<double, 3>;
  using Line = Eigen::ParametrizedLine<double, 3>;

  // The particle line/projection
  auto pPos = dInput.hit.get().position();
  auto pDir = dInput.hit.get().unitDirectionBefore();
  Line sLine(pPos, pDir);

  // The surface properties
  const auto& surface = dInput.surface;
  const auto normal = surface->normal(dInput.geoContext);
  const auto center = surface->center(dInput.geoContext);

  // The module surface planes
  Plane plane0(normal, center - 0.5 * thickness * normal);
  Plane plane1(normal, center + 0.5 * thickness * normal);
  auto edge0 = pPos + sLine.intersection(plane0) * pDir;
  auto edge1 = pPos + sLine.intersection(plane1) * pDir;
  // Projection direction
  double proj = dDir.norm() > Acts::s_epsilon ? dDir.dot(normal) : 0.;

  // Drift the edges to the readout surface
  Acts::Vector3D drifted0 = edge0;
  double driftLength0 = 0.5 * thickness;
  Acts::Vector3D drifted1 = edge1;
  double driftLength1 = -0.5 * thickness;
  if (proj < 0.) {
    driftLength0 = 0.;
    Line dLine(edge1, dDir);
    driftLength1 = dLine.intersection(plane0);
    drifted1 = drifted1 + driftLength1 * pDir;
  } else if (proj > 0.) {
    Line dLine(edge0, dDir);
    driftLength0 = dLine.intersection(plane1);
    drifted0 = drifted0 + driftLength0 * pDir;
    driftLength1 = 0.;
  }

  // Bring to local and apply masking/clipping if necessary
  auto toLocal = surface->transform(dInput.geoContext).inverse();
  auto sResult = surfaceMask.apply(
      *surface,
      {(toLocal * drifted0).segment<2>(0), (toLocal * drifted1).segment<2>(0)});

  if (sResult.ok()) {
    auto segment = sResult.value();
    ProjectedSegment dSegment = {
        ProjectedPosition(segment.first, driftLength0),
        ProjectedPosition(segment.second, driftLength1)};

    return Acts::Result<ProjectedSegment>{dSegment};
  }

  return DigitizationError::MaskingError;
}