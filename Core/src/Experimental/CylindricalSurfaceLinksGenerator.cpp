// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Experimental/CylindricalSurfaceLinksGenerator.hpp"

#include "Acts/Experimental/SurfaceLinks.hpp"
#include "Acts/Experimental/Enumerate.hpp"
#include "Acts/Experimental/GeometricExtent.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include <exception>

Acts::InternalSurfaceLinks Acts::BarrelGridSurfaceLinksGenerator::operator()(
    const GeometryContext& gctx,
    const std::vector<std::shared_ptr<Surface>>& surfaces,
    const VolumeBounds* volumeBounds) const {
  // Bail out if the provided volume does not have cylindrical shape
  if (volumeBounds->type() != VolumeBounds::eCylinder) {
    throw std::invalid_argument(
        "\n *** BarrelGridSurfaceLinkGenerator: volume bounds are not "
        "cylindrical.");
  }

  using namespace detail;

  const auto& values = volumeBounds->values();
  ActsScalar halfZ = values[CylinderVolumeBounds::eHalfLengthZ];

  Axis<AxisType::Equidistant, AxisBoundaryType::Closed> phiAxis(-M_PI, M_PI,
                                                                binsPhi);

  Acts::GeometricExtent zPhiExtent;
  for (const auto& s : surfaces){
    const auto& polyHedron = s->polyhedronRepresentation(gctx,1);
    zPhiExtent.extend(polyHedron.vertices.begin(), polyHedron.vertices.end(), {binZ, binPhi});
  }

  ActsScalar minZ = zPhiExtent.min(binZ);
  ActsScalar maxZ = zPhiExtent.max(binZ);
  EquidistantAxis zAxis(minZ, maxZ, binsZ);

  Grid<int, decltype(phiAxis), decltype(zAxis)> phizGrid(
      {std::move(phiAxis), std::move(zAxis)});

  // Initialize the grid to -1
  for (unsigned int gbin = 0; gbin < binsPhi * binsZ; ++gbin) {
    phizGrid.at(gbin) = -1;
  }
  // Fill the surface indices 
  for (auto [i, s] : enumerate(surfaces)) {
    Vector3 center = s->center(gctx);
    decltype(phizGrid)::point_t phiz = {VectorHelpers::phi(center), center.z()};
    phizGrid.atPosition(phiz) = i;
  }

  // Create a single grid surface
  SingleGridSurfaces barrelSurfaceLinks(std::move(phizGrid), {binPhi, binZ});
  auto portalSurfaceLinks = std::vector<SurfaceLinks>(4u, AllSurfaces{});
  portalSurfaceLinks[2] = barrelSurfaceLinks;
  portalSurfaceLinks[3] = barrelSurfaceLinks;

  return {AllSurfaces{}, portalSurfaceLinks};
}
