// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Detector/detail/ReferenceGenerators.hpp"

const std::vector<Acts::Vector3>
Acts::Experimental::detail::PolyhedronReferenceGenerator::references(
    const GeometryContext& gctx, const Surface& surface) const {
  // Create the return  vector
  std::vector<Vector3> rPositions;
  auto pHedron = surface.polyhedronRepresentation(gctx, nSegements);
  rPositions.insert(rPositions.end(), pHedron.vertices.begin(),
                    pHedron.vertices.end());
  // Add the barycenter if configured
  if (addBarycenter) {
    Vector3 bc(0., 0., 0.);
    std::ranges::for_each(rPositions, [&](const auto& p) { bc += p; });
    bc *= 1. / rPositions.size();
    rPositions.push_back(bc);
  }
  return rPositions;
}

