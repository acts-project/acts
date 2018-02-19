// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Layers/ProtoLayer.hpp"
#include <algorithm>

namespace Acts {

ProtoLayer::ProtoLayer(std::vector<const Surface*> surfaces)
{

  minR   = std::numeric_limits<double>::max();
  maxR   = std::numeric_limits<double>::lowest();
  minX   = std::numeric_limits<double>::max();
  maxX   = std::numeric_limits<double>::lowest();
  minY   = std::numeric_limits<double>::max();
  maxY   = std::numeric_limits<double>::lowest();
  minZ   = std::numeric_limits<double>::max();
  maxZ   = std::numeric_limits<double>::lowest();
  minPhi = std::numeric_limits<double>::max();
  maxPhi = std::numeric_limits<double>::lowest();

  for (const auto& sf : surfaces) {
    // if the associated detector element exists, use
    // it for thickness
    double                     thickness = 0;
    const DetectorElementBase* element   = sf->associatedDetectorElement();
    if (element) thickness               = element->thickness();

    // check the shape
    const PlanarBounds* pBounds
        = dynamic_cast<const PlanarBounds*>(&(sf->bounds()));
    if (pBounds) {

      // get the vertices
      std::vector<Vector2D> vertices  = pBounds->vertices();
      size_t                nVertices = vertices.size();
      // loop over the two sides of the module
      // we only need to run once if no element i.e. no thickness
      for (int side = 0; side < (element ? 2 : 1); ++side) {
        // loop over the vertex combinations
        for (size_t iv = 0; iv < nVertices; ++iv) {
          size_t ivp = iv ? iv - 1 : nVertices - 1;
          // thickness
          double locz = side ? 0.5 * thickness : -0.5 * thickness;
          // p1 & p2 vectors
          Vector3D p2(sf->transform() * Vector3D(vertices.at(iv).x(),
                                                 vertices.at(iv).y(),
                                                 locz));
          Vector3D p1(sf->transform() * Vector3D(vertices.at(ivp).x(),
                                                 vertices.at(ivp).y(),
                                                 locz));

          maxX = std::max(maxX, p2.x());
          minX = std::min(minX, p2.x());

          maxY = std::max(maxY, p2.y());
          minY = std::min(minY, p2.y());

          maxZ = std::max(maxZ, p2.z());
          minZ = std::min(minZ, p2.z());

          maxR = std::max(maxR, p2.perp());
          minR = std::min(minR, radialDistance(p1, p2));

          maxPhi = std::max(maxPhi, p2.phi());
          minPhi = std::min(minPhi, p2.phi());
        }
      }
    } else {
      throw std::domain_error("Not implemented yet for Non-planar bounds");
    }
  }
}

double
ProtoLayer::radialDistance(const Vector3D& pos1, const Vector3D& pos2) const
{
  Vector2D p1(pos1.x(), pos1.y());
  Vector2D p2(pos2.x(), pos2.y());

  Vector2D O(0, 0);
  Vector2D p1p2 = (p2 - p1);
  double   L    = p1p2.norm();
  Vector2D p1O  = (O - p1);

  // don't do division if L is very small
  if (L < 1e-7) return std::numeric_limits<double>::max();
  double f = p1p2.dot(p1O) / L;

  // clamp to [0, |p1p2|]
  f = std::min(L, std::max(0., f));

  Vector2D closest = f * p1p2.unit() + p1;
  double   dist    = (closest - O).norm();

  return dist;
}

std::ostream&
ProtoLayer::dump(std::ostream& sl) const
{
  sl << "ProtoLayer with dimensions (min/max)" << std::endl;
  sl << " - r : " << minR << "/" << maxR << std::endl;
  sl << " - z : " << minZ << "/" << maxZ << std::endl;
  sl << " - phi : " << minPhi << "/" << maxPhi << std::endl;

  return sl;
}

}  // namespace Acts
