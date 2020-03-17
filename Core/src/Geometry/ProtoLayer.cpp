// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>

#include <algorithm>
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Utilities/Helpers.hpp"

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

namespace Acts {

ProtoLayer::ProtoLayer(const GeometryContext& gctx,
                       const std::vector<const Surface*>& surfaces) {
  measure(gctx, surfaces);
}

ProtoLayer::ProtoLayer(
    const GeometryContext& gctx,
    const std::vector<std::shared_ptr<const Surface>>& surfaces) {
  measure(gctx, unpack_shared_vector(surfaces));
}

double ProtoLayer::radialDistance(const Vector3D& pos1,
                                  const Vector3D& pos2) const {
  Vector2D p1(pos1.x(), pos1.y());
  Vector2D p2(pos2.x(), pos2.y());

  Vector2D O(0, 0);
  Vector2D p1p2 = (p2 - p1);
  double L = p1p2.norm();
  Vector2D p1O = (O - p1);

  // don't do division if L is very small
  if (L < 1e-7) {
    return std::numeric_limits<double>::max();
  }
  double f = p1p2.dot(p1O) / L;

  // clamp to [0, |p1p2|]
  f = std::min(L, std::max(0., f));

  Vector2D closest = f * p1p2.normalized() + p1;
  double dist = (closest - O).norm();

  return dist;
}

std::ostream& ProtoLayer::toStream(std::ostream& sl) const {
  sl << "ProtoLayer with dimensions (min/max)" << std::endl;
  sl << " - r : " << minR << " - " << envR.first << " / " << maxR << " + "
     << envR.second << std::endl;
  sl << " - z : " << minZ << " - " << envZ.first << " / " << maxZ << " + "
     << envZ.second << std::endl;
  sl << " - phi : " << minPhi << " - " << envPhi.first << " / " << maxPhi
     << " + " << envPhi.second << std::endl;

  return sl;
}

void ProtoLayer::measure(const GeometryContext& gctx,
                         const std::vector<const Surface*>& surfaces) {
  minR = std::numeric_limits<double>::max();
  maxR = std::numeric_limits<double>::lowest();
  minX = std::numeric_limits<double>::max();
  maxX = std::numeric_limits<double>::lowest();
  minY = std::numeric_limits<double>::max();
  maxY = std::numeric_limits<double>::lowest();
  minZ = std::numeric_limits<double>::max();
  maxZ = std::numeric_limits<double>::lowest();
  minPhi = std::numeric_limits<double>::max();
  maxPhi = std::numeric_limits<double>::lowest();

  for (const auto& sf : surfaces) {
    // if the associated detector element exists, use
    // it for thickness
    double thickness = 0;
    const DetectorElementBase* element = sf->associatedDetectorElement();
    if (element != nullptr) {
      thickness = element->thickness();
    }

    // check the shape
    const PlanarBounds* pBounds =
        dynamic_cast<const PlanarBounds*>(&(sf->bounds()));

    const CylinderSurface* cylSurface =
        dynamic_cast<const CylinderSurface*>(sf);

    if (pBounds != nullptr) {
      const auto& sTransform = sf->transform(gctx);

      // get the vertices
      std::vector<Vector2D> vertices = pBounds->vertices();
      size_t nVertices = vertices.size();
      // loop over the two sides of the module
      // we only need to run once if no element i.e. no thickness
      for (int side = 0; side < (element != nullptr ? 2 : 1); ++side) {
        // loop over the vertex combinations
        for (size_t iv = 0; iv < nVertices; ++iv) {
          size_t ivp = iv != 0u ? iv - 1 : nVertices - 1;
          // thickness
          double locz = side != 0 ? 0.5 * thickness : -0.5 * thickness;
          // p1 & p2 vectors
          Vector3D p2(sTransform *
                      Vector3D(vertices.at(iv).x(), vertices.at(iv).y(), locz));
          Vector3D p1(sTransform * Vector3D(vertices.at(ivp).x(),
                                            vertices.at(ivp).y(), locz));

          maxX = std::max(maxX, p2.x());
          minX = std::min(minX, p2.x());

          maxY = std::max(maxY, p2.y());
          minY = std::min(minY, p2.y());

          maxZ = std::max(maxZ, p2.z());
          minZ = std::min(minZ, p2.z());

          maxR = std::max(maxR, perp(p2));
          minR = std::min(minR, radialDistance(p1, p2));

          maxPhi = std::max(maxPhi, phi(p2));
          minPhi = std::min(minPhi, phi(p2));
        }
      }
    } else if (cylSurface != nullptr) {
      // this is an explicit cast and if right now.
      // It should work with all Polyhedrons
      // @TODO: Remove the cast and if as soon as ::polyhedronRepresentation()
      //        makes it into the Surface base class
      //        The envelopes might need special treatments though

      Polyhedron ph = cylSurface->polyhedronRepresentation(gctx, 1);
      // evaluate at all vertices
      for (const auto& vtx : ph.vertices) {
        maxX = std::max(maxX, vtx.x());
        minX = std::min(minX, vtx.x());

        maxY = std::max(maxY, vtx.y());
        minY = std::min(minY, vtx.y());

        maxZ = std::max(maxZ, vtx.z());
        minZ = std::min(minZ, vtx.z());

        maxR = std::max(maxR, perp(vtx));

        maxPhi = std::max(maxPhi, phi(vtx));
        minPhi = std::min(minPhi, phi(vtx));
      }

      // trace all face connections to possibly catch min-r approach
      for (const auto& face : ph.faces) {
        for (size_t i = 0; i < face.size(); i++) {
          Vector3D p1 = ph.vertices.at(face.at(i));
          Vector3D p2 = ph.vertices.at(face.at((i + 1) % face.size()));
          minR = std::min(minR, radialDistance(p1, p2));
        }
      }

      // set envelopes to half radius
      double cylBoundsR = cylSurface->bounds().r();
      double env = cylBoundsR / 2.;
      envX = {env, env};
      envY = {env, env};
      envZ = {env, env};
      envR = {env, env};

      // evaluate impact of r shift on phi
      double cylPosR = perp(cylSurface->center(gctx));
      double dPhi = std::atan((cylBoundsR + env) / cylPosR) -
                    std::atan(cylBoundsR / cylPosR);

      // use this as phi envelope
      envPhi = {dPhi, dPhi};

    } else {
      const CylinderBounds* cBounds =
          dynamic_cast<const CylinderBounds*>(&(sf->bounds()));

      const AnnulusBounds* aBounds =
          dynamic_cast<const AnnulusBounds*>(&(sf->bounds()));

      if (cBounds != nullptr) {
        double r = cBounds->r();
        double z = sf->center(gctx).z();
        double hZ = cBounds->halflengthZ();
        double phi = cBounds->averagePhi();
        double hPhi = cBounds->halfPhiSector();

        maxX = r;
        minX = -r;

        maxY = r;
        minY = -r;

        maxZ = z + hZ;
        minZ = z - hZ;

        maxR = r;
        minR = r;

        maxPhi = phi + hPhi;
        minPhi = phi - hPhi;

      } else if (aBounds != nullptr) {
        minR = aBounds->rMin();
        maxR = aBounds->rMax();
        double z = sf->center(gctx).z();

        maxX = maxR;
        minX = -maxR;

        maxY = maxR;
        minY = -maxR;

        maxZ = z + 0.5 * thickness;
        minZ = z - 0.5 * thickness;

        maxPhi = aBounds->phiMin();
        minPhi = aBounds->phiMax();

      } else {
        throw std::domain_error(
            "Not implemented yet for the given bounds type.");
      }
    }
  }
}

}  // namespace Acts
