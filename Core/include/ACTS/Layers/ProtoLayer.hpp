// This file is part of the ACTS project.
//
// Copyright (C) 2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_LAYERS_PROTOLAYER_H
#define ACTS_LAYERS_PROTOLAYER_H 1

#include <iostream>
#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {

/// @struct ProtoLayer
///
/// Encapsulates min/max boundaries that will be turned into a layer.
/// The struct allows this information to be obtained in a consistent
/// way, or be caller provided.

struct ProtoLayer
{
  double maxX;
  double minX;

  double maxY;
  double minY;

  double maxZ;
  double minZ;

  double maxR;
  double minR;

  double maxPhi;
  double minPhi;

  /// Constructor
  ///
  /// Loops over a provided vector of surface and calculates the various
  /// min/max values in one go. Also takes into account the thickness
  /// of an associated DetectorElement, if it exists.
  /// @param surfaces The vector of surfaces to consider
  ProtoLayer(std::vector<const Surface*> surfaces)
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

  /// Calculates the closest radial distance of a line
  ///
  /// @param pos1 is the first position on the line
  /// @param pos2 is the second position on the line
  ///
  /// @return is the closest distance
  static double
  radialDistance(const Vector3D& pos1, const Vector3D& pos2)
  {
    // following nominclature found in header file and doxygen documentation
    // line one is the straight track
    const Vector3D& ma = pos1;
    const Vector3D  ea = (pos2 - pos1).unit();
    // line two is the line surface
    Vector3D mb(0., 0., 0);
    Vector3D eb(0., 0., 1.);
    // now go ahead and solve for the closest approach
    Vector3D mab(mb - ma);
    double   eaTeb = ea.dot(eb);
    double   denom = 1 - eaTeb * eaTeb;
    if (std::abs(denom) > 10e-7) {
      double lambda0 = (mab.dot(ea) - mab.dot(eb) * eaTeb) / denom;
      // evaluate validaty in terms of bounds
      if (lambda0 < 1. && lambda0 > 0.) return (ma + lambda0 * ea).perp();
      return lambda0 < 0. ? pos1.perp() : pos2.perp();
    }
    return 10e101;
  }
  
  std::ostream&
  dump(std::ostream& sl) const
  {
    sl << "ProtoLayer with dimensions (min/max)" << std::endl;
    sl << " - r : " << minR << "/" << maxR << std::endl;
    sl << " - z : " << minZ << "/" << maxZ << std::endl;
    sl << " - phi : " << minPhi << "/" << maxPhi << std::endl;

    return sl;

  }
};
}

#endif  // ACTS_LAYERS_PROTOLAYER_H
