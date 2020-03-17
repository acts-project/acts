// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/detail/VerticesHelper.hpp"
#include "Acts/Utilities/BinningType.hpp"

void Acts::Polyhedron::merge(const Acts::Polyhedron& other) {
  size_t cvert = vertices.size();
  vertices.insert(vertices.end(), other.vertices.begin(), other.vertices.end());
  /// Add the new faces with offsets
  auto join = [&](std::vector<Face>& existing,
                  const std::vector<Face>& additional) -> void {
    for (const auto& aface : additional) {
      Face nface = aface;
      std::transform(nface.begin(), nface.end(), nface.begin(),
                     [&](size_t x) { return (x + cvert); });
      existing.push_back(nface);
    }
  };
  // For faces and triangular mesh
  join(faces, other.faces);
  join(triangularMesh, other.triangularMesh);
}

Acts::Extent Acts::Polyhedron::extent(const Transform3D& transform) const {
  Extent extent;
  for (const auto& vtx : vertices) {
    extent.check(transform * vtx);
  }
  // Special checks for binR:
  if (std::abs(extent.range(binZ)) < s_onSurfaceTolerance) {
    // Check inclusion of origin (i.e. convex around origin)
    Vector3D origin = transform * Vector3D(0., 0., extent.medium(binZ));
    for (const auto& face : faces) {
      std::vector<Vector3D> tface;
      tface.reserve(face.size());
      for (auto f : face) {
        tface.push_back(transform * vertices[f]);
      }
      if (detail::VerticesHelper::isInsidePolygon(origin, tface)) {
        extent.ranges[binR].first = 0.;
        break;
      }
    }
    // Check for radial extent
    auto radialDistance = [&](const Vector3D& pos1,
                              const Vector3D& pos2) -> double {
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
    };

    for (size_t iv = 1; iv < vertices.size() + 1; ++iv) {
      size_t fpoint = iv < vertices.size() ? iv : 0;
      double testR = radialDistance(transform * vertices[fpoint],
                                    transform * vertices[iv - 1]);
      extent.ranges[binR].first = std::min(extent.ranges[binR].first, testR);
    }
  }
  return extent;
}
