// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/StrawSurface.hpp"

#include "Acts/Geometry/GeometryObject.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/detail/FacesHelper.hpp"
#include "Acts/Surfaces/detail/VerticesHelper.hpp"

#include <numbers>
#include <utility>
#include <vector>

namespace Acts {

Polyhedron StrawSurface::polyhedronRepresentation(
    const GeometryContext& gctx, unsigned int quarterSegments) const {
  // Prepare vertices and faces
  std::vector<Vector3> vertices;
  std::vector<Polyhedron::FaceType> faces;
  std::vector<Polyhedron::FaceType> triangularMesh;

  const Transform3& ctransform = localToGlobalTransform(gctx);
  // Draw the bounds if more than one segment are chosen
  if (quarterSegments > 0u) {
    double r = m_bounds->get(LineBounds::eR);
    // Write the two bows/circles on either side
    std::vector<int> sides = {-1, 1};
    for (auto& side : sides) {
      // Helper method to create the segment
      auto svertices = detail::VerticesHelper::segmentVertices(
          {r, r}, -std::numbers::pi, std::numbers::pi, {}, quarterSegments,
          Vector3(0., 0., side * m_bounds->get(LineBounds::eHalfLengthZ)),
          ctransform);
      vertices.insert(vertices.end(), svertices.begin(), svertices.end());
    }
    auto facesMesh = detail::FacesHelper::cylindricalFaceMesh(vertices);
    faces = facesMesh.first;
    triangularMesh = facesMesh.second;
  }

  std::size_t bvertices = vertices.size();
  Vector3 left(0, 0, -m_bounds->get(LineBounds::eHalfLengthZ));
  Vector3 right(0, 0, m_bounds->get(LineBounds::eHalfLengthZ));
  // The central wire/straw
  vertices.push_back(ctransform * left);
  vertices.push_back(ctransform * right);
  faces.push_back({bvertices, bvertices + 1});
  vertices.push_back(ctransform * Vector3(0., 0., 0.));
  triangularMesh.push_back({bvertices, bvertices + 2, bvertices + 1});

  return Polyhedron(vertices, faces, triangularMesh, false);
}

}  // namespace Acts
