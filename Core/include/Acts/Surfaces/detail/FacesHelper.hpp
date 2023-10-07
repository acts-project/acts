// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/Polyhedron.hpp"

#include <numeric>
#include <utility>
#include <vector>

namespace Acts {

namespace detail {

/// @brief Helper for writing out faces for polyhedron representation
struct FacesHelper {
  using FaceVector = std::vector<Polyhedron::FaceType>;

  /// @brief This method words for all convex type surface setups
  /// It includes:
  ///
  /// Rectangle / Triangle / Polygon
  /// Full disc, ellipse, half ring
  /// @param vertices The vector of vertices
  /// @param centerLast Boolean indicator if the center is given for
  /// a better triangulation method as last element of the vector
  static std::pair<FaceVector, FaceVector> convexFaceMesh(
      const std::vector<Vector3>& vertices, bool centerLast = false) {
    FaceVector faces;
    FaceVector triangularMesh;
    // Write the face
    unsigned int offset = centerLast ? 1 : 0;
    std::vector<size_t> face(vertices.size() - offset);
    std::iota(face.begin(), face.end(), 0);
    faces.push_back(face);
    /// Triangular mesh construction
    unsigned int anker = centerLast ? vertices.size() - 1 : 0;
    for (unsigned int it = 2 - offset; it < vertices.size() - offset; ++it) {
      triangularMesh.push_back({anker, it - 1, it});
    }
    // Close for centered reference point
    if (centerLast) {
      triangularMesh.push_back({anker, vertices.size() - 2, 0});
    }
    return {faces, triangularMesh};
  }

  /// @brief This method works for all concentric type surface setups
  /// It includes :
  ///
  /// - Cylinder (concentric bows on each side, separated by z)
  /// - Cut-off cone
  ///
  /// The single requirement is that the #vertices are equal and the input
  /// vector is splittable in half into the two bows.
  ///
  /// @param vertices The vector of vertices
  /// @param fullTwoPi The indicator if the concentric face is closed
  static std::pair<FaceVector, FaceVector> cylindricalFaceMesh(
      const std::vector<Vector3>& vertices, bool fullTwoPi = true) {
    FaceVector faces;
    FaceVector triangularMesh;
    size_t nqfaces = static_cast<size_t>(0.5 * vertices.size());
    size_t reduce = (not fullTwoPi) ? 1 : 0;
    for (size_t iface = 0; iface < nqfaces - reduce; ++iface) {
      size_t p2 = (iface + 1 == nqfaces) ? 0 : iface + 1;
      std::vector<size_t> face = {iface, p2, p2 + nqfaces, nqfaces + iface};
      faces.push_back(face);
      std::vector<size_t> triA = {iface, p2, p2 + nqfaces};
      triangularMesh.push_back(triA);
      std::vector<size_t> triB = {p2 + nqfaces, nqfaces + iface, iface};
      triangularMesh.push_back(triB);
    }
    return {faces, triangularMesh};
  }
};

}  // namespace detail

}  // namespace Acts
