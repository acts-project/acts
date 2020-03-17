// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <limits>
#include <vector>

namespace Acts {

/// @class Polyhedron
///
/// Struct which contains a cartesian approximation for any surface type.
/// It contains a list of cartesian vertices in the global frame, and
/// additionally
/// a list of lists of indices which indicate which vertices form a face.
/// Each entry in @c faces is a face, which is in turn a list of vertices
/// that need to be connected to form a face.
/// This allows the @c objString method to produce a ready-to-go obj output.
struct Polyhedron {
  using Face = std::vector<size_t>;

  /// Default constructor
  Polyhedron() = default;

  /// Default constructor from a vector of vertices and a vector of faces
  /// @param verticesIn The 3D global vertices that make up the object
  /// @param facesIn List of lists of indices for faces.
  /// @param triangularMeshIn List of lists of indices for a triangular mesh
  /// @param radialCheck A dedicated check for radial extent done
  /// @note This creates copies of the input vectors
  Polyhedron(const std::vector<Vector3D>& verticesIn,
             const std::vector<Face>& facesIn,
             const std::vector<Face>& triangularMeshIn,
             bool radialCheck = false)
      : vertices(verticesIn),
        faces(facesIn),
        triangularMesh(triangularMeshIn) {}

  /// List of 3D vertices as vectors
  std::vector<Vector3D> vertices;

  /// List of faces connecting the vertices.
  /// each face is a list of vertices v
  /// corresponding to the vertex vector above
  std::vector<Face> faces;

  /// List of faces connecting the vertices.
  /// each face is a list of vertices v
  /// - in this case restricted to a triangular representation
  std::vector<Face> triangularMesh;

  /// Merge another Polyhedron into this one
  ///
  /// @param other is the source representation
  void merge(const Polyhedron& other);

  /// Draw method for polyhedrons
  ///
  /// @tparam helper_t The draw helper
  ///
  /// @param helper The draw helper object (visitor pattern)
  /// @param triangulate Force the faces to be a triangular mesh
  /// @param decompose Boolean that forces a decomposition into
  /// individual faces with unique vertices.
  template <typename helper_t>
  void draw(helper_t& helper, bool triangulate = false,
            bool decompose = false) const {
    // vertices and faces are
    if (not decompose and not triangulate) {
      helper.faces(vertices, faces);
    } else if (triangulate) {
      helper.faces(vertices, triangularMesh);
    } else {
      for (const auto& face : faces) {
        std::vector<Vector3D> face_vtx;
        for (size_t i : face) {
          face_vtx.push_back(vertices[i]);
        }
        helper.face(face_vtx);
      }
    }
  }

  /// Maximum extent of the polyhedron in space
  ///
  /// @param transform An (optional) transform
  /// to apply to the vertices for estimation the extent
  /// with respect to a given coordinate frame
  ///
  /// @return ranges that describe the space taken by this surface
  Extent extent(const Transform3D& transform = Transform3D::Identity()) const;
};
}  // namespace Acts
