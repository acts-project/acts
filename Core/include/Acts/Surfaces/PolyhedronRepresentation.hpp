// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>

#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @class PolyhedronRepresentation
///
/// Struct which contains a cartesian approximation for any surface type.
/// It contains a list of cartesian vertices in the global frame, and
/// additionally
/// a list of lists of indices which indicate which vertices form a face.
/// Each entry in @c faces is a face, which is in turn a list of vertices
/// that need to be connected to form a face.
/// This allows the @c objString method to produce a ready-to-go obj output.
struct PolyhedronRepresentation
{

  /// Default constructor from a vector of vertices and a vector of faces
  /// @param verticesIn The 3D global vertices that make up the object
  /// @param facesIn List of lists of indices for faces.
  /// @note This creates copies of the input vectors
  PolyhedronRepresentation(const std::vector<Vector3D>&            verticesIn,
                           const std::vector<std::vector<size_t>>& facesIn)
    : vertices(verticesIn), faces(facesIn)
  {
  }

  /// list of 3D vertices as vectors
  std::vector<Vector3D> vertices;

  /// list of faces connecting the vertices.
  /// each face is a list of vertices v
  /// corresponding to the vertex vector above
  std::vector<std::vector<size_t>> faces;

  /// Method to return an obj string.
  /// @param vtxOffset Optional obj vertex enumeration offset
  /// @note Vertices in obj are enumerated globally. The offset is required
  ///       to allow for multiple output objects in one obj file.
  std::string
  objString(size_t vtxOffset = 0) const;
};
}
