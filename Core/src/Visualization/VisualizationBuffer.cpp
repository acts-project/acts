// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Visualization/VisualizationBuffer.hpp"

namespace Acts {

void VisualizationBuffer::vertex(const Vector3& vtx, Color color) {
  /// vertex in 2D, either rz- or xy-projection
  m_vertexcolors.push_back(color);
  m_vertices.push_back(Vector3{vtx[0], vtx[1], vtx[2]});
}

void VisualizationBuffer::line(const Vector3& a, const Vector3& b,
                               Color color) {
  m_linecolors.push_back(color);
  std::vector<Vector3> segment;
  segment.push_back(Vector3{a[0], a[1], a[2]});
  segment.push_back(Vector3{b[0], b[1], b[2]});
  m_segments.push_back(segment);
}

// overload to allow for lineThickness argument to be accepted
void VisualizationBuffer::line(const Vector3& a, const Vector3& b,
                               const double lineThickness, Color color) {
  m_linecolors.push_back(color);
  m_linethickness.push_back(lineThickness);
  std::vector<Vector3> segment;
  segment.push_back(Vector3{a[0], a[1], a[2]});
  segment.push_back(Vector3{b[0], b[1], b[2]});
  m_segments.push_back(segment);
}

/// @copydoc Acts::IVisualization3D::faces()
void VisualizationBuffer::faces(const std::vector<Vector3>& vtxs,
                                const std::vector<FaceType>&, Color color) {
  m_facecolors.push_back(color);
  std::vector<Vector3> surface;
  for (const auto& v : vtxs) {
    surface.push_back(Vector3{v[0], v[1], v[2]});
  }
  m_surfaces.push_back(surface);
}

}  // namespace Acts
