// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Visualization/GenericVisualization3D.hpp"

#include <fstream>
#include <ostream>

namespace Acts {

void GenericVisualization3D::vertex(const Vector3& vtx, Color color) {
  m_vertices.push_back({vtx, color});
}

void GenericVisualization3D::line(const Vector3& a, const Vector3& b,
                                  Color color) {
  m_lines.push_back({a, b, color});
}

void GenericVisualization3D::face(const std::vector<Vector3>& vtxs,
                                  Color color) {
  Face fc;
  fc.color = color;
  fc.indices.reserve(vtxs.size());
  for (const auto& vtx : vtxs) {
    vertex(vtx, color);
    fc.indices.push_back(m_vertices.size() - 1);
  }
  m_faces.push_back(std::move(fc));
}

void GenericVisualization3D::faces(const std::vector<Vector3>& vtxs,
                                   const std::vector<FaceType>& faces,
                                   Color color) {
  // No faces given - call the face() method
  if (faces.empty()) {
    face(vtxs, color);
    return;
  }
  const std::size_t vtxOffset = m_vertices.size();
  for (const auto& vtx : vtxs) {
    vertex(vtx, color);
  }
  for (const auto& fc : faces) {
    // Two-index faces describe lines, matching ObjVisualization3D
    if (fc.size() == 2) {
      line(vtxs[fc[0]], vtxs[fc[1]], color);
    } else {
      Face storedFace;
      storedFace.color = color;
      storedFace.indices.reserve(fc.size());
      for (std::size_t idx : fc) {
        storedFace.indices.push_back(idx + vtxOffset);
      }
      m_faces.push_back(std::move(storedFace));
    }
  }
}

void GenericVisualization3D::write(std::ostream& os) const {
  os << "GenericVisualization3D: " << m_vertices.size() << " vertices, "
     << m_faces.size() << " faces, " << m_lines.size() << " lines\n";
}

void GenericVisualization3D::write(const std::filesystem::path& path) const {
  std::ofstream os(path);
  write(os);
}

void GenericVisualization3D::clear() {
  m_vertices.clear();
  m_faces.clear();
  m_lines.clear();
}

}  // namespace Acts
