// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/IVisualization.hpp"

namespace Acts {

template <typename T = double>

/// @brief Helper to write out Ply visualization format
class PlyHelper : public IVisualization {
 public:
  static_assert(std::is_same_v<T, double> || std::is_same_v<T, float>,
                "Use either double or float");

  /// Stored value type, should be double or float
  using value_type = T;

  /// Type of a vertex based on the value type
  using vertex_type = ActsVector<value_type, 3>;

  /// @copydoc Acts::IVisualization::vertex()
  void vertex(const Vector3D& vtx,
              IVisualization::color_type color = {120, 120, 120}) override {
    m_vertices.emplace_back(vtx.template cast<value_type>(), color);
  }

  /// @copydoc Acts::IVisualization::line()
  void face(const std::vector<Vector3D>& vtxs,
            IVisualization::color_type color = {120, 120, 120}) override {
    face_type idxs;
    idxs.reserve(vtxs.size());
    for (const auto& vtx : vtxs) {
      vertex(vtx, color);
      idxs.push_back(m_vertices.size() - 1);
    }
    m_faces.push_back(std::move(idxs));
  }

  /// @copydoc Acts::IVisualization::faces()
  void faces(const std::vector<Vector3D>& vtxs, const std::vector<face_type>&,
             color_type color = {120, 120, 120}) {
    face(vtxs, color);
  }

  /// @copydoc Acts::IVisualization::face()
  void line(const Vector3D& a, const Vector3D& b,
            IVisualization::color_type color = {120, 120, 120}) override {
    vertex(a, color);
    size_t idx_a = m_vertices.size() - 1;
    vertex(b, color);
    size_t idx_b = m_vertices.size() - 1;
    m_edges.emplace_back(std::make_pair(std::make_pair(idx_a, idx_b), color));
  }

  /// @copydoc Acts::IVisualization::write()
  void write(std::ostream& os) const override {
    os << "ply\n";
    os << "format ascii 1.0\n";
    os << "element vertex " << m_vertices.size() << "\n";
    os << "property float x\n";
    os << "property float y\n";
    os << "property float z\n";
    os << "property uchar red\n";
    os << "property uchar green\n";
    os << "property uchar blue\n";
    os << "element face " << m_faces.size() << "\n";
    os << "property list uchar int vertex_index\n";
    os << "element edge " << m_edges.size() << "\n";
    os << "property int vertex1\n";
    os << "property int vertex2\n";
    os << "property uchar red\n";
    os << "property uchar green\n";
    os << "property uchar blue\n";
    os << "end_header\n";

    for (const std::pair<vertex_type, IVisualization::color_type>& vtx :
         m_vertices) {
      os << vtx.first.x() << " " << vtx.first.y() << " " << vtx.first.z()
         << " ";
      os << vtx.second[0] << " " << vtx.second[1] << " " << vtx.second[2]
         << "\n";
    }

    for (const face_type& fc : m_faces) {
      os << fc.size();
      for (size_t i = 0; i < fc.size(); i++) {
        os << " " << fc[i];
      }
      os << "\n";
    }

    for (const std::pair<std::pair<size_t, size_t>, IVisualization::color_type>&
             edge : m_edges) {
      std::pair<size_t, size_t> idxs = edge.first;
      os << idxs.first << " " << idxs.second << " ";
      os << edge.second[0] << " " << edge.second[1] << " " << edge.second[2]
         << "\n";
    }
  }

  /// @copydoc Acts::IVisualization::clear()
  void clear() override {
    m_vertices.clear();
    m_faces.clear();
    m_edges.clear();
  }

 private:
  std::vector<std::pair<vertex_type, IVisualization::color_type>> m_vertices;
  std::vector<face_type> m_faces;
  std::vector<std::pair<std::pair<size_t, size_t>, IVisualization::color_type>>
      m_edges;
};
}  // namespace Acts
