// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Visualization/ProjectedVisualization.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>

namespace Acts {

ProjectedVisualization::ProjectedVisualization(
    std::vector<Projection> projections) {
  for (const auto& projection : projections) {
    m_projections[projection.first] = projection.second;
  }
}

void ProjectedVisualization::vertex(const Vector3& vtx, Color color) {
  for (const auto& projection : m_projections) {
    const auto projected = projection.second(vtx);
    if (m_projectedVertices.contains(projection.first)) {
      m_projectedVertices[projection.first].push_back({projected, color});
    } else {
      m_projectedVertices[projection.first] = {{projected, color}};
    }
  }
}

void ProjectedVisualization::face(const std::vector<Vector3>& vtxs,
                                  Color color) {
  for (const auto& projection : m_projections) {
    std::vector<Vector2> projectedVertices;
    projectedVertices.reserve(vtxs.size());
    for (const auto& vtx : vtxs) {
      projectedVertices.push_back(projection.second(vtx));
    }
    if (m_projectedFaces.contains(projection.first)) {
      m_projectedFaces[projection.first].push_back({projectedVertices, color});
    } else {
      m_projectedFaces[projection.first] = {{projectedVertices, color}};
    }
  }
}

void ProjectedVisualization::faces(const std::vector<Vector3>& vtxs,
                                   const std::vector<FaceType>& faces,
                                   Color color) {
  if (faces.empty()) {
    face(vtxs, color);
  } else {
    for (const auto& faceIdxs : faces) {
      std::vector<Vector3> faceVertices;
      faceVertices.reserve(faceIdxs.size());
      std::transform(faceIdxs.begin(), faceIdxs.end(),
                     std::back_inserter(faceVertices),
                     [&](auto vtxIdx) { return vtxs[vtxIdx]; });
      face(faceVertices, color);
    }
  }
}

void ProjectedVisualization::line(const Vector3& a, const Vector3& b,
                                  Color color) {
  for (const auto& projection : m_projections) {
    const auto projectedA = projection.second(a);
    const auto projectedB = projection.second(b);
    if (m_projectedLines.contains(projection.first)) {
      m_projectedLines[projection.first].push_back(
          {{projectedA, projectedB}, color});
    } else {
      m_projectedLines[projection.first] = {{{projectedA, projectedB}, color}};
    }
  }
}

void ProjectedVisualization::write(std::ostream& os) const {
  for (const auto& projection : m_projections) {
    os << "Projection: " << projection.first << "\n";
    const auto verticesIt = m_projectedVertices.find(projection.first);
    if (verticesIt != m_projectedVertices.end()) {
      for (const auto& [projectedVertices, color] : verticesIt->second) {
        os << "Vertex: " << projectedVertices.x() << " "
           << projectedVertices.y() << " " << color[0] << " " << color[1] << " "
           << color[2] << "\n";
      }
    }

    const auto facesIt = m_projectedFaces.find(projection.first);
    if (facesIt != m_projectedFaces.end()) {
      for (const auto& [projectedFaces, color] : facesIt->second) {
        os << "Face: " << projectedFaces.size() << " " << color[0] << " "
           << color[1] << " " << color[2] << "\n";
      }
    }

    const auto linesIt = m_projectedLines.find(projection.first);
    if (linesIt != m_projectedLines.end()) {
      for (const auto& [projectedLines, color] : linesIt->second) {
        os << "Line: " << projectedLines[0].x() << " " << projectedLines[0].y()
           << " " << projectedLines[1].x() << " " << projectedLines[1].y()
           << " " << color[0] << " " << color[1] << " " << color[2] << "\n";
      }
    }
  }
}

const std::map<std::string, ProjectedVisualization::ProjectionFunction>&
ProjectedVisualization::projections() const {
  return m_projections;
}

const std::map<std::string,
               std::vector<std::tuple<std::vector<Vector2>, Color>>>&
ProjectedVisualization::projectedFaces() const {
  return m_projectedFaces;
}

const std::map<std::string,
               std::vector<std::tuple<std::array<Vector2, 2>, Color>>>&
ProjectedVisualization::projectedLines() const {
  return m_projectedLines;
}

const std::map<std::string, std::vector<std::tuple<Vector2, Color>>>&
ProjectedVisualization::projectedVertices() const {
  return m_projectedVertices;
}

void ProjectedVisualization::write(const std::filesystem::path& path) const {
  std::ofstream os(path.string());
  write(os);
  os.close();
}

void ProjectedVisualization::clear() {
  m_projectedVertices.clear();
  m_projectedFaces.clear();
  m_projectedLines.clear();
}

void ProjectedVisualization::object(const std::string&) {}

Vector2 projectToXY(const Vector3& position) {
  return Vector2(position.x(), position.y());
}

Vector2 projectToZPhi(const Vector3& position) {
  return Vector2(position.z(), std::atan2(position.y(), position.x()));
}

Vector2 projectToZR(const Vector3& position) {
  return Vector2(position.z(), std::sqrt(position.x() * position.x() +
                                         position.y() * position.y()));
}

}  // namespace Acts
