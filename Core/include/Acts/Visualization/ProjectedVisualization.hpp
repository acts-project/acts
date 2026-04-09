// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Visualization/IVisualization3D.hpp"

#include <array>
#include <filesystem>
#include <functional>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace Acts {

/// Stub visualization backend for projected output.
/// Method bodies are intentionally empty and will be implemented later.
class ProjectedVisualization final : public IVisualization3D {
 public:
  using ProjectionFunction = std::function<Vector2(const Vector3& position)>;

  using Projection = std::pair<std::string, ProjectionFunction>;

  /// Constructor that allows to set the projections
  /// @param projections The projections to be used
  ProjectedVisualization(std::vector<Projection> projections);

  /// Destructor
  ~ProjectedVisualization() override = default;

  /// Copy constructor
  ProjectedVisualization(const ProjectedVisualization& other) = default;

  /// Copy assignment operator
  ProjectedVisualization& operator=(const ProjectedVisualization& other) =
      default;

  /// @copydoc Acts::IVisualization3D::vertex()
  void vertex(const Vector3& vtx, Color color = s_defaultColor) final;

  /// @copydoc Acts::IVisualization3D::face()
  void face(const std::vector<Vector3>& vtxs,
            Color color = s_defaultColor) final;

  /// @copydoc Acts::IVisualization3D::faces()
  void faces(const std::vector<Vector3>& vtxs,
             const std::vector<FaceType>& faces,
             Color color = s_defaultColor) final;

  /// @copydoc Acts::IVisualization3D::line()
  void line(const Vector3& a, const Vector3& b,
            Color color = s_defaultColor) final;

  /// @copydoc Acts::IVisualization3D::write(std::ostream&) const
  void write(std::ostream& os) const final;

  /// @copydoc Acts::IVisualization3D::write(const std::filesystem::path&) const
  void write(const std::filesystem::path& path) const final;

  /// @copydoc Acts::IVisualization3D::clear()
  void clear() final;

  /// @copydoc Acts::IVisualization3D::object()
  void object(const std::string& name) final;

  /// Get the projections
  const std::map<std::string, ProjectionFunction>& projections() const;
  /// Get the projected faces
  const std::map<std::string, std::vector<std::tuple<std::vector<Vector2>, Color>>>& projectedFaces() const;
  /// Get the projected lines
  const std::map<std::string, std::vector<std::tuple<std::array<Vector2, 2>, Color>>>& projectedLines() const;
  /// Get the projected vertices
  const std::map<std::string, std::vector<std::tuple<Vector2, Color>>>& projectedVertices() const;

 private:
  std::map<std::string, ProjectionFunction> m_projections;

  std::map<std::string, std::vector<std::tuple<std::vector<Vector2>, Color>>>
      m_projectedFaces;  // <projected vertices, color>
  std::map<std::string, std::vector<std::tuple<std::array<Vector2, 2>, Color>>>
      m_projectedLines;
  std::map<std::string, std::vector<std::tuple<Vector2, Color>>> m_projectedVertices;
};

// Projection functions
Vector2 projectToXY(const Vector3& position);
Vector2 projectToZPhi(const Vector3& position);
Vector2 projectToZR(const Vector3& position);

}  // namespace Acts
