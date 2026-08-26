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
#include "Acts/Visualization/ViewConfig.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace Acts {

/// This helper collects the visualization primitives in memory instead of
/// writing them to a file in a specific format. The collected vertices, faces
/// and lines can be retrieved through accessors, e.g. to hand them over to an
/// external plotting backend.
class GenericVisualization3D : public IVisualization3D {
 public:
  /// A collected vertex: position and color
  struct Vertex {
    /// The vertex position
    Vector3 position = Vector3::Zero();
    /// The vertex color
    Color color = s_defaultColor;
  };

  /// A collected face: indices into the vertex collection and a color
  struct Face {
    /// The indices of the vertices that make up the face
    FaceType indices{};
    /// The face color
    Color color = s_defaultColor;
  };

  /// A collected line: the two end points and a color
  struct Line {
    /// The start point of the line
    Vector3 a = Vector3::Zero();
    /// The end point of the line
    Vector3 b = Vector3::Zero();
    /// The line color
    Color color = s_defaultColor;
  };

  /// @copydoc Acts::IVisualization3D::vertex()
  void vertex(const Vector3& vtx, Color color = s_defaultColor) final;

  /// @copydoc Acts::IVisualization3D::line()
  /// @note The end points are stored with the line and are not appended to
  /// the vertex collection.
  void line(const Vector3& a, const Vector3& b,
            Color color = s_defaultColor) final;

  /// @copydoc Acts::IVisualization3D::face()
  void face(const std::vector<Vector3>& vtxs,
            Color color = s_defaultColor) final;

  /// @copydoc Acts::IVisualization3D::faces()
  void faces(const std::vector<Vector3>& vtxs,
             const std::vector<FaceType>& faces,
             Color color = s_defaultColor) final;

  /// Write a short summary of the collected data to an output stream
  /// @note This is meant for debugging only, it is not a serialization
  /// format. Use ObjVisualization3D or PlyVisualization3D to write files.
  /// @param os The output stream
  void write(std::ostream& os) const final;

  /// Write a short summary of the collected data to a file
  /// @note This is meant for debugging only, it is not a serialization
  /// format. Use ObjVisualization3D or PlyVisualization3D to write files.
  /// @param path is the file system path for writing the file
  void write(const std::filesystem::path& path) const final;

  /// @copydoc Acts::IVisualization3D::clear()
  void clear() final;

  /// @copydoc Acts::IVisualization3D::object()
  /// @note Object contexts are currently not stored by this helper
  void object(const std::string& /*name*/) final {
    // Unimplemented
  }

  /// Access the collected vertices
  /// @return The vertices collected so far
  const std::vector<Vertex>& vertices() const { return m_vertices; }

  /// Access the collected faces
  /// @note This overloads the 3-argument virtual IVisualization3D::faces();
  /// call on a concrete type to avoid ambiguity.
  /// @return The faces collected so far
  const std::vector<Face>& faces() const { return m_faces; }

  /// Access the collected lines
  /// @return The lines collected so far
  const std::vector<Line>& lines() const { return m_lines; }

 private:
  std::vector<Vertex> m_vertices;
  std::vector<Face> m_faces;
  std::vector<Line> m_lines;
};

}  // namespace Acts
