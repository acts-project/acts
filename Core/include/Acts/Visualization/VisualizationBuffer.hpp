// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Visualization/EventDataView3D.hpp"
#include "Acts/Visualization/IVisualization3D.hpp"
#include "Acts/Visualization/ViewConfig.hpp"

#include <array>
#include <filesystem>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace Acts {

/// This helper produces output for Python visualization. Note that colors are
/// not supported in this implementation.
///
class VisualizationBuffer : public IVisualization3D {
 public:
  /// Stored value type, should be double or float
  using ValueType = double;

  /// Type of a vertex based on the value type
  using VertexType = Eigen::Matrix<ValueType, 3, 1>;

  /// Type of a line
  using LineType = std::pair<std::size_t, std::size_t>;

  /// Constructor
  /// @param prec The output precision with std::setprecision
  /// @param scale An (optional) scaling for the writing out
  explicit VisualizationBuffer(unsigned int /*prec*/, double /*scale*/) {}

  /// @copydoc Acts::IVisualization3D::vertex()
  void vertex(const Vector3& vtx, Color color = s_defaultColor) override;

  /// @copydoc Acts::IVisualization3D::line()
  void line(const Vector3& a, const Vector3& b,
            Color color = s_defaultColor) override;

  // overload to allow for lineThickness argument to be accepted
  /// @copydoc Acts::IVisualization3D::line()
  /// @param lineThickness optional parameter to define the thickness of the line
  void line(const Vector3& a, const Vector3& b, const double lineThickness,
            Color color = s_defaultColor);

  /// @copydoc Acts::IVisualization3D::face()
  void face(const std::vector<Vector3>& /*vtxs*/, Color /*color*/) override {
    throw std::logic_error("face() is not supported for this type");
  }

  /// @copydoc Acts::IVisualization3D::faces()
  void faces(const std::vector<Vector3>& vtxs,
             const std::vector<FaceType>& /*faces*/,
             Color color = s_defaultColor) override;

  /// @copydoc Acts::IVisualization3D::write(const std::filesystem::path&) const
  void write(const std::filesystem::path& /*path*/) const override {
    throw std::logic_error("write() is not supported for this type");
  }

  /// @copydoc Acts::IVisualization3D::write(std::ostream&) const
  void write(std::ostream& /*os*/) const override {
    throw std::logic_error("write() is not supported for this type");
  }

  /// Write the object and the material file
  /// @param os the output stream for the object
  /// @param mos the output stream for the auxiliary material file
  void write(std::ostream& /*os*/, std::ostream& /*mos*/) const {
    throw std::logic_error("write() is not supported for this type");
  }

  ///  @copydoc Acts::IVisualization3D::clear()
  void clear() override {}

  /// Start a new object context with a name
  /// @param name The name of the object
  void object(const std::string& /*name*/) override {}

  /// @return vector of 3D vertices
  const std::vector<Vector3>& vertices3D() const { return m_vertices; }

  /// All vertices that make up a surface are collected in a vector
  /// All of these vectors are again stored as a vector which can be understood
  /// as a vector of surfaces
  /// @return vector of surfaces
  const std::vector<std::vector<Vector3>>& surfaces() const {
    return m_surfaces;
  }

  /// Two vertices that form a line segment are collected in a vector
  /// All of these vectors are again stored as a vector which can be understood
  /// as the whole line
  /// @return vector of line segments
  const std::vector<std::vector<Vector3>>& segments() const {
    return m_segments;
  }

  /// Vector that contains the colors of a set of faces forming a surface
  /// @return vector of Color objects belonging to a face
  const std::vector<Color>& faceColors() const { return m_facecolors; }

  /// Vector that contains the colors of a set of line segments forming a line
  /// @return vector of Color objects belonging to a line segment
  const std::vector<Color>& lineColors() const { return m_linecolors; }

  /// Vector that contains the line thickness of a set of line segments forming
  /// a line
  /// @return vector of doubles that define the line thckness of a line segment
  const std::vector<double>& lineThickness() const { return m_linethickness; }

  /// Vector that contains the colors of vertices
  /// @return vector of Color objects belonging to a vertex
  const std::vector<Color>& vertexColors() const { return m_vertexcolors; }

 private:
  std::vector<Vector3> m_vertices;
  std::vector<std::vector<Vector3>> m_surfaces;
  std::vector<std::vector<Vector3>> m_segments;
  std::vector<Color> m_facecolors;
  std::vector<Color> m_vertexcolors;
  std::vector<Color> m_linecolors;
  std::vector<double> m_linethickness;
};

}  // namespace Acts
