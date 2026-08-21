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

  /// Constructor that allows to set scalor and precision
  /// @param prec The output precision with std::setprecision
  /// @param scale An (optional) scaling for the writing out
  explicit VisualizationBuffer(unsigned int = 4, double = 1.) {}

  /// @copydoc Acts::IVisualization3D::vertex()
  void vertex(const Vector3& vtx, Color color = s_defaultColor);

  void line(const Vector3& a, const Vector3& b, Color color = s_defaultColor);

  // overload to allow for lineThickness argument to be accepted
  void line(const Vector3& a, const Vector3& b, const double lineThickness,
            Color color = s_defaultColor);

  /// @copydoc Acts::IVisualization3D::face()
  void face(const std::vector<Vector3>&, Color) {
    throw std::logic_error("face() is not supported for this type");
  }

  /// @copydoc Acts::IVisualization3D::faces()
  void faces(const std::vector<Vector3>& vtxs, const std::vector<FaceType>&,
             Color color = s_defaultColor);

  /// @copydoc Acts::IVisualization3D::write(const std::filesystem::path&) const
  void write(const std::filesystem::path&) const {
    throw std::logic_error("write() is not supported for this type");
  }

  /// @copydoc Acts::IVisualization3D::write(std::ostream&) const
  void write(std::ostream&) const {
    throw std::logic_error("write() is not supported for this type");
  }

  /// Write the object and the material file
  /// @param os the output stream for the object
  /// @param mos the output stream for the auxiliary material file
  void write(std::ostream&, std::ostream&) const {
    throw std::logic_error("write() is not supported for this type");
  }

  ///  @copydoc Acts::IVisualization3D::clear()
  void clear() {}

  /// Start a new object context with a name
  /// @param name The name of the object
  void object(const std::string&) {}

  const std::vector<Vector3>& vertices3D() const { return m_vertices; }

  const std::vector<std::vector<Vector3>>& surfaces() const {
    return m_surfaces;
  }

  const std::vector<std::vector<Vector3>>& segments() const {
    return m_segments;
  }

  const std::vector<Color>& faceColors() const { return m_facecolors; }

  const std::vector<Color>& lineColors() const { return m_linecolors; }

  const std::vector<double>& lineThickness() const { return m_linethickness; }

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
