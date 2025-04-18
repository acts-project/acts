// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Visualization/IVisualization3D.hpp"
#include "Acts/Visualization/ViewConfig.hpp"

#include <array>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

namespace Acts {

template <typename T = double>

/// @brief Helper to write out PlyVisualization3D visualization format
class PlyVisualization3D : public IVisualization3D {
 public:
  static_assert(std::is_same_v<T, double> || std::is_same_v<T, float>,
                "Use either double or float");

  /// Stored value type, should be double or float
  using ValueType = T;

  /// Type of a vertex based on the value type
  using VertexType = Eigen::Matrix<ValueType, 3, 1>;

  /// @copydoc Acts::IVisualization3D::vertex()
  void vertex(const Vector3& vtx, ColorRGB color = {120, 120, 120}) final;

  /// @copydoc Acts::IVisualization3D::face()
  void face(const std::vector<Vector3>& vtxs,
            ColorRGB color = {120, 120, 120}) final;

  /// @copydoc Acts::IVisualization3D::faces()
  void faces(const std::vector<Vector3>& vtxs,
             const std::vector<FaceType>& faces,
             ColorRGB color = {120, 120, 120}) final;

  /// @copydoc Acts::IVisualization3D::line()
  void line(const Vector3& a, const Vector3& b,
            ColorRGB color = {120, 120, 120}) final;

  /// @copydoc Acts::IVisualization3D::write(const std::string&) const
  void write(const std::string& path) const final;

  /// @copydoc Acts::IVisualization3D::write(std::ostream&) const
  void write(std::ostream& os) const final;

  /// @copydoc Acts::IVisualization3D::clear()
  void clear() final;

 private:
  std::vector<std::pair<VertexType, ColorRGB>> m_vertices;
  std::vector<FaceType> m_faces;
  std::vector<std::pair<std::pair<std::size_t, std::size_t>, ColorRGB>> m_edges;
};

#ifndef DOXYGEN
#include "detail/PlyVisualization3D.ipp"
#endif

}  // namespace Acts
