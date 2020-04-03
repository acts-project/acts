// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Visualization/IVisualization.hpp"

#include <array>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

namespace Acts {

template <typename T = double>

/// @brief Helper to write out PlyVisualization visualization format
class PlyVisualization : public IVisualization {
 public:
  static_assert(std::is_same_v<T, double> || std::is_same_v<T, float>,
                "Use either double or float");

  /// Stored value type, should be double or float
  using ValueType = T;

  /// Type of a vertex based on the value type
  using VertexType = ActsVector<ValueType, 3>;

  /// @copydoc Acts::IVisualization::vertex()
  void vertex(const Vector3D& vtx,
              IVisualization::ColorType color = {120, 120, 120}) final;

  /// @copydoc Acts::IVisualization::line()
  void face(const std::vector<Vector3D>& vtxs,
            IVisualization::ColorType color = {120, 120, 120}) final;

  /// @copydoc Acts::IVisualization::faces()
  void faces(const std::vector<Vector3D>& vtxs, const std::vector<FaceType>&,
             ColorType color = {120, 120, 120}) final;

  /// @copydoc Acts::IVisualization::face()
  void line(const Vector3D& a, const Vector3D& b,
            IVisualization::ColorType color = {120, 120, 120}) final;

  /// @copydoc Acts::IVisualization::write()
  void write(const std::string& path) const final;

  /// @copydoc Acts::IVisualization::write()
  void write(std::ostream& os) const final;

  /// @copydoc Acts::IVisualization::clear()
  void clear() final;

 private:
  std::vector<std::pair<VertexType, IVisualization::ColorType>> m_vertices;
  std::vector<FaceType> m_faces;
  std::vector<std::pair<std::pair<size_t, size_t>, IVisualization::ColorType>>
      m_edges;
};

#include "detail/PlyVisualization.ipp"

}  // namespace Acts
