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
#include "Acts/Visualization/ViewConfig.hpp"

#include <array>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace Acts {

/// This helper produces output in the OBJ format. Note that colors are not
/// supported in this implementation.
///
template <typename T = double>
class ObjVisualization : public IVisualization {
 public:
  static_assert(std::is_same_v<T, double> || std::is_same_v<T, float>,
                "Use either double or float");

  /// Stored value type, should be double or float
  using ValueType = T;

  /// Type of a vertex based on the value type
  using VertexType = ActsVector<ValueType, 3>;

  /// Type of a line
  using LineType = std::pair<size_t, size_t>;

  /// Constructor that allows to set scalor and precision
  /// @param prec The output precission with std::setprecision
  /// @param scale An (optional) scaling for the writing out
  ObjVisualization(unsigned int prec = 4, double scale = 1.)
      : m_outputPrecision(prec), m_outputScalor(scale) {}

  /// @copydoc Acts::IVisualization::vertex()
  void vertex(const Vector3D& vtx, ColorRGB color = {0, 0, 0}) final;

  /// @copydoc Acts::IVisualization::line()
  void line(const Vector3D& a, const Vector3D& b,
            ColorRGB color = {0, 0, 0}) final;

  /// @copydoc Acts::IVisualization::face()
  void face(const std::vector<Vector3D>& vtxs,
            ColorRGB color = {0, 0, 0}) final;

  /// @copydoc Acts::IVisualization::faces()
  void faces(const std::vector<Vector3D>& vtxs,
             const std::vector<FaceType>& faces,
             ColorRGB color = {0, 0, 0}) final;

  /// @copydoc Acts::IVisualization::write()
  void write(const std::string& path) const final;

  /// @copydoc Acts::IVisualization::write()
  void write(std::ostream& os) const final;

  /// Write the object and the material file
  /// @param os the output stream for the object
  /// @param mos the output stream for the auxiliary material file
  void write(std::ostream& os, std::ostream& mos) const;

  ///  @copydoc Acts::IVisualization::clear()
  void clear() final;

 private:
  /// The output parameters
  unsigned int m_outputPrecision = 4;
  double m_outputScalor = 1.;
  /// The object data to be written
  std::vector<VertexType> m_vertices;
  std::vector<FaceType> m_faces;
  std::vector<LineType> m_lines;
  /// Map of colors to be written at given index position
  std::map<size_t, ColorRGB> m_lineColors;
  std::map<size_t, ColorRGB> m_vertexColors;
  std::map<size_t, ColorRGB> m_faceColors;
};

#include "detail/ObjVisualization.ipp"

}  // namespace Acts
