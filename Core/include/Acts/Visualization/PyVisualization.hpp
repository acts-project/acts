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
#include <map>
#include <string>
#include <vector>
#include <array>
#include <iostream>

namespace Acts {

/// This helper produces output for Python visualization. Note that colors are not
/// supported in this implementation.
///
class PyVisualization : public IVisualization3D {
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
  explicit PyVisualization(unsigned int prec = 4, double scale = 1., std::string projection = "rz")
      : m_outputPrecision(prec), m_outputScalor(scale), m_projection(projection) {}

  /// @copydoc Acts::IVisualization3D::vertex()
  void vertex(const Vector3& vtx, Color color = s_defaultColor) {}

  /// @copydoc Acts::IVisualization3D::line()
  void line(const Vector3& a, const Vector3& b,
            Color color = s_defaultColor) {}

  /// @copydoc Acts::IVisualization3D::face()
  void face(const std::vector<Vector3>& vtxs,
            Color color = s_defaultColor) {}

  /// @copydoc Acts::IVisualization3D::faces()
  void faces(const std::vector<Vector3>& vtxs,
             const std::vector<FaceType>& faces,
             Color color = s_defaultColor )
    {
      if (faces.size() > 1){
        std::cout << "Surface consists of more than one face. Do nothing.\n";
      }

      // rz projection
      if (m_projection == "rz"){
        std::vector<std::array<double,2>> surface;
        for(const auto& v : vtxs){
          float r = std::sqrt(v[0]*v[0] + v[1]*v[1]);
          surface.push_back({v[2], r});
          //std::cout << "(x,y,z)" << v[0] << v[1] << v[2] << "r=" << r << std::endl; 
          //std::cout << "(x,y,z)" << v[2] << ',' << r << std::endl; 
        }
        m_surfaces.push_back(surface);
      }

      if (m_projection == "xy"){
        std::vector<std::array<double,2>> surface;
        for(const auto& v : vtxs){
          surface.push_back({v[0], v[1]});
          //std::cout << "(x,y,z)" << v[0] << v[1] << v[2]<< std::endl; 
        }
        m_surfaces.push_back(surface);
      }
    }

  const std::vector<std::vector<std::array<double, 2>>> &surfaces() const {
        return m_surfaces;
      }

  const Color &color() const {
        return s_defaultColor;
  }

    /// @copydoc Acts::IVisualization3D::write(const std::filesystem::path&) const
  void write(const std::filesystem::path& path) const {};

  /// @copydoc Acts::IVisualization3D::write(std::ostream&) const
  void write(std::ostream& os) const {};

  /// Write the object and the material file
  /// @param os the output stream for the object
  /// @param mos the output stream for the auxiliary material file
  void write(std::ostream& os, std::ostream& mos) const {};


  ///  @copydoc Acts::IVisualization3D::clear()
  void clear() {}

  /// Start a new object context with a name
  /// @param name The name of the object
  void object(const std::string& name) {}

 private:
  struct Object {
    std::string name;
    std::vector<VertexType> vertices{};
    std::vector<FaceType> faces{};
    std::vector<LineType> lines{};

    /// The object data to be written
    /// Map of colors to be written at given index position
    std::map<std::size_t, Color> lineColors{};
    std::map<std::size_t, Color> vertexColors{};
    std::map<std::size_t, Color> faceColors{};
  };

  Object& object();
  const Object& object() const;

  /// The output parameters
  unsigned int m_outputPrecision = 4;
  double m_outputScalor = 1.;

  std::vector<Object> m_objects;

  std::string m_projection;
  std::vector<std::vector<std::array<double,2>>> m_surfaces;
};

}  // namespace Acts
