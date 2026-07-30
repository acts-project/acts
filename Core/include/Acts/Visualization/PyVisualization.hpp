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
#include "Acts/Visualization/EventDataView3D.hpp"

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
      : m_outputPrecision(prec), m_outputScalor(scale){}

  /// @copydoc Acts::IVisualization3D::vertex()
  void vertex(const Vector3& vtx, Color color = s_defaultColor) {
    /// vertex in 2D, either rz- or xy-projection
    m_vertexcolor = color;
    m_vertex = Vector3{vtx[0], vtx[1], vtx[2]};
  }



  void line(const Vector3& a, const Vector3& b,
            Color color = s_defaultColor)
    {
      m_linecolors.push_back(color);
      std::vector<std::array<double,3>> segment;
      segment.push_back({a[0], a[1], a[2]});
      segment.push_back({b[0], b[1], b[2]});
      m_segments.push_back(segment);   
    }  

    // overload to allow for lineThickness argument to be accepted
    void line(const Vector3& a, const Vector3& b,  const double lineThickness,
            Color color = s_defaultColor)
    {
      m_linecolors.push_back(color);
      m_linethickness.push_back(lineThickness);
      std::vector<std::array<double,3>> segment;
      segment.push_back({a[0], a[1], a[2]});
      segment.push_back({b[0], b[1], b[2]});
      m_segments.push_back(segment);      
    } 

    

  /// @copydoc Acts::IVisualization3D::face()
  void face(const std::vector<Vector3>& vtxs,
            Color color = s_defaultColor) {}

  /// @copydoc Acts::IVisualization3D::faces()
  void faces(const std::vector<Vector3>& vtxs,
             const std::vector<FaceType>& faces,
             Color color = s_defaultColor )
    {
      m_facecolors.push_back(color);
      std::vector<std::array<double,3>> surface;
      for(const auto& v : vtxs){
          surface.push_back({v[0], v[1], v[2]});
        }
        m_surfaces.push_back(surface);
    }


  const Vector3 &vertex3D() const{
        return m_vertex;
  }

  const std::vector<std::vector<std::array<double, 3>>> &surfaces() const {
        return m_surfaces;
      }

  const std::vector<std::vector<std::array<double,3>>> &segments() const {
        return m_segments;
      }

  const std::vector<Color> &faceColors() const {
        return m_facecolors;
  }

  const std::vector<Color> &lineColors() const {
        return m_linecolors;
  }

  const std::vector<double> &lineThickness() const {
        return m_linethickness;
  }

  const Color &vertexColor() const {
        return m_vertexcolor;
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

  Vector3 m_vertex;
  std::string m_projection;
  std::vector<std::vector<std::array<double,3>>> m_surfaces;
  std::vector<std::vector<std::array<double,3>>> m_segments;
  std::vector<Color> m_facecolors;
  Color m_vertexcolor;
  std::vector<Color> m_linecolors;
  std::vector<double> m_linethickness;
};

}  // namespace Acts
