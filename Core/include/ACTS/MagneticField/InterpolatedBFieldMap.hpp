// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <sstream>
#include <string>
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace Acts {

namespace Test {
  struct InterpolatedBFieldTester;
}

/// @brief interpolate magnetic field read from map
///
/// This class implements a magnetic field service which is initialized by a
/// field map. The field map needs to be provided upon construction and only
/// the overall normalization can be adjusted afterwards. Magnetic field vectors
/// at points inside the provided field map region are calculated from the
/// magnetic field at the eight surrounding grid points using tri-linear
/// interpolation.
///
/// You can find example field maps in the "Additional data" section on the <a
/// href="http://acts.web.cern.ch/ACTS/">ACTS webpage</a> (e.g. an ATLAS
/// magnetic field map).
class InterpolatedBFieldMap final
{
  /// unit test helper class
  friend struct Test::InterpolatedBFieldTester;

public:
  /// @brief configuration object for InterpolatedBFieldMap class
  struct Config
  {
    /// @brief global B-field scaling factor
    ///
    /// @note Negative values for @p scale are accepted and will invert the
    ///       direction of the magnetic field.
    double scale = 1.;

    /// @brief path to file containing field map
    ///
    /// @note For the description of the format of such a file, please refer to
    /// the documentation for InterpolatedBFieldMap::InterpolatedBFieldMap.
    std::string fieldMapFile = "";

    /// length unit used in field map
    double lengthUnit = units::_mm;

    /// magnetic field unit used in field map
    double BFieldUnit = units::_T;

    /// @brief number of grid points in field map
    ///
    /// @note This information is only used as a hint for the required size of
    ///       the internal vectors. A correct value is not needed, but will help
    ///       to speed up the field map initialization process.
    size_t nPoints = 1000;
  };

  /// @brief construct from magnetic field map
  ///
  /// @param [in] config configuration object
  ///
  /// Initialize the object with information read from the field map file.
  ///
  /// @note The format of the field map file specified by Config::fieldMapFile
  ///       must be the following:
  /// - Every row contains the coordinates of a point and the magnetic field
  ///   vector at this point.
  /// - The delimiter between numbers is a whitespace and the dot is used as
  ///   decimal point.
  /// - The order in every row is: <tt>x y z Bx By Bz</tt>
  /// - Comment lines starting with '#' as first character and empty lines are
  ///   skipped.
  /// - The points must form an equidistant grid in 3D.
  /// - Rows must be ordered according to ascending values first in x, then in y
  ///   and finally in z, e.g.
  /// @code
  /// 1 1 1 ...
  /// 1 1 2 ...
  /// 1 1 3 ...
  /// 1 2 1 ...
  /// 1 2 2 ...
  /// 1 2 3 ...
  /// 1 3 1 ...
  /// 1 3 2 ...
  /// 1 3 3 ...
  /// 2 1 1 ...
  /// 2 1 2 ...
  /// ...
  /// @endcode
  InterpolatedBFieldMap(Config config) : m_config(std::move(config))
  {
    initialize(config);
  }

  /// @brief get configuration object
  ///
  /// @return copy of the internal configuration object
  Config
  getConfiguration() const
  {
    return m_config;
  }

  /// @brief retrieve magnetic field value
  ///
  /// @param [in]  xyz   global position
  /// @param [out] bxyz  magnetic field vector at given position
  void
  getField(const double* xyz, double* bxyz) const
  {
    Vector3D    pos(xyz[0], xyz[1], xyz[2]);
    const auto& B = getField(pos);
    bxyz[0]       = B.x();
    bxyz[1]       = B.y();
    bxyz[2]       = B.z();
  }

  /// @brief retrieve magnetic field value
  ///
  /// @param [in] pos global position
  ///
  /// @return magnetic field vector at given position
  Vector3D
  getField(const Vector3D& pos) const
  {
    if (!insideGrid(pos)) return Vector3D(0, 0, 0);

    return interpolate(pos);
  }

  /// @brief get global scaling factor for magnetic field
  ///
  /// @return global factor for scaling the magnetic field
  double
  getScale() const
  {
    return m_config.scale;
  }

  /// @brief update configuration
  ///
  /// @param config new configuration object
  ///
  /// Reset the configuration of this object.
  ///
  /// @note Calling this function will trigger a re-initialization of the
  ///       magnetic field map.
  void
  setConfiguration(const Config& config)
  {
    m_config = config;
  }

private:
  /// @brief determine bin numbers for given position along each axis
  ///
  /// @param [in]  pos  global position
  /// @param [out] xBin bin number along global x-axis
  /// @param [out] yBin bin number along global y-axis
  /// @param [out] zBin bin number along global z-axis
  ///
  /// @pre #insideGrid(@p pos) must be @c true
  ///
  /// @post @p 0 <= xBin < #m_NPointsX - 1
  /// @post @p 0 <= yBin < #m_NPointsY - 1
  /// @post @p 0 <= zBin < #m_NPointsZ - 1
  ///
  /// This function determines in which bin along each global coordinate axis
  /// the given point is contained. Bin numbers start with 0 for the first bin.
  /// In case the given point is below the minimum/above the maximum value
  /// covered by field grid along one axis, the first/last bin along
  /// this axis is returned.
  void
  getBinNumbers(const Vector3D& pos,
                size_t&         xBin,
                size_t&         yBin,
                size_t&         zBin) const
  {
    const auto x_it
        = std::upper_bound(m_xPoints.begin(), m_xPoints.end(), pos.x());
    const auto y_it
        = std::upper_bound(m_yPoints.begin(), m_yPoints.end(), pos.y());
    const auto z_it
        = std::upper_bound(m_zPoints.begin(), m_zPoints.end(), pos.z());

    xBin = std::distance(m_xPoints.begin(), x_it) - 1;
    yBin = std::distance(m_yPoints.begin(), y_it) - 1;
    zBin = std::distance(m_zPoints.begin(), z_it) - 1;
  }

  /// @brief calculate global index for accessing internal storage
  ///
  /// @param [in] xBin bin along global x-axis within given field map
  /// @param [in] yBin bin along global y-axis within given field map
  /// @param [in] zBin bin along global z-axis within given field map
  ///
  /// @pre 0 <= @p xBin < #m_NPointsX
  /// @pre 0 <= @p yBin < #m_NPointsY
  /// @pre 0 <= @p zBin < #m_NPointsZ
  ///
  /// The magnetic field values are stored internally in a vector. For reasons
  /// of data locality, the magnetic field components at the same one point are
  /// stored consecutively in this vector, that is
  /// @code
  /// B = {Bx(0), By(0), Bz(0), Bx(1), By(1), Bz(1), Bx(2), By(2), Bz(2), ... }
  /// @endcode
  /// The triple of bin numbers along each global axis (@p xBin,@p yBin, @p
  /// zBin) can be mapped onto a global index. The mapping is done such that the
  /// returned global index points to the x-component of the magnetic field at
  /// the point \f$(x_\mathrm{low},y_\mathrm{low},z_\mathrm{low})\f$.
  ///
  /// @return global index
  size_t
  globalIndex(size_t xBin, size_t yBin, size_t zBin) const
  {
    return 3 * (xBin * (m_NPointsY * m_NPointsZ) + yBin * m_NPointsZ + zBin);
  }

  /// @brief initialize internal field map
  ///
  /// @param [in] config configuration object
  ///
  /// The internal storage for the field map is cleared and re-initialized with
  /// the information from the configuration object provided. Please refer to
  /// InterpolatedBFieldMap::InterpolatedBFieldMap for a documentation of the
  /// format for the file specified in Config::fieldMapFile.
  void
  initialize(const Config& config)
  {
    m_xPoints.clear();
    m_yPoints.clear();
    m_zPoints.clear();
    m_BField.clear();

    m_xPoints.reserve(m_config.nPoints);
    m_yPoints.reserve(m_config.nPoints);
    m_zPoints.reserve(m_config.nPoints);
    m_BField.reserve(3 * m_config.nPoints);

    std::ifstream map_file(m_config.fieldMapFile.c_str(), std::ios::in);
    std::string   line;
    double        x = 0, y = 0, z = 0;
    double        bx, by, bz;
    while (std::getline(map_file, line)) {
      if (line.empty() || line[0] == '#'
          || line.find_first_not_of(' ') == std::string::npos)
        continue;

      std::istringstream tmp(line);
      tmp >> x >> y >> z >> bx >> by >> bz;
      m_xPoints.push_back(x * m_config.lengthUnit);
      m_yPoints.push_back(y * m_config.lengthUnit);
      m_zPoints.push_back(z * m_config.lengthUnit);
      m_BField.push_back(bx * m_config.BFieldUnit);
      m_BField.push_back(by * m_config.BFieldUnit);
      m_BField.push_back(bz * m_config.BFieldUnit);
    }
    map_file.close();

    std::sort(m_xPoints.begin(), m_xPoints.end());
    std::sort(m_yPoints.begin(), m_yPoints.end());
    std::sort(m_zPoints.begin(), m_zPoints.end());

    m_xPoints.erase(std::unique(m_xPoints.begin(), m_xPoints.end()),
                    m_xPoints.end());
    m_yPoints.erase(std::unique(m_yPoints.begin(), m_yPoints.end()),
                    m_yPoints.end());
    m_zPoints.erase(std::unique(m_zPoints.begin(), m_zPoints.end()),
                    m_zPoints.end());

    m_xPoints.shrink_to_fit();
    m_yPoints.shrink_to_fit();
    m_zPoints.shrink_to_fit();

    m_NPointsX = m_xPoints.size();
    m_NPointsY = m_yPoints.size();
    m_NPointsZ = m_zPoints.size();
  }

  /// @brief check that point is inside the 3D field map grid
  ///
  /// @param [in] pos global position
  ///
  /// This function checks whether the given position is inside the 3D region
  /// covered by the magnetic field map. A point is considered inside if
  /// \f$x_\mathrm{min} \le x < x_\mathrm{max}\f$ and respectively for the y-
  /// and z-coordinate.
  ///
  /// @return @c true if point is inside the field map region, otherwise @c
  ///         false
  bool
  insideGrid(const Vector3D& pos) const
  {
    return (m_xPoints.at(0) <= pos.x()
            && pos.x() < m_xPoints.at(m_NPointsX - 1))
        && (m_yPoints.at(0) <= pos.y()
            && pos.y() < m_yPoints.at(m_NPointsY - 1))
        && (m_zPoints.at(0) <= pos.z()
            && pos.z() < m_zPoints.at(m_NPointsZ - 1));
  }

  /// @brief perform trilinear interpolation of magnetic field vectors
  ///
  /// @param [in] pos global position
  ///
  /// @pre #insideGrid(@p pos) must be @c true
  ///
  /// This function performs the trilinear interpolation using information from
  /// the magnetic field map to obtain an estimate of the magnetic field vector
  /// at the given position. In a first step, the corresponding @a bin in the
  /// field map is determined. As a second step, the magnetic field vectors of
  /// the eight corner points of this point are used for the trilinear
  /// interpolation.
  ///
  /// @return magnetic field vector at given position
  Vector3D
  interpolate(const Vector3D& pos) const
  {
    typedef Eigen::Map<const Vector3D> MappedVector3D;

    size_t xBin = 0, yBin = 0, zBin = 0;
    getBinNumbers(pos, xBin, yBin, zBin);

    const double xd = (pos.x() - m_xPoints.at(xBin))
        / (m_xPoints.at(xBin + 1) - m_xPoints.at(xBin));
    const double yd = (pos.y() - m_yPoints.at(yBin))
        / (m_yPoints.at(yBin + 1) - m_yPoints.at(yBin));
    const double zd = (pos.z() - m_zPoints.at(zBin))
        / (m_zPoints.at(zBin + 1) - m_zPoints.at(zBin));

    MappedVector3D B1(&m_BField.at(globalIndex(xBin, yBin, zBin)));
    MappedVector3D B2(&m_BField.at(globalIndex(xBin, yBin, zBin + 1)));
    MappedVector3D B3(&m_BField.at(globalIndex(xBin, yBin + 1, zBin)));
    MappedVector3D B4(&m_BField.at(globalIndex(xBin, yBin + 1, zBin + 1)));
    MappedVector3D B5(&m_BField.at(globalIndex(xBin + 1, yBin, zBin)));
    MappedVector3D B6(&m_BField.at(globalIndex(xBin + 1, yBin, zBin + 1)));
    MappedVector3D B7(&m_BField.at(globalIndex(xBin + 1, yBin + 1, zBin)));
    MappedVector3D B8(&m_BField.at(globalIndex(xBin + 1, yBin + 1, zBin + 1)));

    const auto V1 = (1 - zd) * B1 + zd * B2;
    const auto V2 = (1 - zd) * B3 + zd * B4;
    const auto V3 = (1 - zd) * B5 + zd * B6;
    const auto V4 = (1 - zd) * B7 + zd * B8;

    const auto W1 = (1 - yd) * V1 + yd * V2;
    const auto W2 = (1 - yd) * V3 + yd * V4;

    return m_config.scale * ((1 - xd) * W1 + xd * W2);
  }

  /// @brief configuration object
  Config m_config;
  /// number of grid points along global x-axis
  size_t m_NPointsX;
  /// number of grid points along global y-axis
  size_t m_NPointsY;
  /// number of grid points along global z-axis
  size_t m_NPointsZ;
  /// position of grid points on global x-axis
  std::vector<double> m_xPoints;
  /// position of grid points on global y-axis
  std::vector<double> m_yPoints;
  /// position of grid points on global z-axis
  std::vector<double> m_zPoints;
  /// components of magnetic field on grid points
  std::vector<double> m_BField;
};

}  // namespace Acts
