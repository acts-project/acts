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
class InterpolatedBFieldMap
{
  friend struct Test::InterpolatedBFieldTester;

public:
  /// @brief construct from magnetic field map
  ///
  /// @param [in] sFieldMapFile path to file containing field map
  /// @param [in] scale         global scaling factor
  /// @param [in] lengthUnit    length unit used in field map file
  /// @param [in] BFieldUnit    magnetic field unit used in field map file
  ///
  /// Initialize the object with information read from the field map file.
  ///
  /// @note The format of the field map file must be the following:
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
  InterpolatedBFieldMap(const std::string& sFieldMapFile,
                        double             scale      = 1.,
                        double             lengthUnit = units::_mm,
                        double             BFieldUnit = units::_T)
    : m_dScale(scale)
  {
    m_xPoints.reserve(1000);
    m_yPoints.reserve(1000);
    m_zPoints.reserve(1000);

    std::ifstream map_file(sFieldMapFile.c_str(), std::ios::in);
    std::string   line;
    double        x = 0, y = 0, z = 0;
    double        bx, by, bz;
    while (std::getline(map_file, line)) {
      if (line.empty() || line[0] == '#') continue;

      std::istringstream tmp(line);
      tmp >> x >> y >> z >> bx >> by >> bz;
      m_xPoints.push_back(x * lengthUnit);
      m_yPoints.push_back(y * lengthUnit);
      m_zPoints.push_back(z * lengthUnit);
      m_BField.push_back(bx * BFieldUnit);
      m_BField.push_back(by * BFieldUnit);
      m_BField.push_back(bz * BFieldUnit);
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

    m_NPointsX = m_xPoints.size();
    m_NPointsY = m_yPoints.size();
    m_NPointsZ = m_zPoints.size();
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
    bxyz[0]       = B.y();
    bxyz[0]       = B.z();
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
    return m_dScale;
  }

  /// @brief update global scaling factor for magnetic field
  ///
  /// @param scale new global scaling factor for magnetic field
  ///
  /// @note Negative values for @p scale are accepted and will invert the
  ///       direction of the magnetic field.
  void
  setScale(double scale)
  {
    m_dScale = scale;
  }

private:
  /// @brief determine bin numbers along each axis
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
  /// @param xBin bin along global x-axis within given field map (starts at 0)
  /// @param yBin bin along global y-axis within given field map (starts at 0)
  /// @param zBin bin along global z-axis within given field map (starts at 0)
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

  /// @brief check that point is inside the 3d field map grid
  ///
  /// @param pos global position
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

    return m_dScale * ((1 - xd) * W1 + xd * W2);
  }

  /// global B-field scaling factor
  double m_dScale;
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
