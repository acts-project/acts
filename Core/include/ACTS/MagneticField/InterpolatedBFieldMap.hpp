// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include "ACTS/Utilities/Units.hpp"

namespace Acts {

/// @brief magnetic field read from map
///
/// This class implements a magnetic field service which is initialized by a
/// field map. The field map needs to be provided upon construction and only
/// the overall normalization can be adjusted afterwards.
class MappedBField
{
public:
  /// @brief construct from magnetic field map
  ///
  /// Initialize the object with information read from the field map file.
  ///
  /// @note The format of the field map file must be the following:
  /// - Every row contains the coordinates of a point and the magnetic field
  ///   vector at this point.
  /// - The delimiter between numbers is a whitespace and the dot is used as
  ///   decimal point.
  /// - The order in every row is: <tt>x y z Bx By Bz</tt>
  /// - No comments nor empty rows are allowed.
  /// - The points must form an equidistant grid in 3D.
  ///
  /// @param [in] sFieldMapFile path to file containing field map
  /// @param [in] scale         global scaling factor
  /// @param [in] lengthUnit    length unit used in field map file
  /// @param [in] BFieldUnit    magnetic field unit used in field map file
  MappedBField(const std::string& sFieldMapFile,
               double             scale      = 1.,
               double             lengthUnit = units::_mm,
               double             BFieldUnit = units::_T)
    : m_dScale(scale)
  {
    std::ifstream map_file(sFieldMapFile.c_str(), std::ios::in);

    m_xPoints.reserve(1000);
    m_yPoints.reserve(1000);
    m_zPoints.reserve(1000);

    double x, y, z;
    double bx, by, bz;
    while (map_file.good()) {
      map_file >> x >> y >> z >> bx >> by >> bz;
      m_xPoints.push_back(x * lengthUnit);
      m_yPoints.push_back(y * lengthUnit);
      m_zPoints.push_back(z * lengthUnit);
      m_BField.push_back(bx * BFieldUnit);
      m_BField.push_back(by * BFieldUnit);
      m_BField.push_back(bz * BFieldUnit);
    }

    m_NPointsX = m_xPoints.size();
    m_NPointsY = m_yPoints.size();
    m_NPointsZ = m_zPoints.size();
  }

  /// @brief retrieve magnetic field value
  ///
  /// @param [in]  xyz   global position
  /// @param [out] bxyz  magnetic field vector
  void
  getField(const double* xyz, double* bxyz) const
  {
    const auto x_it
        = std::lower_bound(m_xPoints.begin(), m_xPoints.end(), xyz[0]);
    const auto y_it
        = std::lower_bound(m_yPoints.begin(), m_yPoints.end(), xyz[1]);
    const auto z_it
        = std::lower_bound(m_zPoints.begin(), m_zPoints.end(), xyz[2]);

    const size_t xBin = std::distance(m_xPoints.begin(), x_it);
    const size_t yBin = std::distance(m_yPoints.begin(), y_it);
    const size_t zBin = std::distance(m_zPoints.begin(), z_it);

    const size_t globalBin = globalBin(xBin, yBin, zBin);
    bxyz[0]                = m_dScale * m_BField.at(globalBin);
    bxyz[1]                = m_dScale * m_BField.at(globalBin + 1);
    bxyz[2]                = m_dScale * m_BField.at(globalBin + 2);
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
  /// @note Negative values for @p scale are accepted and will invert the
  ///       direction of the magnetic field.
  ///
  /// @param scale new global scaling factor for magnetic field
  void
  setScale(double scale)
  {
    m_dScale = scale;
  }

private:
  /// @brief calculate global bin number for accessing internal storage
  ///
  /// The magnetic field values are stored internally in a vector. For reasons
  /// of data locality, the magnetic field components at one point are stored
  /// consecutively in this vector, that is
  /// @code
  /// B = {Bx(0), By(0), Bz(0), Bx(1), By(1), Bz(1), Bx(2), By(2), Bz(2), ... }
  /// @endcode
  /// The triple of bin numbers along each global axis (@p xBin,@p yBin, @p
  /// zBin) can be mapped onto a global bin number. The mapping is done such
  /// that the returned global bin number points to the x-component of magnetic
  /// field at the point \f$(x_{low},y_{low},z_{low})\f$. From the total
  /// number of bins along each axis, one can then compute the indices for the
  /// magnetic field components at all eight corner points.
  ///
  /// @param xBin bin along global x-axis within given field map (starts at 0)
  /// @param yBin bin along global y-axis within given field map (starts at 0)
  /// @param zBin bin along global z-axis within given field map (starts at 0)
  ///
  /// @return global bin number
  size_t
  globalBin(size_t xBin, size_t yBin, size_t zBin) const
  {
    return xBin * (m_NPointsY * m_NPointsZ) + yBin * m_NPointsZ + zBin;
  }

  /// global B-field scaling factor
  double m_dScale;
  /// number of grid points on global x-axis
  size_t m_NPointsX;
  /// number of grid points on global y-axis
  size_t m_NPointsY;
  /// number of grid points on global z-axis
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
