// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <functional>
#include "ACTS/MagneticField/concept/AnyFieldLookup.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/concept/AnyGrid.hpp"

namespace Acts {

/// @ingroup MagneticField
/// @brief interpolate magnetic field value from field values on a given grid
///
/// This class implements a magnetic field service which is initialized by a
/// field map defined by:
/// - a list of field values on a regular grid in some n-dimensional space,
/// - a transformation of global 3D coordinates onto this n-dimensional space.
///
/// The magnetic field value for a given global position is then determined by:
/// - mapping the position onto the grid,
/// - looking up the magnetic field values on the closest grid points,
/// - doing a linear interpolation of these magnetic field values.
class InterpolatedBFieldMap final
{
public:
  /// @brief struct for mapping global 3D positions to field values
  ///
  /// @tparam DIM dimensionality of magnetic field map
  ///
  /// Global 3D positions are transformed into a @c DIM dimensional vector which
  /// is used to look up the magnetic field value in the underlying field map.
  template <unsigned int DIM>
  struct FieldMapper
  {
  public:
    /// @brief default constructor
    ///
    /// @param transform transformation to map global 3D coordinates onto grid
    /// @param grid      grid storing magnetic field values
    FieldMapper(std::function<ActsVectorD<DIM>(const Vector3D&)> transform,
                concept::AnyGrid<Vector3D, ActsVectorD<DIM>> grid)
      : m_transform(std::move(transform)), m_grid(std::move(grid))
    {
    }

    /// @brief retrieve field at given position
    ///
    /// @param [in] position global 3D position
    /// @return magnetic field value at the given position
    ///
    /// @pre The given @c position must lie within the range of the underlying
    ///      magnetic field map.
    Vector3D
    getField(const Vector3D& position) const
    {
      return m_grid.interpolate(m_transform(position));
    }

    /// @brief check whether given 3D position is inside look-up domain
    ///
    /// @param [in] position global 3D position
    /// @return @c true if position is inside the defined look-up grid,
    ///         otherwise @c false
    bool
    isInside(const Vector3D& position) const
    {
      return m_grid.isInside(m_transform(position));
    }

  private:
    /// geometric transformation applied to global 3D positions
    std::function<ActsVectorD<DIM>(const Vector3D&)> m_transform;

    /// grid storing magnetic field values
    concept::AnyGrid<Vector3D, ActsVectorD<DIM>> m_grid;
  };

  /// @brief configuration object for magnetic field interpolation
  struct Config
  {
    /// @brief global B-field scaling factor
    ///
    /// @note Negative values for @p scale are accepted and will invert the
    ///       direction of the magnetic field.
    double scale = 1.;

    /// @brief object for global coordinate transformation and interpolation
    ///
    /// This object performs the mapping of the global 3D coordinates onto the
    /// field grid and the interpolation of the field values on close-by grid
    /// points.
    concept::AnyFieldLookup<> mapper;
  };

  /// @brief create interpolated magnetic field map
  ///
  /// @param [in] config configuration object
  InterpolatedBFieldMap(Config config) : m_config(std::move(config)) {}

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
  /// @param [in] position global 3D position
  ///
  /// @return magnetic field vector at given position
  Vector3D
  getField(const Vector3D& position) const
  {
    return m_config.mapper.getField(position);
  }

  /// @brief retrieve magnetic field value
  ///
  /// @param [in]  position   global 3D position
  /// @param [out] derivative gradient of magnetic field vector as (3x3) matrix
  /// @return magnetic field vector
  ///
  /// @note currently the derivative is not calculated
  /// @todo return derivative
  Vector3D
  getFieldGradient(const Vector3D& position,
                   ActsMatrixD<3, 3>& derivative) const
  {
    return m_config.mapper.getField(position);
  }

  /// @brief get global scaling factor for magnetic field
  ///
  /// @return global factor for scaling the magnetic field
  double
  getScale() const
  {
    return m_config.scale;
  }

  bool
  isInside(const Vector3D& position) const
  {
    return m_config.mapper.isInside(position);
  }

  /// @brief update configuration
  ///
  /// @param [in] config new configuration object
  void
  setConfiguration(const Config& config)
  {
    m_config = config;
  }

private:
  /// @brief configuration object
  Config m_config;
};

}  // namespace Acts
