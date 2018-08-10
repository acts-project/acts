// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialStep.hpp, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <memory>
#include "Acts/Material/MaterialProperties.hpp"

#if !defined(__CLING__)
#include "Acts/Utilities/Definitions.hpp"
#endif

namespace Acts {

/// @class MaterialStep
///
/// @brief class holding the material properties of a certain step
///
/// The MaterialStep class is needed to store the material properties (material
/// and step length) at a given global position for the material mapping
/// process.
///
/// @todo Currently a specific Position struct is used instead of the
/// Acts::Vector3D
/// to simplify the creation of a ROOT dictionary for this class, since
/// Acts::Vector3D is an Eigen class. In future the Vector3D should be used
/// to guarantee consitency and avoid conversions.

class MaterialStep
{
public:
  /// @struct Position
  /// the global three dimensional position of the material step
  /// @todo replace by Acts::Vector3D
  struct Position
  {
    /// X Coordinate of the material step
    double x;
    /// Y Coordinate of the material step
    double y;
    /// Z Coordinate of the material step
    double z;

    /// Default constructor creating a position at the origin
    Position() : x(0.), y(0.), z(0.) {}

    /// Constructor to set the three coordinates
    Position(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    /// Copy Constructor
    Position(const Position& pos) : x(pos.x), y(pos.y), z(pos.z) {}

#if !defined(__CLING__)
    /// Constructor from Vector3D
    Position(const Vector3D& pos) : x(pos.x()), y(pos.y()), z(pos.z()) {}

    /// assignment operator from Vector3D
    Position&
    operator=(const Vector3D& pos)
    {
      x = pos.x();
      y = pos.y();
      z = pos.z();
      return (*this);
    }
#endif
  };

  /// Default constructor
  /// setting the position to the origin and making default material properties
  MaterialStep();

  /// Constructor to set th material properties at a certain position
  /// @param mat the material (+ step length) at the given position
  /// @param pos three dimensional global position of the step
  /// @param geoID is the geoId value (optional)
  MaterialStep(const MaterialProperties& mat,
               const Position&           pos,
               uint64_t                  geoId = 0);

  /// Copy Constructor
  MaterialStep(const MaterialStep& mstep);

  /// Default Destructor
  ~MaterialStep() = default;

  /// Assignment operator
  MaterialStep&
  operator=(const MaterialStep& mstep);

#if !defined(__CLING__)
  const Vector3D
  position() const;
#else
  /// return method for the position of the step
  const Position
  position() const;
#endif

  /// return method for the material properties
  const MaterialProperties&
  materialProperties() const;

  // return the value of the geometry id
  uint64_t
  geoID() const;

private:
  /// the global three dimensional position of the material step
  Position m_position;

  /// the accumulated material of the step
  /// containing the material and the step length
  MaterialProperties m_material;

  /// the geometry id
  uint64_t m_geoID;
};

}  /// namespace

#if !defined(__CLING__)
inline const Acts::Vector3D
Acts::MaterialStep::position() const
{
  return Acts::Vector3D(m_position.x, m_position.y, m_position.z);
}
#else
inline const Acts::MaterialStep::Position
Acts::MaterialStep::position() const
{
  return m_position;
}
#endif

/// return method for the material properties
inline const Acts::MaterialProperties&
Acts::MaterialStep::materialProperties() const
{
  return m_material;
}

inline uint64_t
Acts::MaterialStep::geoID() const
{
  return m_geoID;
}