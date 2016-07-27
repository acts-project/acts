// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialStep.hpp, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_MATERIAL_MATERIALSTEP_H
#define ACTS_MATERIAL_MATERIALSTEP_H

#include <memory>
#include "ACTS/Material/MaterialProperties.hpp"

namespace Acts {

/// @class MaterialStep
///
/// @brief class holding the material properties at a certain point
///
/// The MaterialStep class is needed to store the material properties (material
/// + step length) at a given
/// global position for the material mapping process.
/// For this class a ROOT dictionary is created in order to store it in a ROOT
/// tree.

class MaterialStep
{
public:
  /// @struct Position
  /// the global three dimensional position of the material step

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
    Position(double x, double y, double z) : x(x), y(y), z(z) {}
  };
  /// Default constructor
  /// setting the position to the origin and making default material properties
  MaterialStep();
  /// Constructor to set th material properties at a certain position
  /// @param mat the material properties (material + step length) at the given
  /// position
  /// @param pos three dimensional global position of the step
  MaterialStep(const MaterialProperties& mat, const Position& pos);
  /// Copy Constructor
  MaterialStep(const MaterialStep& mstep);
  /// Implicit contructor
  /// - uses the copy constructor
  MaterialStep*
  clone() const;
  /// Default Destructor
  ~MaterialStep() = default;
  /// Assignment operator
  MaterialStep&
  operator=(const Acts::MaterialStep& mstep);
  /// return method for the position of the step
  const Position
  position();
  /// return method for the material properties
  const MaterialProperties
  material();

private:
  /// the global three dimensional position of the material step
  Position m_position;
  /// the accumulated material of the step containing the material and the step
  /// length
  MaterialProperties m_material;
};
}

inline const Acts::MaterialStep::Position
Acts::MaterialStep::position()
{
  return m_position;
}
/// return method for the material properties
inline const Acts::MaterialProperties
Acts::MaterialStep::material()
{
  return m_material;
}

#endif  // ACTS_MATERIAL_MATERIALSTEP_H