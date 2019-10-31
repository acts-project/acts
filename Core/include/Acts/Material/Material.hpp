// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <climits>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// Material definition for interactions with matter.
///
/// The following parameters are used to specify the material and its
/// interactions with traversing particles:
///
///     * X0:  radiation length (native length units)
///     * L0:  nuclear interaction length (native length units)
///     * A:   relative atomic mass (unitless number)
///     * Z:   atomic number (unitless number)
///     * rho: density (native mass unit / (native length unit)³)
///
/// The parameters can be effective or average parameters when e.g. a mixture
/// of materials is described.
class Material {
 public:
  /// Index of the parameters in the encoded parameters vector.
  enum Param {
    eX0 = 0,
    eL0 = 1,
    eA = 2,
    eZ = 3,
    eRho = 4,
  };

  /// Construct a vacuum representation.
  Material() = default;
  /// Construct from material parameters.
  ///
  /// @param iX0  is the radiation length parameter
  /// @param iL0  is the nuclear interaction length
  /// @param iA   is the relative atomic mass
  /// @param iZ   is the atomic number
  /// @param iRho is the density
  Material(float iX0, float iL0, float iA, float iZ, float iRho);
  /// Construct from an encoded parameters vector.
  Material(const ActsVectorF<5>& parameters);
  ~Material() = default;

  Material(Material&& mat) = default;
  Material(const Material& mat) = default;
  Material& operator=(Material&& mat) = default;
  Material& operator=(const Material& mat) = default;

  /// Check if the material is valid, i.e. it is not vacuum.
  constexpr operator bool() const { return m_a != 0.0f; }

  /// Return the radition length. Infinity in case of vacuum.
  constexpr float X0() const { return m_x0; }
  /// Return the nuclear interaction length. Infinity in case of vacuum.
  constexpr float L0() const { return m_l0; }
  /// Return the relative atomic mass.
  constexpr float A() const { return m_a; }
  /// Return the atomic number.
  constexpr float Z() const { return m_z; }
  /// Return the density.
  constexpr float rho() const { return m_rho; }

  /// Return the electron density in mol / (native length unit)³.
  ///
  /// Use mol instead of the real number of electrons to avoid large numbers
  /// which could result in numerical instabilities somewhere else.
  constexpr float zOverAtimesRho() const { return m_ne; }

  /// Encode the properties into a parameter vector.
  ActsVectorF<5> classificationNumbers();

  /// spit out as a string
  std::string toString() const;

 private:
  float m_x0 = std::numeric_limits<float>::infinity();
  float m_l0 = std::numeric_limits<float>::infinity();
  float m_a = 0.0f;
  float m_z = 0.0f;
  float m_rho = 0.0f;
  float m_ne = 0.0f;

  friend constexpr bool operator==(const Material& lhs, const Material& rhs) {
    return (lhs.m_x0 == rhs.m_x0) and (lhs.m_l0 == rhs.m_l0) and
           (lhs.m_a == rhs.m_a) and (lhs.m_z == rhs.m_z) and
           (lhs.m_rho == rhs.m_rho) and (lhs.m_ne == rhs.m_ne);
  }
  friend constexpr bool operator!=(const Material& lhs, const Material& rhs) {
    return !(lhs == rhs);
  }
};

}  // namespace Acts
