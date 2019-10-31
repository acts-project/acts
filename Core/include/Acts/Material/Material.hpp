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
    matX0 = 0,
    matL0 = 1,
    matA = 2,
    matZ = 3,
    matrho = 4,
    matZ_AR = 5,
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
  Material(float iX0, float iL0, float iA, float iZ, float iRho)
      : m_vacuum(false), m_store({iX0, iL0, iA, iZ, iRho, 0.}) {
    float zOaTr = (iA > 0. ? iZ / iA * iRho : 0.);
    m_store[5] = zOaTr;
  }
  /// Construct from an encoded parameters vector.
  Material(const ActsVectorF<5>& parameters)
      : Material(parameters[0], parameters[1], parameters[2], parameters[3],
                 parameters[4]) {}
  ~Material() = default;

  Material(Material&& mat) = default;
  Material(const Material& mat) = default;
  Material& operator=(Material&& mat) = default;
  Material& operator=(const Material& mat) = default;

  /// Check if the material is valid, i.e. it is not vacuum.
  operator bool() const { return (!m_vacuum); }

  /// Return the radition length. Infinity in case of vacuum.
  float X0() const { return m_store[matX0]; }
  /// Return the nuclear interaction length. Infinity in case of vacuum.
  float L0() const { return m_store[matL0]; }
  /// Return the relative atomic mass.
  float A() const { return m_store[matA]; }
  /// Return the atomic number.
  float Z() const { return m_store[matZ]; }
  /// Return the density.
  float rho() const { return m_store[matrho]; }

  /// Return the electron density in mol / (native length unit)³.
  ///
  /// Use mol instead of the real number of electrons to avoid large numbers
  /// which could result in numerical instabilities somewhere else.
  float zOverAtimesRho() const { return m_store[matZ_AR]; }

  /// Encode the properties into a parameter vector.
  ActsVectorF<5> classificationNumbers() {
    ActsVectorF<5> numbers;
    numbers << X0(), L0(), A(), Z(), rho();
    return numbers;
  }

  /// spit out as a string
  std::string toString() const {
    std::ostringstream sout;
    sout << std::setiosflags(std::ios::fixed) << std::setprecision(4);
    sout << " | ";
    if (m_vacuum) {
      sout << " vacuum | ";
    } else {
      for (auto& mat : m_store) {
        sout << mat << " | ";
      }
    }
    return sout.str();
  }

 private:
  /// define it is vacuum or not
  bool m_vacuum = true;
  /// standard x0, l0, A, Z, rho description
  std::array<float, 6> m_store = {std::numeric_limits<float>::infinity(),
                                  std::numeric_limits<float>::infinity(),
                                  0.,
                                  0.,
                                  0.,
                                  0.};

  friend bool operator==(const Material& lhs, const Material& rhs) {
    return (lhs.m_store == rhs.m_store);
  }
  friend bool operator!=(const Material& lhs, const Material& rhs) {
    return !(lhs == rhs);
  }
};

}  // namespace Acts
