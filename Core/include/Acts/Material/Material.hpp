// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// Material.hpp, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include <climits>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include "Acts/Material/MaterialComposition.hpp"

namespace Acts {

/// @class Material
///
/// A common object to be contained by
/// - MaterialStep ( for mapping)
/// - MaterialProperties ( for reconstruction )
/// - it is optimized for T/P split
///
///  X0 - radiation length
///  L0 - nuclear interaction length
///  A - nuclear mass
///  Z - nuclear number
///  rho - density
///  Z/A*rho
class Material
{
public:
  /// @brief Accessor enums
  enum Param {
    matX0   = 0,  ///< X0 - this is in mm
    matL0   = 1,  ///< L0 - this is in mm
    matA    = 2,  ///< A - in au
    matZ    = 3,  ///< Z - in e
    matrho  = 4,  ///< rho
    matZ_AR = 5   ///< Z/A*rho
  };

  /// @brief Default Constructor - vacuum material
  Material() = default;

  /// @brief Constructor with arguments
  ///
  /// @param iX0 is the radiation length parameter
  /// @param iL0 is the nuclear interaction length
  /// @param iA is the average atomic weight
  /// @param iZ is the average atomic number
  /// @param iRho is the average density
  /// @param imc is the material composition
  Material(float               iX0,
           float               iL0,
           float               iA,
           float               iZ,
           float               iRho,
           MaterialComposition imc = {})
    : m_store({iX0, iL0, iA, iZ, iRho}), m_composition(imc)
  {
    float zOaTr = (iA > 0. ? iZ / iA * iRho : 0.);
    m_store.push_back(zOaTr);
  }

  /// @brief Copy Constructor
  ///
  /// @param mat copy constructor
  Material(const Material& mat) = default;

  /// @brief Copy Move constructor
  ///
  /// @param  mat copy constructor
  Material(Material&& mat) = default;

  /// @brief Assignment operator
  ///
  /// @param mat is the source material
  Material&
  operator=(const Material& mat)
      = default;

  /// @brief Assignment Move operator
  ///
  /// @param mat is the source material
  Material&
  operator=(Material&& mat)
      = default;

  /// @brief Desctructor
  ~Material() = default;

  /// @brief Equality operator
  ///
  /// @param mat is the source material
  bool
  operator==(const Material& mat) const;

  /// @brief Inequality operator
  ///
  /// @param mat is the source material
  bool
  operator!=(const Material& mat) const;

  /// @brief Boolean operator to check if this is
  /// vacuum has 0 zero size and will indicate false
  operator bool() const { return !m_store.empty(); }

  /// @brief Access to X0
  /// if it's vacuum, infinity
  float
  X0() const
  {
    if (!m_store.empty()) {
      return m_store[matX0];
    } else {
      return std::numeric_limits<float>::infinity();
    }
  }

  /// @brief Access to L0
  /// if it's vacuum, infinity
  float
  L0() const
  {
    if (!m_store.empty()) {
      return m_store[matL0];
    } else {
      return std::numeric_limits<float>::infinity();
    }
  }

  /// @brief Access to A
  float
  A() const
  {
    if (!m_store.empty()) {
      return m_store[matA];
    } else {
      return 0.;
    }
  }

  /// @brief Access to Z
  float
  Z() const
  {
    if (!m_store.empty()) {
      return m_store[matZ];
    } else {
      return 0.;
    }
  }
  /// @brief Access to rho
  float
  rho() const
  {
    if (!m_store.empty()) {
      return m_store[matrho];
    } else {
      return 0.;
    }
  }

  ///  @brief Access to z/A*tho
  float
  zOverAtimesRho() const
  {
    if (m_store.size() > 4) {
      return m_store[matZ_AR];
    }
    return 0.;
  }

  /// spit out as a string
  std::string
  toString() const
  {
    std::ostringstream sout;
    sout << std::setiosflags(std::ios::fixed) << std::setprecision(4);
    sout << " | ";
    for (auto& mat : m_store) {
      sout << mat << " | ";
    }
    return sout.str();
  }

private:
  /// standard x0, l0, A, Z, rho description
  std::vector<float> m_store = {};

  /// optional composition parameter
  MaterialComposition m_composition = MaterialComposition();
};

inline bool
Material::operator==(const Material& mat) const
{
  return (m_store == mat.m_store && m_composition == mat.m_composition);
}

inline bool
Material::operator!=(const Material& mat) const
{
  return !operator==(mat);
}
}