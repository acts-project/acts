// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// Material.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_MATERIAL_MATERIAL_H
#define ACTS_MATERIAL_MATERIAL_H

#include <algorithm>
#include <climits>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace Acts {

static double s_oneOverUcharMax = 1. / double(UCHAR_MAX);

/// @class ElementFraction
class ElementFraction : public std::pair<unsigned char, unsigned char>
{
public:
  /// Default Constructor
  ElementFraction() : std::pair<unsigned char, unsigned char>(0, 0) {}

  /// Copy Constructor from base class
  ///
  /// @param ef is the element fraction source object
  ElementFraction(const std::pair<unsigned char, unsigned char>& ef)
    : std::pair<unsigned char, unsigned char>(ef)
  {
  }

  /// Constructor from arguments
  ///
  /// @param iz is the z value of the element as an unsigned int
  /// @param ifrac is the associated fraction of that element
  ElementFraction(unsigned int iz, float ifrac)
    : std::pair<unsigned char, unsigned char>(
          (unsigned char)iz,
          (unsigned char)(ifrac * double(UCHAR_MAX)))
  {
  }

  /// Constructor from arguments
  ///
  /// @param iz is the z value of the element
  /// @param ifrac is the associated fraction of that element
  ElementFraction(unsigned char iz, unsigned char ifrac)
    : std::pair<unsigned char, unsigned char>(iz, ifrac)
  {
  }

  /// assignment operator from base class
  ///
  /// @param ef is the element fraction source object
  ElementFraction&
  operator=(const std::pair<unsigned char, unsigned char>& ef)
  {
    if (this != &ef) {
      std::pair<unsigned char, unsigned char>::operator=(ef);
    }
    return (*this);
  }

  /// Return in a nice format
  /// - casts back to an unsigned integer
  unsigned int
  element() const
  {
    return static_cast<unsigned int>((*this).first);
  }

  /// Return in a nice format
  /// - casts char to an unsigned int and then into double
  double
  fraction() const
  {
    return (static_cast<unsigned int>((*this).second)) * s_oneOverUcharMax;
  }

  /// Define the equality operator
  bool
  operator==(const ElementFraction& ef) const
  {
    return ((*this).first == ef.first && (*this).second == ef.second);
  }

  /// Define smaller operator for sorting
  /// we always sort by fraction for fastest access to the
  /// most probable fraction
  bool
  operator<(const ElementFraction& ef) const
  {
    return ((*this).second < ef.second);
  }
};

/// @struct MaterialComposition
///
/// This helper struct allows to create a material composition
/// as terms of element fraction objects
class MaterialComposition : public std::vector<ElementFraction>
{
public:
  /// Default constructor
  MaterialComposition() : std::vector<ElementFraction>() {}

  /// Destructor
  ~MaterialComposition() {}

  /// Constructor from vector of pairs
  ///
  /// @param efracs are the element fractions
  MaterialComposition(
      const std::vector<std::pair<unsigned char, unsigned char>>& efracs)
  {
    reserve(efracs.size());
    for (auto& efracIt : efracs) push_back(efracIt);
    std::sort(begin(), end());
  }

  /// Constructor from ElementFraction
  ///
  /// @param mc is the elment fraction vector
  MaterialComposition(const std::vector<ElementFraction>& mc)
    : std::vector<ElementFraction>(mc)
  {
    std::sort(begin(), end());
  }

  /// Assignment operator from base class
  ///
  /// @param mc is the source object
  MaterialComposition&
  operator=(const std::vector<ElementFraction>& mc)
  {
    if (this != &mc) {
      std::vector<ElementFraction>::operator=(mc);
    }
    std::sort(begin(), end());
    return (*this);
  }

  /// Euality operator
  bool
  operator==(const std::vector<ElementFraction>& mc) const
  {
    if (mc.size() != size()) return false;
    for (size_t ef = 0; ef < mc.size(); ++ef) {
      if (!(mc[ef] == (*this)[ef])) return false;
    }
    return true;
  }
};

/// @class Material
///
/// A common object to be contained by
/// - MaterialStep ( for mapping)
/// - MaterialProperties ( for reconstruction )
/// - it is optimized for T/P split
///
class Material
{
public:
  /// Accessor enums
  enum Param {
    matX0   = 0,  ///< Z0
    matL0   = 1,  ///< L0
    matA    = 2,  ///< A
    matZ    = 3,  ///< Z
    matrho  = 4,  ///< rho
    matZ_AR = 5   ///< Z/A*rho
  };

  /// Default Constructor - vacuum material
  Material() : m_store(), m_composition() {}

  /// Constructor with arguments
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
    float zOaTr = (iA > 0 ? iZ / iA * iRho : 0.);
    m_store.push_back(zOaTr);
  }

  /// Copy Constructor
  ///
  /// @param material copy constructor
  Material(const Material& mat)
    : m_store(mat.m_store), m_composition(mat.m_composition)
  {
  }

  /// Desctructor
  ~Material() {}

  /// boolean operator to check if this is
  /// vacuum has 0 zero size and will indicate false
  operator bool() const { return m_store.size(); }

  /// Assignment operator
  ///
  /// @param mat is the source material
  Material&
  operator=(const Material& amc)
  {
    if (this != &amc) {
      m_store       = amc.m_store;
      m_composition = amc.m_composition;
    }
    return (*this);
  }

  /// access to X0
  /// if it's vacuum, infinity
  float
  X0() const
  {
    if (m_store.size())
      return m_store[matX0];
    else
      return std::numeric_limits<float>::infinity();
  }

  /// access to l0
  /// if it's vacuum, infinity
  float
  L0() const
  {
    if (m_store.size())
      return m_store[matL0];
    else
      return std::numeric_limits<float>::infinity();
  }

  /// access to A
  float
  A() const
  {
    if (m_store.size())
      return m_store[matA];
    else
      return 0.;
  }

  /// access to Z
  float
  Z() const
  {
    if (m_store.size())
      return m_store[matZ];
    else
      return 0.;
  }
  /// access to rho
  float
  rho() const
  {
    if (m_store.size())
      return m_store[matrho];
    else
      return 0.;
  }

  /// access to z/A*tho
  float
  zOverAtimesRho() const
  {
    if (m_store.size() > 4) return m_store[matZ_AR];
    return 0.;
  }

  /// spit out as a string
  std::string
  toString() const
  {
    std::ostringstream sout;
    sout << std::setiosflags(std::ios::fixed) << std::setprecision(4);
    sout << " | ";
    for (auto& mat : m_store) sout << mat << " | ";
    return sout.str();
  }

private:
  /// standard x0, l0, A, Z, rho description
  std::vector<float> m_store;

  /// optional composition parameter
  MaterialComposition m_composition;
};

}

#endif
