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

#include <climits>
#include <iomanip>
#include <iostream>
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
  ~MaterialComposition() {}

  /// Constructor from vector of pairs
  ///
  /// @param efracs are the element fractions
  MaterialComposition(
      const std::vector<std::pair<unsigned char, unsigned char>>& efracs)
  {
    reserve(efracs.size());
    for (auto& efracIt : efracs) push_back(efracIt);
  }

  /// Copy constructor from base class 
  /// 
  /// @param mc is the elment fraction vector
  MaterialComposition(const std::vector<ElementFraction>& mc)
    : std::vector<ElementFraction>(mc)
  {
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
    return (*this);
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
  // standard x0, l0, A, Z, rho description
  float X0;
  float L0;
  float A;
  float Z;
  float rho;
  float dEdX;
  float zOaTr;
  MaterialComposition*
      composition;  //! transient member to ROOT (for the moment)

  /// Default Constructor - vacuum material
  Material()
    : X0(10e10)
    , L0(10e10)
    , A(0.)
    , Z(0.)
    , rho(0.)
    , dEdX(0.)
    , zOaTr(0.)
    , composition(0)
  {
  }

  /// Constructor with arguments 
  ///
  /// @param iX0 is the radiation length parameter
  /// @param iL0 is the nuclear interaction length
  /// @param iA is the average atomic weight
  /// @param iZ is the average atomic number
  /// @param iRho is the average density
  /// @param imc is the material composition
  Material(float                iX0,
           float                iL0,
           float                iA,
           float                iZ,
           float                iRho,
           MaterialComposition  imc = {})
    : X0(iX0)
    , L0(iL0)
    , A(iA)
    , Z(iZ)
    , rho(iRho)
    , zOaTr(iA > 0 ? iZ/iA * iRho : 0.)
    , composition(mc)
  {
  }

  /// Copy Constructor 
  ///
  /// @param material copy constructor
  Material(const Material& mat)
    : X0(amc.X0)
    , L0(amc.L0)
    , A(amc.A)
    , Z(amc.Z)
    , rho(amc.rho)
    , zOaTr(amc.zOaTr)
    , composition(amc.composition)
  {
  }

  /// Desctructor 
  ~Material() {}
  
  /// Assignment operator
  ///
  /// @param mat is the source material
  Material&
  operator=(const Material& amc)
  {
    if (this != &amc) {
      X0    = amc.X0;
      L0    = amc.L0;
      A     = amc.A;
      Z     = amc.Z;
      rho   = amc.rho;
      zOaTr = amc.zOaTr;
      composition = amc.composition;
    }
    return (*this);
  }

  /// access to methods z/A*tho
  float
  zOverAtimesRho() const
  {
    return (*this).zOaTr;
  }
  
  /// access to methods x0
  float
  x0() const
  {
    return (*this).X0;
  }

  /// access to methods x0
  float
  averageZ() const
  {
    return (*this).Z;
  }

  /// spit out as a string 
  std::string
  toString() const
  {
    std::ostringstream sout;
    sout << std::setiosflags(std::ios::fixed) << std::setprecision(4);
    sout << "(" << X0 << " | " << L0 << " | " << A << " | " << Z << " | " << rho
         << ")";
    return sout.str();
  }
};

}

#endif
