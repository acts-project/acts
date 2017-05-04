// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceMaterial.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_MATERIAL_SURFACEMATERIAL_H
#define ACTS_MATERIAL_SURFACEMATERIAL_H

#include "ACTS/Material/MaterialProperties.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include <memory>
#include <vector>

namespace Acts {

class BinUtility;

/// @class SurfaceMaterial
///
/// MaterialProperties that are associated to a surface,
/// extended by certain special representations (binned, homogenous)
///
/// The SurfaceMaterial class inherits from GeometryID,
/// in order to allow storing the material in a file and assigning it uniquely.
///
class SurfaceMaterial
{
public:
  /// Constructor
  SurfaceMaterial() : m_splitFactor(1.) {}

  /// Constructor
  ///
  /// @param splitFactor is the splitting ratio between pre/post update
  SurfaceMaterial(double splitFactor) : m_splitFactor(splitFactor) {}

  /// Destructor
  virtual ~SurfaceMaterial() {}

  /// Pseudo-Constructor clone()
  virtual SurfaceMaterial*
  clone() const = 0;

  /// Scale operator
  ///
  /// @param scale is the scale factor applied
  virtual SurfaceMaterial&
  operator*=(double scale)
      = 0;

  /// Return method for full material description of the Surface
  /// - from local coordinate on the surface
  ///
  /// @param lp is the local position used for the (eventual) lookup
  ///
  /// @retun const MaterialProperties, nullptr indicates no material
  virtual const MaterialProperties*
  material(const Vector2D& lp) const = 0;

  /// Return method for full material description of the Surface
  /// - from the global coordinates
  ///
  /// @param gp is the global position used for the (eventual) lookup
  ///
  /// @retun const MaterialProperties, nullptr indicates no material
  virtual const MaterialProperties*
  material(const Vector3D& gp) const = 0;

  /// Direct access via bins to the MaterialProperties
  ///
  /// @param ib0 is the material bin in dimension 0
  /// @param ib1 is the material bin in dimension 1
  virtual const MaterialProperties*
  material(size_t ib0, size_t ib1) const = 0;

  /// Update pre factor
  double
  factor(PropDirection pDir, MaterialUpdateStage mStage) const;

  /// Output Method for std::ostream, to be overloaded by child classes
  virtual std::ostream&
  dump(std::ostream& sl) const = 0;

protected:
  double m_splitFactor;  //!< the split factor in favour of oppositePre
};

/// inline return methods for the pre/post factors
inline double
SurfaceMaterial::factor(PropDirection pDir, MaterialUpdateStage mStage) const
{
  if (mStage == Acts::fullUpdate) return 1.;
  return (pDir * mStage > 0 ? m_splitFactor : 1. - m_splitFactor);
  }

  /// Overload of << operator for std::ostream for debug output
  std::ostream&
  operator<<(std::ostream& sl, const SurfaceMaterial& sm);

}  // end of namespace

#endif  // ACTS_MATERIAL_SURFACEMATERIAL_H
