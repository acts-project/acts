// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BinnedSurfaceMaterial.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Material/SurfaceMaterial.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @class BinnedSurfaceMaterial
///
/// It extends the SurfaceMaterial base class and is an array pf
/// MaterialProperties. This is not memory optimised as every bin
/// holds one material property object.

class BinnedSurfaceMaterial : public SurfaceMaterial
{
public:
  /// Default Constructor - deleted
  BinnedSurfaceMaterial() = delete;

  /// Explicit constructor with only full MaterialProperties,
  /// for one-dimensional binning.
  ///
  /// The split factors:
  ///    - 1. : oppositePre
  ///    - 0. : alongPre
  ///  ===> 1 Dimensional array
  ///
  /// @param binUtility defines the binning structure on the surface
  /// @param fullProperties is the vector of properties as recorded
  /// @param splitFactor is the pre/post splitting directive
  /// @param entries is the (optional) number of mapping entries
  BinnedSurfaceMaterial(const BinUtility&               binUtility,
                        const MaterialPropertiesVector& fullProperties,
                        double                          splitFactor = 0.);

  /// Explicit constructor with only full MaterialProperties,
  /// for two-dimensional binning.
  ///
  /// The split factors:
  ///    - 1. : oppositePre
  ///    - 0. : alongPre
  ///  ===> 1 Dimensional array
  ///
  /// @param binUtility defines the binning structure on the surface
  /// @param fullProperties is the vector of properties as recorded
  /// @param splitFactor is the pre/post splitting directive
  /// @param entries is the (optional) number of mapping entries
  BinnedSurfaceMaterial(const BinUtility&               binUtility,
                        const MaterialPropertiesMatrix& fullProperties,
                        double                          splitFactor = 0.);

  /// Copy Move Constructor
  ///
  /// @param lmp is the source object to be copied
  BinnedSurfaceMaterial(BinnedSurfaceMaterial&& lmp) = default;

  /// Copy Constructor
  ///
  /// @param lmp is the source object to be copied
  BinnedSurfaceMaterial(const BinnedSurfaceMaterial& lmp) = default;

  /// Assignment Move operator
  BinnedSurfaceMaterial&
  operator=(BinnedSurfaceMaterial&& lmp)
      = default;

  /// Assignment operator
  BinnedSurfaceMaterial&
  operator=(const BinnedSurfaceMaterial& lmp)
      = default;

  /// Destructor
  ~BinnedSurfaceMaterial() override = default;

  /// Scale operator
  ///
  /// @param scale is the scale factor for the full material
  BinnedSurfaceMaterial&
  operator*=(double scale) final;

  /// Return the BinUtility
  const BinUtility&
  binUtility() const;

  /// @copydoc SurfaceMaterial::fullMaterial
  const MaterialPropertiesMatrix&
  fullMaterial() const;

  /// @copydoc SurfaceMaterial::materialProperties(const Vector2D&)
  const MaterialProperties&
  materialProperties(const Vector2D& lp) const final;

  /// @copydoc SurfaceMaterial::materialProperties(const Vector3D&)
  const MaterialProperties&
  materialProperties(const Vector3D& gp) const final;

  /// @copydoc SurfaceMaterial::materialProperties(size_t, size_t)
  const MaterialProperties&
  materialProperties(size_t bin0, size_t bin1) const final;

  /// Output Method for std::ostream, to be overloaded by child classes
  std::ostream&
  dump(std::ostream& sl) const final;

private:
  /// The helper for the bin finding
  BinUtility m_binUtility;

  /// The five different MaterialProperties
  MaterialPropertiesMatrix m_fullMaterial;

};

inline const BinUtility&
BinnedSurfaceMaterial::binUtility() const
{
  return (m_binUtility);
}

inline const MaterialPropertiesMatrix&
BinnedSurfaceMaterial::fullMaterial() const
{
  return m_fullMaterial;
}

inline const MaterialProperties&
BinnedSurfaceMaterial::materialProperties(size_t bin0, size_t bin1) const
{
  return m_fullMaterial[bin1][bin0];
}

}