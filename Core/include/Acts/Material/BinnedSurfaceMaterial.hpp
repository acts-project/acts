// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @class BinnedSurfaceMaterial
///
/// It extends the SurfaceMaterial base class and is an array pf
/// MaterialProperties. This is not memory optimised as every bin
/// holds one material property object.

class BinnedSurfaceMaterial : public ISurfaceMaterial {
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
  /// @param binUtility defines the binning structure on the surface (copied)
  /// @param fullProperties is the vector of properties as recorded (moved)
  /// @param splitFactor is the pre/post splitting directive
  BinnedSurfaceMaterial(const BinUtility& binUtility,
                        MaterialPropertiesVector fullProperties,
                        double splitFactor = 0.);

  /// Explicit constructor with only full MaterialProperties,
  /// for two-dimensional binning.
  ///
  /// The split factors:
  ///    - 1. : oppositePre
  ///    - 0. : alongPre
  ///  ===> 1 Dimensional array
  ///
  /// @param binUtility defines the binning structure on the surface (copied)
  /// @param fullProperties is the vector of properties as recorded (moved)
  /// @param splitFactor is the pre/post splitting directive
  BinnedSurfaceMaterial(const BinUtility& binUtility,
                        MaterialPropertiesMatrix fullProperties,
                        double splitFactor = 0.);

  /// Copy Move Constructor
  ///
  /// @param bsm is the source object to be copied
  BinnedSurfaceMaterial(BinnedSurfaceMaterial&& bsm) = default;

  /// Copy Constructor
  ///
  /// @param bsm is the source object to be copied
  BinnedSurfaceMaterial(const BinnedSurfaceMaterial& bsm) = default;

  /// Assignment Move operator
  BinnedSurfaceMaterial& operator=(BinnedSurfaceMaterial&& bsm) = default;

  /// Assignment operator
  BinnedSurfaceMaterial& operator=(const BinnedSurfaceMaterial& bsm) = default;

  /// Destructor
  ~BinnedSurfaceMaterial() override = default;

  /// Scale operator
  ///
  /// @param scale is the scale factor for the full material
  BinnedSurfaceMaterial& operator*=(double scale) final;

  /// Return the BinUtility
  const BinUtility& binUtility() const;

  /// @copydoc SurfaceMaterial::fullMaterial
  const MaterialPropertiesMatrix& fullMaterial() const;

  /// @copydoc SurfaceMaterial::materialProperties(const Vector2D&)
  const MaterialProperties& materialProperties(const Vector2D& lp) const final;

  /// @copydoc SurfaceMaterial::materialProperties(const Vector3D&)
  const MaterialProperties& materialProperties(const Vector3D& gp) const final;

  /// @copydoc SurfaceMaterial::materialProperties(size_t, size_t)
  const MaterialProperties& materialProperties(size_t bin0,
                                               size_t bin1) const final;

  /// Output Method for std::ostream, to be overloaded by child classes
  std::ostream& toStream(std::ostream& sl) const final;

 private:
  /// The helper for the bin finding
  BinUtility m_binUtility;

  /// The five different MaterialProperties
  MaterialPropertiesMatrix m_fullMaterial;
};

inline const BinUtility& BinnedSurfaceMaterial::binUtility() const {
  return (m_binUtility);
}

inline const MaterialPropertiesMatrix& BinnedSurfaceMaterial::fullMaterial()
    const {
  return m_fullMaterial;
}

inline const MaterialProperties& BinnedSurfaceMaterial::materialProperties(
    size_t bin0, size_t bin1) const {
  return m_fullMaterial[bin1][bin0];
}
}  // namespace Acts
