// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/BinUtility.hpp"

#include <iosfwd>

namespace Acts {

/// @ingroup material
///
/// It extends the @ref ISurfaceMaterial base class and is an array pf
/// MaterialSlab. This is not memory optimised as every bin
/// holds one material property object.
///
/// The split factors:
///    - 1. : oppositePre
///    - 0. : alongPre
class BinnedSurfaceMaterial : public ISurfaceMaterial {
 public:
  /// Explicit constructor with only full MaterialSlab,
  /// for one-dimensional binning.
  ///
  /// @param binUtility defines the binning structure on the surface (copied)
  /// @param materialVector is the vector of material slabs as recorded (moved)
  /// @param splitFactor is the pre/post splitting directive
  /// @param mappingType is the type of surface mapping associated to the surface
  BinnedSurfaceMaterial(const BinUtility& binUtility,
                        MaterialSlabVector materialVector,
                        double splitFactor = 0.,
                        MappingType mappingType = MappingType::Default);

  /// Explicit constructor with only full MaterialSlab,
  /// for two-dimensional binning.
  ///
  /// @param binUtility defines the binning structure on the surface (copied)
  /// @param materialMatrix is the matrix of material slabs as recorded (moved)
  /// @param splitFactor is the pre/post splitting directive
  /// @param mappingType is the type of surface mapping associated to the surface
  BinnedSurfaceMaterial(const BinUtility& binUtility,
                        MaterialSlabMatrix materialMatrix,
                        double splitFactor = 0.,
                        MappingType mappingType = MappingType::Default);

  /// Scale operation
  ///
  /// @param factor is the scale factor for the full material
  /// @return Reference to this object after scaling
  BinnedSurfaceMaterial& scale(double factor) final;

  /// Return the BinUtility
  /// @return Reference to the bin utility used for material binning
  const BinUtility& binUtility() const { return m_binUtility; }

  /// @brief Retrieve the entire material slab matrix
  /// @return Reference to the complete matrix of material slabs
  const MaterialSlabMatrix& fullMaterial() const { return m_fullMaterial; }

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector2&) const
  const MaterialSlab& materialSlab(const Vector2& lp) const final;

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector3&) const
  [[deprecated(
      "Use materialSlab(const Vector2& lp) with a prior "
      "Surface::globalToLocal() call instead")]] const MaterialSlab&
  materialSlab(const Vector3& gp) const final;

  using ISurfaceMaterial::materialSlab;

  /// @copydoc ISurfaceMaterial::localAxisDirections() const
  std::vector<AxisDirection> localAxisDirections() const final;

  /// Output Method for std::ostream, to be overloaded by child classes
  /// @param sl The output stream to write to
  /// @return Reference to the output stream after writing
  std::ostream& toStream(std::ostream& sl) const final;

 private:
  /// The helper for the bin finding
  BinUtility m_binUtility;

  /// The five different MaterialSlab
  MaterialSlabMatrix m_fullMaterial;
};

}  // namespace Acts
