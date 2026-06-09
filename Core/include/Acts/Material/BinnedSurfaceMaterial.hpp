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
#include "Acts/Utilities/ProtoAxis.hpp"

#include <array>
#include <iosfwd>

namespace Acts {

/// @ingroup material
///
/// It extends the @ref ISurfaceMaterial base class and is an array of
/// MaterialSlab. This is not memory optimised as every bin
/// holds one material property object.
///
/// The binning is described by two @ref DirectedProtoAxis objects. A
/// conceptually 1D surface uses a second axis with a single bin.
///
/// The split factors:
///    - 1. : oppositePre
///    - 0. : alongPre
class BinnedSurfaceMaterial : public ISurfaceMaterial {
 public:
  /// Primary constructor using two directed proto axes (always 2D).
  ///
  /// @param axes defines the 2D binning structure on the surface; for a
  ///             conceptually 1D surface use a second axis with a single bin
  /// @param materialMatrix is the matrix of material slabs (moved)
  /// @param splitFactor is the pre/post splitting directive
  /// @param mappingType is the type of surface mapping associated to the surface
  BinnedSurfaceMaterial(std::array<DirectedProtoAxis, 2> axes,
                        MaterialSlabMatrix materialMatrix,
                        double splitFactor = 0.,
                        MappingType mappingType = MappingType::Default);

  /// @deprecated Use the DirectedProtoAxis-based constructor instead.
  ///
  /// Explicit constructor with only full MaterialSlab,
  /// for one-dimensional binning.
  ///
  /// @param binUtility defines the binning structure on the surface (copied)
  /// @param materialVector is the vector of material slabs as recorded (moved)
  /// @param splitFactor is the pre/post splitting directive
  /// @param mappingType is the type of surface mapping associated to the surface
  [[deprecated(
      "Use BinnedSurfaceMaterial(std::array<DirectedProtoAxis, 2>, "
      "MaterialSlabMatrix) instead")]]
  BinnedSurfaceMaterial(const BinUtility& binUtility,
                        MaterialSlabVector materialVector,
                        double splitFactor = 0.,
                        MappingType mappingType = MappingType::Default);

  /// @deprecated Use the DirectedProtoAxis-based constructor instead.
  ///
  /// Explicit constructor with only full MaterialSlab,
  /// for two-dimensional binning.
  ///
  /// @param binUtility defines the binning structure on the surface (copied)
  /// @param materialMatrix is the matrix of material slabs as recorded (moved)
  /// @param splitFactor is the pre/post splitting directive
  /// @param mappingType is the type of surface mapping associated to the surface
  [[deprecated(
      "Use BinnedSurfaceMaterial(std::array<DirectedProtoAxis, 2>, "
      "MaterialSlabMatrix) instead")]]
  BinnedSurfaceMaterial(const BinUtility& binUtility,
                        MaterialSlabMatrix materialMatrix,
                        double splitFactor = 0.,
                        MappingType mappingType = MappingType::Default);

  /// Scale operation
  ///
  /// @param factor is the scale factor for the full material
  /// @return Reference to this object after scaling
  BinnedSurfaceMaterial& scale(double factor) final;

  /// Return the two directed proto axes describing the binning.
  /// @return const reference to the axes array
  const std::array<DirectedProtoAxis, 2>& axes() const { return m_axes; }

  /// @deprecated Use axes() instead.
  ///
  /// Return a BinUtility reconstructed from the internal axes.
  /// Single-bin dummy axes are omitted so the result matches what was
  /// originally provided to the deprecated BinUtility constructors.
  ///
  /// @return reconstructed BinUtility
  [[deprecated("Use axes() instead")]] BinUtility binUtility() const;

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
  /// The two directed axes defining the binning (axis 0 → lp[0], axis 1 →
  /// lp[1])
  std::array<DirectedProtoAxis, 2> m_axes;

  /// The MaterialSlab matrix (row = axis-1 bin, col = axis-0 bin)
  MaterialSlabMatrix m_fullMaterial;
};

}  // namespace Acts
