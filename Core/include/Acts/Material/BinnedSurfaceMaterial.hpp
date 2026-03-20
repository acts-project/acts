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
class BinnedSurfaceMaterial : public ISurfaceMaterial {
 public:
  /// Default Constructor - deleted
  BinnedSurfaceMaterial() = delete;


  ///Alias for Grid Level values of elementZ vector
  /// [bin0][bin1] -> Z values per bin
  using ElementZMatrix = 
      std::vector<std::vector<std::vector<unsigned int>>>;

  ///Alias for Grid Level values of ElementFrac vetcor
  using ElementFracMatrix = 
      std::vector<std::vector<std::vector<float>>>;


  /// Explicit constructor with only full MaterialSlab,
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
  /// @param mappingType is the type of surface mapping associated to the surface
  BinnedSurfaceMaterial(const BinUtility& binUtility,
                        MaterialSlabVector fullProperties,
                        double splitFactor = 0.,
                        MappingType mappingType = MappingType::Default,
			ElementZMatrix elementZ = {},
			ElementFracMatrix elementFrac = {});

  /// Explicit constructor with only full MaterialSlab,
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
  /// @param mappingType is the type of surface mapping associated to the surface
  BinnedSurfaceMaterial(const BinUtility& binUtility,
                        MaterialSlabMatrix fullProperties,
                        double splitFactor = 0.,
                        MappingType mappingType = MappingType::Default,
			ElementZMatrix elementZ = {},
			ElementFracMatrix elementFrac = {});

  /// Copy Move Constructor
  ///
  /// @param bsm is the source object to be copied
  BinnedSurfaceMaterial(BinnedSurfaceMaterial&& bsm) = default;

  /// Copy Constructor
  ///
  /// @param bsm is the source object to be copied
  BinnedSurfaceMaterial(const BinnedSurfaceMaterial& bsm) = default;

  /// Assignment Move operator
  /// @param bsm The source object to move from
  /// @return Reference to this object after move assignment
  BinnedSurfaceMaterial& operator=(BinnedSurfaceMaterial&& bsm) = default;

  /// Assignment operator
  /// @param bsm The source object to copy from
  /// @return Reference to this object after copy assignment
  BinnedSurfaceMaterial& operator=(const BinnedSurfaceMaterial& bsm) = default;

  /// Destructor
  ~BinnedSurfaceMaterial() override = default;

  /// Scale operation
  ///
  /// @param factor is the scale factor for the full material
  /// @return Reference to this object after scaling
  BinnedSurfaceMaterial& scale(double factor) final;

  /// Return the BinUtility
  /// @return Reference to the bin utility used for material binning
  const BinUtility& binUtility() const;

  /// @brief Retrieve the entire material slab matrix
  /// @return Reference to the complete matrix of material slabs
  const MaterialSlabMatrix& fullMaterial() const;

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector2&) const
  const MaterialSlab& materialSlab(const Vector2& lp) const final;

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector3&) const
  const MaterialSlab& materialSlab(const Vector3& gp) const final;

  using ISurfaceMaterial::materialSlab;

  ///Return 3D Grid of element Z values
  ///@return Reference to the ElementZMatrix
  const ElementZMatrix& elementZMatrix() const;

  ///Return 3D Grid of fractions for each element Z
  ///@return Reference to the ElementFracMatrix
  const ElementFracMatrix& elementFracMatrix() const;

  ///Return elementZ for a specific bin
  ///@return Reference to the elementZ vector for that bin
  const std::vector<unsigned int>& elementZAt(std::size_t bin0,
		 std::size_t bin1) const;

  ///Return elementFrac for a specific bin
  ///@return Reference to the elementFrac vector for that bin
  const std::vector<float>& elementFracAt(std::size_t bin0,
		std::size_t bin1) const; 


  /// Output Method for std::ostream, to be overloaded by child classes
  /// @param sl The output stream to write to
  /// @return Reference to the output stream after writing
  std::ostream& toStream(std::ostream& sl) const final;

 private:
  /// The helper for the bin finding
  BinUtility m_binUtility;

  /// The five different MaterialSlab
  MaterialSlabMatrix m_fullMaterial;

  ///The Per-bin element list - indexed[bin1][bin0]
  ///empty if using older material recording (element list was not recorded)
  ElementZMatrix m_elementZ;

  ///The Per-bin list of fractions weighted by thickness- sum to 1 per bin
  /// Empty if using older material recording (fraction was not recorded)
  ElementFracMatrix m_elementFrac;
};

inline const BinUtility& BinnedSurfaceMaterial::binUtility() const {
  return m_binUtility;
}

inline const MaterialSlabMatrix& BinnedSurfaceMaterial::fullMaterial() const {
  return m_fullMaterial;
}

inline const BinnedSurfaceMaterial::ElementZMatrix&
BinnedSurfaceMaterial::elementZMatrix() const {
  return m_elementZ;
}

inline const BinnedSurfaceMaterial::ElementFracMatrix&
BinnedSurfaceMaterial::elementFracMatrix() const {
  return m_elementFrac;
}

inline const std::vector<unsigned int>& BinnedSurfaceMaterial::elementZAt(
    std::size_t bin0, std::size_t bin1) const{
  return m_elementZ[bin0][bin1];
}

inline const std::vector<float>& BinnedSurfaceMaterial::elementFracAt(
    std::size_t bin0, std::size_t bin1) const{
  return m_elementFrac[bin0][bin1];
}

}  // namespace Acts
