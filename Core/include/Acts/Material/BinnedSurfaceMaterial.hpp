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
#include "Acts/Utilities/GridAccessHelpers.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"

#include <array>
#include <iosfwd>
#include <stdexcept>
#include <utility>
#include <vector>

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
  [[deprecated(
      "BinnedSurfaceMaterial BinUtility constructor is deprecated. "
      "Use DirectedProtoAxis constructor instead.")]]
  BinnedSurfaceMaterial(const BinUtility& binUtility,
                        MaterialSlabVector fullProperties,
                        double splitFactor = 0.,
                        MappingType mappingType = MappingType::Default);

  /// Explicit constructor with only full MaterialSlab,
  /// for one-dimensional binning.
  ///
  /// @param directedProtoAxis defines axis direction and axis binning
  /// @param fullProperties is the vector of properties as recorded (moved)
  /// @param splitFactor is the pre/post splitting directive
  /// @param mappingType is the type of surface mapping associated to the surface
  BinnedSurfaceMaterial(DirectedProtoAxis directedProtoAxis,
                        MaterialSlabVector fullProperties,
                        double splitFactor = 0.,
                        MappingType mappingType = MappingType::Default);

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
  [[deprecated(
      "BinnedSurfaceMaterial BinUtility constructor is deprecated. "
      "Use std::array<DirectedProtoAxis, 2> constructor instead.")]]
  BinnedSurfaceMaterial(const BinUtility& binUtility,
                        MaterialSlabMatrix fullProperties,
                        double splitFactor = 0.,
                        MappingType mappingType = MappingType::Default);

  /// Explicit constructor with only full MaterialSlab,
  /// for two-dimensional binning.
  ///
  /// @param directedProtoAxes defines axis directions and axis binnings
  /// @param fullProperties is the matrix of properties as recorded (moved)
  /// @param splitFactor is the pre/post splitting directive
  /// @param mappingType is the type of surface mapping associated to the surface
  BinnedSurfaceMaterial(std::array<DirectedProtoAxis, 2> directedProtoAxes,
                        MaterialSlabMatrix fullProperties,
                        double splitFactor = 0.,
                        MappingType mappingType = MappingType::Default);

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
  /// @return BinUtility used for material binning (created on-the-fly)
  [[deprecated(
      "BinnedSurfaceMaterial::binUtility() is deprecated. "
      "Use directedProtoAxes() instead.")]]
  BinUtility binUtility() const;

  /// Return the DirectedProtoAxis descriptors
  /// @return Reference to the directed proto axes used for material binning
  const std::vector<DirectedProtoAxis>& directedProtoAxes() const;

  /// @brief Retrieve the entire material slab matrix
  /// @return Reference to the complete matrix of material slabs
  const MaterialSlabMatrix& fullMaterial() const;

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector2&) const
  const MaterialSlab& materialSlab(const Vector2& lp) const final;

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector3&) const
  ///
  /// For this class the provided 3D position is expected to already be in the
  /// local surface frame.
  const MaterialSlab& materialSlab(const Vector3& lp3D) const final;

  using ISurfaceMaterial::materialSlab;

  /// Output Method for std::ostream, to be overloaded by child classes
  /// @param sl The output stream to write to
  /// @return Reference to the output stream after writing
  std::ostream& toStream(std::ostream& sl) const final;

 private:
  /// Convert legacy BinUtility setup to DirectedProtoAxis.
  static std::vector<DirectedProtoAxis> convertBinUtility(
      const BinUtility& binUtility);

  /// Correct potentially underflow/overflow IAxis bins into matrix bins.
  static std::size_t correctedBinIndex(const IAxis& axis, double value);

  /// The helper for axis-based bin finding
  std::vector<DirectedProtoAxis> m_directedProtoAxes;

  /// The five different MaterialSlab
  MaterialSlabMatrix m_fullMaterial;
};

inline BinUtility BinnedSurfaceMaterial::binUtility() const {
  BinUtility converted;
  for (const auto& directedProtoAxis : m_directedProtoAxes) {
    converted += BinUtility(BinningData(directedProtoAxis));
  }
  return converted;
}

inline const std::vector<DirectedProtoAxis>&
BinnedSurfaceMaterial::directedProtoAxes() const {
  return m_directedProtoAxes;
}

inline std::vector<DirectedProtoAxis> BinnedSurfaceMaterial::convertBinUtility(
    const BinUtility& binUtility) {
  const auto& binningData = binUtility.binningData();
  if (binningData.empty() || binningData.size() > 2u) {
    throw std::invalid_argument(
        "BinnedSurfaceMaterial supports only 1D/2D binning.");
  }

  std::vector<DirectedProtoAxis> directedProtoAxes;
  directedProtoAxes.reserve(binningData.size());
  for (const auto& bData : binningData) {
    const AxisBoundaryType boundaryType = bData.option == closed
                                              ? AxisBoundaryType::Closed
                                              : AxisBoundaryType::Bound;
    if (bData.type == equidistant) {
      directedProtoAxes.emplace_back(
          bData.binvalue, boundaryType, static_cast<double>(bData.min),
          static_cast<double>(bData.max), bData.bins());
      continue;
    }
    std::vector<double> edges;
    edges.reserve(bData.boundaries().size());
    for (const auto edge : bData.boundaries()) {
      edges.push_back(static_cast<double>(edge));
    }
    directedProtoAxes.emplace_back(bData.binvalue, boundaryType, edges);
  }

  return directedProtoAxes;
}

inline std::size_t BinnedSurfaceMaterial::correctedBinIndex(const IAxis& axis,
                                                            double value) {
  const std::size_t rawBin = axis.getBin(value);
  const std::size_t nBins = axis.getNBins();
  if (nBins == 0u) {
    return 0u;
  }
  if (rawBin == 0u) {
    return 0u;
  }
  if (rawBin > nBins) {
    return nBins - 1u;
  }
  return rawBin - 1u;
}

inline const MaterialSlabMatrix& BinnedSurfaceMaterial::fullMaterial() const {
  return m_fullMaterial;
}

}  // namespace Acts
