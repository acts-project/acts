// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/AccumulatedMaterialSlab.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

namespace Acts {

class ISurfaceMaterial;

/// @class AccumulatedSurfaceMaterial
///
/// @ingroup material_mapping
///
/// This class is used by the SurfaceMaterialMapper in order to
/// accumulate/collect material information during the mapping process.
///
/// It performs event- and run-average when called, and returns
/// a new SurfaceMaterial object as a unique_ptr after finalisation
class AccumulatedSurfaceMaterial {
 public:
  /// Type alias for vector of accumulated material slabs
  using AccumulatedVector = std::vector<AccumulatedMaterialSlab>;
  /// Type alias for matrix (vector of vectors) of accumulated material slabs
  using AccumulatedMatrix = std::vector<AccumulatedVector>;

  /// Default Constructor - for homogeneous material
  ///
  /// @param splitFactor is the pre/post splitting directive
  explicit AccumulatedSurfaceMaterial(double splitFactor = 0.);

  /// Explicit constructor with only full MaterialSlab,
  /// for one-dimensional binning.
  ///
  /// The split factors:
  ///    - 1. : oppositePre
  ///    - 0. : alongPre
  ///  ===> 1 Dimensional array
  ///
  /// @param directedProtoAxes defines the binning structure on the surface
  /// @param globalToLocalTransform transform from global to local 3D frame
  /// @param splitFactor is the pre/post splitting directive
  explicit AccumulatedSurfaceMaterial(
      std::vector<DirectedProtoAxis> directedProtoAxes,
      Transform3 globalToLocalTransform = Transform3::Identity(),
      double splitFactor = 0.);

  /// Copy Constructor
  ///
  /// @param asma is the source object to be copied
  AccumulatedSurfaceMaterial(const AccumulatedSurfaceMaterial& asma) = default;

  /// Copy Move Constructor
  ///
  /// @param asma is the source object to be copied
  AccumulatedSurfaceMaterial(AccumulatedSurfaceMaterial&& asma) = default;

  /// Assignment Move operator
  ///
  /// @param asma is the source object to be copied
  /// @return Reference to this object after move assignment
  AccumulatedSurfaceMaterial& operator=(AccumulatedSurfaceMaterial&& asma) =
      default;

  /// Assignment operator
  ///
  /// @param asma is the source object to be copied
  /// @return Reference to this object after copy assignment
  AccumulatedSurfaceMaterial& operator=(
      const AccumulatedSurfaceMaterial& asma) = default;

  /// Destructor
  ~AccumulatedSurfaceMaterial() = default;

  /// Return the DirectedProtoAxis descriptors
  /// @return Reference to the directed proto axes used for material binning
  const std::vector<DirectedProtoAxis>& directedProtoAxes() const;

  /// Return the transform from global to local 3D frame
  /// @return Reference to global-to-local transform
  const Transform3& globalToLocalTransform() const;

  /// Assign a material properties object
  ///
  /// @param lp local position for the bin assignment
  /// @param mp material properties to be assigned
  /// @param pathCorrection Correction factor for the effective path length
  ///
  /// @return the bin triple to which the material was assigned
  std::array<std::size_t, 3> accumulate(const Vector2& lp,
                                        const MaterialSlab& mp,
                                        double pathCorrection = 1.);

  /// Assign a material properties object
  ///
  /// @param gp global position for the bin assignment
  /// @param mp material properties to be assigned
  /// @param pathCorrection Correction factor for the effective path length
  ///
  /// @return the bin triple to which the material was assigned
  std::array<std::size_t, 3> accumulate(const Vector3& gp,
                                        const MaterialSlab& mp,
                                        double pathCorrection = 1.);

  /// Use the accumulated material to update the material variance
  ///
  /// @param trackBins The bins that were touched by this event
  /// @param emptyHit indicator if this is an empty assignment
  /// @param slabReference reference slab (from the map) used to compute the variance
  /// If none is given, the average runs over all bins in the surface map
  void trackVariance(const std::vector<std::array<std::size_t, 3>>& trackBins,
                     MaterialSlab slabReference, bool emptyHit = false);

  /// Use the accumulated material to update the material variance
  ///
  /// @param gp global position for the bin assignment
  /// @param emptyHit indicator if this is an empty assignment
  /// @param slabReference indicator if this is an empty assignment
  void trackVariance(const Vector3& gp, MaterialSlab slabReference,
                     bool emptyHit = false);

  /// Average the information accumulated from one mapped track
  ///
  /// @param trackBins The bins that were touched by this event
  /// @param emptyHit indicator if this is an empty assignment
  /// If none is given, the average runs over all bins in the surface map
  void trackAverage(
      const std::vector<std::array<std::size_t, 3>>& trackBins = {},
      bool emptyHit = false);

  /// Average the information accumulated from one mapped track
  ///
  /// @param gp global position for the bin assignment
  /// @param emptyHit indicator if this is an empty assignment
  void trackAverage(const Vector3& gp, bool emptyHit = false);

  /// Total average creates SurfaceMaterial
  /// @return Unique pointer to the averaged surface material
  std::unique_ptr<const ISurfaceMaterial> totalAverage();

  /// Access to the accumulated material
  /// @return Reference to the matrix of accumulated material data
  const AccumulatedMatrix& accumulatedMaterial() const;

  /// Access to the split factor
  /// @return The split factor used for material averaging
  double splitFactor() const;

 private:
  /// Correct potentially underflow/overflow IAxis bins into matrix bins.
  static std::size_t correctedBinIndex(const IAxis& axis, double value);

  /// Return material bins corresponding to a local 3D position.
  std::array<std::size_t, 3> lookupBins(const Vector3& localPosition) const;

  /// The helper for axis-based bin finding
  std::vector<DirectedProtoAxis> m_directedProtoAxes{};

  /// Global to local transform for bin lookups.
  Transform3 m_globalToLocalTransform = Transform3::Identity();

  /// the split factor
  double m_splitFactor{0.};

  /// The stored accumulated material matrix
  AccumulatedMatrix m_accumulatedMaterial;
};


inline const std::vector<DirectedProtoAxis>&
AccumulatedSurfaceMaterial::directedProtoAxes() const {
  return m_directedProtoAxes;
}

inline const Transform3& AccumulatedSurfaceMaterial::globalToLocalTransform()
    const {
  return m_globalToLocalTransform;
}

inline std::size_t AccumulatedSurfaceMaterial::correctedBinIndex(
    const IAxis& axis, double value) {
  const std::size_t rawBin = axis.getBin(value);
  const std::size_t nBins = axis.getNBins();
  if (nBins == 0u || rawBin == 0u) {
    return 0u;
  }
  if (rawBin > nBins) {
    return nBins - 1u;
  }
  return rawBin - 1u;
}

inline std::array<std::size_t, 3> AccumulatedSurfaceMaterial::lookupBins(
    const Vector3& localPosition) const {
  if (m_directedProtoAxes.empty()) {
    return {0u, 0u, 0u};
  }
  std::array<std::size_t, 3> bins = {0u, 0u, 0u};
  bins[0] = correctedBinIndex(
      m_directedProtoAxes[0u].getAxis(),
      VectorHelpers::cast(localPosition,
                          m_directedProtoAxes[0u].getAxisDirection()));
  if (m_directedProtoAxes.size() > 1u) {
    bins[1] = correctedBinIndex(
        m_directedProtoAxes[1u].getAxis(),
        VectorHelpers::cast(localPosition,
                            m_directedProtoAxes[1u].getAxisDirection()));
  }
  return bins;
}

inline const AccumulatedSurfaceMaterial::AccumulatedMatrix&
AccumulatedSurfaceMaterial::accumulatedMaterial() const {
  return m_accumulatedMaterial;
}

inline double AccumulatedSurfaceMaterial::splitFactor() const {
  return m_splitFactor;
}
}  // namespace Acts
