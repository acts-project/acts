// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <vector>

#include "Acts/Material/AccumulatedMaterialProperties.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @class AccumulatedSurfaceMaterial
///
/// This class is used by the SurfaceMaterialMapper in order to
/// accumulate/collect material information during the mapping process.
///
/// It performs event- and run-average when called, and returns
/// a new SurfaceMaterial object as a unique_ptr after finalisation
class AccumulatedSurfaceMaterial {
 public:
  using AccumulatedVector = std::vector<AccumulatedMaterialProperties>;
  using AccumulatedMatrix = std::vector<AccumulatedVector>;

  /// Default Constructor - for homogeneous material
  ///
  /// @param splitFactor is the pre/post splitting directive
  AccumulatedSurfaceMaterial(double splitFactor = 0.);

  /// Explicit constructor with only full MaterialProperties,
  /// for one-dimensional binning.
  ///
  /// The split factors:
  ///    - 1. : oppositePre
  ///    - 0. : alongPre
  ///  ===> 1 Dimensional array
  ///
  /// @param binUtility defines the binning structure on the surface
  /// @param splitFactor is the pre/post splitting directive
  AccumulatedSurfaceMaterial(const BinUtility& binUtility,
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
  AccumulatedSurfaceMaterial& operator=(AccumulatedSurfaceMaterial&& asma) =
      default;

  /// Assignment operator
  ///
  /// @param asma is the source object to be copied
  AccumulatedSurfaceMaterial& operator=(
      const AccumulatedSurfaceMaterial& asma) = default;

  /// Destructor
  ~AccumulatedSurfaceMaterial() = default;

  /// Return the BinUtility
  const BinUtility& binUtility() const;

  /// Assign a material properites object
  ///
  /// @param lp local position for the bin assignment
  /// @param mp material properties to be assigned
  ///
  /// @return the bin triple to which the material was assigned
  std::array<size_t, 3> accumulate(const Vector2D& lp,
                                   const MaterialProperties& mp,
                                   double pathCorrection = 1.);

  /// Assign a material properites object
  ///
  /// @param gp global position for the bin assignment
  /// @param mp material properties to be assigned
  ///
  /// @return the bin triple to which the material was assigned
  std::array<size_t, 3> accumulate(const Vector3D& gp,
                                   const MaterialProperties& mp,
                                   double pathCorrection = 1.);

  /// Average the information accumulated from one mapped track
  ///
  /// @param trackBins The bins that were touched by this event
  /// @param emptyHit indicator if this is an empty assignment
  /// If none is given, the average runs over all bins in the surface map
  void trackAverage(const std::vector<std::array<size_t, 3>>& trackBins = {},
                    bool emptyHit = false);

  /// Average the information accumulated from one mapped track
  ///
  /// @param gp global position for the bin assignment
  /// @param emptyHit indicator if this is an empty assignment
  void trackAverage(const Vector3D& gp, bool emptyHit = false);

  /// Total average creates SurfaceMaterial
  std::unique_ptr<const ISurfaceMaterial> totalAverage();

  /// Access to the accumulated material
  const AccumulatedMatrix& accumulatedMaterial() const;

  /// Access to the split factor
  double splitFactor() const;

 private:
  /// The helper for the bin finding
  BinUtility m_binUtility{};

  /// the split factor
  double m_splitFactor{0.};

  /// The stored accumulated material matrix
  AccumulatedMatrix m_accumulatedMaterial;
};

inline const BinUtility& AccumulatedSurfaceMaterial::binUtility() const {
  return (m_binUtility);
}

inline const AccumulatedSurfaceMaterial::AccumulatedMatrix&
AccumulatedSurfaceMaterial::accumulatedMaterial() const {
  return (m_accumulatedMaterial);
}

inline double AccumulatedSurfaceMaterial::splitFactor() const {
  return m_splitFactor;
}
}  // namespace Acts
