// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialProperties.hpp"

namespace Acts {

/// Accumulate material properties from multiple hits/track and multiple tracks.
///
/// This is a helper class for the `SurfaceMaterialMapper` to handle material
/// accumulating and averaging for one surface bin. The accumulation procedure
/// is done in two steps:
///
/// 1.  The per-track store accumulates material steps from one track/particle.
///     Multiple material steps can be assigned to the same bin time by
///     one particle, e.g. if the simulation has created more than one step in
///     the material or if several components are compressed into one
///     description. Multiple steps are treated as if they are passed one after
///     the other.
/// 2.  The total store averages the accumulated material properties over all
///     tracks. Each track contributes equally.
class AccumulatedMaterialProperties {
 public:
  // this class does not have a custom default constructor and thus should not
  // provide any custom default cstors, dstor, or assignment. see ISOCPP C.20.

  /// Add the material to the current per-track store.
  ///
  /// Vacuum steps with a non-zero thickness can be added to account for holes
  /// in material structures.
  ///
  /// @param slab Equivalent material slab for this step
  /// @param pathCorreciton Correction factor due to non-perpendicular incident
  void accumulate(MaterialProperties slab, float pathCorrection = 1);

  /// Add the accumulated material for the current track to the total average.
  ///
  /// This finishes the material accumulation for the current track and resets
  /// the per-track store. Subsequent calls to `.accumulate(...)` will start
  /// accumulating material for a new track.
  ///
  /// Each track contributes equally to the total average regardless of its
  /// measured path within the material. An empty per-track store, i.e.
  /// vanishing per-track material thickness, does not contribute to the total
  /// unless explicitely requested.
  ///
  /// @param useEmptyTrack indicate whether to consider an empty track store
  void trackAverage(bool useEmptyTrack = false);

  /// Return the average material properties from all accumulated tracks.
  ///
  /// Only contains the information up to the last `.trackAverag(...)` call. If
  /// there have been additional calls to `.accumulate(...)` afterwards, the
  /// information is not yet part of the total average.
  ///
  /// @returns Average material properties and the number of contributing tracks
  /// @note The averaged material properties are **always** given for unit
  ///   thickness.
  std::pair<MaterialProperties, unsigned int> totalAverage() const;

 private:
  /// Averaged properties for a single track.
  MaterialProperties m_trackAverage;
  /// Averaged properties over multiple tracks.
  MaterialProperties m_totalAverage;
  // Number of tracks contributing to the total average.
  unsigned int m_totalCount = 0u;
};

}  // namespace Acts
