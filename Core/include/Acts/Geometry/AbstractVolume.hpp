// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/Volume.hpp"

#include <memory>
#include <vector>

namespace Acts {

class AbstractVolume;
template <class volume_t>
class BoundarySurfaceT;

using BoundarySurfacePtr =
    std::shared_ptr<const BoundarySurfaceT<AbstractVolume>>;

class VolumeBounds;

using VolumeBoundsPtr = std::shared_ptr<const VolumeBounds>;

/// @class AbstractVolume
///
/// AbstractVolume description inside the tracking realm. This is the purely
/// geometrical object volume.
///
/// The Acts::AbstractVolume is constructed by giving a pointer to a Transform3
/// and a pointer to Acts::VolumeBounds, this implies that the ownership of the
/// objects pointed to is passed as well. For memory optimisation, the
/// AbstractVolume can also be constructed with shared_ptr objects.
///
/// A Acts::AbstractVolume is at first a collection class of
/// Acts::BoundarySurface, the vector of Acts::BoundarySurface is returned by
/// the Acts::VolumeBounds that carry a decompose method.
///
/// Boundary surfaces can be shared between AbstractVolumes to enhance automatic
/// navigation between AbstractVolumes, therefor they are reference counted by a
/// std::shared_ptr holder class.

class AbstractVolume : public Volume {
 public:
  /// Constructor with shared Transform3*, VolumeBounds*
  ///
  /// @param transform is the gobal 3d transformation into the volume frame
  /// @param volbounds is the boundary definition
  AbstractVolume(const Transform3& transform, VolumeBoundsPtr volbounds);

  AbstractVolume(const AbstractVolume& vol) = default;

  AbstractVolume() = delete;
  ~AbstractVolume() override = default;
  AbstractVolume& operator=(const AbstractVolume& vol) = delete;

  /// Method to return the BoundarySurfaces
  ///
  /// @return the vector of boundary surfaces
  const std::vector<BoundarySurfacePtr>& boundarySurfaces() const;

 private:
  /// Private method to create BoundarySurfaces
  void createBoundarySurfaces();

  /// boundary Surfaces for this volume
  std::vector<BoundarySurfacePtr> m_boundarySurfaces;
};

}  // namespace Acts