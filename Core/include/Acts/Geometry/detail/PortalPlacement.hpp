// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Surfaces/SurfacePlacementBase.hpp"

namespace Acts {
class VolumePlacementBase;

namespace detail {
/// Implementation of the `SurfacePlacementBase` to synchronize the alignment
/// of an alignable volume with the alignment of the boundary surfaces
/// associated with the Volume.
class PortalPlacement final : public SurfacePlacementBase {
 public:
  ///  Allow the VolumePlacementBase as only class to construct the portal
  ///  placement
  friend class Acts::VolumePlacementBase;

  /// Returns the transform to switch from the portal reference
  /// frame into the experiment's global frame taking the alignment
  /// correction into account
  /// @param gctx The current geometry context object, e.g. alignment
  const Transform3& localToGlobalTransform(
      const GeometryContext& gctx) const override;

  ///  Returns the const reference to portal surface connected with this
  ///  placement instance
  const RegularSurface& surface() const override;

  ///  Returns the mutable reference to portal surface connected with this
  ///  placement instance
  RegularSurface& surface() override;

  /// Returns the pointer to the hold surface
  const std::shared_ptr<RegularSurface>& surfacePtr();

  /// Returns the pointer to the hold surface
  std::shared_ptr<const RegularSurface> surfacePtr() const;

  /// Declares the surface object to be non-sensitive
  bool isSensitive() const override;

  ///  Returns the face index of the associated portal surface
  std::size_t index() const;

  /// Returns the transform from the portal to the associated volume reference
  /// frame
  const Transform3& portalToVolumeCenter() const;

  /// Delete the copy constructor
  PortalPlacement(const PortalPlacement& other) = delete;

  /// Delete the copy assignment operator
  PortalPlacement& operator=(const PortalPlacement& other) = delete;

  /// Constructor to instantiate a PortalPlacement
  /// @param portalIdx: Internal index to associated the portal with the
  ///                   i-th boundary surface of the volume
  /// @param portalTrf: Transform from the portal's frame into the
  ///                   volume (I.e. the  orientation of the portal w.r.t.
  ///                   volume)
  /// @param parent: Pointer to the parent which is hosting the placement and also
  ///                providing the override transforms from the portal ->
  ///                experiment's frame
  /// @param surface: Pointer to the portal surface itself which is becoming alignable
  ///                 with the construction of this PortalPlacement
  PortalPlacement(const std::size_t portalIdx, const Transform3& portalTrf,
                  const VolumePlacementBase* parent,
                  std::shared_ptr<RegularSurface> surface);

 private:
  /// Assembles the transform to switch from the portal's frame
  /// to the experiment's global frame which is essentially the
  /// same as a call of `localToGlobalTransform` but this method
  /// is solely intended to be used during the population stage
  /// of the geometry context
  /// @param gctx The current geometry context object, e.g. alignment
  Transform3 assembleFullTransform(const GeometryContext& gctx) const;

  /// Orientation of the portal surface w.r.t the volume
  Transform3 m_portalToVolumeCenter{Transform3::Identity()};
  /// Pointer to the surface held by the placement
  std::shared_ptr<RegularSurface> m_surface{};
  /// Pointer to the parent managing this instance
  const VolumePlacementBase* m_parent{nullptr};
  /// Internal index of the portal
  std::size_t m_portalIdx{0ul};
};
}  // namespace detail
}  // namespace Acts
