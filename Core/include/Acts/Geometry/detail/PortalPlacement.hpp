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
/// @brief Implementation of the `SurfacePlacementBase` to synchronize the alignment
///        of an alignable volume with the alignment of the boundary surfaces
///        associated with the Volume.
class PortalPlacement : public SurfacePlacementBase {
 public:
  /// @brief Allow the VolumePlacementBase as only class to construct
  ///        the portal placement
  friend class Acts::VolumePlacementBase;

  /// @brief Returns the transform to switch from the portal reference
  ///       frame into the experiment's global frame taking the alignment
  ///       correction into account
  /// @param gctx The current geometry context object, e.g. alignment
  const Transform3& localToGlobalTransform(
      const GeometryContext& gctx) const final;

  /// @brief Returns the const reference to portal surface connected with
  ///        this placement instance
  const Surface& surface() const final;

  /// @brief Returns the mutable reference to portal surface connected with
  ///        this placement instance
  Surface& surface() final;

  /// @brief Returns the pointer to the hold surface
  std::shared_ptr<RegularSurface> surfacePtr();

  /// @brief Returns the pointer to the hold surface
  std::shared_ptr<const RegularSurface> surfacePtr() const;

  /// @brief Declares the surface object to be non-sensitive
  bool isSensitive() const final;

  /// @brief Returns the face index of the associated portal surface
  std::size_t index() const;

  /// @brief Returns the transform from the portal to the
  ///        associated volume reference frame
  const Transform3& portalToVolumeCenter() const;

  /// @brief Assembles the transform to switch from the portal's frame
  ///        to the experiment's global frame which is essentially the
  ///        same as a call of `localToGlobalTransform` but this method
  ///        is solely intended to be used during the population stage
  ///        of the geometry context
  /// @param gctx The current geometry context object, e.g. alignment
  Transform3 assembleFullTransform(const GeometryContext& gctx) const;

  /// @brief Delete the copy constructor
  PortalPlacement(const PortalPlacement& other) = delete;

  /// @brief Delete the copy assignment operator
  PortalPlacement& operator=(const PortalPlacement& other) = delete;

 protected:
  /// @brief Constructor to instantiate a PortalPlacement
  /// @param portalIdx:
  /// @param portalTrf:
  /// @param parent:
  /// @param surface:
  PortalPlacement(const std::size_t portalIdx, const Transform3& portalTrf,
                  VolumePlacementBase* parent,
                  std::shared_ptr<RegularSurface> surface);

 private:
  Transform3 m_interalTrf{Transform3::Identity()};
  std::shared_ptr<RegularSurface> m_surface{};
  VolumePlacementBase* m_parent{nullptr};

  std::size_t m_portalIdx{0ul};
};
}  // namespace detail
}  // namespace Acts
