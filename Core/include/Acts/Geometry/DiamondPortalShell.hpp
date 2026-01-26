// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/DiamondVolumeBounds.hpp"
#include "Acts/Geometry/PortalShell.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

/// Base class for diamond  shaped portal shells, e.g
/// single volumes with polygon shape or stacked (multiple) volumes (TODO)
class DiamondPortalShell : public PortalShellBase {
 public:
  using Face = DiamondVolumeBounds::Face;

  using enum DiamondVolumeBounds::Face;

  /// Retrieve a shared_ptr for the portal associated to the given face. Can be
  /// nullptr if unset.
  /// @param face The face to retrieve the portal for
  /// @return The portal associated to the face
  virtual std::shared_ptr<Portal> portalPtr(Face face) = 0;

  /// Set the portal associated to the given face.
  /// @param portal The portal to set
  /// @param face The face to set the portal
  virtual void setPortal(std::shared_ptr<Portal> portal, Face face) = 0;

  /// @copydoc PortalShellBase::fill
  void fill(TrackingVolume& volume) override;
};
// Output stream operator for the CuboidPortalShell::Face enum
/// @param os The output stream
/// @param face The face to output
/// @return The output stream
std::ostream& operator<<(std::ostream& os, DiamondPortalShell::Face face);

/// Implementation of a portal shell class for a single convex polygon volume
class SingleDiamondPortalShell : public DiamondPortalShell {
 public:
  /// Constructor of a convex polygon shape portal shell for the given volume
  /// @param volume The tracking volume this portal shell is associated with
  explicit SingleDiamondPortalShell(TrackingVolume& volume);

  /// @copydoc DiamondPortalShell::portalPtr
  std::shared_ptr<Portal> portalPtr(Face face) override;

  /// @copydoc DiamondPortalShell::setPortal
  void setPortal(std::shared_ptr<Portal> portal, Face face) override;

  /// @copydoc PortalShellBase::applyToVolume
  void applyToVolume() override;

  /// @copydoc PortalShellBase::isValid
  bool isValid() const override;

  /// @copydoc PortalShellBase::size
  std::size_t size() const override;

  /// @copydoc PortalShellBase::label
  std::string label() const override;

 private:
  std::array<std::shared_ptr<Portal>, 8> m_portals;

  TrackingVolume* m_volume{nullptr};
};
}  // namespace Acts
