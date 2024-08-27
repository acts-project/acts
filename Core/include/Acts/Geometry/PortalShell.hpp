// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BoundarySurfaceFace.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <memory>

namespace Acts {

class Portal;
class TrackingVolume;
class GeometryContext;

class PortalShellBase {
 public:
  // virtual void apply(TrackingVolume& volume) const = 0;

  virtual ~PortalShellBase() = default;

 private:
};

class CylinderPortalShellBase : public PortalShellBase {
 public:
  // These values are synchronized with the BoundarySurfaceFace enum.
  // Once Gen1 is removed, this can be changed.
  enum class Face : unsigned int {
    PositiveDisc = BoundarySurfaceFace::positiveFaceXY,
    NegativeDisc = BoundarySurfaceFace::negativeFaceXY,
    OuterCylinder = BoundarySurfaceFace::tubeOuterCover,
    InnerCylinder = BoundarySurfaceFace::tubeInnerCover,
    NegativePhiPlane = BoundarySurfaceFace::tubeSectorNegativePhi,
    PositivePhiPlane = BoundarySurfaceFace::tubeSectorPositivePhi
  };

  using enum Face;

  virtual std::size_t size() const = 0;

  virtual Portal* portal(Face face) = 0;
  virtual const std::shared_ptr<Portal>& portalPtr(Face face) = 0;

  virtual void setPortal(std::shared_ptr<Portal> portal, Face face) = 0;
};

class SingleCylinderPortalShell : public CylinderPortalShellBase {
 public:
  explicit SingleCylinderPortalShell(TrackingVolume& volume);

  std::size_t size() const final;

  Portal* portal(Face face) final;
  const std::shared_ptr<Portal>& portalPtr(Face face) final;
  void setPortal(std::shared_ptr<Portal> portal, Face face) final;

 private:
  std::array<std::shared_ptr<Portal>, 6> m_portals{};
};

class CylinderStackPortalShell : public CylinderPortalShellBase {
 public:
  /// @note The shells must be ordered in the given direction
  CylinderStackPortalShell(const GeometryContext& gctx,
                           std::vector<CylinderPortalShellBase*> shells,
                           BinningValue direction,
                           const Logger& logger = getDummyLogger());

  std::size_t size() const final;
  Portal* portal(Face face) final;
  const std::shared_ptr<Portal>& portalPtr(Face face) final;
  void setPortal(std::shared_ptr<Portal> portal, Face face) final;

 private:
  BinningValue m_direction;
  std::vector<CylinderPortalShellBase*> m_shells;
  bool m_hasInnerCylinder{true};
};

}  // namespace Acts
