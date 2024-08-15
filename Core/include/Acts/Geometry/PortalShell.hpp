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

#include <array>
#include <memory>

namespace Acts {

class Portal;
class TrackingVolume;

class PortalShell {
 public:
  // virtual void apply(TrackingVolume& volume) const = 0;

  virtual ~PortalShell() = default;

 private:
};

class CylinderPortalShell : public PortalShell {
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

  virtual Portal* portal(Face face) = 0;
};

class SingleCylinderPortalShell : public CylinderPortalShell {
 public:
  SingleCylinderPortalShell(TrackingVolume& volume);

  std::size_t size() const {
    std::size_t count = 0;
    std::ranges::for_each(
        m_portals, [&count](const auto& portal) { count += portal ? 1 : 0; });
    return count;
  }

  Portal* portal(Face face) final;

 private:
  std::array<std::shared_ptr<Portal>, 6> m_portals{};
};

class CylinderStackPortalShell : public CylinderPortalShell {
 public:
  /// @note The shells must be ordered in the given direction
  CylinderStackPortalShell(std::vector<CylinderPortalShell*> shells,
                           BinningValue direction);

  Portal* portal(Face face) final;

 private:
  BinningValue m_direction;
  std::vector<CylinderPortalShell*> m_shells;
};

}  // namespace Acts
