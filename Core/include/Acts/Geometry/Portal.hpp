// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <memory>

namespace Acts {

class RegularSurface;
class GeometryContext;
class TrackingVolume;
class CylinderSurface;
class PlaneSurface;
class DiscSurface;
class Surface;

class PortalLinkBase;

class PortalMergingException : public std::exception {
  const char* what() const noexcept override {
    return "Failure to merge portals";
  }
};

class PortalFusingException : public std::exception {
  const char* what() const noexcept override {
    return "Failure to fuse portals";
  }
};

class Portal {
 public:
  Portal(const GeometryContext& gctx, Direction direction,
         std::unique_ptr<PortalLinkBase> link);

  Portal(const GeometryContext& gctx, Direction direction,
         std::shared_ptr<RegularSurface> surface, TrackingVolume& volume);

  Portal(const GeometryContext& gctx,
         std::unique_ptr<PortalLinkBase> alongNormal,
         std::unique_ptr<PortalLinkBase> oppositeNormal);

  struct Config {
    struct Link {
      Link() = default;
      Link(std::shared_ptr<RegularSurface> surface, TrackingVolume& volume)
          : m_surface(std::move(surface)), m_volume(&volume) {}

      std::shared_ptr<RegularSurface> m_surface = nullptr;
      TrackingVolume* m_volume = nullptr;
    };

    Link alongNormal;
    Link oppositeNormal;
  };

  Portal(const GeometryContext& gctx, Config&& config);

  ///    portal1   portal2
  ///      +---+   +---+
  ///      |   |   |   |
  ///      |   |   |   |
  /// <----+   | + |   +---->
  ///      |   |   |   |
  ///      |   |   |   |
  ///      +---+   +---+
  /// @note This is a destructive operaion on the portals involved
  /// @TODO: Handle material on portal surfaces
  static Portal fuse(const GeometryContext& gctx, Portal& aPortal,
                     Portal& bPortal, const Logger& logger = getDummyLogger());

  ///         ^                     ^
  ///         |                     |
  ///  portal1|              portal2|
  /// +-------+-------+     +-------+-------+
  /// |               |  +  |               |
  /// +-------+-------+     +-------+-------+
  ///         |                     |
  ///         |                     |
  ///         v                     v
  /// @note This is a destructive operation on both portals, their
  ///       links will be moved to produce merged links, which can fail
  ///       if the portal links are not compatible
  static Portal merge(const GeometryContext& gctx, Portal& aPortal,
                      Portal& bPortal, BinningValue direction,
                      const Logger& logger = getDummyLogger());

  Result<const TrackingVolume*> resolveVolume(const GeometryContext& gctx,
                                              const Vector3& position,
                                              const Vector3& direction) const;

  void setLink(const GeometryContext& gctx, Direction direction,
               std::unique_ptr<PortalLinkBase> link);
  void setLink(const GeometryContext& gctx, Direction direction,
               std::shared_ptr<RegularSurface> surface, TrackingVolume& volume);

  const PortalLinkBase* getLink(Direction direction) const;

  bool isValid() const;

  const RegularSurface& surface() const;

 private:
  static bool isSameSurface(const GeometryContext& gctx, const Surface& a,
                            const Surface& b);

  // @TODO: Potentially short circuit the virtual call
  // using VolumeResolver = Delegate<const TrackingVolume*(const Vector3&
  // position)>;

  std::shared_ptr<RegularSurface> m_surface;

  std::unique_ptr<PortalLinkBase> m_alongNormal;
  std::unique_ptr<PortalLinkBase> m_oppositeNormal;
};

}  // namespace Acts
