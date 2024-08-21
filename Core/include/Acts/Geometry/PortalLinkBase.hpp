// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>

namespace Acts {

class RegularSurface;
class TrackingVolume;
class GeometryContext;

class PortalLinkBase {
 public:
  PortalLinkBase(std::shared_ptr<RegularSurface> surface)
      : m_surface(std::move(surface)) {}

  virtual ~PortalLinkBase() = default;

  // @TODO: Does this need boundary tolerance?
  virtual const TrackingVolume* resolveVolume(
      const GeometryContext& gctx, const Vector3& position) const = 0;

  virtual const TrackingVolume* resolveVolume(
      const GeometryContext& gctx, const Vector2& position) const = 0;

  static std::unique_ptr<PortalLinkBase> merge(
      std::unique_ptr<PortalLinkBase> a, std::unique_ptr<PortalLinkBase> b,
      BinningValue direction, const Logger& logger = getDummyLogger());

  virtual void toStream(std::ostream& os) const = 0;

  friend std::ostream& operator<<(std::ostream& os,
                                  const PortalLinkBase& link) {
    link.toStream(os);
    return os;
  }

  const RegularSurface& surface() const { return *m_surface; }
  void setSurface(std::shared_ptr<RegularSurface> surface) {
    m_surface = std::move(surface);
  }

  std::shared_ptr<RegularSurface> surfacePtr() const { return m_surface; }

 protected:
  static void checkMergePreconditions(const PortalLinkBase& a,
                                      const PortalLinkBase& b,
                                      BinningValue direction);

  std::shared_ptr<RegularSurface> m_surface;
};

}  // namespace Acts
