// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/PortalLinkBase.hpp"

namespace Acts {

class TrackingVolume;
class GridPortalLink;

class TrivialPortalLink final : public PortalLinkBase {
 public:
  TrivialPortalLink(std::shared_ptr<RegularSurface> surface,
                    TrackingVolume& volume)
      : PortalLinkBase(std::move(surface)), m_volume{&volume} {}

  std::unique_ptr<GridPortalLink> makeGrid(BinningValue direction) const;

  void toStream(std::ostream& os) const final {
    os << "TrivialPortalLink<vol=" << m_volume << ">";
  }

  const TrackingVolume* resolveVolume(const GeometryContext& gctx,
                                      const Vector2& position) const final;

  const TrackingVolume* resolveVolume(const GeometryContext& gctx,
                                      const Vector3& position) const final;

 private:
  TrackingVolume* m_volume;
};

}  // namespace Acts
