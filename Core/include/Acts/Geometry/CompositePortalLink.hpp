// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Portal.hpp"

namespace Acts {

class CompositePortalLink final : public PortalLinkBase {
 public:
  void toStream(std::ostream& os) const final { os << "CompositePortalLink"; }

  const TrackingVolume* resolveVolume(const GeometryContext& gctx,
                                      const Vector2& position) const final;

  std::unique_ptr<PortalLinkBase> mergeImpl(
      const PortalLinkBase& other, const RegularSurface& surfaceA,
      const RegularSurface& surfaceB, BinningValue direction,
      const Logger& logger = getDummyLogger()) const final;

  std::unique_ptr<PortalLinkBase> mergeImpl(
      const CompositePortalLink& other, const RegularSurface& surfaceA,
      const RegularSurface& surfaceB, BinningValue direction,
      const Logger& logger = getDummyLogger()) const final;

  std::unique_ptr<PortalLinkBase> mergeImpl(
      const TrivialPortalLink& other, const RegularSurface& surfaceA,
      const RegularSurface& surfaceB, BinningValue direction,
      const Logger& logger = getDummyLogger()) const final;

  std::unique_ptr<PortalLinkBase> mergeImpl(
      const GridPortalLink& other, const RegularSurface& surfaceA,
      const RegularSurface& surfaceB, BinningValue direction,
      const Logger& logger = getDummyLogger()) const final;
};

}  // namespace Acts
