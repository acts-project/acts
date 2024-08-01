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

  CompositePortalLink(const std::shared_ptr<const PortalLinkBase>& a,
                      const std::shared_ptr<const PortalLinkBase>& b)
      : PortalLinkBase(nullptr) {
    (void)a;
    (void)b;
    // @TODO: Implement
    throw std::logic_error{"Not implemented"};
  };

  // @TODO: Implement optimization in case on or both of the arguments are already composites
};

}  // namespace Acts
