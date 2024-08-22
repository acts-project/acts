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

#include <iosfwd>

#include <boost/container/small_vector.hpp>

namespace Acts {

class CompositePortalLink final : public PortalLinkBase {
 public:
  CompositePortalLink(std::unique_ptr<PortalLinkBase> a,
                      std::unique_ptr<PortalLinkBase> b, BinningValue direction,
                      bool flatten = true);

  void toStream(std::ostream& os) const final;

  const TrackingVolume* resolveVolume(const GeometryContext& gctx,
                                      const Vector2& position) const final;

  const TrackingVolume* resolveVolume(const GeometryContext& gctx,
                                      const Vector3& position) const final;

  std::size_t depth() const;

  std::size_t size() const;

 private:
  static std::shared_ptr<RegularSurface> mergedSurface(const PortalLinkBase* a,
                                                       const PortalLinkBase* b,
                                                       BinningValue direction);

  boost::container::small_vector<std::unique_ptr<PortalLinkBase>, 4>
      m_children{};
};

}  // namespace Acts
