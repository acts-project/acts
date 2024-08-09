// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Portal.hpp"

#include <boost/container/small_vector.hpp>

namespace Acts {

class CompositePortalLink final : public PortalLinkBase {
 public:
  void toStream(std::ostream& os) const final { os << "CompositePortalLink"; }

  const TrackingVolume* resolveVolume(const GeometryContext& gctx,
                                      const Vector2& position) const final;

  CompositePortalLink(const std::shared_ptr<PortalLinkBase>& a,
                      const std::shared_ptr<PortalLinkBase>& b,
                      BinningValue direction, bool flatten = true);

  std::size_t depth() const;

  std::size_t size() const;

 private:
  static std::shared_ptr<RegularSurface> mergedSurface(const PortalLinkBase& a,
                                                       const PortalLinkBase& b,
                                                       BinningValue direction);

  boost::container::small_vector<std::shared_ptr<PortalLinkBase>, 4>
      m_children{};
};

}  // namespace Acts
