// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CompositePortalLink.hpp"

#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"

#include <stdexcept>

namespace Acts {

CompositePortalLink::CompositePortalLink(
    const std::shared_ptr<const PortalLinkBase>& a,
    const std::shared_ptr<const PortalLinkBase>& b, BinningValue direction,
    bool flatten)
    : PortalLinkBase(mergedSurface(*a, *b, direction)) {
  if (!flatten) {
    m_children = {a, b};
  } else {
    auto handle = [&](const std::shared_ptr<const PortalLinkBase>& link) {
      if (const auto* composite =
              dynamic_cast<const CompositePortalLink*>(link.get());
          composite != nullptr) {
        m_children.insert(m_children.end(), composite->m_children.begin(),
                          composite->m_children.end());
      } else {
        m_children.push_back(link);
      }
    };

    handle(a);
    handle(b);
  }
}

const TrackingVolume* CompositePortalLink::resolveVolume(
    const GeometryContext& gctx, const Vector2& position) const {
  Vector3 global = m_surface->localToGlobal(gctx, position);

  for (const auto& child : m_children) {
    if (child->surface().isOnSurface(gctx, global)) {
      return child->resolveVolume(gctx, position);
    }
  }

  throw std::runtime_error{"CompositePortalLink: Neither portal is on surface"};
}

std::shared_ptr<RegularSurface> CompositePortalLink::mergedSurface(
    const PortalLinkBase& a, const PortalLinkBase& b, BinningValue direction) {
  if (a.surface().type() != b.surface().type()) {
    throw std::invalid_argument{"Cannot merge surfaces of different types"};
  }

  if (const auto* cylinderA =
          dynamic_cast<const CylinderSurface*>(&a.surface());
      cylinderA != nullptr) {
    const auto& cylinderB = dynamic_cast<const CylinderSurface&>(b.surface());

    auto [merged, reversed] = cylinderA->mergedWith(cylinderB, direction, true);
    return merged;
  } else if (const auto* discA = dynamic_cast<const DiscSurface*>(&a.surface());
             discA != nullptr) {
    const auto& discB = dynamic_cast<const DiscSurface&>(b.surface());
    auto [merged, reversed] = discA->mergedWith(discB, direction, true);
    return merged;
  } else if (const auto* planeA =
                 dynamic_cast<const PlaneSurface*>(&a.surface());
             planeA != nullptr) {
    throw std::logic_error{"Plane surfaces not implemented yet"};
  } else {
    throw std::invalid_argument{"Unsupported surface type"};
  }
}

std::size_t CompositePortalLink::depth() const {
  std::size_t result = 1;

  for (const auto& child : m_children) {
    if (const auto* composite =
            dynamic_cast<const CompositePortalLink*>(child.get())) {
      result = std::max(result, composite->depth() + 1);
    }
  }

  return result;
}

std::size_t CompositePortalLink::size() const {
  return m_children.size();
}

}  // namespace Acts
