// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Portal.hpp"

#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/AxisFwd.hpp"

#include <type_traits>

namespace Acts {

Portal::Portal(std::shared_ptr<RegularSurface> surface)
    : m_surface(std::move(surface)) {
  throw_assert(m_surface, "Portal surface is nullptr");
}

const Volume* Portal::resolveVolume(const GeometryContext& gctx,
                                    const Vector3& position,
                                    const Vector3& direction) const {
  const Vector3 normal = m_surface->normal(gctx, position);
  Direction side = Direction::fromScalar(normal.dot(direction));

  const std::unique_ptr<PortalLinkBase>& link =
      side == Direction::AlongNormal ? m_alongNormal : m_oppositeNormal;

  return nullptr;
  if (link == nullptr) {
    // no link is attached in this direction => this is the end of the world as
    // we know it. (i feel fine)
    return nullptr;
  } else {
    // return link->resolveVolume(position);
  }
}

// MARK: - PortalLinkBase

std::ostream& operator<<(std::ostream& os, const PortalLinkBase& link) {
  link.toStream(os);
  return os;
}

// MARK: - GridPortalLink

// std::unique_ptr<PortalLinkBase> GridPortalLink::merge(
//     const PortalLinkBase& other) const {
//   std::cout << "Merge GridPortalLink + PortalLinkBase" << std::endl;
//   return other.merge(*this);
// }
//
// std::unique_ptr<PortalLinkBase> GridPortalLink::merge(
//     const GridPortalLink& other) const {
//   std::cout << "Merge GridPortalLink + GridPortalLink" << std::endl;
//   return nullptr;
// }

namespace {
std::unique_ptr<GridPortalLink1> merge(const GridPortalLink1& a,
                                       const GridPortalLink1& b) {}

std::unique_ptr<GridPortalLink2> merge(const GridPortalLink2& a,
                                       const GridPortalLink2& b) {}
}  // namespace

// MARK: - GridPortalLink1

std::unique_ptr<PortalLinkBase> GridPortalLink1::merge(
    const PortalLinkBase& other) const {
  std::cout << "Merge GridPortalLink1 + PortalLinkBase" << std::endl;
  return other.merge(*this);
}
std::unique_ptr<PortalLinkBase> GridPortalLink1::merge(
    const GridPortalLink1& other) const {
  std::cout << "Merge GridPortalLink1 + GridPortalLink1" << std::endl;
  return nullptr;
}

std::unique_ptr<PortalLinkBase> GridPortalLink1::merge(
    const GridPortalLink2& other) const {
  std::cout << "Merge GridPortalLink1 + GridPortalLink2" << std::endl;
  return nullptr;
}

// MARK: - GridPortalLink2

std::unique_ptr<PortalLinkBase> GridPortalLink2::merge(
    const PortalLinkBase& other) const {
  std::cout << "Merge GridPortalLink2 + PortalLinkBase" << std::endl;
  return other.merge(*this);
}

std::unique_ptr<PortalLinkBase> GridPortalLink2::merge(
    const GridPortalLink2& other) const {
  std::cout << "Merge GridPortalLink2 + GridPortalLink2" << std::endl;
  return nullptr;
}

std::unique_ptr<PortalLinkBase> GridPortalLink2::merge(
    const GridPortalLink1& other) const {
  std::cout << "Merge GridPortalLink2 + GridPortalLink1" << std::endl;
  return nullptr;
}

// template <typename axis_0_t, typename axis_1_t>
// class GridPortalLink2 : public GridPortalLink {
//   void toStream(std::ostream& os) const override {
//     os << "GridPortalLink2<" << axis_0_t::type << ", " << axis_1_t::type <<
//     ">";
//   }
// };

// std::unique_ptr<GridPortalLink> GridPortalLink::make(ActsScalar min,
//                                                      ActsScalar max,
//                                                      std::size_t bins) {
//   return std::make_unique<GridPortalLinkNT<detail::Axis<
//       detail::AxisType::Equidistant, detail::AxisBoundaryType::Closed>>>();
// }

}  // namespace Acts
