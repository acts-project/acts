// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
// #include "Acts/Utilities/Delegate.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Surfaces/SurfaceConcept.hpp"
#include "Acts/Utilities/Concepts.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <concepts>
#include <memory>
#include <type_traits>

namespace Acts {

class PortalLinkBase;

// template <std::size_t N>
// class GridPortalLink;

class Portal {
 public:
  Portal(std::shared_ptr<RegularSurface> surface);

  static std::shared_ptr<Portal> fuse(const std::shared_ptr<Portal>& aPortal,
                                      const std::shared_ptr<Portal>& bPortal);

  static std::shared_ptr<Portal> mergeAdjacent(
      const std::shared_ptr<Portal>& aPortal,
      const std::shared_ptr<Portal>& bPortal);

  const Volume* resolveVolume(const GeometryContext& gctx,
                              const Vector3& position,
                              const Vector3& direction) const;

 private:
  // @TODO: Potentially short circuit the virtual call
  // using VolumeResolver = Delegate<TrackingVolume*(const Vector3& position)>;

  std::shared_ptr<RegularSurface> m_surface;

  std::unique_ptr<PortalLinkBase> m_alongNormal;
  std::unique_ptr<PortalLinkBase> m_oppositeNormal;
};

class GridPortalLink;
class GridPortalLink1;
class GridPortalLink2;

// template <std::size_t N>
// class GridPortalLinkN;

class PortalLinkBase {
 public:
  PortalLinkBase(const RegularSurface& surface) : m_surface(&surface) {}

  virtual ~PortalLinkBase() = default;
  // virtual const Volume* resolveVolume(const Vector3& position) const = 0;

  // virtual bool inside(const Vector2& position) const = 0;

  std::unique_ptr<PortalLinkBase> merge(
      const PortalLinkBase& other,
      const Logger& logger = getDummyLogger()) const;

  virtual void toStream(std::ostream& os) const = 0;

  friend std::ostream& operator<<(std::ostream& os, const PortalLinkBase& link);

  const RegularSurface& surface() const { return *m_surface; }

 protected:
  const RegularSurface* m_surface;
};

class BinaryPortalLink : public PortalLinkBase {};

namespace detail {
#if defined(__cpp_concepts)
template <typename S>
concept RegularSurfaceConcept =
    std::is_same_v<S, CylinderSurface> || std::is_same_v<S, DiscSurface> ||
    std::is_same_v<S, PlaneSurface>;

#endif
}  // namespace detail

template <ACTS_CONCEPT(detail::RegularSurfaceConcept) surface_t,
          typename... Axes>
class GridPortalLinkT;

class GridPortalLink : public PortalLinkBase {
 public:
  enum class Direction { loc0, loc1 };

  GridPortalLink(const RegularSurface& surface, Direction direction)
      : PortalLinkBase(surface), m_direction(direction) {}

  template <ACTS_CONCEPT(detail::RegularSurfaceConcept) surface_t,
            typename... Axes, typename = std::enable_if_t<sizeof...(Axes) == 1>>
  static auto make(const surface_t& surface, Direction direction,
                   Axes&&... axes) {
    return std::make_unique<GridPortalLinkT<surface_t, Axes...>>(
        surface, direction, std::forward<Axes>(axes)...);
  }

  template <ACTS_CONCEPT(detail::RegularSurfaceConcept) surface_t,
            typename... Axes, typename = std::enable_if_t<sizeof...(Axes) == 2>>
  static auto make(const surface_t& surface, Axes&&... axes) {
    return std::make_unique<GridPortalLinkT<surface_t, Axes...>>(
        surface, Direction::loc0, std::forward<Axes>(axes)...);
  }

  virtual const IGrid& grid() const = 0;

  Direction direction() const { return m_direction; }

 private:
  Direction m_direction;
};

template <ACTS_CONCEPT(detail::RegularSurfaceConcept) surface_t,
          typename... Axes>
class GridPortalLinkT : public GridPortalLink {
 public:
  using GridType = Grid<const Volume*, Axes...>;
  static constexpr std::size_t DIM = sizeof...(Axes);

  static_assert(DIM == 1 || DIM == 2,
                "GridPortalLinks only support 1D or 2D grids");

  GridPortalLinkT(const RegularSurface& surface, Direction direction,
                  Axes&&... axes)
      : GridPortalLink(surface, direction),
        m_grid(std::tuple{std::move(axes)...}) {
    checkConsistency(this->surface());
  }

  const GridType& grid() const override { return m_grid; }

  void toStream(std::ostream& os) const override { os << m_grid; }

 private:
  const surface_t& surface() const {
    return static_cast<const surface_t&>(*m_surface);
  }

  void checkConsistency(const CylinderSurface& cyl) const {
    auto checkZ = [&](const auto& axis) {
      ActsScalar hlZ = cyl.bounds().get(CylinderBounds::eHalfLengthZ);
      if (axis.getMin() != -hlZ || axis.getMax() != hlZ) {
        throw std::invalid_argument(
            "GridPortalLink: CylinderBounds: invalid length setup.");
      }
    };
    auto checkRPhi = [&](const auto& axis) {
      ActsScalar hlRPhi = cyl.bounds().get(CylinderBounds::eR) *
                          cyl.bounds().get(CylinderBounds::eHalfPhiSector);

      if (axis.getMin() != -hlRPhi || axis.getMax() != hlRPhi) {
        throw std::invalid_argument(
            "GridPortalLink: CylinderBounds: invalid phi sector setup.");
      }
    };

    if constexpr (DIM == 1) {
      auto [axisLoc0] = m_grid.axesTuple();
      if (direction() == Direction::loc0) {
        checkRPhi(axisLoc0);
      } else {
        checkZ(axisLoc0);
      }
    } else {  // DIM == 2
      auto [axisLoc0, axisLoc1] = m_grid.axesTuple();
      checkRPhi(axisLoc0);
      checkZ(axisLoc1);
    }
  }

  GridType m_grid;
};

std::ostream& operator<<(std::ostream& os, GridPortalLink::Direction direction);

}  // namespace Acts
