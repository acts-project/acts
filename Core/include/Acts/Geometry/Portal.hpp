// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/RegularSurface.hpp"
// #include "Acts/Utilities/Delegate.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Utilities/Logger.hpp"

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

  // std::unique_ptr<PortalLinkBase> merge(
  //     const PortalLinkBase& other, const Vector2& offset,
  //     const Logger& logger = getDummyLogger()) const;

  virtual void toStream(std::ostream& os) const = 0;

  friend std::ostream& operator<<(std::ostream& os, const PortalLinkBase& link);

  const RegularSurface& surface() const { return *m_surface; }

 protected:
  const RegularSurface* m_surface;
};

class BinaryPortalLink : public PortalLinkBase {};

template <typename... Axes>
class GridPortalLinkT;

class GridPortalLink : public PortalLinkBase {
 public:
  enum class Direction { loc0, loc1 };

  GridPortalLink(const RegularSurface& surface, Direction direction)
      : PortalLinkBase(surface), m_direction(direction) {}

  template <typename... Axes, typename = std::enable_if_t<sizeof...(Axes) == 1>>
  static std::unique_ptr<GridPortalLinkT<Axes...>> make(
      const RegularSurface& surface, Direction direction, Axes&&... axes) {
    return std::make_unique<GridPortalLinkT<Axes...>>(
        surface, direction, std::forward<Axes>(axes)...);
  }

  template <typename... Axes, typename = std::enable_if_t<sizeof...(Axes) == 2>>
  static std::unique_ptr<GridPortalLinkT<Axes...>> make(
      const RegularSurface& surface, Axes&&... axes) {
    return std::make_unique<GridPortalLinkT<Axes...>>(
        surface, Direction::loc0, std::forward<Axes>(axes)...);
  }

  virtual const IGrid& grid() const = 0;

  Direction direction() const { return m_direction; }

 private:
  Direction m_direction;
};

template <typename... Axes>
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
    checkConsistency();
  }

  const GridType& grid() const override { return m_grid; }

  void toStream(std::ostream& os) const override { os << m_grid; }

 private:
  void checkConsistency() const {
    // if (DIM == 1) {
    // } else {
    // }
  }

  GridType m_grid;
};

std::ostream& operator<<(std::ostream& os, GridPortalLink::Direction direction);

}  // namespace Acts
