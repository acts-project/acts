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

enum class PortalDirection { loc0, loc1 };

std::ostream& operator<<(std::ostream& os, PortalDirection direction);

class PortalLinkBase {
 public:
  virtual ~PortalLinkBase() = default;
  // virtual const Volume* resolveVolume(const Vector3& position) const = 0;

  // virtual bool inside(const Vector2& position) const = 0;

  // std::unique_ptr<PortalLinkBase> merge(
  //     const PortalLinkBase& other, const Vector2& offset,
  //     const Logger& logger = getDummyLogger()) const;

  virtual std::unique_ptr<PortalLinkBase> merge(
      const PortalLinkBase& other, const Vector2& offset,
      const Logger& logger = getDummyLogger()) const = 0;

  virtual std::unique_ptr<PortalLinkBase> merge(
      const GridPortalLink1& other, const Vector2& offset,
      const Logger& logger = getDummyLogger()) const = 0;

  virtual std::unique_ptr<PortalLinkBase> merge(
      const GridPortalLink2& other, const Vector2& offset,
      const Logger& logger = getDummyLogger()) const = 0;

  static std::unique_ptr<GridPortalLink1> merge1d(const GridPortalLink1& a,
                                                  const GridPortalLink1& b,
                                                  const Vector2& offset,
                                                  const Logger& logger);

  virtual void toStream(std::ostream& os) const = 0;

  friend std::ostream& operator<<(std::ostream& os, const PortalLinkBase& link);
};

class BinaryPortalLink : public PortalLinkBase {};

template <typename... Axes>
class GridPortalLinkN;

class GridPortalLink : public PortalLinkBase {
 public:
  template <typename... Axes>
  static std::unique_ptr<GridPortalLinkN<Axes...>> make(Axes&&... axes);

  virtual const IGrid& grid() const = 0;
};

class GridPortalLink1 : public GridPortalLink {
 public:
  std::unique_ptr<PortalLinkBase> merge(
      const PortalLinkBase& other, const Vector2& offset,
      const Logger& logger = getDummyLogger()) const override;

  std::unique_ptr<PortalLinkBase> merge(
      const GridPortalLink1& other, const Vector2& offset,
      const Logger& logger = getDummyLogger()) const override;

  std::unique_ptr<PortalLinkBase> merge(
      const GridPortalLink2& other, const Vector2& offset,
      const Logger& logger = getDummyLogger()) const override;

  PortalDirection direction() const { return m_direction; }

 protected:
  GridPortalLink1(PortalDirection direction) : m_direction(direction) {}

 private:
  /// Local direction of the 1D grid:
  /// - 0 means the axis is along
  PortalDirection m_direction;
};

class GridPortalLink2 : public GridPortalLink {
 public:
  std::unique_ptr<PortalLinkBase> merge(
      const PortalLinkBase& other, const Vector2& offset,
      const Logger& logger = getDummyLogger()) const override;

  std::unique_ptr<PortalLinkBase> merge(
      const GridPortalLink1& other, const Vector2& offset,
      const Logger& logger = getDummyLogger()) const override;

  std::unique_ptr<PortalLinkBase> merge(
      const GridPortalLink2& other, const Vector2& offset,
      const Logger& logger = getDummyLogger()) const override;
};

template <typename... Axes>
class GridPortalLinkN
    : public std::conditional_t<sizeof...(Axes) == 1, GridPortalLink1,
                                GridPortalLink2> {
 public:
  static constexpr std::size_t dim = sizeof...(Axes);
  static_assert(dim == 1 || dim == 2, "Only 1D and 2D grids are supported");

  using GridType = Grid<int, Axes...>;

  template <std::size_t D = dim, typename = std::enable_if_t<D == 1>>
  GridPortalLinkN(GridType grid,
                  PortalDirection direction = PortalDirection::loc0)
      : GridPortalLink1(direction), m_grid(std::move(grid)) {}

  template <std::size_t D = dim, typename = std::enable_if_t<D == 2>>
  GridPortalLinkN(GridType grid) : m_grid(std::move(grid)) {}

  void toStream(std::ostream& os) const override {
    os << "GridPortalLinkN<" << dim << ">";
  }

  const IGrid& grid() const override { return m_grid; }

 private:
  GridType m_grid;
};

template <typename... Axes>
std::unique_ptr<GridPortalLinkN<Axes...>> GridPortalLink::make(Axes&&... axes) {
  using portal_link_t = GridPortalLinkN<Axes...>;
  using grid_t = typename portal_link_t::GridType;
  return std::make_unique<portal_link_t>(
      grid_t{std::tuple{std::forward<Axes>(axes)...}});
}

}  // namespace Acts
