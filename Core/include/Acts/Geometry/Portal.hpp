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

class PortalLinkBase {
 public:
  virtual ~PortalLinkBase() = default;
  // virtual const Volume* resolveVolume(const Vector3& position) const = 0;

  // virtual bool inside(const Vector2& position) const = 0;

  virtual std::unique_ptr<PortalLinkBase> merge(
      const PortalLinkBase& other) const = 0;

  // virtual std::unique_ptr<PortalLinkBase> merge(
  //     const GridPortalLink& other) const = 0;
  //
  virtual std::unique_ptr<PortalLinkBase> merge(
      const GridPortalLink1& other) const = 0;

  virtual std::unique_ptr<PortalLinkBase> merge(
      const GridPortalLink2& other) const = 0;

  // virtual std::unique_ptr<PortalLinkBase> merge(
  //     const GridPortalLink<1>& other) const = 0;
  //
  // virtual std::unique_ptr<PortalLinkBase> merge(
  //     const GridPortalLink<2>& other) const = 0;

  virtual void toStream(std::ostream& os) const = 0;

  friend std::ostream& operator<<(std::ostream& os, const PortalLinkBase& link);
};

class BinaryPortalLink : public PortalLinkBase {};

class GridPortalLink : public PortalLinkBase {
 public:
  template <typename... Axes>
  static std::unique_ptr<GridPortalLink> make(Axes&&... axes);
  // static std::unique_ptr<GridPortalLink> make(ActsScalar min, ActsScalar max,
  //                                             std::size_t bins);
  //
  // static std::unique_ptr<GridPortalLink> make(std::vector<ActsScalar> bins);

  // static std::unique_ptr<GridPortalLink> make(ActsScalar min0, ActsScalar
  // max0,
  //                                             std::size_t bins0,
  //                                             ActsScalar min1, ActsScalar
  //                                             max1, std::size_t bins1);
  //
  // static std::unique_ptr<GridPortalLink> make(std::vector<ActsScalar> bins0,
  //                                             std::vector<ActsScalar> bins1);

  // std::unique_ptr<PortalLinkBase> merge(
  //     const PortalLinkBase& other) const override;

  // std::unique_ptr<PortalLinkBase> merge(
  //     const GridPortalLink& other) const override;
};

// template <std::size_t N>
// class GridPortalLinkN : public GridPortalLink {
//   std::unique_ptr<PortalLinkBase> merge(
//       const PortalLinkBase& other) const override {
//     std::cout << "Merge GridPortalLinkN<" << N << "> + PortalLinkBase"
//               << std::endl;
//     return other.merge(*this);
//   }
//
//   std::unique_ptr<PortalLinkBase> merge(
//       const GridPortalLinkN<N>& other) const override {
//     std::cout << "Merge GridPortalLinkN<" << N << "> + GridPortalLinkN<" << N
//               << ">" << std::endl;
//     return nullptr;
//   }
// };

class GridPortalLink1 : public GridPortalLink {
  std::unique_ptr<PortalLinkBase> merge(
      const PortalLinkBase& other) const override;

  std::unique_ptr<PortalLinkBase> merge(
      const GridPortalLink1& other) const override;

  std::unique_ptr<PortalLinkBase> merge(
      const GridPortalLink2& other) const override;
};

class GridPortalLink2 : public GridPortalLink {
  std::unique_ptr<PortalLinkBase> merge(
      const PortalLinkBase& other) const override;

  std::unique_ptr<PortalLinkBase> merge(
      const GridPortalLink1& other) const override;

  std::unique_ptr<PortalLinkBase> merge(
      const GridPortalLink2& other) const override;
};

template <typename... Axes>
class GridPortalLinkN
    : public std::conditional_t<sizeof...(Axes) == 1, GridPortalLink1,
                                GridPortalLink2> {
 public:
  static constexpr std::size_t dim = sizeof...(Axes);
  static_assert(dim == 1 || dim == 2, "Only 1D and 2D grids are supported");

  using GridType = Grid<int, Axes...>;

  GridPortalLinkN(Axes... axes) : m_grid(axes...) {}

  void toStream(std::ostream& os) const override {
    os << "GridPortalLinkN<" << dim << ">";
  }

  // std::unique_ptr<PortalLinkBase> merge(
  //     const PortalLinkBase& other) const override {
  //   std::cout << "Merge GridPortalLink1 + PortalLinkBase" << std::endl;
  //   return other.merge(*this);
  // }

 private:
  GridType m_grid;
};

template <typename... Axes>
std::unique_ptr<GridPortalLink> GridPortalLink::make(Axes&&... axes) {
  return std::make_unique<GridPortalLinkN<Axes...>>(
      std::forward<Axes>(axes)...);
}

}  // namespace Acts
