// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
// #include "Acts/Utilities/Delegate.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Surfaces/SurfaceConcept.hpp"
#include "Acts/Utilities/BinningType.hpp"
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

// template <std::size_t N>
// class GridPortalLinkN;

class PortalLinkBase {
 public:
  PortalLinkBase(const RegularSurface& surface) : m_surface(&surface) {}

  virtual ~PortalLinkBase() = default;
  // virtual const Volume* resolveVolume(const Vector3& position) const = 0;

  // virtual bool inside(const Vector2& position) const = 0;

  std::unique_ptr<PortalLinkBase> merge(
      const GeometryContext& gctx, const PortalLinkBase& other,
      BinningValue direction, const Logger& logger = getDummyLogger()) const;

  virtual void toStream(std::ostream& os) const = 0;

  friend std::ostream& operator<<(std::ostream& os, const PortalLinkBase& link);

  const RegularSurface& surface() const { return *m_surface; }

 protected:
  virtual std::unique_ptr<PortalLinkBase> mergeImpl(
      const PortalLinkBase& other, const RegularSurface& surfaceA,
      const RegularSurface& surfaceB, BinningValue direction,
      const Logger& logger = getDummyLogger()) const;

  const RegularSurface* m_surface;
};

class BinaryPortalLink : public PortalLinkBase {};

class TrivialPortalLink : public PortalLinkBase {};

namespace detail {
#if defined(__cpp_concepts)
template <typename S>
concept RegularSurfaceConcept =
    std::is_same_v<S, CylinderSurface> || std::is_same_v<S, DiscSurface> ||
    std::is_same_v<S, PlaneSurface>;

#endif
}  // namespace detail

template <typename... Axes>
class GridPortalLinkT;

class GridPortalLink : public PortalLinkBase {
 protected:
  GridPortalLink(const RegularSurface& surface, BinningValue direction)
      : PortalLinkBase(surface), m_direction(direction) {}

 public:
  template <typename... Axes, typename = std::enable_if_t<sizeof...(Axes) == 1>>
  static auto make(const CylinderSurface& surface, BinningValue direction,
                   Axes&&... axes) {
    if (direction != BinningValue::binZ && direction != BinningValue::binRPhi) {
      throw std::invalid_argument{"Invalid binning direction"};
    }

    return std::make_unique<GridPortalLinkT<Axes...>>(
        surface, direction, std::forward<Axes>(axes)...);
  }

  template <typename... Axes, typename = std::enable_if_t<sizeof...(Axes) == 1>>
  static auto make(const DiscSurface& surface, BinningValue direction,
                   Axes&&... axes) {
    if (direction != BinningValue::binR && direction != BinningValue::binPhi) {
      throw std::invalid_argument{"Invalid binning direction"};
    }

    return std::make_unique<GridPortalLinkT<Axes...>>(
        surface, direction, std::forward<Axes>(axes)...);
  }

  template <typename... Axes, typename = std::enable_if_t<sizeof...(Axes) == 2>>
  static auto make(const CylinderSurface& surface, Axes&&... axes) {
    return std::make_unique<GridPortalLinkT<Axes...>>(
        surface, BinningValue::binRPhi, std::forward<Axes>(axes)...);
  }

  template <typename... Axes, typename = std::enable_if_t<sizeof...(Axes) == 2>>
  static auto make(const DiscSurface& surface, Axes&&... axes) {
    return std::make_unique<GridPortalLinkT<Axes...>>(
        surface, BinningValue::binR, std::forward<Axes>(axes)...);
  }

  std::unique_ptr<PortalLinkBase> mergeImpl(
      const PortalLinkBase& other, const RegularSurface& surfaceA,
      const RegularSurface& surfaceB, BinningValue direction,
      const Logger& logger = getDummyLogger()) const override;

  virtual const IGrid& grid() const = 0;
  virtual unsigned int dim() const = 0;
  virtual std::unique_ptr<GridPortalLink> make2DGrid() const = 0;

  BinningValue direction() const { return m_direction; }

 private:
  BinningValue m_direction;
};

template <typename... Axes>
class GridPortalLinkT : public GridPortalLink {
 public:
  using GridType = Grid<const Volume*, Axes...>;
  static constexpr std::size_t DIM = sizeof...(Axes);

  static_assert(DIM == 1 || DIM == 2,
                "GridPortalLinks only support 1D or 2D grids");

  GridPortalLinkT(const RegularSurface& surface, BinningValue direction,
                  Axes&&... axes)
      : GridPortalLink(surface, direction),
        m_grid(std::tuple{std::move(axes)...}) {
    if (const auto* cylinder = dynamic_cast<const CylinderSurface*>(&surface)) {
      checkConsistency(*cylinder);
    } else {
      throw std::logic_error{"Surface type is not supported"};
    }
  }

  const GridType& grid() const override { return m_grid; }
  unsigned int dim() const override { return DIM; }

  void toStream(std::ostream& os) const override {
    os << "GridPortalLink<dim=" << dim() << ">";
  }

  std::unique_ptr<GridPortalLink> make2DGrid() const override {
    if constexpr (DIM == 2) {
      return std::make_unique<GridPortalLinkT<Axes...>>(*this);
    } else {
      if (const auto* cylinder =
              dynamic_cast<const CylinderSurface*>(m_surface)) {
        return extendTo2D(*cylinder);
      } else {
        throw std::logic_error{
            "Surface type is not supported (this should not happen)"};
      }
    }
  }

 private:
  void checkConsistency(const CylinderSurface& cyl) const {
    constexpr auto tolerance = s_onSurfaceTolerance;
    auto same = [](auto a, auto b) { return std::abs(a - b) < tolerance; };

    auto checkZ = [&](const auto& axis) {
      ActsScalar hlZ = cyl.bounds().get(CylinderBounds::eHalfLengthZ);
      if (!same(axis.getMin(), -hlZ) || !same(axis.getMax(), hlZ)) {
        throw std::invalid_argument(
            "GridPortalLink: CylinderBounds: invalid length setup.");
      }
    };
    auto checkRPhi = [&](const auto& axis) {
      if (cyl.bounds().get(CylinderBounds::eAveragePhi) != 0) {
        throw std::invalid_argument(
            "GridPortalLink: CylinderBounds: only average phi == 0 is "
            "supported. Rotate the cylinder surface.");
      };

      ActsScalar hlPhi = cyl.bounds().get(CylinderBounds::eHalfPhiSector);
      ActsScalar r = cyl.bounds().get(CylinderBounds::eR);
      ActsScalar hlRPhi = r * hlPhi;

      if (!same(axis.getMin(), -hlRPhi) || !same(axis.getMax(), hlRPhi)) {
        throw std::invalid_argument(
            "GridPortalLink: CylinderBounds: invalid phi sector setup: axes "
            "don't match bounds");
      }

      // If full cylinder, make sure axis wraps around
      if (same(hlPhi, M_PI) &&
          axis.getBoundaryType() != Acts::AxisBoundaryType::Closed) {
        throw std::invalid_argument(
            "GridPortalLink: CylinderBounds: invalid phi sector setup: axis is "
            "not closed.");
      }
    };

    if constexpr (DIM == 1) {
      auto [axisLoc0] = m_grid.axesTuple();
      if (direction() == BinningValue::binRPhi) {
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

  std::unique_ptr<GridPortalLink> extendTo2D(
      const CylinderSurface& surface) const {
    if (direction() == BinningValue::binRPhi) {
      auto [axisRPhi] = m_grid.axesTuple();
      // 1D direction is binRPhi, so add a Z axis
      ActsScalar hlZ = surface.bounds().get(CylinderBounds::eHalfLengthZ);

      Axis axisZ{AxisBound, -hlZ, hlZ, 1};
      return GridPortalLink::make(surface, std::move(axisRPhi),
                                  std::move(axisZ));
    } else {
      auto [axisZ] = m_grid.axesTuple();
      // 1D direction is binZ, so add an rPhi axis
      ActsScalar r = surface.bounds().get(CylinderBounds::eR);
      ActsScalar hlPhi = surface.bounds().get(CylinderBounds::eHalfPhiSector);
      ActsScalar hlRPhi = r * hlPhi;

      auto axis = [&](auto bdt) {
        return GridPortalLink::make(surface, Axis{bdt, -hlRPhi, hlRPhi, 1},
                                    std::move(axisZ));
      };

      if (surface.bounds().coversFullAzimuth()) {
        return axis(AxisClosed);
      } else {
        return axis(AxisBound);
      }
    }
  }

  GridType m_grid;
};

}  // namespace Acts
