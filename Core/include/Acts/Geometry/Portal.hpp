// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Volume.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Grid.hpp"
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

template <typename S>
concept PortalSurfaceConcept =
    std::is_same_v<S, CylinderSurface> || std::is_same_v<S, DiscSurface> ||
    std::is_same_v<S, PlaneSurface>;

class PortalLinkBase {
 public:
  PortalLinkBase(std::shared_ptr<RegularSurface> surface)
      : m_surface(std::move(surface)) {}

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

  std::shared_ptr<RegularSurface> m_surface;
};

class BinaryPortalLink : public PortalLinkBase {};

class TrivialPortalLink : public PortalLinkBase {};

template <typename... Axes>
class GridPortalLinkT;

class GridPortalLink : public PortalLinkBase {
 protected:
  GridPortalLink(std::shared_ptr<RegularSurface> surface,
                 BinningValue direction)
      : PortalLinkBase(std::move(surface)), m_direction(direction) {}

 public:
  template <PortalSurfaceConcept surface_t, typename axis_t>
  static std::unique_ptr<GridPortalLink> make(
      std::shared_ptr<surface_t> surface, BinningValue direction,
      axis_t&& axis) {
    if constexpr (std::is_same_v<surface_t, CylinderSurface>) {
      if (direction != BinningValue::binZ &&
          direction != BinningValue::binRPhi) {
        throw std::invalid_argument{"Invalid binning direction"};
      }
    } else if constexpr (std::is_same_v<surface_t, DiscSurface>) {
      if (direction != BinningValue::binR &&
          direction != BinningValue::binPhi) {
        throw std::invalid_argument{"Invalid binning direction"};
      }
    }

    return std::make_unique<GridPortalLinkT<axis_t>>(surface, direction,
                                                     std::move(axis));
  }

  template <PortalSurfaceConcept surface_t, typename axis_1_t,
            typename axis_2_t>
  static std::unique_ptr<GridPortalLink> make(
      std::shared_ptr<surface_t> surface, axis_1_t axis1, axis_2_t axis2) {
    std::optional<BinningValue> direction;
    if constexpr (std::is_same_v<surface_t, CylinderSurface>) {
      direction = BinningValue::binRPhi;
    } else if constexpr (std::is_same_v<surface_t, DiscSurface>) {
      direction = BinningValue::binR;
    }

    return std::make_unique<GridPortalLinkT<axis_1_t, axis_2_t>>(
        surface, direction.value(), std::move(axis1), std::move(axis2));
  }

  std::unique_ptr<PortalLinkBase> mergeImpl(
      const PortalLinkBase& other, const RegularSurface& surfaceA,
      const RegularSurface& surfaceB, BinningValue direction,
      const Logger& logger = getDummyLogger()) const override;

  virtual const IGrid& grid() const = 0;
  virtual unsigned int dim() const = 0;
  virtual std::unique_ptr<GridPortalLink> make2DGrid() const = 0;

  BinningValue direction() const { return m_direction; }

 protected:
  void checkConsistency(const CylinderSurface& disc) const;
  void checkConsistency(const DiscSurface& disc) const;

  std::unique_ptr<GridPortalLink> extendTo2D(
      const std::shared_ptr<CylinderSurface>& surface) const;
  std::unique_ptr<GridPortalLink> extendTo2D(
      const std::shared_ptr<DiscSurface>& surface) const;

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

  GridPortalLinkT(std::shared_ptr<RegularSurface> surface,
                  BinningValue direction, Axes&&... axes)
      : GridPortalLink(std::move(surface), direction),
        m_grid(std::tuple{std::move(axes)...}) {
    if (const auto* cylinder =
            dynamic_cast<const CylinderSurface*>(m_surface.get())) {
      checkConsistency(*cylinder);
    } else if (const auto* disc =
                   dynamic_cast<const DiscSurface*>(m_surface.get())) {
      checkConsistency(*disc);
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
      if (auto cylinder =
              std::dynamic_pointer_cast<CylinderSurface>(m_surface)) {
        return extendTo2D(cylinder);
      } else if (auto disc =
                     std::dynamic_pointer_cast<DiscSurface>(m_surface)) {
        return extendTo2D(disc);
      } else {
        throw std::logic_error{
            "Surface type is not supported (this should not happen)"};
      }
    }
  }

  GridType m_grid;
};

}  // namespace Acts
