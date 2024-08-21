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
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

class IGrid;

template <typename... Axes>
  requires(sizeof...(Axes) <= 2)
class GridPortalLinkT;

class GridPortalLink : public PortalLinkBase {
 protected:
  GridPortalLink(std::shared_ptr<RegularSurface> surface,
                 BinningValue direction)
      : PortalLinkBase(std::move(surface)), m_direction(direction) {}

 public:
  template <AxisConcept axis_t>
  static std::unique_ptr<GridPortalLinkT<axis_t>> make(
      std::shared_ptr<RegularSurface> surface, BinningValue direction,
      axis_t&& axis) {
    if (dynamic_cast<const CylinderSurface*>(surface.get()) != nullptr) {
      if (direction != BinningValue::binZ &&
          direction != BinningValue::binRPhi) {
        throw std::invalid_argument{"Invalid binning direction"};
      }
    } else if (dynamic_cast<const DiscSurface*>(surface.get()) != nullptr) {
      if (direction != BinningValue::binR &&
          direction != BinningValue::binPhi) {
        throw std::invalid_argument{"Invalid binning direction"};
      }
    }

    return std::make_unique<GridPortalLinkT<axis_t>>(surface, direction,
                                                     std::move(axis));
  }

  template <AxisConcept axis_1_t, AxisConcept axis_2_t>
  static std::unique_ptr<GridPortalLinkT<axis_1_t, axis_2_t>> make(
      std::shared_ptr<RegularSurface> surface, axis_1_t axis1, axis_2_t axis2) {
    std::optional<BinningValue> direction;
    if (dynamic_cast<const CylinderSurface*>(surface.get()) != nullptr) {
      direction = BinningValue::binRPhi;
    } else if (dynamic_cast<const DiscSurface*>(surface.get()) != nullptr) {
      direction = BinningValue::binR;
    }

    return std::make_unique<GridPortalLinkT<axis_1_t, axis_2_t>>(
        surface, direction.value(), std::move(axis1), std::move(axis2));
  }

  static std::unique_ptr<GridPortalLink> make(
      const std::shared_ptr<RegularSurface>& surface, TrackingVolume& volume,
      BinningValue direction);

  /// +----------------------------------------------------------+
  /// |Colinear                                                  |
  /// |                                                          |
  /// | +-------+-------+-------+     +-------+-------+-------+  |
  /// | |       |       |       |     |       |       |       |  |
  /// | |       |       |       |  +  |       |       |       |  |
  /// | |       |       |       |     |       |       |       |  |
  /// | +-------+-------+-------+     +-------+-------+-------+  |
  /// |       +-------+-------+-------+-------+-------+-------+  |
  /// |       |       |       |       |       |       |       |  |
  /// |    =  |       |       |       |       |       |       |  |
  /// |       |       |       |       |       |       |       |  |
  /// |       +-------+-------+-------+-------+-------+-------+  |
  /// |                                                          |
  /// +----------------------------------------------------------+
  ///
  /// +-------------------------------------------------+
  /// |Parallel                                         |
  /// |                                                 |
  /// | +-------+      +-------+     +-------+-------+  |
  /// | |       |      |       |     |       |       |  |
  /// | |       |      |       |     |       |       |  |
  /// | |       |      |       |     |       |       |  |
  /// | +-------+      +-------+     +-------+-------+  |
  /// | |       |      |       |     |       |       |  |
  /// | |       |  +   |       |  =  |       |       |  |
  /// | |       |      |       |     |       |       |  |
  /// | +-------+      +-------+     +-------+-------+  |
  /// | |       |      |       |     |       |       |  |
  /// | |       |      |       |     |       |       |  |
  /// | |       |      |       |     |       |       |  |
  /// | +-------+      +-------+     +-------+-------+  |
  /// |                                                 |
  /// +-------------------------------------------------+
  ///
  /// +-------------------------------------------+
  /// |Perpendicular                              |
  /// |                                           |
  /// |  +-------+                                |
  /// |  |       |                                |
  /// |  |       |                                |
  /// |  |       |                                |
  /// |  +-------+     +-------+-------+-------+  |
  /// |  |       |     |       |       |       |  |
  /// |  |       |  +  |       |       |       |  |
  /// |  |       |     |       |       |       |  |
  /// |  +-------+     +-------+-------+-------+  |
  /// |  |       |                                |
  /// |  |       |                                |
  /// |  |       |     +-------+-------+-------+  |
  /// |  +-------+     |       |       |       |  |
  /// |                |       |       |       |  |
  /// |                |       |       |       |  |
  /// |                +-------+-------+-------+  |
  /// |                |       |       |       |  |
  /// |             =  |       |       |       |  |
  /// |                |       |       |       |  |
  /// |                +-------+-------+-------+  |
  /// |                |       |       |       |  |
  /// |                |       |       |       |  |
  /// |                |       |       |       |  |
  /// |                +-------+-------+-------+  |
  /// |                                           |
  /// +-------------------------------------------+
  static std::unique_ptr<PortalLinkBase> merge(
      const GridPortalLink& a, const GridPortalLink& b, BinningValue direction,
      const Logger& logger = getDummyLogger());

  virtual const IGrid& grid() const = 0;
  virtual void setVolume(TrackingVolume* volume) = 0;
  virtual unsigned int dim() const = 0;
  virtual std::unique_ptr<GridPortalLink> make2DGrid(
      const IAxis* other) const = 0;

  BinningValue direction() const { return m_direction; }

  /// This is primarily for testing / inspection
  virtual void visitBins(
      const std::function<void(const TrackingVolume*)> func) const = 0;

  static void fillMergedGrid(const GridPortalLink& a, const GridPortalLink& b,
                             GridPortalLink& merged, BinningValue direction,
                             const Logger& logger);

  void printContents(std::ostream& os) const;

 protected:
  void checkConsistency(const CylinderSurface& cyl) const;
  void checkConsistency(const DiscSurface& disc) const;

  std::unique_ptr<GridPortalLink> extendTo2D(
      const std::shared_ptr<CylinderSurface>& surface,
      const IAxis* other) const;
  std::unique_ptr<GridPortalLink> extendTo2D(
      const std::shared_ptr<DiscSurface>& surface, const IAxis* other) const;

  enum class FillDirection {
    loc0,
    loc1,
  };
  static void fillGrid1dTo2d(FillDirection dir, const GridPortalLink& grid1d,
                             GridPortalLink& grid2d);

  using IndexType = boost::container::static_vector<std::size_t, 2>;

  // These should only be used to synchronize the bins between grids
  virtual IndexType numLocalBins() const = 0;
  virtual TrackingVolume*& atLocalBins(IndexType indices) = 0;

  virtual TrackingVolume* atLocalBins(IndexType indices) const = 0;

 private:
  BinningValue m_direction;
};

class GridPortalLink1 : public GridPortalLink {};
class GridPortalLink2 : public GridPortalLink {};

template <typename... Axes>
  requires(sizeof...(Axes) <= 2)
class GridPortalLinkT final : public GridPortalLink {
 public:
  using GridType = Grid<TrackingVolume*, Axes...>;
  static constexpr std::size_t DIM = sizeof...(Axes);

  GridPortalLinkT(std::shared_ptr<RegularSurface> surface,
                  BinningValue direction, Axes&&... axes)
      : GridPortalLink(std::move(surface), direction),
        m_grid(std::tuple{std::move(axes)...}) {
    if (const auto* cylinder =
            dynamic_cast<const CylinderSurface*>(m_surface.get())) {
      checkConsistency(*cylinder);

      if (direction == BinningValue::binRPhi) {
        m_projection = &projection<CylinderSurface, BinningValue::binRPhi>;
      } else if (direction == BinningValue::binZ) {
        m_projection = &projection<CylinderSurface, BinningValue::binZ>;
      } else {
        throw std::invalid_argument{"Invalid binning direction"};
      }

    } else if (const auto* disc =
                   dynamic_cast<const DiscSurface*>(m_surface.get())) {
      checkConsistency(*disc);

      if (direction == BinningValue::binR) {
        m_projection = &projection<DiscSurface, BinningValue::binR>;
      } else if (direction == BinningValue::binPhi) {
        m_projection = &projection<DiscSurface, BinningValue::binPhi>;
      } else {
        throw std::invalid_argument{"Invalid binning direction"};
      }

    } else {
      throw std::logic_error{"Surface type is not supported"};
    }
  }

  const GridType& grid() const final { return m_grid; }
  GridType& grid() { return m_grid; }
  unsigned int dim() const final { return DIM; }

  void toStream(std::ostream& os) const final {
    os << "GridPortalLink<dim=" << dim() << ">";
  }

  std::unique_ptr<GridPortalLink> make2DGrid(const IAxis* other) const final {
    if constexpr (DIM == 2) {
      return std::make_unique<GridPortalLinkT<Axes...>>(*this);
    } else {
      if (auto cylinder =
              std::dynamic_pointer_cast<CylinderSurface>(m_surface)) {
        return extendTo2D(cylinder, other);
      } else if (auto disc =
                     std::dynamic_pointer_cast<DiscSurface>(m_surface)) {
        return extendTo2D(disc, other);
      } else {
        throw std::logic_error{
            "Surface type is not supported (this should not happen)"};
      }
    }
  }

  void setVolume(TrackingVolume* volume) final {
    auto loc = m_grid.numLocalBins();
    if constexpr (GridType::DIM == 1) {
      for (std::size_t i = 1; i <= loc[0]; i++) {
        m_grid.atLocalBins({i}) = volume;
      }
    } else {
      for (std::size_t i = 1; i <= loc[0]; i++) {
        for (std::size_t j = 1; j <= loc[1]; j++) {
          m_grid.atLocalBins({i, j}) = volume;
        }
      }
    }
  }

  const TrackingVolume* resolveVolume(const GeometryContext& gctx,
                                      const Vector3& position) const final {
    auto res = m_surface->globalToLocal(gctx, position);
    if (!res.ok()) {
      throw std::invalid_argument{"Cannot resolve position"};
    }

    const Vector2& local = *res;
    return resolveVolume(gctx, local);
  }

  const TrackingVolume* resolveVolume(const GeometryContext& /*gctx*/,
                                      const Vector2& position) const final {
    // @TODO: Should this be an exception or nullptr? Or not even checked at all?
    if (!surface().insideBounds(position, BoundaryTolerance::None())) {
      throw std::invalid_argument{"Position is outside surface bounds"};
    }
    return m_grid.atPosition(m_projection(position));
  }

  void visitBins(
      const std::function<void(const TrackingVolume*)> func) const final {
    auto loc = m_grid.numLocalBins();
    if constexpr (GridType::DIM == 1) {
      for (std::size_t i = 1; i <= loc[0]; i++) {
        func(m_grid.atLocalBins({i}));
      }
    } else {
      for (std::size_t i = 1; i <= loc[0]; i++) {
        for (std::size_t j = 1; j <= loc[1]; j++) {
          func(m_grid.atLocalBins({i, j}));
        }
      }
    }
  }

 protected:
  IndexType numLocalBins() const final {
    typename GridType::index_t idx = m_grid.numLocalBins();
    IndexType result;
    for (std::size_t i = 0; i < DIM; i++) {
      result.push_back(idx[i]);
    }
    return result;
  }

  TrackingVolume*& atLocalBins(IndexType indices) final {
    throw_assert(indices.size() == DIM, "Invalid number of indices");
    typename GridType::index_t idx;
    for (std::size_t i = 0; i < DIM; i++) {
      idx[i] = indices[i];
    }
    return m_grid.atLocalBins(idx);
  }

  TrackingVolume* atLocalBins(IndexType indices) const final {
    throw_assert(indices.size() == DIM, "Invalid number of indices");
    typename GridType::index_t idx;
    for (std::size_t i = 0; i < DIM; i++) {
      idx[i] = indices[i];
    }
    return m_grid.atLocalBins(idx);
  }

 private:
  template <class surface_t, BinningValue direction>
  static ActsVector<DIM> projection(const Vector2& position) {
    if constexpr (DIM == 2) {
      return position;
    } else {
      if constexpr (std::is_same_v<surface_t, CylinderSurface>) {
        if constexpr (direction == BinningValue::binRPhi) {
          return ActsVector<1>{position[0]};
        } else if constexpr (direction == BinningValue::binZ) {
          return ActsVector<1>{position[1]};
        } else {
          []<bool flag = false>() {
            static_assert(flag, "invalid direction");
          }();
        }
      } else if constexpr (std::is_same_v<surface_t, DiscSurface>) {
        if constexpr (direction == BinningValue::binR) {
          return ActsVector<1>{position[0]};
        } else if constexpr (direction == BinningValue::binPhi) {
          return ActsVector<1>{position[1]};
        } else {
          []<bool flag = false>() {
            static_assert(flag, "invalid direction");
          }();
        }
      } else if constexpr (std::is_same_v<surface_t, PlaneSurface>) {
        if constexpr (direction == BinningValue::binX) {
          return ActsVector<1>{position[0]};
        } else if constexpr (direction == BinningValue::binY) {
          return ActsVector<1>{position[1]};
        } else {
          []<bool flag = false>() {
            static_assert(flag, "invalid direction");
          }();
        }
      }
    }
  }

  GridType m_grid;

  ActsVector<DIM> (*m_projection)(const Vector2& position);
};

}  // namespace Acts
