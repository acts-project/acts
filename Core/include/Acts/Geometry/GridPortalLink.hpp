// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Utilities/Grid.hpp"

#include <format>
#include <vector>

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
      const std::shared_ptr<RegularSurface>& surface,
      const TrackingVolume& volume, BinningValue direction);

  virtual const IGrid& grid() const = 0;
  virtual void setVolume(const TrackingVolume& volume) = 0;
  virtual unsigned int dim() const = 0;
  virtual std::unique_ptr<GridPortalLink> make2DGrid(
      const IAxis* other) const = 0;

  BinningValue direction() const { return m_direction; }

  /// This is primarily for testing / inspection
  virtual void visitBins(
      const std::function<void(const TrackingVolume*)> func) const = 0;

  std::unique_ptr<PortalLinkBase> mergeImpl(
      const PortalLinkBase& other, const RegularSurface& surfaceA,
      const RegularSurface& surfaceB, BinningValue direction,
      const Logger& logger = getDummyLogger()) const override;

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

  // These should only be used to synchronize the bins between grids
  virtual std::vector<std::size_t> numLocalBins() const = 0;
  virtual const TrackingVolume*& atLocalBins(
      const std::vector<std::size_t>) = 0;

  virtual const TrackingVolume* atLocalBins(
      const std::vector<std::size_t>) const = 0;

 private:
  BinningValue m_direction;
};

class GridPortalLink1 : public GridPortalLink {};
class GridPortalLink2 : public GridPortalLink {};

template <typename... Axes>
  requires(sizeof...(Axes) <= 2)
class GridPortalLinkT final : public GridPortalLink {
 public:
  using GridType = Grid<const TrackingVolume*, Axes...>;
  static constexpr std::size_t DIM = sizeof...(Axes);

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

  void setVolume(const TrackingVolume& volume) final {
    for (std::size_t i = 0; i < m_grid.size(); i++) {
      m_grid.at(i) = &volume;
    }
  }

  const TrackingVolume* resolveVolume(const GeometryContext& /*gctx*/,
                                      const Vector2& position) const final {
    if (!surface().insideBounds(position, BoundaryTolerance::None())) {
      // @FIXME: Should this throw an exception?
      throw std::invalid_argument{"Position is outside of the surface bounds"};
    }
    return m_grid.atPosition(position);
  }

  void visitBins(
      const std::function<void(const TrackingVolume*)> func) const final {
    for (std::size_t i = 0; i < m_grid.size(); i++) {
      func(m_grid.at(i));
    }
  }

 protected:
  std::vector<std::size_t> numLocalBins() const final {
    typename GridType::index_t idx = m_grid.numLocalBins();
    std::vector<std::size_t> result;
    for (std::size_t i = 0; i < DIM; i++) {
      result.push_back(idx[i]);
    }
    return result;
  }

  const TrackingVolume*& atLocalBins(
      const std::vector<std::size_t> indices) final {
    throw_assert(indices.size() == DIM, "Invalid number of indices");
    typename GridType::index_t idx;
    for (std::size_t i = 0; i < DIM; i++) {
      idx[i] = indices[i];
    }
    return m_grid.atLocalBins(idx);
  }

  const TrackingVolume* atLocalBins(
      const std::vector<std::size_t> indices) const final {
    throw_assert(indices.size() == DIM, "Invalid number of indices");
    typename GridType::index_t idx;
    for (std::size_t i = 0; i < DIM; i++) {
      idx[i] = indices[i];
    }
    return m_grid.atLocalBins(idx);
  }

 private:
  GridType m_grid;
};

}  // namespace Acts
