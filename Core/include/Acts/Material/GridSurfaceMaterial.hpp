// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/GridAccessHelpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <ostream>
#include <stdexcept>
#include <vector>

namespace Acts {

/// @brief  This is an accessor for cases where the material is directly stored
/// in the grid, it simply forwards the grid entry in const and non-const way.
struct GridMaterialAccessor {
  /// @brief  Direct const access to the material slap sorted in the grid
  /// @tparam grid_type the type of the grid, also defines the point type
  /// @param grid the grid
  /// @param point the lookup point (already casted from global, or filled from local)
  ///
  /// @return the material slab from the grid bin associated to the lookup point
  template <typename grid_type>
  inline MaterialSlab& slab(grid_type& grid,
                            const typename grid_type::point_t& point) const {
    return grid.atPosition(point);
  }

  /// @brief Scale the material (by scaling the thickness)
  ///
  /// @param grid the grid (ignored)
  /// @param scale the amount of the scaling
  ///
  /// @note this is not particularly fast
  template <typename grid_type>
  void scale(grid_type& grid, ActsScalar scale) {
    for (auto& m : grid.values()) {
      m.scaleThickness(scale);
    }
  }
};

/// @brief  This is an accessor for cases where the material is filled in a vector
/// and then indexed by the grid
struct IndexedMaterialAccessor {
  /// @brief The internal storage of the material
  std::vector<MaterialSlab> material;
  /// @brief  Direct const access to the material slap sorted in the grid
  /// @tparam grid_type the type of the grid, also defines the point type
  /// @param grid the grid
  /// @param point the lookup point (already casted from global, or filled from local)
  ///
  /// @return the material slab from the grid bin associated to the lookup point
  template <typename grid_type>
  inline const MaterialSlab& slab(
      const grid_type& grid, const typename grid_type::point_t& point) const {
    auto index = grid.atPosition(point);
    return material[index];
  }

  /// @brief Scale the material (by scaling the thickness)
  ///
  /// @param scale the amount of the scaling
  template <typename grid_type>
  void scale(grid_type& /*unused*/, ActsScalar scale) {
    for (auto& m : material) {
      m.scaleThickness(scale);
    }
  }
};

/// @class GridSurfaceMaterialBase
///
/// It extends the @c ISurfaceMaterial base class and allows to create
/// material maps associated to a grid structure
///
/// @tparam grid_type is the type of the grid used here
/// @tparam material_accessor_type is the type of the accessor to the material
///
/// It is templated on the material type and a slab accessor type in order
/// to allow it to be used in the material recording as well.
template <typename grid_t, typename material_accessor_t = GridMaterialAccessor>
class GridSurfaceMaterialBase : public ISurfaceMaterial {
 public:
  /// Broadcast grid type
  using grid_type = grid_t;

  /// Broadcast material accessor type
  using material_accessor_type = material_accessor_t;

  /// @brief Constructor for indexed surface material
  /// @param grid the index grid steering the access to the material vector
  /// @param materialAccessor the material accessor: from grid, from indexed vector
  /// @param gcasts global casts -> casts a grid position from global
  /// @param laccessors local accessors -> accessors to the local grid
  /// @param transform  transform from global frame into map frame
  GridSurfaceMaterialBase(grid_type&& grid,
                          material_accessor_type&& materialAccessor,
                          const std::vector<BinningValue>& gcasts,
                          const std::vector<std::size_t>& laccessors,
                          const Transform3& transform = Transform3::Identity())
      : m_grid(std::move(grid)),
        m_materialAccessor(std::move(materialAccessor)),
        m_globalCasts(gcasts),
        m_localAccessors(laccessors),
        m_transform(transform) {
    if (gcasts.size() != grid_type::DIM) {
      throw std::invalid_argument(
          "The number of casts must match the grid dimension");
    }

    if (gcasts.size() != laccessors.size()) {
      throw std::invalid_argument(
          "The number of casts must match the number of local accessors");
    }
  }

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector2&) const
  const MaterialSlab& materialSlab(const Vector2& lp) const final {
    return m_materialAccessor.slab(
        m_grid,
        GridAccessHelpers::accessLocal<grid_type>(lp, m_localAccessors));
  }

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector3&) const
  const MaterialSlab& materialSlab(const Vector3& gp) const final {
    // Access via (transformed) global position lookup
    return m_materialAccessor.slab(
        m_grid, GridAccessHelpers::castPosition<grid_type>(m_transform * gp,
                                                           m_globalCasts));
  }

  /// @copydoc ISurfaceMaterial::materialSlab(std::size_t bin0, std::size_t bin1) const
  ///
  /// The ISurface material class expects bins in [ 0 - n ] x [ 0, m ], we will
  /// convert into this format, but it is gonna be slow
  const MaterialSlab& materialSlab(std::size_t bin0,
                                   std::size_t bin1) const final {
    // Access via bin0 and bin1
    Acts::Vector2 lposition =
        GridAccessHelpers::toLocal(m_grid, bin0, bin1, m_localAccessors[0u]);
    return materialSlab(lposition);
  }

  /// Scale operator
  ///
  /// @param scale is the scale factor applied
  ISurfaceMaterial& operator*=(ActsScalar scale) final {
    m_materialAccessor.scale(m_grid, scale);
    return (*this);
  }

  /// Output Method for std::ostream, to be overloaded by child classes
  std::ostream& toStream(std::ostream& sl) const final {
    sl << "GridSurfaceMaterial - material access via accessor.";
    return sl;
  }

 private:
  /// @brief The grid
  grid_type m_grid;

  /// @brief The stored material
  material_accessor_type m_materialAccessor;

  /// @brief The casts from global to local
  std::vector<BinningValue> m_globalCasts;

  /// @brief Local accessors - empty assumes direct translation
  std::vector<std::size_t> m_localAccessors;

  /// @brief The transform to the local frame
  Transform3 m_transform;
};

// Indexed Surface material
template <typename grid_type>
using IndexedSurfaceMaterial =
    GridSurfaceMaterialBase<grid_type, IndexedMaterialAccessor>;

// Grid Surface material
template <typename grid_type>
using GridSurfaceMaterial = GridSurfaceMaterialBase<grid_type>;

}  // namespace Acts
