// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/AnyGridView.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/GridAccessHelpers.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"

#include <ostream>
#include <stdexcept>
#include <vector>

namespace Acts {

/// @brief Base class for material accessors, this is needed
/// for the I/O of the different grid material types, in the actual
/// implementation the material accessor is a template parameter.
struct IGridMaterialAccessor {
  virtual ~IGridMaterialAccessor() = default;
};

/// @brief  This is an accessor for cases where the material is directly stored
/// in the grid, it simply forwards the grid entry in const and non-const way.
struct GridMaterialAccessor : public IGridMaterialAccessor {
  /// @brief  Broadcast the type of the material slab
  using grid_value_type = MaterialSlab;
  /// @brief  Direct const access to the material slap sorted in the grid
  /// @tparam grid_type the type of the grid, also defines the point type
  /// @param grid the grid
  /// @param point the lookup point (already casted from global, or filled from local)
  ///
  /// @return the material slab from the grid bin associated to the lookup point
  template <typename grid_type>
  inline const MaterialSlab& slab(
      grid_type& grid, const typename grid_type::point_t& point) const {
    return grid.atPosition(point);
  }
};

/// @brief  This is an accessor for cases where the material is indexed
///
/// It can work for globally or locally indexed materials, depending if the
/// material store is shared or not.
struct IndexedMaterialAccessor : public IGridMaterialAccessor {
  /// Constructor for the indexed material accessor
  /// @param idxMaterial is the indexed material vector
  explicit IndexedMaterialAccessor(
      std::shared_ptr<std::vector<MaterialSlab>> idxMaterial)
      : IGridMaterialAccessor(), indexedMaterial(std::move(idxMaterial)) {}

  /// Broadcast the grid_value_type
  using grid_value_type = std::size_t;

  /// @brief The storage vector of the material
  std::shared_ptr<std::vector<MaterialSlab>> indexedMaterial = nullptr;

  /// @brief  Direct const access to the material slap sorted in the grid
  ///
  /// @tparam grid_type the type of the grid, also defines the point type
  ///
  /// @param grid the grid holding the indices into the global material vector
  /// @param point the lookup point (already casted from global, or filled from local)
  ///
  /// @return the material slab from the grid bin associated to the lookup point
  template <typename grid_type>
  inline const MaterialSlab& slab(
      const grid_type& grid, const typename grid_type::point_t& point) const {
    auto index = grid.atPosition(point);
    return (*indexedMaterial)[index];
  }
};

/// Intermediate interface to the grid surface material given access to the grid
/// and the material accessor.
template <typename grid_value_t>
class IGridSurfaceMaterial : public ISurfaceMaterial {
 public:
  /// @brief Accessor to the grid interface
  virtual const IGrid& grid() const = 0;

  /// @brief Accessor to the material accessor
  virtual const IGridMaterialAccessor& materialAccessor() const = 0;

  /// @brief Accessor to the bound to grid local delegate
  virtual const GridAccess::IBoundToGridLocal& boundToGridLocal() const = 0;

  /// @brief Accessor to the global to grid local delegate
  virtual const GridAccess::IGlobalToGridLocal& globalToGridLocal() const = 0;

  /// Return the type erased grid view
  virtual AnyGridView<grid_value_t> gridView() = 0;

  /// Return the type erased (const) grid view
  virtual AnyGridConstView<grid_value_t> gridConstView() const = 0;
};

/// @brief GridSurfaceMaterialT
///
/// It extends the @c ISurfaceMaterial base class and allows to create
/// material maps associated to a grid structure
///
/// @tparam grid_type is the type of the grid used here
/// @tparam material_accessor_type is the type of the accessor to the material
///
/// It is templated on the material type and a slab accessor type in order
/// to allow it to be used in the material recording as well.
template <typename grid_t, typename material_accessor_t>
class GridSurfaceMaterialT
    : public IGridSurfaceMaterial<
          typename material_accessor_t::grid_value_type> {
 public:
  // Definition of bound (on surface) to grid local representation delegate
  using BoundToGridLocalDelegate =
      OwningDelegate<typename grid_t::point_t(const Vector2&),
                     GridAccess::IBoundToGridLocal>;

  // Definition of global to grid local representation delegate
  using GlobalToGridLocalDelegate =
      OwningDelegate<typename grid_t::point_t(const Vector3&),
                     GridAccess::IGlobalToGridLocal>;

  /// Broadcast grid type
  using grid_type = grid_t;

  /// Broadcast material accessor type
  using material_accessor_type = material_accessor_t;

  /// @brief Constructor for indexed surface material
  ///
  /// @param grid the index grid steering the access to the material vector
  /// @param materialAccessor the material accessor: from grid, from indexed vector
  /// @param boundToGridLocal the delegation from bound to grid local frame
  /// @param globalToGridLocal the delegation from global into grid local frame
  GridSurfaceMaterialT(grid_type&& grid,
                       material_accessor_type&& materialAccessor,
                       BoundToGridLocalDelegate boundToGridLocal,
                       GlobalToGridLocalDelegate globalToGridLocal)
      : m_grid(std::move(grid)),
        m_materialAccessor(std::move(materialAccessor)),
        m_globalToGridLocal(std::move(globalToGridLocal)),
        m_boundToGridLocal(std::move(boundToGridLocal)) {
    if (!m_globalToGridLocal.connected()) {
      throw std::invalid_argument(
          "GridSurfaceMaterialT: GlobalToGridLocalDelegate is not connected.");
    }
    if (!m_boundToGridLocal.connected()) {
      throw std::invalid_argument(
          "GridSurfaceMaterialT: BoundToGridLocalDelegate is not connected.");
    }
  }

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector2&) const
  const MaterialSlab& materialSlab(const Vector2& lp) const final {
    return m_materialAccessor.slab(m_grid, m_boundToGridLocal(lp));
  }

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector3&) const
  const MaterialSlab& materialSlab(const Vector3& gp) const final {
    return m_materialAccessor.slab(m_grid, m_globalToGridLocal(gp));
  }

  /// Output Method for std::ostream, to be overloaded by child classes
  std::ostream& toStream(std::ostream& sl) const final {
    sl << "GridSurfaceMaterial - material access via accessor.";
    return sl;
  }

  /// @brief Accessor to the grid
  const grid_type& grid() const final { return m_grid; }

  // Return a type-erased indexed grid view
  AnyGridView<typename material_accessor_t::grid_value_type> gridView() final {
    return AnyGridView<typename material_accessor_t::grid_value_type>(m_grid);
  }

  // Return a type-erased indexed const grid view
  AnyGridConstView<typename material_accessor_t::grid_value_type>
  gridConstView() const final {
    return AnyGridConstView<typename material_accessor_t::grid_value_type>(
        m_grid);
  }

  /// @brief Accessor to the material accessor
  const material_accessor_type& materialAccessor() const final {
    return m_materialAccessor;
  }

  /// @brief Accessor to the bound to grid local delegate
  const GridAccess::IBoundToGridLocal& boundToGridLocal() const final {
    return *(m_boundToGridLocal.instance());
  }

  /// @brief Accessor to the bound to grid local delegate
  const BoundToGridLocalDelegate& boundToGridLocalDelegate() const {
    return m_boundToGridLocal;
  }

  /// @brief Accessor to the global to grid local delegate
  const GridAccess::IGlobalToGridLocal& globalToGridLocal() const final {
    return *(m_globalToGridLocal.instance());
  }

  /// @brief Accessor to the global to grid local delegate
  const GlobalToGridLocalDelegate& globalToGridLocalDelegate() const {
    return m_globalToGridLocal;
  }

 private:
  /// @brief The grid
  grid_type m_grid;

  /// @brief The stored material accessor
  material_accessor_type m_materialAccessor;

  /// The global to grid local delegate
  GlobalToGridLocalDelegate m_globalToGridLocal;

  /// The bound to grid local delegate
  BoundToGridLocalDelegate m_boundToGridLocal;
};

// Indexed Surface material
template <typename grid_type>
using IndexedSurfaceMaterial =
    GridSurfaceMaterialT<grid_type, IndexedMaterialAccessor>;

// Grid Surface material
template <typename grid_type>
using GridSurfaceMaterial =
    GridSurfaceMaterialT<grid_type, GridMaterialAccessor>;

}  // namespace Acts
