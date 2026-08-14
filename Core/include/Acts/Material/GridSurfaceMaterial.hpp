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
#include "Acts/Utilities/Diagnostics.hpp"
#include "Acts/Utilities/Grid.hpp"

#include <memory>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace Acts {

/// @addtogroup material
/// @{

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

  /// @brief Scale the material (by scaling the thickness)
  ///
  /// @param grid the grid (ignored)
  /// @param scale the amount of the scaling
  ///
  /// @note this is not particularly fast
  template <typename grid_type>
  void scale(grid_type& grid, double scale) {
    // Loop through the grid bins, get the indices and scale the material
    for (std::size_t ib = 0; ib < grid.size(); ++ib) {
      grid.at(ib).scaleThickness(static_cast<float>(scale));
    }
  }
};

/// @brief  This is an accessor for cases where the material is filled in a vector
/// and then indexed by the grid
struct IndexedMaterialAccessor : public IGridMaterialAccessor {
  /// Broadcast the grid_value_type
  using grid_value_type = std::size_t;

  /// @brief The internal storage of the material
  /// @param mmaterial Vector of material slabs to store and access by index
  explicit IndexedMaterialAccessor(std::vector<MaterialSlab>&& mmaterial)
      : IGridMaterialAccessor(), material(std::move(mmaterial)) {}

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
      const grid_type& grid, const typename grid_type::point_t& point) const
    requires(std::is_same_v<typename grid_type::value_type, grid_value_type>)
  {
    std::size_t index = grid.atPosition(point);
    return material[index];
  }

  /// @brief Scale the material (by scaling the thickness)
  ///
  /// @param scale the amount of the scaling
  template <typename grid_type>
  void scale(grid_type& /*grid*/, double scale) {
    for (auto& m : material) {
      m.scaleThickness(static_cast<float>(scale));
    }
  }
};

/// @brief  This is an accessor for cases where the material is filled in a global
/// material vector that is accessed from the different material grids.
struct GloballyIndexedMaterialAccessor : public IGridMaterialAccessor {
  /// Constructor with global material vector
  /// @param gMaterial Shared pointer to global material vector
  /// @param shared Whether material entries are shared across grid points
  explicit GloballyIndexedMaterialAccessor(
      std::shared_ptr<std::vector<MaterialSlab>> gMaterial, bool shared = false)
      : IGridMaterialAccessor(),
        globalMaterial(std::move(gMaterial)),
        sharedEntries(shared) {}

  /// Broadcast the grid_value_type
  using grid_value_type = std::size_t;

  /// @brief The internal storage of the material
  std::shared_ptr<std::vector<MaterialSlab>> globalMaterial = nullptr;

  /// Indicate if you have entries bins across different grids, e.g. by
  /// running a compression/clustering algorithm.
  ///
  /// It is the responsibility of the user to set this flag correctly.
  bool sharedEntries = false;

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
    return (*globalMaterial)[index];
  }

  /// @brief Scale the material (by scaling the thickness)
  ///
  /// @param grid the grid holding the indices into the global material vector
  /// @param scale the amount of the scaling
  ///
  /// @note this will scale only the bins touched by this grid, however,
  /// if there are shared bins, then it will throw an exception as the
  /// outcome is unpredictable.
  ///
  template <typename grid_type>
  void scale(grid_type& grid, double scale) {
    if (sharedEntries) {
      throw std::invalid_argument(
          "GloballyIndexedMaterialAccessor: shared entry scaling is not "
          "supported.");
    }
    // Loop through the grid bins, get the indices and scale the material
    for (std::size_t ib = 0; ib < grid.size(); ++ib) {
      auto index = grid.at(ib);
      (*globalMaterial)[index].scaleThickness(static_cast<float>(scale));
    }
  }
};

/// Base class for the concrete templated grid surface material types.
/// This allows referning to all template instances as the same base class type.
class IGridSurfaceMaterialBase : public ISurfaceMaterial {};

/// Intermediate interface to the grid surface material given access to the grid
/// and the material accessor.
template <typename grid_value_t>
class IGridSurfaceMaterial : public IGridSurfaceMaterialBase {
 public:
  /// @brief Accessor to the grid interface
  /// @return Reference to the grid interface
  virtual const IGrid& grid() const = 0;

  /// @brief Accessor to the material accessor
  /// @return Reference to the material accessor
  virtual const IGridMaterialAccessor& materialAccessor() const = 0;

  /// Return the type erased grid view
  /// @return Type-erased grid view for accessing grid contents
  virtual AnyGridView<grid_value_t> gridView() = 0;

  /// Return the type erased (const) grid view
  /// @return Type-erased const grid view for read-only access to grid contents
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
///
/// The grid is always 2D: local (bound) lookup is assumed to already be
/// expressed in the grid's own 2D coordinate system, i.e. @c lp[0] maps
/// directly onto axis 0 and @c lp[1] onto axis 1. Only local lookup is
/// supported - there is no global (Vector3) lookup, since that would require
/// a global-to-local coordinate transform.
template <typename grid_t, typename material_accessor_t>
class GridSurfaceMaterialT
    : public IGridSurfaceMaterial<
          typename material_accessor_t::grid_value_type> {
  static_assert(
      grid_t::DIM == 2,
      "GridSurfaceMaterialT requires a 2D grid; local lookup is "
      "always expressed as (loc0, loc1) directly onto the grid axes.");

 public:
  /// Broadcast grid type
  using grid_type = grid_t;

  /// Broadcast material accessor type
  using material_accessor_type = material_accessor_t;

  /// @brief Constructor for indexed surface material
  ///
  /// @param grid the index grid steering the access to the material vector
  /// @param materialAccessor the material accessor: from grid, from indexed vector
  GridSurfaceMaterialT(grid_type&& grid,
                       material_accessor_type&& materialAccessor)
      : m_grid(std::move(grid)),
        m_materialAccessor(std::move(materialAccessor)) {}

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector2&) const
  const MaterialSlab& materialSlab(const Vector2& lp) const final {
    return m_materialAccessor.slab(m_grid,
                                   typename grid_type::point_t{lp[0], lp[1]});
  }

  /// @copydoc ISurfaceMaterial::localAxisDirections() const
  std::vector<AxisDirection> localAxisDirections() const final { return {}; }

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector3&) const
  /// @deprecated Global (Vector3) lookup is not supported; use
  ///             materialSlab(const Vector2&) with a prior
  ///             Surface::globalToLocal() call instead.
  [[deprecated(
      "Global lookup is not supported; use materialSlab(const Vector2& lp) "
      "with a prior Surface::globalToLocal() call "
      "instead")]] const MaterialSlab&
  materialSlab(const Vector3& /*gp*/) const final {
    throw std::logic_error(
        "GridSurfaceMaterialT: global (Vector3) material lookup is not "
        "supported, use materialSlab(const Vector2&) instead.");
  }

  using ISurfaceMaterial::materialSlab;

  /// Scale operator
  ///
  /// @param factor is the scale factor applied
  /// @return Reference to this surface material for method chaining
  ISurfaceMaterial& scale(double factor) final {
    m_materialAccessor.scale(m_grid, factor);
    return (*this);
  }

  /// Output Method for std::ostream, to be overloaded by child classes
  /// @param sl Output stream to write to
  /// @return Reference to the output stream
  std::ostream& toStream(std::ostream& sl) const final {
    sl << "GridSurfaceMaterial - material access via accessor.";
    return sl;
  }

  /// @brief Accessor to the grid
  /// @return Reference to the underlying grid
  const grid_type& grid() const final { return m_grid; }

  // Return a type-erased indexed grid view
  /// @return Type-erased grid view for accessing grid contents
  AnyGridView<typename material_accessor_t::grid_value_type> gridView() final {
    return AnyGridView<typename material_accessor_t::grid_value_type>(m_grid);
  }

  // Return a type-erased indexed const grid view
  /// @return Type-erased const grid view for read-only access to grid contents
  AnyGridConstView<typename material_accessor_t::grid_value_type>
  gridConstView() const final {
    return AnyGridConstView<typename material_accessor_t::grid_value_type>(
        m_grid);
  }

  /// @brief Accessor to the material accessor
  /// @return Reference to the material accessor
  const material_accessor_type& materialAccessor() const final {
    return m_materialAccessor;
  }

 private:
  /// @brief The grid
  grid_type m_grid;

  /// @brief The stored material accessor
  material_accessor_type m_materialAccessor;
};

/// @brief Concrete, non-template surface material backed by a grid that
/// stores material slabs directly.
///
/// Wraps a type-erased @c IGridSurfaceMaterial<MaterialSlab> built with a
/// @c GridMaterialAccessor so that callers do not need to name the concrete
/// grid/axis types.
class GridSurfaceMaterial final : public ISurfaceMaterial {
 public:
  /// Construct from a grid that stores material slabs directly.
  ///
  /// @param gridMaterial type-erased grid material implementation
  /// @param splitFactor pre/post splitting directive
  /// @param mappingType type of surface mapping associated to the surface
  /// @throws std::invalid_argument if @p gridMaterial is null
  explicit GridSurfaceMaterial(
      std::unique_ptr<IGridSurfaceMaterial<MaterialSlab>> gridMaterial,
      double splitFactor = 1., MappingType mappingType = MappingType::Default)
      : ISurfaceMaterial(splitFactor, mappingType),
        m_gridMaterial(std::move(gridMaterial)) {
    if (m_gridMaterial == nullptr) {
      throw std::invalid_argument(
          "GridSurfaceMaterial: grid material must not be null.");
    }
  }

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector2&) const
  const MaterialSlab& materialSlab(const Vector2& lp) const final {
    return m_gridMaterial->materialSlab(lp);
  }

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector3&) const
  ACTS_PUSH_IGNORE_DEPRECATED()
  [[deprecated(
      "Use materialSlab(const Vector2& lp) with a prior "
      "Surface::globalToLocal() call instead")]] const MaterialSlab&
  materialSlab(const Vector3& gp) const final {
    return m_gridMaterial->materialSlab(gp);
  }
  ACTS_POP_IGNORE_DEPRECATED()

  using ISurfaceMaterial::materialSlab;

  /// @copydoc ISurfaceMaterial::localAxisDirections() const
  std::vector<AxisDirection> localAxisDirections() const final {
    return m_gridMaterial->localAxisDirections();
  }

  /// @copydoc ISurfaceMaterial::scale(double)
  ISurfaceMaterial& scale(double factor) final {
    m_gridMaterial->scale(factor);
    return *this;
  }

  /// @copydoc ISurfaceMaterial::toStream(std::ostream&) const
  std::ostream& toStream(std::ostream& sl) const final {
    return m_gridMaterial->toStream(sl);
  }

  /// Return the type-erased grid.
  /// @return Reference to the type-erased grid
  const IGrid& grid() const { return m_gridMaterial->grid(); }

  /// Return a type-erased const view of the stored material slabs.
  /// @return Const view onto the stored grid values
  AnyGridConstView<MaterialSlab> gridConstView() const {
    return m_gridMaterial->gridConstView();
  }

  /// Return the material accessor.
  /// @return Reference to the material accessor
  const IGridMaterialAccessor& materialAccessor() const {
    return m_gridMaterial->materialAccessor();
  }

 private:
  std::unique_ptr<IGridSurfaceMaterial<MaterialSlab>> m_gridMaterial;
};

/// @brief Concrete, non-template surface material backed by a locally
/// indexed material grid.
///
/// Wraps a type-erased @c IGridSurfaceMaterial<std::size_t> built with an
/// @c IndexedMaterialAccessor so that callers do not need to name the
/// concrete grid/axis types.
class IndexedGridSurfaceMaterial final : public ISurfaceMaterial {
 public:
  /// Construct from a grid that stores material indices, indexed by a
  /// per-surface (locally owned) material vector.
  ///
  /// @param gridMaterial type-erased indexed grid material implementation
  /// @param splitFactor pre/post splitting directive
  /// @param mappingType type of surface mapping associated to the surface
  /// @throws std::invalid_argument if @p gridMaterial is null or does not use
  ///         an @c IndexedMaterialAccessor
  explicit IndexedGridSurfaceMaterial(
      std::unique_ptr<IGridSurfaceMaterial<std::size_t>> gridMaterial,
      double splitFactor = 1., MappingType mappingType = MappingType::Default)
      : ISurfaceMaterial(splitFactor, mappingType),
        m_gridMaterial(std::move(gridMaterial)) {
    if (m_gridMaterial == nullptr) {
      throw std::invalid_argument(
          "IndexedGridSurfaceMaterial: grid material must not be null.");
    }
    if (dynamic_cast<const IndexedMaterialAccessor*>(
            &m_gridMaterial->materialAccessor()) == nullptr) {
      throw std::invalid_argument(
          "IndexedGridSurfaceMaterial: grid material must use an "
          "IndexedMaterialAccessor.");
    }
  }

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector2&) const
  const MaterialSlab& materialSlab(const Vector2& lp) const final {
    return m_gridMaterial->materialSlab(lp);
  }

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector3&) const
  ACTS_PUSH_IGNORE_DEPRECATED()
  [[deprecated(
      "Use materialSlab(const Vector2& lp) with a prior "
      "Surface::globalToLocal() call instead")]] const MaterialSlab&
  materialSlab(const Vector3& gp) const final {
    return m_gridMaterial->materialSlab(gp);
  }
  ACTS_POP_IGNORE_DEPRECATED()

  using ISurfaceMaterial::materialSlab;

  /// @copydoc ISurfaceMaterial::localAxisDirections() const
  std::vector<AxisDirection> localAxisDirections() const final {
    return m_gridMaterial->localAxisDirections();
  }

  /// @copydoc ISurfaceMaterial::scale(double)
  ISurfaceMaterial& scale(double factor) final {
    m_gridMaterial->scale(factor);
    return *this;
  }

  /// @copydoc ISurfaceMaterial::toStream(std::ostream&) const
  std::ostream& toStream(std::ostream& sl) const final {
    return m_gridMaterial->toStream(sl);
  }

  /// Return the type-erased grid.
  /// @return Reference to the type-erased grid
  const IGrid& grid() const { return m_gridMaterial->grid(); }

  /// Return a type-erased const view of the stored material indices.
  /// @return Const view onto the stored grid values
  AnyGridConstView<std::size_t> gridConstView() const {
    return m_gridMaterial->gridConstView();
  }

  /// Return the material accessor.
  /// @return Reference to the material accessor
  const IGridMaterialAccessor& materialAccessor() const {
    return m_gridMaterial->materialAccessor();
  }

 private:
  std::unique_ptr<IGridSurfaceMaterial<std::size_t>> m_gridMaterial;
};

/// @brief Concrete, non-template surface material backed by a globally
/// indexed material grid.
///
/// Wraps a type-erased @c IGridSurfaceMaterial<std::size_t> built with a
/// @c GloballyIndexedMaterialAccessor so that callers do not need to name the
/// concrete grid/axis types.
class GloballyIndexedGridSurfaceMaterial final : public ISurfaceMaterial {
 public:
  /// Construct from a grid that stores material indices into a globally
  /// shared material vector.
  ///
  /// @param gridMaterial type-erased indexed grid material implementation
  /// @param splitFactor pre/post splitting directive
  /// @param mappingType type of surface mapping associated to the surface
  /// @throws std::invalid_argument if @p gridMaterial is null or does not use
  ///         a @c GloballyIndexedMaterialAccessor
  explicit GloballyIndexedGridSurfaceMaterial(
      std::unique_ptr<IGridSurfaceMaterial<std::size_t>> gridMaterial,
      double splitFactor = 1., MappingType mappingType = MappingType::Default)
      : ISurfaceMaterial(splitFactor, mappingType),
        m_gridMaterial(std::move(gridMaterial)) {
    if (m_gridMaterial == nullptr) {
      throw std::invalid_argument(
          "GloballyIndexedGridSurfaceMaterial: grid material must not be "
          "null.");
    }
    if (dynamic_cast<const GloballyIndexedMaterialAccessor*>(
            &m_gridMaterial->materialAccessor()) == nullptr) {
      throw std::invalid_argument(
          "GloballyIndexedGridSurfaceMaterial: grid material must use a "
          "GloballyIndexedMaterialAccessor.");
    }
  }

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector2&) const
  const MaterialSlab& materialSlab(const Vector2& lp) const final {
    return m_gridMaterial->materialSlab(lp);
  }

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector3&) const
  ACTS_PUSH_IGNORE_DEPRECATED()
  [[deprecated(
      "Use materialSlab(const Vector2& lp) with a prior "
      "Surface::globalToLocal() call instead")]] const MaterialSlab&
  materialSlab(const Vector3& gp) const final {
    return m_gridMaterial->materialSlab(gp);
  }
  ACTS_POP_IGNORE_DEPRECATED()

  using ISurfaceMaterial::materialSlab;

  /// @copydoc ISurfaceMaterial::localAxisDirections() const
  std::vector<AxisDirection> localAxisDirections() const final {
    return m_gridMaterial->localAxisDirections();
  }

  /// @copydoc ISurfaceMaterial::scale(double)
  ISurfaceMaterial& scale(double factor) final {
    m_gridMaterial->scale(factor);
    return *this;
  }

  /// @copydoc ISurfaceMaterial::toStream(std::ostream&) const
  std::ostream& toStream(std::ostream& sl) const final {
    return m_gridMaterial->toStream(sl);
  }

  /// Return the type-erased grid.
  /// @return Reference to the type-erased grid
  const IGrid& grid() const { return m_gridMaterial->grid(); }

  /// Return a type-erased const view of the stored material indices.
  /// @return Const view onto the stored grid values
  AnyGridConstView<std::size_t> gridConstView() const {
    return m_gridMaterial->gridConstView();
  }

  /// Return the material accessor.
  /// @return Reference to the material accessor
  const IGridMaterialAccessor& materialAccessor() const {
    return m_gridMaterial->materialAccessor();
  }

 private:
  std::unique_ptr<IGridSurfaceMaterial<std::size_t>> m_gridMaterial;
};

/// @}

}  // namespace Acts
