// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/GridAccessHelpers.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"

#include <memory>
#include <vector>

namespace Acts::GridSurfaceMaterialFactory {

/// Create and fill from a single proto axis
///
/// @param pAxis the type of the ProtoAxis
/// @param materialAccessor the material accessor
/// @param boundToGridLocal the delegate from bound to grid local frame
/// @param globalToGridLocal the delegate from global into grid local frame
/// @param payload the grid payload (material slab / indices)
///
/// @return a unique pointer to the surface material
template <typename material_accessor_t>
std::unique_ptr<
    IGridSurfaceMaterial<typename material_accessor_t::grid_value_type>>
create1D(
    const ProtoAxis& pAxis, material_accessor_t&& materialAccessor,
    GridAccess::BoundToGridLocal1DimDelegate boundToGridLocal,
    GridAccess::GlobalToGridLocal1DimDelegate globalToGridLocal,
    const std::vector<typename material_accessor_t::grid_value_type>& payload) {
  // Visit the axis type and create the grid surface material
  auto ism = pAxis.getAxis().visit(
      [&]<typename AxisType>(const AxisType& axis)
          -> std::unique_ptr<IGridSurfaceMaterial<
              typename material_accessor_t::grid_value_type>> {
        using GridType =
            Grid<typename material_accessor_t::grid_value_type, AxisType>;
        return std::make_unique<
            GridSurfaceMaterialT<GridType, material_accessor_t>>(
            GridType(axis), std::forward<material_accessor_t>(materialAccessor),
            std::move(boundToGridLocal), std::move(globalToGridLocal));
      });
  // Fill it via the grid view
  AnyGridView<typename material_accessor_t::grid_value_type> gv =
      ism->gridView();
  auto indices = gv.numLocalBins();
  for (std::size_t i0 = 0; i0 < indices[0]; ++i0) {
    // Offset comes from overflow/underflow bin
    gv.atLocalBins({i0 + 1u}) = payload[i0];
  }
  return ism;
}

/// Static creation method for the with ProtoAxis objects
///
/// @param pAxis0 proto axis in direction 0
/// @param pAxis1 proto axis in direction 1
/// @param materialAccessor the material accessor
/// @param boundToGridLocal the delegate from bound to grid local frame
/// @param globalToGridLocal the delegate from global into grid local frame
/// @param payload the grid payload in 2D (material slab / indices)
/// the payload has to be column major, i.e. [i0][i1]
///
/// @return a unique pointer to the surface material
template <typename material_accessor_t>
std::unique_ptr<
    IGridSurfaceMaterial<typename material_accessor_t::grid_value_type>>
create2D(
    const ProtoAxis& pAxis0, const ProtoAxis& pAxis1,
    material_accessor_t&& materialAccessor,
    GridAccess::BoundToGridLocal2DimDelegate boundToGridLocal,
    GridAccess::GlobalToGridLocal2DimDelegate globalToGridLocal,
    const std::vector<
        std::vector<typename material_accessor_t::grid_value_type>>& payload) {
  // Validate axis compatibility
  if (pAxis0.getAxisDirection() == pAxis1.getAxisDirection()) {
    throw std::invalid_argument(
        "createGridSurfaceMaterial: ProtoAxes must have different directions");
  }
  auto ism = pAxis0.getAxis().visit(
      [&]<typename AxisTypeA>(const AxisTypeA& axisA)
          -> std::unique_ptr<IGridSurfaceMaterial<
              typename material_accessor_t::grid_value_type>> {
        return pAxis1.getAxis().visit(
            [&]<typename AxisTypeB>(const AxisTypeB& axisB)
                -> std::unique_ptr<IGridSurfaceMaterial<
                    typename material_accessor_t::grid_value_type>> {
              using GridType =
                  Grid<typename material_accessor_t::grid_value_type, AxisTypeA,
                       AxisTypeB>;
              return std::make_unique<
                  GridSurfaceMaterialT<GridType, material_accessor_t>>(
                  GridType(axisA, axisB),
                  std::forward<material_accessor_t>(materialAccessor),
                  std::move(boundToGridLocal), std::move(globalToGridLocal));
            });
      });

  // Fill it via the grid view
  AnyGridView<typename material_accessor_t::grid_value_type> gv =
      ism->gridView();
  auto indices = gv.numLocalBins();
  for (std::size_t i0 = 0; i0 < indices[0]; ++i0) {
    for (std::size_t i1 = 0; i1 < indices[1]; ++i1) {
      // Offset comes from overflow/underflow bin
      gv.atLocalBins({i0 + 1, i1 + 1}) = payload[i0][i1];
    }
  }

  return ism;
}

/// The resolved functions to reduce compile time template bloat
///  - GridMaterial 1D
/// @param pAxis the proto axis
/// @param materialAccessor the material accessor
/// @param boundToGridLocal the delegate from bound to grid local frame
/// @param globalToGridLocal the delegate from global into grid local frame
/// @param payload the grid payload (material slab / indices)
std::unique_ptr<IGridSurfaceMaterial<MaterialSlab>> create(
    const ProtoAxis& pAxis, GridMaterialAccessor&& materialAccessor,
    GridAccess::BoundToGridLocal1DimDelegate boundToGridLocal,
    GridAccess::GlobalToGridLocal1DimDelegate globalToGridLocal,
    const std::vector<MaterialSlab>& payload);

/// The resolved functions to reduce compile time template bloat
/// - IndexedMaterial 1D
/// @param pAxis the proto axis
/// @param materialAccessor the material accessor
/// @param boundToGridLocal the delegate from bound to grid local frame
/// @param globalToGridLocal the delegate from global into grid local frame
/// @param payload the grid payload (material slab / indices)
std::unique_ptr<IGridSurfaceMaterial<IndexedMaterialAccessor::grid_value_type>>
create(const ProtoAxis& pAxis, IndexedMaterialAccessor&& materialAccessor,
       GridAccess::BoundToGridLocal1DimDelegate boundToGridLocal,
       GridAccess::GlobalToGridLocal1DimDelegate globalToGridLocal,
       const std::vector<IndexedMaterialAccessor::grid_value_type>& payload);

/// The resolved functions to reduce compile time template bloat
/// - GloballyIndexedMaterial 1D
/// @param pAxis the proto axis
/// @param materialAccessor the material accessor
/// @param boundToGridLocal the delegate from bound to grid local frame
/// @param globalToGridLocal the delegate from global into grid local frame
/// @param payload the grid payload (material slab / indices)
std::unique_ptr<
    IGridSurfaceMaterial<GloballyIndexedMaterialAccessor::grid_value_type>>
create(const ProtoAxis& pAxis,
       GloballyIndexedMaterialAccessor&& materialAccessor,
       GridAccess::BoundToGridLocal1DimDelegate boundToGridLocal,
       GridAccess::GlobalToGridLocal1DimDelegate globalToGridLocal,
       const std::vector<GloballyIndexedMaterialAccessor::grid_value_type>&
           payload);

/// The resolved functions to reduce compile time template bloat
/// - GridMaterial 2D
/// @param pAxis0 the proto axis in direction 0
/// @param pAxis1 the proto axis in direction 1
/// @param materialAccessor the material accessor
/// @param boundToGridLocal the delegate from bound to grid local frame
/// @param globalToGridLocal the delegate from global into grid local frame
/// @param payload the grid payload (material slab / indices)
std::unique_ptr<IGridSurfaceMaterial<MaterialSlab>> create(
    const ProtoAxis& pAxis0, const ProtoAxis& pAxis1,
    GridMaterialAccessor&& materialAccessor,
    GridAccess::BoundToGridLocal2DimDelegate boundToGridLocal,
    GridAccess::GlobalToGridLocal2DimDelegate globalToGridLocal,
    const std::vector<std::vector<MaterialSlab>>& payload);

/// The resolved functions to reduce compile time template bloat
/// - IndexedMaterial 2D
/// @param pAxis0 the proto axis in direction 0
/// @param pAxis1 the proto axis in direction 1
/// @param materialAccessor the material accessor
/// @param boundToGridLocal the delegate from bound to grid local frame
/// @param globalToGridLocal the delegate from global into grid local frame
/// @param payload the grid payload (material slab / indices)
std::unique_ptr<IGridSurfaceMaterial<IndexedMaterialAccessor::grid_value_type>>
create(const ProtoAxis& pAxis0, const ProtoAxis& pAxis1,
       IndexedMaterialAccessor&& materialAccessor,
       GridAccess::BoundToGridLocal2DimDelegate boundToGridLocal,
       GridAccess::GlobalToGridLocal2DimDelegate globalToGridLocal,
       const std::vector<std::vector<IndexedMaterialAccessor::grid_value_type>>&
           payload);

/// The resolved functions to reduce compile time template bloat
/// - GloballyIndexedMaterial 2D
/// @param pAxis0 the proto axis in direction 0
/// @param pAxis1 the proto axis in direction 1
/// @param materialAccessor the material accessor
/// @param boundToGridLocal the delegate from bound to grid local frame
/// @param globalToGridLocal the delegate from global into grid local frame
/// @param payload the grid payload (material slab / indices)
std::unique_ptr<
    IGridSurfaceMaterial<GloballyIndexedMaterialAccessor::grid_value_type>>
create(const ProtoAxis& pAxis0, const ProtoAxis& pAxis1,
       GloballyIndexedMaterialAccessor&& materialAccessor,
       GridAccess::BoundToGridLocal2DimDelegate boundToGridLocal,
       GridAccess::GlobalToGridLocal2DimDelegate globalToGridLocal,
       const std::vector<
           std::vector<GloballyIndexedMaterialAccessor::grid_value_type>>&
           payload);
}  // namespace Acts::GridSurfaceMaterialFactory
