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

#include <memory>
#include <vector>

namespace Acts::GridSurfaceMaterialFactory {

/// Create and fill from two axes
///
/// @param axis0 the type-erased axis in direction 0
/// @param axis1 the type-erased axis in direction 1
/// @param materialAccessor the material accessor
/// @param payload the grid payload in 2D (material slab / indices)
/// the payload has to be column major, i.e. [i0][i1]
///
/// @return a unique pointer to the surface material
template <typename material_accessor_t>
std::unique_ptr<
    IGridSurfaceMaterial<typename material_accessor_t::grid_value_type>>
create2D(
    const IAxis& axis0, const IAxis& axis1,
    material_accessor_t&& materialAccessor,
    const std::vector<
        std::vector<typename material_accessor_t::grid_value_type>>& payload) {
  // Validate axis compatibility
  auto ism =
      axis0.visit([&]<typename AxisTypeA>(const AxisTypeA& axisA)
                      -> std::unique_ptr<IGridSurfaceMaterial<
                          typename material_accessor_t::grid_value_type>> {
        return axis1.visit(
            [&]<typename AxisTypeB>(const AxisTypeB& axisB)
                -> std::unique_ptr<IGridSurfaceMaterial<
                    typename material_accessor_t::grid_value_type>> {
              using GridType =
                  Grid<typename material_accessor_t::grid_value_type, AxisTypeA,
                       AxisTypeB>;
              return std::make_unique<
                  GridSurfaceMaterialT<GridType, material_accessor_t>>(
                  GridType(axisA, axisB),
                  std::forward<material_accessor_t>(materialAccessor));
            });
      });

  // Fill it via the grid view
  AnyGridView<typename material_accessor_t::grid_value_type> gv =
      ism->gridView();
  auto indices = gv.multiAxisAny().getNBinsAny();
  for (std::size_t i0 = 0; i0 < indices[0]; ++i0) {
    for (std::size_t i1 = 0; i1 < indices[1]; ++i1) {
      // Offset comes from overflow/underflow bin
      gv.atLocalBins({i0 + 1, i1 + 1}) = payload[i0][i1];
    }
  }

  return ism;
}

/// The resolved functions to reduce compile time template bloat
/// - GridMaterial 2D
/// @param axis0 the axis in direction 0
/// @param axis1 the axis in direction 1
/// @param materialAccessor the material accessor
/// @param payload the grid payload (material slab / indices)
std::unique_ptr<GridSurfaceMaterial> create(
    const IAxis& axis0, const IAxis& axis1,
    GridMaterialAccessor&& materialAccessor,
    const std::vector<std::vector<MaterialSlab>>& payload);

/// The resolved functions to reduce compile time template bloat
/// - IndexedMaterial 2D
/// @param axis0 the axis in direction 0
/// @param axis1 the axis in direction 1
/// @param materialAccessor the material accessor
/// @param payload the grid payload (material slab / indices)
std::unique_ptr<IndexedGridSurfaceMaterial> create(
    const IAxis& axis0, const IAxis& axis1,
    IndexedMaterialAccessor&& materialAccessor,
    const std::vector<std::vector<IndexedMaterialAccessor::grid_value_type>>&
        payload);

/// The resolved functions to reduce compile time template bloat
/// - GloballyIndexedMaterial 2D
/// @param axis0 the axis in direction 0
/// @param axis1 the axis in direction 1
/// @param materialAccessor the material accessor
/// @param payload the grid payload (material slab / indices)
std::unique_ptr<GloballyIndexedGridSurfaceMaterial> create(
    const IAxis& axis0, const IAxis& axis1,
    GloballyIndexedMaterialAccessor&& materialAccessor,
    const std::vector<
        std::vector<GloballyIndexedMaterialAccessor::grid_value_type>>&
        payload);

}  // namespace Acts::GridSurfaceMaterialFactory
