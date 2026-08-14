// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/GridSurfaceMaterialFactory.hpp"

#include "Acts/Utilities/Diagnostics.hpp"

namespace Acts {

std::unique_ptr<GridSurfaceMaterial> GridSurfaceMaterialFactory::create(
    const IAxis& axis0, const IAxis& axis1,
    GridMaterialAccessor&& materialAccessor,
    const std::vector<std::vector<MaterialSlab>>& payload) {
  return std::make_unique<GridSurfaceMaterial>(create2D<GridMaterialAccessor>(
      axis0, axis1, std::move(materialAccessor), payload));
}

std::unique_ptr<IndexedGridSurfaceMaterial> GridSurfaceMaterialFactory::create(
    const IAxis& axis0, const IAxis& axis1,
    IndexedMaterialAccessor&& materialAccessor,
    const std::vector<std::vector<IndexedMaterialAccessor::grid_value_type>>&
        payload) {
  return std::make_unique<IndexedGridSurfaceMaterial>(
      create2D<IndexedMaterialAccessor>(axis0, axis1,
                                        std::move(materialAccessor), payload));
}

std::unique_ptr<GloballyIndexedGridSurfaceMaterial>
GridSurfaceMaterialFactory::create(
    const IAxis& axis0, const IAxis& axis1,
    GloballyIndexedMaterialAccessor&& materialAccessor,
    const std::vector<
        std::vector<GloballyIndexedMaterialAccessor::grid_value_type>>&
        payload) {
  return std::make_unique<GloballyIndexedGridSurfaceMaterial>(
      create2D<GloballyIndexedMaterialAccessor>(
          axis0, axis1, std::move(materialAccessor), payload));
}

// Unlike the declarations, the definitions of the deprecated overloads do not
// carry the deprecation attribute themselves, so they need the suppression
ACTS_PUSH_IGNORE_DEPRECATED()

std::unique_ptr<GridSurfaceMaterial> GridSurfaceMaterialFactory::create(
    const ProtoAxis& pAxis0, const ProtoAxis& pAxis1,
    GridMaterialAccessor&& materialAccessor,
    const std::vector<std::vector<MaterialSlab>>& payload) {
  return std::make_unique<GridSurfaceMaterial>(
      create2D<GridMaterialAccessor>(pAxis0.getAxis(), pAxis1.getAxis(),
                                     std::move(materialAccessor), payload));
}

std::unique_ptr<IndexedGridSurfaceMaterial> GridSurfaceMaterialFactory::create(
    const ProtoAxis& pAxis0, const ProtoAxis& pAxis1,
    IndexedMaterialAccessor&& materialAccessor,
    const std::vector<std::vector<IndexedMaterialAccessor::grid_value_type>>&
        payload) {
  return std::make_unique<IndexedGridSurfaceMaterial>(
      create2D<IndexedMaterialAccessor>(pAxis0.getAxis(), pAxis1.getAxis(),
                                        std::move(materialAccessor), payload));
}

std::unique_ptr<GloballyIndexedGridSurfaceMaterial>
GridSurfaceMaterialFactory::create(
    const ProtoAxis& pAxis0, const ProtoAxis& pAxis1,
    GloballyIndexedMaterialAccessor&& materialAccessor,
    const std::vector<
        std::vector<GloballyIndexedMaterialAccessor::grid_value_type>>&
        payload) {
  return std::make_unique<GloballyIndexedGridSurfaceMaterial>(
      create2D<GloballyIndexedMaterialAccessor>(
          pAxis0.getAxis(), pAxis1.getAxis(), std::move(materialAccessor),
          payload));
}

ACTS_POP_IGNORE_DEPRECATED()

}  // namespace Acts
