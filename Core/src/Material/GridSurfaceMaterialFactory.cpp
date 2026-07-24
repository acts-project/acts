// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/GridSurfaceMaterialFactory.hpp"

#include "Acts/Utilities/Diagnostics.hpp"

std::unique_ptr<Acts::IGridSurfaceMaterial<Acts::MaterialSlab>>
Acts::GridSurfaceMaterialFactory::create(
    const IAxis& axis, GridMaterialAccessor&& materialAccessor,
    GridAccess::BoundToGridLocal1DimDelegate boundToGridLocal,
    GridAccess::GlobalToGridLocal1DimDelegate globalToGridLocal,
    const std::vector<MaterialSlab>& payload) {
  return create1D<GridMaterialAccessor>(axis, std::move(materialAccessor),
                                        std::move(boundToGridLocal),
                                        std::move(globalToGridLocal), payload);
}

std::unique_ptr<
    Acts::IGridSurfaceMaterial<Acts::IndexedMaterialAccessor::grid_value_type>>
Acts::GridSurfaceMaterialFactory::create(
    const IAxis& axis, IndexedMaterialAccessor&& materialAccessor,
    GridAccess::BoundToGridLocal1DimDelegate boundToGridLocal,
    GridAccess::GlobalToGridLocal1DimDelegate globalToGridLocal,
    const std::vector<IndexedMaterialAccessor::grid_value_type>& payload) {
  return create1D<IndexedMaterialAccessor>(
      axis, std::move(materialAccessor), std::move(boundToGridLocal),
      std::move(globalToGridLocal), payload);
}

std::unique_ptr<Acts::IGridSurfaceMaterial<
    Acts::GloballyIndexedMaterialAccessor::grid_value_type>>
Acts::GridSurfaceMaterialFactory::create(
    const IAxis& axis, GloballyIndexedMaterialAccessor&& materialAccessor,
    GridAccess::BoundToGridLocal1DimDelegate boundToGridLocal,
    GridAccess::GlobalToGridLocal1DimDelegate globalToGridLocal,
    const std::vector<GloballyIndexedMaterialAccessor::grid_value_type>&
        payload) {
  return create1D<GloballyIndexedMaterialAccessor>(
      axis, std::move(materialAccessor), std::move(boundToGridLocal),
      std::move(globalToGridLocal), payload);
}

std::unique_ptr<Acts::IGridSurfaceMaterial<Acts::MaterialSlab>>
Acts::GridSurfaceMaterialFactory::create(
    const IAxis& axis0, const IAxis& axis1,
    GridMaterialAccessor&& materialAccessor,
    GridAccess::BoundToGridLocal2DimDelegate boundToGridLocal,
    GridAccess::GlobalToGridLocal2DimDelegate globalToGridLocal,
    const std::vector<std::vector<MaterialSlab>>& payload) {
  return create2D<GridMaterialAccessor>(
      axis0, axis1, std::move(materialAccessor), std::move(boundToGridLocal),
      std::move(globalToGridLocal), payload);
}

std::unique_ptr<
    Acts::IGridSurfaceMaterial<Acts::IndexedMaterialAccessor::grid_value_type>>
Acts::GridSurfaceMaterialFactory::create(
    const IAxis& axis0, const IAxis& axis1,
    IndexedMaterialAccessor&& materialAccessor,
    GridAccess::BoundToGridLocal2DimDelegate boundToGridLocal,
    GridAccess::GlobalToGridLocal2DimDelegate globalToGridLocal,
    const std::vector<std::vector<IndexedMaterialAccessor::grid_value_type>>&
        payload) {
  return create2D<IndexedMaterialAccessor>(
      axis0, axis1, std::move(materialAccessor), std::move(boundToGridLocal),
      std::move(globalToGridLocal), payload);
}

std::unique_ptr<Acts::IGridSurfaceMaterial<
    Acts::GloballyIndexedMaterialAccessor::grid_value_type>>
Acts::GridSurfaceMaterialFactory::create(
    const IAxis& axis0, const IAxis& axis1,
    GloballyIndexedMaterialAccessor&& materialAccessor,
    GridAccess::BoundToGridLocal2DimDelegate boundToGridLocal,
    GridAccess::GlobalToGridLocal2DimDelegate globalToGridLocal,
    const std::vector<
        std::vector<GloballyIndexedMaterialAccessor::grid_value_type>>&
        payload) {
  return create2D<GloballyIndexedMaterialAccessor>(
      axis0, axis1, std::move(materialAccessor), std::move(boundToGridLocal),
      std::move(globalToGridLocal), payload);
}

ACTS_PUSH_IGNORE_DEPRECATED()

std::unique_ptr<Acts::IGridSurfaceMaterial<Acts::MaterialSlab>>
Acts::GridSurfaceMaterialFactory::create(
    const ProtoAxis& pAxis, GridMaterialAccessor&& materialAccessor,
    GridAccess::BoundToGridLocal1DimDelegate boundToGridLocal,
    GridAccess::GlobalToGridLocal1DimDelegate globalToGridLocal,
    const std::vector<MaterialSlab>& payload) {
  return create1D<GridMaterialAccessor>(
      pAxis.getAxis(), std::move(materialAccessor), std::move(boundToGridLocal),
      std::move(globalToGridLocal), payload);
}

std::unique_ptr<
    Acts::IGridSurfaceMaterial<Acts::IndexedMaterialAccessor::grid_value_type>>
Acts::GridSurfaceMaterialFactory::create(
    const ProtoAxis& pAxis, IndexedMaterialAccessor&& materialAccessor,
    GridAccess::BoundToGridLocal1DimDelegate boundToGridLocal,
    GridAccess::GlobalToGridLocal1DimDelegate globalToGridLocal,
    const std::vector<IndexedMaterialAccessor::grid_value_type>& payload) {
  return create1D<IndexedMaterialAccessor>(
      pAxis.getAxis(), std::move(materialAccessor), std::move(boundToGridLocal),
      std::move(globalToGridLocal), payload);
}

std::unique_ptr<Acts::IGridSurfaceMaterial<
    Acts::GloballyIndexedMaterialAccessor::grid_value_type>>
Acts::GridSurfaceMaterialFactory::create(
    const ProtoAxis& pAxis, GloballyIndexedMaterialAccessor&& materialAccessor,
    GridAccess::BoundToGridLocal1DimDelegate boundToGridLocal,
    GridAccess::GlobalToGridLocal1DimDelegate globalToGridLocal,
    const std::vector<GloballyIndexedMaterialAccessor::grid_value_type>&
        payload) {
  return create1D<GloballyIndexedMaterialAccessor>(
      pAxis.getAxis(), std::move(materialAccessor), std::move(boundToGridLocal),
      std::move(globalToGridLocal), payload);
}

std::unique_ptr<Acts::IGridSurfaceMaterial<Acts::MaterialSlab>>
Acts::GridSurfaceMaterialFactory::create(
    const ProtoAxis& pAxis0, const ProtoAxis& pAxis1,
    GridMaterialAccessor&& materialAccessor,
    GridAccess::BoundToGridLocal2DimDelegate boundToGridLocal,
    GridAccess::GlobalToGridLocal2DimDelegate globalToGridLocal,
    const std::vector<std::vector<MaterialSlab>>& payload) {
  return create2D<GridMaterialAccessor>(
      pAxis0.getAxis(), pAxis1.getAxis(), std::move(materialAccessor),
      std::move(boundToGridLocal), std::move(globalToGridLocal), payload);
}

std::unique_ptr<
    Acts::IGridSurfaceMaterial<Acts::IndexedMaterialAccessor::grid_value_type>>
Acts::GridSurfaceMaterialFactory::create(
    const ProtoAxis& pAxis0, const ProtoAxis& pAxis1,
    IndexedMaterialAccessor&& materialAccessor,
    GridAccess::BoundToGridLocal2DimDelegate boundToGridLocal,
    GridAccess::GlobalToGridLocal2DimDelegate globalToGridLocal,
    const std::vector<std::vector<IndexedMaterialAccessor::grid_value_type>>&
        payload) {
  return create2D<IndexedMaterialAccessor>(
      pAxis0.getAxis(), pAxis1.getAxis(), std::move(materialAccessor),
      std::move(boundToGridLocal), std::move(globalToGridLocal), payload);
}

std::unique_ptr<Acts::IGridSurfaceMaterial<
    Acts::GloballyIndexedMaterialAccessor::grid_value_type>>
Acts::GridSurfaceMaterialFactory::create(
    const ProtoAxis& pAxis0, const ProtoAxis& pAxis1,
    GloballyIndexedMaterialAccessor&& materialAccessor,
    GridAccess::BoundToGridLocal2DimDelegate boundToGridLocal,
    GridAccess::GlobalToGridLocal2DimDelegate globalToGridLocal,
    const std::vector<
        std::vector<GloballyIndexedMaterialAccessor::grid_value_type>>&
        payload) {
  return create2D<GloballyIndexedMaterialAccessor>(
      pAxis0.getAxis(), pAxis1.getAxis(), std::move(materialAccessor),
      std::move(boundToGridLocal), std::move(globalToGridLocal), payload);
}

ACTS_POP_IGNORE_DEPRECATED()
