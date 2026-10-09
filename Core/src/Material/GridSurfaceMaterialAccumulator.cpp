// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/GridSurfaceMaterialAccumulator.hpp"

#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AxisSpec.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"

#include <array>
#include <set>
#include <stdexcept>
#include <utility>

namespace Acts {

namespace {

GridSurfaceMaterialAccumulator::State& checkedState(
    ISurfaceMaterialAccumulator::State& state) {
  auto* concrete = dynamic_cast<GridSurfaceMaterialAccumulator::State*>(&state);
  if (concrete == nullptr || !concrete->materialSurfaceRegistry) {
    throw std::invalid_argument(
        "State was not created by a grid surface material accumulator");
  }
  return *concrete;
}

AccumulatedMaterialSlab& findBin(GridSurfaceMaterialAccumulator::State& state,
                                 const GeometryContext& gctx,
                                 const Surface* surface,
                                 const Vector3& position,
                                 const Vector3& direction) {
  if (surface == nullptr) {
    throw std::invalid_argument("Null material assignment surface");
  }
  auto found = state.accumulatedMaterial.find(surface->geometryId());
  if (found == state.accumulatedMaterial.end() ||
      found->second.surface != surface) {
    throw std::invalid_argument(
        "Surface material is not found, inconsistent configuration.");
  }
  auto local = surface->globalToLocal(gctx, position, direction);
  if (!local.ok()) {
    throw std::invalid_argument(
        "Material assignment cannot be converted to surface-local coordinates");
  }
  auto& grid = found->second;
  Vector2 point = local.value();
  if (grid.swapAxes) {
    std::swap(point[0], point[1]);
  }
  auto bin = grid.multiAxis->getGlobalBinFromPoint({point[0], point[1]});
  return grid.material.at(bin);
}

MultiAxisSpec2D binningFromAxes(const IMultiAxis2D& axes) {
  return MultiAxisSpec2D(std::array{AxisSpec::FromAxis(axes.getAxis(0)),
                                    AxisSpec::FromAxis(axes.getAxis(1))});
}

}  // namespace

GridSurfaceMaterialAccumulator::GridSurfaceMaterialAccumulator(
    const Config& cfg, std::unique_ptr<const Logger> mlogger)
    : m_cfg(cfg), m_logger(std::move(mlogger)) {}

std::unique_ptr<ISurfaceMaterialAccumulator::State>
GridSurfaceMaterialAccumulator::createState(
    const GeometryContext& /*gctx*/) const {
  auto state = std::make_unique<State>();
  state->materialSurfaceRegistry.emplace(m_cfg.materialSurfaces);
  for (const auto& [id, surface] : state->materialSurfaceRegistry->surfaces()) {
    const ISurfaceMaterial* material = surface->surfaceMaterial();
    std::unique_ptr<IMultiAxis2D> axes;
    bool swapAxes = false;
    if (const auto* proto =
            dynamic_cast<const ProtoSurfaceMaterial*>(material)) {
      axes = resolveMultiAxis(proto->binning(), *surface);
    } else if (const auto* grid =
                   dynamic_cast<const GridSurfaceMaterial*>(material)) {
      axes = grid->binning().buildMultiAxis();
      const auto directions = grid->localAxisDirections();
      swapAxes =
          !directions.empty() && directions[0] != surface->localAxes()[0];
    } else if (dynamic_cast<const HomogeneousSurfaceMaterial*>(material) !=
               nullptr) {
      axes = resolveMultiAxis(ProtoSurfaceMaterial{}.binning(), *surface);
    } else {
      throw std::invalid_argument(
          "Grid surface material accumulation requires prototype, grid or "
          "homogeneous surface material");
    }
    ACTS_DEBUG("Resolved grid binning for Surface " << id << ": "
                                                    << binningFromAxes(*axes));
    auto nBins = axes->getNTotalBins(true);
    state->accumulatedMaterial.emplace(
        id, State::AccumulatedGrid{surface, std::move(axes),
                                   std::vector<AccumulatedMaterialSlab>(nBins),
                                   swapAxes});
  }
  return state;
}

void GridSurfaceMaterialAccumulator::accumulate(
    ISurfaceMaterialAccumulator::State& state, const GeometryContext& gctx,
    const std::vector<MaterialInteraction>& interactions,
    const std::vector<IAssignmentFinder::SurfaceAssignment>&
        surfacesWithoutAssignment) const {
  auto& concrete = checkedState(state);
  std::set<AccumulatedMaterialSlab*> touchedBins;
  for (const auto& interaction : interactions) {
    auto& bin = findBin(concrete, gctx, interaction.surface,
                        interaction.intersection, interaction.direction);
    bin.accumulate(interaction.materialSlab, interaction.pathCorrection);
    touchedBins.insert(&bin);
  }
  if (m_cfg.emptyBinCorrection) {
    for (const auto& assignment : surfacesWithoutAssignment) {
      touchedBins.insert(&findBin(concrete, gctx, assignment.surface,
                                  assignment.position, assignment.direction));
    }
  }
  for (auto* bin : touchedBins) {
    bin->trackAverage(true);
  }
}

std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>
GridSurfaceMaterialAccumulator::finalizeMaterial(
    ISurfaceMaterialAccumulator::State& state,
    const GeometryContext& /*gctx*/) const {
  auto& concrete = checkedState(state);
  std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>
      materials;
  for (const auto& [id, grid] : concrete.accumulatedMaterial) {
    ACTS_DEBUG("Finalizing grid material for Surface " << id);
    GridSurfaceMaterial::Direct slabs;
    slabs.reserve(grid.material.size());
    for (const auto& bin : grid.material) {
      slabs.push_back(bin.totalAverage().first);
    }
    materials.emplace(
        id, std::make_shared<GridSurfaceMaterial>(
                binningFromAxes(*grid.multiAxis), std::move(slabs), 0.));
  }
  return materials;
}

TrackingGeometryMaterial GridSurfaceMaterialAccumulator::finalizeMaps(
    ISurfaceMaterialAccumulator::State& state,
    const GeometryContext& gctx) const {
  auto& concrete = checkedState(state);
  return concrete.materialSurfaceRegistry->materialMaps(
      finalizeMaterial(state, gctx));
}

}  // namespace Acts
