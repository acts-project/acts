// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/GridSurfaceMaterialAccumulator.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/GridSurfaceMaterialFactory.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"

#include <set>
#include <sstream>
#include <stdexcept>

Acts::GridSurfaceMaterialAccumulator::GridSurfaceMaterialAccumulator(
    const Config& cfg, std::unique_ptr<const Logger> mlogger)
    : m_cfg(cfg), m_logger(std::move(mlogger)) {}

std::unique_ptr<Acts::ISurfaceMaterialAccumulator::State>
Acts::GridSurfaceMaterialAccumulator::createState(
    const GeometryContext& /*gctx*/) const {
  auto state = std::make_unique<State>();

  for (const auto& surface : m_cfg.materialSurfaces) {
    GeometryIdentifier geoID = surface->geometryId();

    const ISurfaceMaterial* surfaceMaterial = surface->surfaceMaterial();
    if (surfaceMaterial == nullptr) {
      throw std::invalid_argument(
          "Surface material is not set, inconsistent configuration.");
    }

    // Legacy BinUtility-based proto/binned material is not supported here -
    // BinnedSurfaceMaterialAccumulator remains the accumulator for that
    if (dynamic_cast<const ProtoSurfaceMaterial*>(surfaceMaterial) != nullptr ||
        dynamic_cast<const BinnedSurfaceMaterial*>(surfaceMaterial) !=
            nullptr) {
      throw std::invalid_argument(
          "GridSurfaceMaterialAccumulator: legacy BinUtility-based "
          "ProtoSurfaceMaterial/BinnedSurfaceMaterial is not supported, use "
          "ProtoGridSurfaceMaterial (or BinnedSurfaceMaterialAccumulator for "
          "the legacy pipeline) instead.");
    }

    // First attempt: (possibly deferred) ProtoGridSurfaceMaterial
    auto psgm = dynamic_cast<const ProtoGridSurfaceMaterial*>(surfaceMaterial);
    if (psgm != nullptr) {
      ACTS_DEBUG("       - (proto) binning from ProtoGridSurfaceMaterial is "
                 << psgm->binning());
      auto multiAxis = resolveMultiAxis(psgm->binning(), *surface);
      ACTS_DEBUG("       - resolved binning has "
                 << multiAxis->getNTotalBins(true) << " bins (incl. "
                 << "under-/overflow)");
      std::size_t nBins = multiAxis->getNTotalBins(true);
      State::GridAccumulation grid{std::move(multiAxis),
                                   std::vector<AccumulatedMaterialSlab>(nBins)};
      state->gridMaterial.emplace(geoID, std::move(grid));
      continue;
    }

    // Second attempt: an already concrete GridSurfaceMaterial, e.g. for
    // re-accumulation/refinement
    auto gsm = dynamic_cast<const GridSurfaceMaterial*>(surfaceMaterial);
    if (gsm != nullptr) {
      ACTS_DEBUG("       - binning from GridSurfaceMaterial is "
                 << gsm->binning());
      auto multiAxis = gsm->binning().buildMultiAxis();
      std::size_t nBins = multiAxis->getNTotalBins(true);
      State::GridAccumulation grid{std::move(multiAxis),
                                   std::vector<AccumulatedMaterialSlab>(nBins)};
      state->gridMaterial.emplace(geoID, std::move(grid));
      continue;
    }

    // Fall back: homogeneous material
    ACTS_DEBUG("       - this is homogeneous material.");
    state->homogeneousMaterial.emplace(geoID, AccumulatedMaterialSlab());
  }
  return state;
}

Acts::Vector2 Acts::GridSurfaceMaterialAccumulator::resolveLocalPosition(
    const GeometryContext& gctx, const Surface& surface,
    const Vector3& position, const Vector3& direction) {
  auto lpResult = surface.globalToLocal(gctx, position, direction);
  if (!lpResult.ok()) {
    std::stringstream ss;
    ss << "GridSurfaceMaterialAccumulator: failed to resolve local position "
          "on surface "
       << surface.geometryId();
    throw std::invalid_argument(ss.str());
  }
  return lpResult.value();
}

void Acts::GridSurfaceMaterialAccumulator::accumulate(
    ISurfaceMaterialAccumulator::State& state, const GeometryContext& gctx,
    const std::vector<MaterialInteraction>& interactions,
    const std::vector<IAssignmentFinder::SurfaceAssignment>&
        surfacesWithoutAssignment) const {
  // Cast into the right state object (guaranteed by upstream algorithm)
  State* cState = static_cast<State*>(&state);
  if (cState == nullptr) {
    throw std::invalid_argument(
        "Invalid state object provided, something is seriously wrong.");
  }

  // Every touched bin is tracked individually, so a track that touches
  // several distinct bins on the same surface still contributes to all of
  // them (unlike a per-surface touched-bin cache, which only remembers one).
  std::set<AccumulatedMaterialSlab*> touchedSlabs;

  // Assign the hits
  for (const auto& mi : interactions) {
    const Surface* surface = mi.surface;
    GeometryIdentifier geoID = surface->geometryId();

    auto gridIt = cState->gridMaterial.find(geoID);
    if (gridIt != cState->gridMaterial.end()) {
      auto& grid = gridIt->second;
      Vector2 lp =
          resolveLocalPosition(gctx, *surface, mi.intersection, mi.direction);
      std::size_t bin = grid.multiAxis->getGlobalBinFromPoint({lp[0], lp[1]});
      AccumulatedMaterialSlab& slab = grid.accumulatedMaterial.at(bin);
      slab.accumulate(mi.materialSlab, mi.pathCorrection);
      touchedSlabs.insert(&slab);
      continue;
    }

    auto homIt = cState->homogeneousMaterial.find(geoID);
    if (homIt != cState->homogeneousMaterial.end()) {
      homIt->second.accumulate(mi.materialSlab, mi.pathCorrection);
      touchedSlabs.insert(&(homIt->second));
      continue;
    }

    throw std::invalid_argument(
        "Surface material is not found, inconsistent configuration.");
  }

  // After mapping this track, average the touched bins
  for (auto* slab : touchedSlabs) {
    slab->trackAverage(true);
  }

  // Empty bin correction
  if (m_cfg.emptyBinCorrection) {
    for (const auto& [surface, position, direction] :
         surfacesWithoutAssignment) {
      GeometryIdentifier geoID = surface->geometryId();

      auto gridIt = cState->gridMaterial.find(geoID);
      if (gridIt != cState->gridMaterial.end()) {
        auto& grid = gridIt->second;
        Vector2 lp = resolveLocalPosition(gctx, *surface, position, direction);
        std::size_t bin = grid.multiAxis->getGlobalBinFromPoint({lp[0], lp[1]});
        grid.accumulatedMaterial.at(bin).trackAverage(true);
        continue;
      }

      auto homIt = cState->homogeneousMaterial.find(geoID);
      if (homIt != cState->homogeneousMaterial.end()) {
        homIt->second.trackAverage(true);
        continue;
      }

      throw std::invalid_argument(
          "Surface material is not found, inconsistent configuration.");
    }
  }
}

std::map<Acts::GeometryIdentifier,
         std::shared_ptr<const Acts::ISurfaceMaterial>>
Acts::GridSurfaceMaterialAccumulator::finalizeMaterial(
    ISurfaceMaterialAccumulator::State& state,
    const GeometryContext& /*gctx*/) const {
  std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>
      sMaterials;

  // Cast into the right state object (guaranteed by upstream algorithm)
  State* cState = static_cast<State*>(&state);
  if (cState == nullptr) {
    throw std::invalid_argument(
        "Invalid state object provided, something is seriously wrong.");
  }

  for (auto& [geoID, grid] : cState->gridMaterial) {
    ACTS_DEBUG("Finalizing grid map for Surface " << geoID);
    const IAxis& axis0 = grid.multiAxis->getAxis(0);
    const IAxis& axis1 = grid.multiAxis->getAxis(1);
    IMultiAxis2D::LocalBins nBins = grid.multiAxis->getNBins();

    // Regular-bin-only payload, matching GridSurfaceMaterialFactory's
    // convention - under-/overflow bins are not persisted/reconstructed
    std::vector<std::vector<MaterialSlab>> payload(
        nBins[0], std::vector<MaterialSlab>(nBins[1]));
    for (std::size_t i0 = 0; i0 < nBins[0]; ++i0) {
      for (std::size_t i1 = 0; i1 < nBins[1]; ++i1) {
        IMultiAxis2D::LocalBins lbin{i0 + 1, i1 + 1};
        std::size_t bin = grid.multiAxis->getGlobalBinFromLocalBins(lbin);
        payload[i0][i1] = grid.accumulatedMaterial.at(bin).totalAverage().first;
      }
    }

    if (m_cfg.storageKind == StorageKind::Direct) {
      sMaterials[geoID] =
          GridSurfaceMaterialFactory::createDirect(axis0, axis1, payload);
      continue;
    }

    // Indexed: one (unique) entry per bin, in bin order - no merging of
    // equal/similar slabs here, that is left to a postprocessing step, same
    // as GloballyIndexed storage is
    std::vector<MaterialSlab> material;
    material.reserve(nBins[0] * nBins[1]);
    std::vector<std::vector<std::size_t>> indices(
        nBins[0], std::vector<std::size_t>(nBins[1]));
    for (std::size_t i0 = 0; i0 < nBins[0]; ++i0) {
      for (std::size_t i1 = 0; i1 < nBins[1]; ++i1) {
        indices[i0][i1] = material.size();
        material.push_back(payload[i0][i1]);
      }
    }
    sMaterials[geoID] = GridSurfaceMaterialFactory::createIndexed(
        axis0, axis1, std::move(material), indices);
  }

  for (auto& [geoID, slab] : cState->homogeneousMaterial) {
    ACTS_DEBUG("Finalizing homogeneous map for Surface " << geoID);
    sMaterials[geoID] = std::make_shared<const HomogeneousSurfaceMaterial>(
        slab.totalAverage().first);
  }

  return sMaterials;
}
