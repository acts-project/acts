// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/BinnedSurfaceMaterialAccumulater.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinAdjustment.hpp"
#include "Acts/Utilities/BinUtility.hpp"

Acts::BinnedSurfaceMaterialAccumulater::BinnedSurfaceMaterialAccumulater(
    const Config& cfg, std::unique_ptr<const Logger> mlogger)
    : m_cfg(cfg), m_logger(std::move(mlogger)) {}

std::unique_ptr<Acts::ISurfaceMaterialAccumulater::State>
Acts::BinnedSurfaceMaterialAccumulater::createState() const {
  auto state = std::make_unique<State>();

  /// Create the surface accumulation
  for (const auto& surface : m_cfg.materialSurfaces) {
    GeometryIdentifier geoID = surface->geometryId();

    // Get the Surface Material
    const ISurfaceMaterial* surfaceMaterial = surface->surfaceMaterial();
    if (surfaceMaterial == nullptr) {
      throw std::invalid_argument(
          "Surface material is not set, inconsistent configuration.");
    }

    // We need a dynamic_cast to either a surface material proxy or
    // proper surface material
    auto psm = dynamic_cast<const ProtoSurfaceMaterial*>(surfaceMaterial);

    // Get the bin utility: try proxy material first
    const BinUtility* bu = (psm != nullptr) ? (&psm->binning()) : nullptr;
    if (bu != nullptr) {
      // Screen output for Binned Surface material
      ACTS_DEBUG("       - (proto) binning is " << *bu);
      // Now update
      BinUtility buAdjusted = adjustBinUtility(*bu, *surface, m_cfg.geoContext);
      // Screen output for Binned Surface material
      ACTS_DEBUG("       - adjusted binning is " << buAdjusted);
      state->accumulatedMaterial[geoID] =
          AccumulatedSurfaceMaterial(buAdjusted);
      // Material accumulation  is created for this
      continue;
    }

    // Second attempt: binned material
    auto bmp = dynamic_cast<const BinnedSurfaceMaterial*>(surfaceMaterial);
    bu = (bmp != nullptr) ? (&bmp->binUtility()) : nullptr;
    // Creaete a binned type of material
    if (bu != nullptr) {
      // Screen output for Binned Surface material
      ACTS_DEBUG("       - binning is " << *bu);
      state->accumulatedMaterial[geoID] = AccumulatedSurfaceMaterial(*bu);
      // Material accumulation  is created for this
      continue;
    } else {
      // Create a homogeneous type of material
      ACTS_DEBUG("       - this is homogeneous material.");
      state->accumulatedMaterial[geoID] = AccumulatedSurfaceMaterial();
      // Material accumulation  is created for this
      continue;
    }
  }
  return state;
}

void Acts::BinnedSurfaceMaterialAccumulater::accumulate(
    ISurfaceMaterialAccumulater::State& state,
    const std::vector<MaterialInteraction>& interactions,
    const std::vector<SurfaceAssignment>& surfacesWithoutAssignment) const {
  // Cast into the right state object (guaranteed by upstream algroithm)
  State* cState = static_cast<State*>(&state);
  if (cState == nullptr) {
    throw std::invalid_argument(
        "Invalid state object provided, something is seriously wrong.");
  }

  using MapBin =
      std::pair<AccumulatedSurfaceMaterial*, std::array<std::size_t, 3>>;
  std::map<AccumulatedSurfaceMaterial*, std::array<std::size_t, 3>>
      touchedMapBins;

  // Assign the hits
  for (const auto& mi : interactions) {
    // Get the surface
    const Surface* surface = mi.surface;
    GeometryIdentifier geoID = surface->geometryId();
    // Get the accumulated material
    auto accMaterial = cState->accumulatedMaterial.find(geoID);
    if (accMaterial == cState->accumulatedMaterial.end()) {
      throw std::invalid_argument(
          "Surface material is not found, inconsistent configuration.");
    }
    // Accumulate the material - remember the touched bin
    auto tBin = accMaterial->second.accumulate(mi.position, mi.materialSlab,
                                               mi.pathCorrection);
    touchedMapBins.insert(MapBin(&(accMaterial->second), tBin));
  }

  // After mapping this track, average the touched bins
  for (auto tmapBin : touchedMapBins) {
    std::vector<std::array<std::size_t, 3>> trackBins = {tmapBin.second};
    tmapBin.first->trackAverage(trackBins, true);
  }

  // Empty bin correction
  if (m_cfg.emptyBinCorrection) {
    for (auto [surface, position, direction] : surfacesWithoutAssignment) {
      // Get the accumulated material
      auto missedMaterial =
          cState->accumulatedMaterial.find(surface->geometryId());
      if (missedMaterial == cState->accumulatedMaterial.end()) {
        throw std::invalid_argument(
            "Surface material is not found, inconsistent configuration.");
      }
      // Apply empty hit correction
      missedMaterial->second.trackAverage(position, true);
    }
  }
}

std::map<Acts::GeometryIdentifier,
         std::shared_ptr<const Acts::ISurfaceMaterial>>
Acts::BinnedSurfaceMaterialAccumulater::finalizeMaterial(
    ISurfaceMaterialAccumulater::State& state) const {
  std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>
      sMaterials;

  // Cast into the right state object (guaranteed by upstream algroithm)
  State* cState = static_cast<State*>(&state);
  if (cState == nullptr) {
    throw std::invalid_argument(
        "Invalid state object provided, something is seriously wrong.");
  }

  // iterate over the map to call the total average
  for (auto& accMaterial : cState->accumulatedMaterial) {
    ACTS_DEBUG("Finalizing map for Surface " << accMaterial.first);
    auto sMaterial = accMaterial.second.totalAverage();
    sMaterials[accMaterial.first] = std::move(sMaterial);
  }

  return sMaterials;
}
