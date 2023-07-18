// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/DetectorComponentsExtractors.hpp"

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Navigation/NavigationState.hpp"

#include <stdexcept>

const std::vector<const Acts::Experimental::Portal*>
Acts::Experimental::AllPortalsExtractor::extract(
    [[maybe_unused]] const GeometryContext& gctx,
    const NavigationState& nState) {
  if (nState.currentVolume == nullptr) {
    throw std::runtime_error("AllPortalsExtractor: no detector volume given.");
  }
  return nState.currentVolume->portals();
}

const std::vector<const Acts::Surface*>
Acts::Experimental::AllSurfacesExtractor::extract(
    [[maybe_unused]] const GeometryContext& gctx, const NavigationState& nState,
    [[maybe_unused]] const std::vector<size_t>& indices) {
  if (nState.currentVolume == nullptr) {
    throw std::runtime_error("AllSurfacesExtractor: no detector volume given.");
  }
  return nState.currentVolume->surfaces();
}

const std::vector<const Acts::Surface*>
Acts::Experimental::IndexedSurfacesExtractor::extract(
    [[maybe_unused]] const GeometryContext& gctx, const NavigationState& nState,
    const std::vector<size_t>& indices) {
  if (nState.currentVolume == nullptr) {
    throw std::runtime_error(
        "IndexedSurfacesExtractor: no detector volume given.");
  }
  // Get the surface container
  const auto& surfaces = nState.currentVolume->surfaces();
  // The extracted surfaces
  std::vector<const Surface*> eSurfaces;
  eSurfaces.reserve(indices.size());
  std::for_each(indices.begin(), indices.end(),
                [&](const auto& i) { eSurfaces.push_back(surfaces[i]); });
  return eSurfaces;
}

const std::vector<const Acts::Experimental::DetectorVolume*>
Acts::Experimental::AllSubVolumesExtractor::extract(
    [[maybe_unused]] const GeometryContext& gctx, const NavigationState& nState,
    [[maybe_unused]] const std::vector<size_t>& indices) {
  if (nState.currentVolume == nullptr) {
    throw std::runtime_error(
        "AllSubVolumesExtractor: no detector volume given.");
  }
  return nState.currentVolume->volumes();
}

const std::vector<const Acts::Experimental::DetectorVolume*>
Acts::Experimental::IndexedSubVolumesExtractor::extract(
    [[maybe_unused]] const GeometryContext& gctx, const NavigationState& nState,
    const std::vector<size_t>& indices) {
  if (nState.currentVolume == nullptr) {
    throw std::runtime_error(
        "AllSubVolumesExtractor: no detector volume given.");
  }
  // Get the sub volumes container
  const auto& volumes = nState.currentVolume->volumes();
  // The extracted volumes
  std::vector<const DetectorVolume*> eVolumes;
  eVolumes.reserve(indices.size());
  std::for_each(indices.begin(), indices.end(),
                [&](const auto& i) { eVolumes.push_back(volumes[i]); });
  return eVolumes;
}
