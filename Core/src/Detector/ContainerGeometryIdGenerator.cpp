// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/ContainerGeometryIdGenerator.hpp"

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Surfaces/Surface.hpp"

IGeometryIdGenerator::GeoIdCache
Acts::Experimental::ContainerGeometryIdGenerator::generateCache() const {
  return Cache{};
}

void Acts::Experimental::ContainerGeometryIdGenerator::assignGeometryId(
    IGeometryIdGenerator::GeoIdCache& cache, DetectorVolume& dVolume) const {
  auto ccache = std::any_cast<Cache&>(cache);
}

void Acts::Experimental::ContainerGeometryIdGenerator::assignGeometryId(
    IGeometryIdGenerator::GeoIdCache& cache, Portal& portal) const {
  auto ccache = std::any_cast<Cache&>(cache);
}

void Acts::Experimental::ContainerGeometryIdGenerator::assignGeometryId(
    IGeometryIdGenerator::GeoIdCache& cache, Surface& surface) const {
  auto ccache = std::any_cast<Cache&>(cache);
}
