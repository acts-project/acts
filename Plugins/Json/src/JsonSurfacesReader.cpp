// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/JsonSurfacesReader.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"
#include "ActsPlugins/Json/GeometryIdentifierJsonConverter.hpp"
#include "ActsPlugins/Json/SurfaceJsonConverter.hpp"

#include <fstream>
#include <iostream>

namespace Acts {

Acts::GeometryHierarchyMap<std::shared_ptr<Acts::Surface>>
JsonSurfacesReader::readHierarchyMap(
    const JsonSurfacesReader::Options& options) {
  // Read the json file into a json object
  nlohmann::json j;
  std::ifstream in(options.inputFile);
  in >> j;
  in.close();

  using SurfaceHierachyMap =
      Acts::GeometryHierarchyMap<std::shared_ptr<Acts::Surface>>;
  std::vector<SurfaceHierachyMap::InputElement> surfaceElements;

  // Walk down the path to the surface entries
  nlohmann::json jSurfaces = j;
  for (const auto& jep : options.jsonEntryPath) {
    jSurfaces = jSurfaces[jep];
  }

  // Loop over the surfaces
  surfaceElements.reserve(jSurfaces.size());
  for (const auto& jSurface : jSurfaces) {
    // Decode the surface identifier
    GeometryIdentifier geoId =
        GeometryIdentifierJsonConverter::decodeIdentifier(jSurface);
    auto surface = Acts::SurfaceJsonConverter::fromJson(jSurface["value"]);
    surfaceElements.emplace_back(geoId, surface);
  }
  return SurfaceHierachyMap(std::move(surfaceElements));
}

std::vector<std::shared_ptr<Acts::Surface>> JsonSurfacesReader::readVector(
    const Options& options) {
  // Read the json file into a json object
  nlohmann::json j;
  std::ifstream in(options.inputFile);
  in >> j;
  in.close();

  // Walk down the path to the surface entries
  nlohmann::json jSurfaces = j;
  for (const auto& jep : options.jsonEntryPath) {
    jSurfaces = jSurfaces[jep];
  }

  std::vector<std::shared_ptr<Acts::Surface>> surfaces;
  for (const auto& jSurface : jSurfaces) {
    auto surface = Acts::SurfaceJsonConverter::fromJson(jSurface);
    surfaces.push_back(surface);
  }
  return surfaces;
}

std::vector<std::shared_ptr<Acts::JsonDetectorElement>>
JsonSurfacesReader::readDetectorElements(const Options& options,
                                         double thickness = 0.0) {
  nlohmann::json j;
  {
    std::ifstream in(options.inputFile);
    in >> j;
  }

  // Walk down the path to the surface entries
  nlohmann::json jSurfaces = j;
  for (const auto& jep : options.jsonEntryPath) {
    jSurfaces = jSurfaces[jep];
  }

  std::vector<std::shared_ptr<Acts::JsonDetectorElement>> elements;
  for (const auto& jSurface : jSurfaces) {
    elements.emplace_back(
        std::make_shared<Acts::JsonDetectorElement>(jSurface, thickness));
  }
  return elements;
}

}  // namespace Acts
