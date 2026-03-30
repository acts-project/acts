// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "ActsPlugins/Json/JsonDetectorElement.hpp"

#include <memory>
#include <string>
#include <vector>

namespace Acts {

class Surface;

/// @addtogroup json_plugin
/// @{

namespace JsonSurfacesReader {

/// @brief Options specification for surface reading
/// The file should contain an array of json surfaces
/// as produced by the SurfaceJsonConverter tools
struct Options {
  /// @brief  Which input file to read from
  std::string inputFile = "";
  /// The entry path until you reach the surface entries
  std::vector<std::string> jsonEntryPath = {"Surfaces", "entries"};
};

/// @brief Read the surfaces from the input file
/// For details on the file format see the options struct
///
/// @param options specifies which file and what to read
///
/// @return  a vector of surfaces
Acts::GeometryHierarchyMap<std::shared_ptr<Acts::Surface>> readHierarchyMap(
    const Options& options);

/// @brief Read the flat surfaces from the input file
/// For details on the file format see the options struct
///
/// @param options options for surface reading
///
/// @return  a vector of surfaces
std::vector<std::shared_ptr<Acts::Surface>> readVector(const Options& options);

/// @brief Read the surfaces from the input file and create
/// detector elements
/// For details on the file format see the options struct
///
/// @param options options for surface reading
/// @param thickness the thickness used to construct the detector element
///
/// @return  a vector of surfaces
std::vector<std::shared_ptr<Acts::JsonDetectorElement>> readDetectorElements(
    const Options& options, double thickness);

}  // namespace JsonSurfacesReader

/// @}

}  // namespace Acts
