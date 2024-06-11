// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Plugins/Json/DetectorVolumeJsonConverter.hpp"

// Custom Json encoder/decoders

namespace Acts {

namespace Experimental {
class Detector;
}

namespace DetrayJsonConverter {

struct Options {
  DetectorVolumeJsonConverter::Options volumeOptions =
      DetectorVolumeJsonConverter::Options{};
};

/// Convert an axis to json - detray style
///
/// @param ia the axis
///
/// @return a json object to represent the axis in detray json
nlohmann::json toJsonDetray(const IAxis& ia);

/// @brief Templated grid conversion to json
///
/// @tparam grid_type the type of the grid
/// @param grid the grid object
/// @param swapAxis - this is needed for detray
///
/// @note detray has a different offset for the
/// local indices, it starts at 0
///
/// @return a json object to represent the grid
template <typename grid_type>
nlohmann::json toJsonDetray(const grid_type& grid, bool swapAxis = false) {
  nlohmann::json jGrid;
  // Get the grid axes & potentially swap them
  std::array<const Acts::IAxis*, grid_type::DIM> axes = grid.axes();
  if (swapAxis && grid_type::DIM == 2u) {
    std::swap(axes[0u], axes[1u]);
  }

  nlohmann::json jAxes;

  // Fill the axes in the order they are
  for (unsigned int ia = 0u; ia < grid_type::DIM; ++ia) {
    auto jAxis = AxisJsonConverter::toJsonDetray(*axes[ia]);
    jAxis["label"] = ia;
    jAxes.push_back(jAxis);
  }
  jGrid["axes"] = jAxes;

  nlohmann::json jData;
  // 1D connections
  if constexpr (grid_type::DIM == 1u) {
    for (unsigned int ib0 = 1u; ib0 <= axes[0u]->getNBins(); ++ib0) {
      // Lookup bin
      typename grid_type::index_t lbin;
      lbin[0u] = ib0;
      nlohmann::json jBin;
      jBin["content"] = grid.atLocalBins(lbin);
      // Corrected bin for detray
      lbin[0u] = ib0 - 1u;
      jBin["loc_index"] = lbin;
      jData.push_back(jBin);
    }
  }

  // 2D connections
  if constexpr (grid_type::DIM == 2u) {
    for (unsigned int ib0 = 1u; ib0 <= axes[0u]->getNBins(); ++ib0) {
      for (unsigned int ib1 = 1u; ib1 <= axes[1u]->getNBins(); ++ib1) {
        typename grid_type::index_t lbin;
        // Lookup bin - respect swap (if it happened) for the lookup
        lbin[0u] = swapAxis ? ib1 : ib0;
        lbin[1u] = swapAxis ? ib0 : ib1;
        nlohmann::json jBin;
        jBin["content"] = grid.atLocalBins(lbin);
        // Corrected bin for detray
        lbin[0u] = ib0 - 1u;
        lbin[1u] = ib1 - 1u;
        jBin["loc_index"] = lbin;
        jData.push_back(jBin);
      }
    }
  }
  jGrid["bins"] = jData;

  return jGrid;
}


/// Interface with detray conversion option
/// @param bounds is the bounds object
/// @param portal is the flag for conversion into detray portal format
///
/// @note reading back detray json format to Acts is not supported
///
/// @return the json object
nlohmann::json toJsonDetray(const SurfaceBounds& bounds, bool portal = false);


/// Contextual conversion of a surface - Detray export
///
/// @param gctx the geometry context for this
/// @param surface the surface to be converted
/// @param options the writing options for the surfaces
///
/// @note reading back detray json is not supported and will fail
///
/// @return a json object representing the surface
nlohmann::json toJsonDetray(const GeometryContext& gctx, const Surface& surface,
                            const Options& options = Options{});
                            
/// @brief Convert to detray json format
///
/// @param gctx the geometry context
/// @param detector the detector instance
/// @param options the writing options that propagate
///        to the downstream converters
///
/// @return a json object in detray format
nlohmann::json toJsonDetray(const GeometryContext& gctx,
                            const Experimental::Detector& detector,
                            const Options& options = Options{});


/// @brief Convert to json format - dedicated Detray function
///
/// @param gctx the geometry context
/// @param portal the portal instance
/// @param ip is the portal index that could be used to pick the oriented surface
/// @param volume is the detector volume to which these surfaces belong to
/// @param orientedSurfaces are the bounding surfaces (may still need to be split)
/// @param detectorVolumes is the list of all detector voluems for portal pointing
/// @param options how to write this thing out
///
/// @note that detray relies on singly connected masks, hence a portal from ACTS
/// with multi volume link needs to be split into the multiple volumes
///
/// @note detray also only has outside pointing links
///
/// @return a tuple of json object
std::tuple<std::vector<nlohmann::json>, std::vector<std::shared_ptr<Surface>>>
toJsonDetray(
    const GeometryContext& gctx, const Experimental::Portal& portal,
    std::size_t ip, const Experimental::DetectorVolume& volume,
    const std::vector<OrientedSurface>& orientedSurfaces,
    const std::vector<const Experimental::DetectorVolume*>& detectorVolumes,
    const Options& options = Options{});


/// @brief Convert a surface material to json - detray format
///
/// @param surfaceMaterial is the surface material to be converted
/// @param surface is the surface the material is attached to
/// @param surfaceIndex is the index of the surface
/// @param gridLink [in, out] is the grid index in the volume
///
/// @note the surface is needed to shift the z boundaries for concentric cylinders
///
/// @return a json object representing the surface material in detray format
nlohmann::json toJsonDetray(const ISurfaceMaterial& material,
                            const Acts::Surface& surface,
                            std::size_t surfaceIndex,
                            std::map<std::size_t, std::size_t>& gridLink);

/// @brief Convert a bin utility to json - detray format
///
/// @param binUtility is the bin utility to be converted
/// @param surface is the surface the material is attached to
///
/// @note the surface is needed to shift the z boundaries for concentric cylinders
///
/// @return a json object representing the bin utility in detray format
nlohmann::json toJsonDetray(const Acts::BinUtility& binUtility,
                            const Acts::Surface& surface);



} // namespace DetectorJsonDetrayConverter

}  // namespace Acts

