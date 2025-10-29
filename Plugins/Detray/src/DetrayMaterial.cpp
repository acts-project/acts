// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Detray/DetrayPayloadConverter.hpp"
//
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "ActsPlugins/Detray/DetrayConversionUtils.hpp"

#include <stdexcept>

#include <detray/io/frontend/payloads.hpp>

using namespace Acts;

namespace ActsPlugins {

using DetraySurfaceMaterial = DetrayPayloadConverter::DetraySurfaceMaterial;
using DetraySurfaceGrid = DetrayPayloadConverter::DetraySurfaceGrid;

std::optional<DetraySurfaceMaterial>
DetrayPayloadConverter::convertBinnedSurfaceMaterial(
    const BinnedSurfaceMaterial& material) {
  using enum AxisDirection;
  // Get the bin utility and convert to 2D if needed
  // Detray expects 2-dimensional grid, currently supported are
  // x-y, r-phi, phi-z
  auto [bUtility, swapped] =
      DetrayConversionUtils::convertBinUtilityTo2D(material.binUtility());

  AxisDirection bVal0 = bUtility.binningData()[0u].binvalue;
  AxisDirection bVal1 = bUtility.binningData()[1u].binvalue;

  // Translate into grid index type
  detray::io::material_id gridIndexType = detray::io::material_id::unknown;
  if (bVal0 == AxisR && bVal1 == AxisPhi) {
    gridIndexType = detray::io::material_id::ring2_map;
  } else if (bVal0 == AxisPhi && bVal1 == AxisZ) {
    gridIndexType = detray::io::material_id::concentric_cylinder2_map;
  } else if (bVal0 == AxisX && bVal1 == AxisY) {
    gridIndexType = detray::io::material_id::rectangle2_map;
  } else {
    std::runtime_error(
        "DetrayMaterialConverter: Unsupported binning for Detray");
  }

  detray::io::grid_payload<detray::io::material_slab_payload,
                           detray::io::material_id>
      materialGrid;

  detray::io::typed_link_payload<detray::io::material_id> linkPayload{
      gridIndexType, 0u};
  materialGrid.grid_link = linkPayload;

  // Now convert the (modified) bin utility
  for (const auto& bData : bUtility.binningData()) {
    auto axisPayload = DetrayConversionUtils::convertBinningData(bData);
    materialGrid.axes.push_back(axisPayload);
  }

  // Convert the material slabs from the material matrix
  auto materialMatrix = material.fullMaterial();
  for (std::size_t ib1 = 0; ib1 < materialMatrix.size(); ++ib1) {
    for (std::size_t ib0 = 0; ib0 < materialMatrix[0u].size(); ++ib0) {
      // Translate into a local bin
      std::size_t lb0 = swapped ? ib1 : ib0;
      std::size_t lb1 = swapped ? ib0 : ib1;
      detray::io::material_slab_payload slab =
          DetrayConversionUtils::convertMaterialSlab(materialMatrix[ib1][ib0]);
      detray::io::grid_bin_payload<detray::io::material_slab_payload> slabBin{
          {static_cast<unsigned int>(lb0), static_cast<unsigned int>(lb1)},
          {slab}};
      // Fill into the grid
      materialGrid.bins.push_back(slabBin);
    }
  }
  return materialGrid;
}

std::optional<DetraySurfaceMaterial>
DetrayPayloadConverter::convertGridSurfaceMaterial(
    const IGridSurfaceMaterialBase& /*material*/) {
  throw DetrayUnsupportedMaterialException("detail::IGridSurfaceMaterialBase");
}

std::optional<DetraySurfaceMaterial>
DetrayPayloadConverter::convertHomogeneousSurfaceMaterial(
    const HomogeneousSurfaceMaterial& material) {
  return DetrayConversionUtils::convertMaterialSlab(material.materialSlab());
}

std::optional<DetraySurfaceMaterial>
DetrayPayloadConverter::convertProtoSurfaceMaterialBinUtility(
    const ProtoSurfaceMaterialT<Acts::BinUtility>& /*material*/) {
  return std::nullopt;
}

std::optional<DetraySurfaceMaterial>
DetrayPayloadConverter::convertProtoSurfaceMaterialProtoAxes(
    const ProtoSurfaceMaterialT<std::vector<DirectedProtoAxis>>& /*material*/) {
  return std::nullopt;
}

DetrayUnsupportedMaterialException::DetrayUnsupportedMaterialException(
    std::string_view name)
    : std::runtime_error(std::string("Material type ") + std::string(name) +
                         " not supported by detray") {}

}  // namespace ActsPlugins
