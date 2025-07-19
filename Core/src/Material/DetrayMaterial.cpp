// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/DetrayExceptions.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"

#include <numbers>

#include <detray/io/frontend/payloads.hpp>

namespace Acts {

namespace {
detray::io::material_slab_payload convertMaterialSlab(
    const MaterialSlab& slab) {
  detray::io::material_slab_payload payload;
  // Fill the material parameters and the thickness
  const auto& material = slab.material();
  payload.thickness = slab.thickness();
  payload.mat = detray::io::material_payload{
      {material.X0(), material.L0(), material.Ar(), material.Z(),
       material.massDensity(), material.molarDensity(), 0.}};
  payload.type = detray::io::material_id::slab;
  return payload;
}

detray::axis::label convertAxisDirection(AxisDirection bValue) {
  switch (bValue) {
    case AxisDirection::AxisX:
      return detray::axis::label::e_x;
    case AxisDirection::AxisY:
      return detray::axis::label::e_y;
    case AxisDirection::AxisZ:
      return detray::axis::label::e_z;
    case AxisDirection::AxisR:
      return detray::axis::label::e_r;
    case AxisDirection::AxisPhi:
      return detray::axis::label::e_phi;
    case AxisDirection::AxisRPhi:
      return detray::axis::label::e_rphi;
    default:
      throw std::invalid_argument(
          "DetrayMaterialConverter: unknown binning value detected.");
  }
}

detray::axis::bounds convertBinningOption(BinningOption bOption) {
  // That's a bit of a mind bender, but the conversion is correct
  // closed -> axis are closed, i.e. circular
  // open -> axis are not closed, but the range is closed (no overflow bin) ->
  // closed
  switch (bOption) {
    case BinningOption::closed:
      return detray::axis::bounds::e_circular;
    case BinningOption::open:
      return detray::axis::bounds::e_closed;
    default:
      throw std::invalid_argument(
          "DetrayMaterialConverter: unknown binning option detected.");
  }
}

detray::axis::binning convertBinningType(BinningType bType) {
  switch (bType) {
    case BinningType::equidistant:
      return detray::axis::binning::e_regular;
    case BinningType::arbitrary:
      return detray::axis::binning::e_irregular;
    default:
      throw std::invalid_argument(
          "DetrayMaterialConverter: unknown binning type detected.");
  }
}

detray::io::axis_payload convertBinningData(const BinningData& bData) {
  detray::io::axis_payload axis;

  axis.bins = bData.bins();
  // Set the binning type
  axis.binning = convertBinningType(bData.type);
  // Set the binning option
  axis.bounds = convertBinningOption(bData.option);
  // Set the binning value
  axis.label = convertAxisDirection(bData.binvalue);
  // Set the binning range
  axis.edges = {};
  if (bData.type == BinningType::equidistant) {
    axis.edges = {bData.min, bData.max};
  } else {
    axis.edges.insert(axis.edges.end(), bData.boundaries().begin(),
                      bData.boundaries().end());
  }
  return axis;
}
}  // namespace

std::unique_ptr<DetraySurfaceMaterial> BinnedSurfaceMaterial::toDetrayPayload()
    const {
  // BinUtility modifications
  bool swapped = false;
  // Get the bin utility (make a copy as we may modify it)
  // Detray expects 2-dimensional grid, currently supported are
  // x-y, r-phi, phi-z
  BinUtility bUtility = binUtility();
  // Turn the bin value into a 2D grid
  if (bUtility.dimensions() == 1u) {
    if (bUtility.binningData()[0u].binvalue == AxisDirection::AxisX) {
      // Turn to X-Y
      bUtility += BinUtility(1u, std::numeric_limits<float>::lowest(),
                             std::numeric_limits<float>::max(),
                             BinningOption::closed, AxisDirection::AxisY);
    } else if (bUtility.binningData()[0u].binvalue == AxisDirection::AxisY) {
      // Turn to X-Y
      BinUtility nbUtility(1u, std::numeric_limits<float>::lowest(),
                           std::numeric_limits<float>::max(),
                           BinningOption::closed, AxisDirection::AxisX);
      nbUtility += bUtility;
      bUtility = std::move(nbUtility);
      swapped = true;
    } else if (bUtility.binningData()[0u].binvalue == AxisDirection::AxisR) {
      // Turn to R-Phi
      bUtility += BinUtility(1u, -std::numbers::pi, std::numbers::pi, closed,
                             AxisDirection::AxisPhi);
    } else if (bUtility.binningData()[0u].binvalue == AxisDirection::AxisZ) {
      // Turn to Phi-Z - swap needed
      BinUtility nbUtility(1u, -std::numbers::pi, std::numbers::pi, closed,
                           AxisDirection::AxisPhi);
      nbUtility += bUtility;
      bUtility = std::move(nbUtility);
      swapped = true;
    } else {
      std::invalid_argument("Unsupported binning for Detray");
    }
  } else if (bUtility.dimensions() == 2u &&
             bUtility.binningData()[0u].binvalue == AxisDirection::AxisZ &&
             bUtility.binningData()[1u].binvalue == AxisDirection::AxisPhi) {
    BinUtility nbUtility(bUtility.binningData()[1u]);
    nbUtility += BinUtility{bUtility.binningData()[0u]};
    bUtility = std::move(nbUtility);
    swapped = true;
  }

  AxisDirection bVal0 = bUtility.binningData()[0u].binvalue;
  AxisDirection bVal1 = bUtility.binningData()[1u].binvalue;

  // Translate into grid index type
  detray::io::material_id gridIndexType = detray::io::material_id::unknown;
  if (bVal0 == AxisDirection::AxisR && bVal1 == AxisDirection::AxisPhi) {
    gridIndexType = detray::io::material_id::ring2_map;
  } else if (bVal0 == AxisDirection::AxisPhi && bVal1 == AxisDirection::AxisZ) {
    gridIndexType = detray::io::material_id::concentric_cylinder2_map;
  } else if (bVal0 == AxisDirection::AxisX && bVal1 == AxisDirection::AxisY) {
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
    auto axisPayload = convertBinningData(bData);
    materialGrid.axes.push_back(axisPayload);
  }

  // Convert the material slabs from the material matrix
  auto materialMatrix = fullMaterial();
  for (std::size_t ib1 = 0; ib1 < materialMatrix.size(); ++ib1) {
    for (std::size_t ib0 = 0; ib0 < materialMatrix[0u].size(); ++ib0) {
      // Translate into a local bin
      std::size_t lb0 = swapped ? ib1 : ib0;
      std::size_t lb1 = swapped ? ib0 : ib1;
      detray::io::material_slab_payload slab =
          convertMaterialSlab(materialMatrix[ib1][ib0]);
      detray::io::grid_bin_payload<detray::io::material_slab_payload> slabBin{
          {static_cast<unsigned int>(lb0), static_cast<unsigned int>(lb1)},
          {slab}};
      // Fill into the grid
      materialGrid.bins.push_back(slabBin);
    }
  }
  return std::make_unique<DetraySurfaceMaterial>(materialGrid);
}

template <>
std::unique_ptr<DetraySurfaceMaterial>
ProtoSurfaceMaterialT<Acts::BinUtility>::toDetrayPayload() const {
  // Does not apply to detray
  return nullptr;
}

template <>
std::unique_ptr<DetraySurfaceMaterial>
ProtoSurfaceMaterialT<std::vector<DirectedProtoAxis>>::toDetrayPayload() const {
  // Does not apply to detray
  return nullptr;
}

std::unique_ptr<DetraySurfaceMaterial>
detail::IGridSurfaceMaterialBase::toDetrayPayload() const {
  throw DetrayUnsupportedMaterialException("detail::IGridSurfaceMaterialBase");
}

std::unique_ptr<DetraySurfaceMaterial>
HomogeneousSurfaceMaterial::toDetrayPayload() const {
  return std::make_unique<DetraySurfaceMaterial>(
      convertMaterialSlab(materialSlab()));
}

}  // namespace Acts
