// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Detray/DetrayConversionUtils.hpp"

#include <limits>
#include <numbers>
#include <stdexcept>

using namespace Acts;

namespace ActsPlugins {

detray::axis::label DetrayConversionUtils::convertAxisDirection(
    AxisDirection bValue) {
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

detray::axis::bounds DetrayConversionUtils::convertBinningOption(
    BinningOption bOption) {
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

detray::axis::binning DetrayConversionUtils::convertBinningType(
    BinningType bType) {
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

detray::io::axis_payload DetrayConversionUtils::convertBinningData(
    const BinningData& bData) {
  detray::io::axis_payload axisPayload;

  // Number of bins
  axisPayload.bins = bData.bins();
  // Set the binning type
  axisPayload.binning = convertBinningType(bData.type);
  // Set the binning option
  axisPayload.bounds = convertBinningOption(bData.option);
  // Set the binning value
  axisPayload.label = convertAxisDirection(bData.binvalue);
  // Set the binning range
  axisPayload.edges = {};
  if (bData.type == BinningType::equidistant) {
    axisPayload.edges = {bData.min, bData.max};
  } else {
    axisPayload.edges.insert(axisPayload.edges.end(),
                             bData.boundaries().begin(),
                             bData.boundaries().end());
  }
  return axisPayload;
}

detray::io::axis_payload DetrayConversionUtils::convertAxis(
    const Acts::IAxis& axis) {
  detray::io::axis_payload axisPayload;

  // Number of bins
  axisPayload.bins = axis.getNBins();
  // Set the binning type and bin edges
  if (axis.isEquidistant()) {
    axisPayload.binning = detray::axis::binning::e_regular;
    axisPayload.edges = {axis.getMin(), axis.getMax()};
  } else {
    axisPayload.binning = detray::axis::binning::e_irregular;
    axisPayload.edges = axis.getBinEdges();
  }
  // Axis boundary behaviour
  switch (axis.getBoundaryType()) {
    using enum Acts::AxisBoundaryType;
    case Open:
      // Open interval: Overflow bins
      axisPayload.bounds = detray::axis::bounds::e_open;
      break;
    case Closed:
      // Periodic boundary conditions
      axisPayload.bounds = detray::axis::bounds::e_circular;
      break;
    case Bound:
      // Closed interval: no overflow bins
      axisPayload.bounds = detray::axis::bounds::e_closed;
      break;
  }

  return axisPayload;
}

detray::io::surface_material_payload DetrayConversionUtils::convertMaterialSlab(
    const Acts::MaterialSlab& slab) {
  detray::io::surface_material_payload matPayload;

  // Fill the material parameters and the thickness
  const auto& material = slab.material();
  matPayload.thickness = slab.thickness();
  matPayload.mat = detray::io::material_param_payload{
      {material.X0(), material.L0(), material.Ar(), material.Z(),
       material.massDensity(), material.molarDensity(), 0.}};
  matPayload.type = detray::io::material_id::slab;

  return matPayload;
}

detray::io::transform_payload DetrayConversionUtils::convertTransform(
    const Acts::Transform3& transform) {
  detray::io::transform_payload trfPayload;

  Eigen::Map<Acts::Vector3> tr{trfPayload.tr.data()};
  tr = transform.translation();

  Eigen::Map<Acts::SquareMatrix3> rot{trfPayload.rot.data()};
  rot = transform.linear();

  return trfPayload;
}

std::tuple<Acts::BinUtility, bool> DetrayConversionUtils::convertBinUtilityTo2D(
    const Acts::BinUtility& bUtility) {
  using enum Acts::AxisDirection;
  using enum Acts::BinningOption;

  // Return as-is if already 2D
  if (bUtility.dimensions() == 2u) {
    // Check if we need to swap for rphi-z -> z-rphi
    if (bUtility.binningData()[0u].binvalue == AxisZ &&
        bUtility.binningData()[1u].binvalue == AxisRPhi) {
      BinUtility nbUtility(bUtility.binningData()[1u]);
      nbUtility += BinUtility{bUtility.binningData()[0u]};
      return {std::move(nbUtility), true};
    }
    return {bUtility, false};
  }

  // Convert 1D to 2D
  if (bUtility.dimensions() == 1u) {
    BinUtility result = bUtility;
    bool swapped = false;

    if (bUtility.binningData()[0u].binvalue == AxisX) {
      // Turn to X-Y
      result += BinUtility(1u, std::numeric_limits<float>::lowest(),
                           std::numeric_limits<float>::max(), closed, AxisY);
    } else if (bUtility.binningData()[0u].binvalue == AxisY) {
      // Turn to X-Y (swap needed)
      BinUtility nbUtility(1u, std::numeric_limits<float>::lowest(),
                           std::numeric_limits<float>::max(), closed, AxisX);
      nbUtility += bUtility;
      result = std::move(nbUtility);
      swapped = true;
    } else if (bUtility.binningData()[0u].binvalue == AxisR) {
      // Turn to R-Phi
      result +=
          BinUtility(1u, -std::numbers::pi, std::numbers::pi, closed, AxisPhi);
    } else if (bUtility.binningData()[0u].binvalue == AxisZ) {
      // Turn to RPhi-Z (swap needed)
      BinUtility nbUtility(1u, -std::numbers::pi, std::numbers::pi, closed,
                           AxisPhi);
      nbUtility += bUtility;
      result = std::move(nbUtility);
      swapped = true;
    } else {
      throw std::invalid_argument("Unsupported binning for Detray");
    }

    return {result, swapped};
  }

  throw std::invalid_argument(
      "DetrayConversionUtils: BinUtility must be 1D or 2D");
}

}  // namespace ActsPlugins
