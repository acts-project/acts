// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/ProtoSurfaceMaterial.hpp"

#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/IAxis.hpp"

#include <ostream>
#include <stdexcept>

namespace Acts {

namespace {

/// Construct a single-bin dummy axis on the given direction.
DirectedProtoAxis makeDummyAxis(AxisDirection dir) {
  return DirectedProtoAxis(dir, AxisBoundaryType::Bound, 0., 1., 1);
}

/// Convert a single BinningData entry to a DirectedProtoAxis.
DirectedProtoAxis binningDataToProtoAxis(const BinningData& bd) {
  AxisBoundaryType bType = (bd.option == closed) ? AxisBoundaryType::Closed
                                                 : AxisBoundaryType::Bound;
  if (bd.type == equidistant) {
    return DirectedProtoAxis(bd.binvalue, bType, static_cast<double>(bd.min),
                             static_cast<double>(bd.max), bd.bins());
  }
  const auto& fbounds = bd.boundaries();
  std::vector<double> edges(fbounds.begin(), fbounds.end());
  return DirectedProtoAxis(bd.binvalue, bType, edges);
}

/// Convert a BinUtility to a pair of DirectedProtoAxis.
std::array<DirectedProtoAxis, 2> binUtilityToAxes(const BinUtility& bu) {
  const auto& bdata = bu.binningData();
  DirectedProtoAxis ax0 = binningDataToProtoAxis(bdata[0]);
  if (bdata.size() >= 2) {
    return {ax0, binningDataToProtoAxis(bdata[1])};
  }
  AxisDirection dummyDir = (bdata[0].binvalue == AxisDirection::AxisX)
                               ? AxisDirection::AxisY
                               : AxisDirection::AxisX;
  return {ax0, makeDummyAxis(dummyDir)};
}

/// Convert a vector of DirectedProtoAxis to a fixed-size pair.
std::array<DirectedProtoAxis, 2> vectorToAxes(
    const std::vector<DirectedProtoAxis>& axes) {
  if (axes.empty()) {
    return {makeDummyAxis(AxisDirection::AxisX),
            makeDummyAxis(AxisDirection::AxisY)};
  }
  if (axes.size() == 1) {
    AxisDirection dummyDir =
        (axes[0].getAxisDirection() == AxisDirection::AxisX)
            ? AxisDirection::AxisY
            : AxisDirection::AxisX;
    return {axes[0], makeDummyAxis(dummyDir)};
  }
  return {axes[0], axes[1]};
}

}  // namespace

std::array<DirectedProtoAxis, 2> protoAxesFromBinUtility(const BinUtility& bu) {
  return binUtilityToAxes(bu);
}

ProtoSurfaceMaterial::ProtoSurfaceMaterial()
    : m_axes{makeDummyAxis(AxisDirection::AxisX),
             makeDummyAxis(AxisDirection::AxisY)} {}

ProtoSurfaceMaterial::ProtoSurfaceMaterial(
    std::array<DirectedProtoAxis, 2> axes, MappingType mappingType)
    : ISurfaceMaterial(1., mappingType), m_axes(std::move(axes)) {}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

ProtoSurfaceMaterial::ProtoSurfaceMaterial(const BinUtility& binUtil,
                                           MappingType mappingType)
    : ISurfaceMaterial(1., mappingType), m_axes(binUtilityToAxes(binUtil)) {}

ProtoSurfaceMaterial::ProtoSurfaceMaterial(
    const std::vector<DirectedProtoAxis>& axes, MappingType mappingType)
    : ISurfaceMaterial(1., mappingType), m_axes(vectorToAxes(axes)) {}

#pragma GCC diagnostic pop

// ---------------------------------------------------------------------------
// Deprecated binning() accessor
// ---------------------------------------------------------------------------

BinUtility ProtoSurfaceMaterial::binning() const {
  std::vector<DirectedProtoAxis> sigAxes;
  for (const auto& ax : m_axes) {
    if (ax.getAxis().getNBins() > 1) {
      sigAxes.push_back(ax);
    }
  }
  return BinUtility(sigAxes);
}

// ---------------------------------------------------------------------------
// localAxisDirections
// ---------------------------------------------------------------------------

std::vector<AxisDirection> ProtoSurfaceMaterial::localAxisDirections() const {
  std::vector<AxisDirection> dirs;
  for (const auto& ax : m_axes) {
    if (ax.getAxis().getNBins() > 1) {
      dirs.push_back(ax.getAxisDirection());
    }
  }
  return dirs;
}

// ---------------------------------------------------------------------------
// toStream
// ---------------------------------------------------------------------------

std::ostream& ProtoSurfaceMaterial::toStream(std::ostream& sl) const {
  sl << "Acts::ProtoSurfaceMaterial : " << std::endl;
  sl << "  axis 0: " << m_axes[0] << std::endl;
  sl << "  axis 1: " << m_axes[1] << std::endl;
  return sl;
}

}  // namespace Acts
