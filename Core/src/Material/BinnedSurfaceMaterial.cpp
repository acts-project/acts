// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/BinnedSurfaceMaterial.hpp"

#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/IAxis.hpp"

#include <algorithm>
#include <ostream>
#include <utility>
#include <vector>

namespace Acts {

namespace {

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
/// For a 1D BinUtility a single-bin dummy axis is added as the second axis.
std::array<DirectedProtoAxis, 2> binUtilityToAxes(const BinUtility& bu) {
  const auto& bdata = bu.binningData();
  DirectedProtoAxis ax0 = binningDataToProtoAxis(bdata[0]);
  if (bdata.size() >= 2) {
    return {ax0, binningDataToProtoAxis(bdata[1])};
  }
  // 1D → pad with a single-bin dummy axis on the orthogonal local direction
  AxisDirection dummyDir = (bdata[0].binvalue == AxisDirection::AxisX)
                               ? AxisDirection::AxisY
                               : AxisDirection::AxisX;
  return {ax0, DirectedProtoAxis(dummyDir, AxisBoundaryType::Bound, 0., 1., 1)};
}

}  // namespace

BinnedSurfaceMaterial::BinnedSurfaceMaterial(
    std::array<DirectedProtoAxis, 2> axes, MaterialSlabMatrix materialMatrix,
    double splitFactor, MappingType mappingType)
    : ISurfaceMaterial(splitFactor, mappingType),
      m_axes(std::move(axes)),
      m_fullMaterial(std::move(materialMatrix)) {
  const std::size_t nBins0 = m_axes[0].getAxis().getNBins();
  const std::size_t nBins1 = m_axes[1].getAxis().getNBins();
  if (m_fullMaterial.size() != nBins1) {
    throw std::invalid_argument(
        "BinnedSurfaceMaterial: number of material rows does not match axis-1 "
        "bin count.");
  }
  for (const auto& row : m_fullMaterial) {
    if (row.size() != nBins0) {
      throw std::invalid_argument(
          "BinnedSurfaceMaterial: number of material columns does not match "
          "axis-0 bin count.");
    }
  }
}
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

BinnedSurfaceMaterial::BinnedSurfaceMaterial(const BinUtility& binUtility,
                                             MaterialSlabVector materialVector,
                                             double splitFactor,
                                             MappingType mappingType)
    : BinnedSurfaceMaterial(binUtilityToAxes(binUtility),
                            MaterialSlabMatrix{std::move(materialVector)},
                            splitFactor, mappingType) {
  if (binUtility.dimensions() != 1) {
    throw std::invalid_argument(
        "BinnedSurfaceMaterial with material vector only supports 1D binning.");
  }
}

BinnedSurfaceMaterial::BinnedSurfaceMaterial(const BinUtility& binUtility,
                                             MaterialSlabMatrix materialMatrix,
                                             double splitFactor,
                                             MappingType mappingType)
    : BinnedSurfaceMaterial(binUtilityToAxes(binUtility),
                            std::move(materialMatrix), splitFactor,
                            mappingType) {
  if (binUtility.dimensions() != 1 && binUtility.dimensions() != 2) {
    throw std::invalid_argument(
        "BinnedSurfaceMaterial with material matrix only supports 1D and 2D "
        "binning.");
  }
}

#pragma GCC diagnostic pop

// ---------------------------------------------------------------------------
// scale
// ---------------------------------------------------------------------------

BinnedSurfaceMaterial& BinnedSurfaceMaterial::scale(double factor) {
  for (auto& row : m_fullMaterial) {
    for (auto& bin : row) {
      bin.scaleThickness(factor);
    }
  }
  return *this;
}

// ---------------------------------------------------------------------------
// Deprecated binUtility() accessor
// ---------------------------------------------------------------------------

BinUtility BinnedSurfaceMaterial::binUtility() const {
  std::vector<DirectedProtoAxis> sigAxes;
  for (const auto& ax : m_axes) {
    if (ax.getAxis().getNBins() > 1) {
      sigAxes.push_back(ax);
    }
  }
  return BinUtility(sigAxes);
}

// ---------------------------------------------------------------------------
// materialSlab(Vector2)
// ---------------------------------------------------------------------------

const MaterialSlab& BinnedSurfaceMaterial::materialSlab(
    const Vector2& lp) const {
  const IAxis& ax0 = m_axes[0].getAxis();
  const IAxis& ax1 = m_axes[1].getAxis();
  // IAxis::getBin is 1-based (0 = underflow, nBins+1 = overflow); clamp to
  // valid range then convert to 0-based index.
  const std::size_t ibin0 =
      std::clamp(ax0.getBin(lp[0]), std::size_t{1}, ax0.getNBins()) - 1;
  const std::size_t ibin1 =
      std::clamp(ax1.getBin(lp[1]), std::size_t{1}, ax1.getNBins()) - 1;
  return m_fullMaterial[ibin1][ibin0];
}

// ---------------------------------------------------------------------------
// Deprecated materialSlab(Vector3)
// ---------------------------------------------------------------------------

const MaterialSlab& BinnedSurfaceMaterial::materialSlab(
    const Vector3& gp) const {
  // Reconstruct a BinUtility for the global-position lookup (deprecated path).
  BinUtility bu(std::vector<DirectedProtoAxis>{m_axes[0], m_axes[1]});
  const std::size_t ibin0 = bu.bin(gp, 0);
  const std::size_t ibin1 = bu.bin(gp, 1);
  return m_fullMaterial[ibin1][ibin0];
}

// ---------------------------------------------------------------------------
// localAxisDirections
// ---------------------------------------------------------------------------

std::vector<AxisDirection> BinnedSurfaceMaterial::localAxisDirections() const {
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

std::ostream& BinnedSurfaceMaterial::toStream(std::ostream& sl) const {
  sl << "BinnedSurfaceMaterial : " << std::endl;
  sl << "   - Number of Material bins [0,1] : "
     << m_axes[0].getAxis().getNBins() << " / "
     << m_axes[1].getAxis().getNBins() << std::endl;
  sl << "   - Parse full update material    : " << std::endl;
  unsigned int imat1 = 0;
  for (const auto& row : m_fullMaterial) {
    unsigned int imat0 = 0;
    for (const auto& bin : row) {
      sl << " Bin [" << imat1 << "][" << imat0 << "] - " << bin;
      ++imat0;
    }
    ++imat1;
  }
  sl << "  - Axes: " << m_axes[0] << " | " << m_axes[1] << std::endl;
  return sl;
}

}  // namespace Acts
