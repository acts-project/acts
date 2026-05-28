// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/AccumulatedSurfaceMaterial.hpp"

#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Utilities/ProtoAxisHelpers.hpp"

#include <utility>

// Default Constructor - for homogeneous material
Acts::AccumulatedSurfaceMaterial::AccumulatedSurfaceMaterial(double splitFactor)
    : m_splitFactor(splitFactor) {
  AccumulatedVector accMat = {{AccumulatedMaterialSlab()}};
  m_accumulatedMaterial = {{accMat}};
}

// Binned Material constructor with split factor
// TODO: Remove this constructor after DirectedProtoAxis migration
Acts::AccumulatedSurfaceMaterial::AccumulatedSurfaceMaterial(
    const BinUtility& binUtility, double splitFactor)
    : m_binUtility(binUtility), m_splitFactor(splitFactor) {
  std::size_t bins0 = m_binUtility.bins(0);
  std::size_t bins1 = m_binUtility.bins(1);
  AccumulatedVector accVec(bins0, AccumulatedMaterialSlab());
  m_accumulatedMaterial = AccumulatedMatrix(bins1, accVec);
}

// Constructor using DirectedProtoAxes for binning
Acts::AccumulatedSurfaceMaterial::AccumulatedSurfaceMaterial(
    const std::vector<DirectedProtoAxis>& axes, double splitFactor)
    : m_axes(axes), m_splitFactor(splitFactor) {
  if (m_axes.size() < 1 || m_axes.size() > 2) {
    throw std::invalid_argument(
        "AccumulatedSurfaceMaterial: At least one, maximum 2  proto axis is "
        "required for binning.");
  }
  if (m_axes.size() == 1) {
    std::size_t bins0 = m_axes[0].getAxis().getNBins();
    AccumulatedVector accVec(bins0, AccumulatedMaterialSlab());
    m_accumulatedMaterial = AccumulatedMatrix(1, accVec);
    return;
  }
  if (m_axes.size() == 2) {
    std::size_t bins0 = m_axes[0].getAxis().getNBins();
    std::size_t bins1 = m_axes[1].getAxis().getNBins();
    AccumulatedVector accVec(bins0, AccumulatedMaterialSlab());
    m_accumulatedMaterial = AccumulatedMatrix(bins1, accVec);
    return;
  }
}

// Assign a material properties object
std::array<std::size_t, 3> Acts::AccumulatedSurfaceMaterial::accumulate(
    const Vector2& lp, const MaterialSlab& mp, double pathCorrection) {
  if (m_axes.empty()) {
    m_accumulatedMaterial[0][0].accumulate(mp, pathCorrection);
    return {0, 0, 0};
  }
  std::size_t bin0 = Acts::ProtoAxisHelpers::binFromProtoAxis(m_axes[0], lp);
  std::size_t bin1 = Acts::ProtoAxisHelpers::binFromProtoAxis(m_axes[1], lp);
  m_accumulatedMaterial[bin1][bin0].accumulate(mp, pathCorrection);
  return {bin0, bin1, 0};
}

// Assign a material properties object
std::array<std::size_t, 3> Acts::AccumulatedSurfaceMaterial::accumulate(
    const Vector3& gp, const MaterialSlab& mp, double pathCorrection) {
  if (m_axes.empty()) {
    m_accumulatedMaterial[0][0].accumulate(mp, pathCorrection);
    return {0, 0, 0};
  }
  std::array<std::size_t, 3> bTriple =
      Acts::ProtoAxisHelpers::binTripleFromProtoAxes(m_axes, gp);
  m_accumulatedMaterial[bTriple[1]][bTriple[0]].accumulate(mp, pathCorrection);
  return bTriple;
}

// Void average for vacuum assignment
void Acts::AccumulatedSurfaceMaterial::trackVariance(const Vector3& gp,
                                                     MaterialSlab slabReference,
                                                     bool emptyHit) {
  if (m_axes.empty()) {
    m_accumulatedMaterial[0][0].trackVariance(slabReference, emptyHit);
    return;
  }
  std::array<std::size_t, 3> bTriple =
      Acts::ProtoAxisHelpers::binTripleFromProtoAxes(m_axes, gp);
  std::vector<std::array<std::size_t, 3>> trackBins = {bTriple};
  trackVariance(trackBins, slabReference);
}

// Average the information accumulated during one event
void Acts::AccumulatedSurfaceMaterial::trackVariance(
    const std::vector<std::array<std::size_t, 3>>& trackBins,
    MaterialSlab slabReference, bool emptyHit) {
  // the homogeneous material case
  if (m_axes.empty()) {
    m_accumulatedMaterial[0][0].trackVariance(slabReference, emptyHit);
    return;
  }
  // The touched bins are known, so you can access them directly
  if (!trackBins.empty()) {
    for (auto bin : trackBins) {
      m_accumulatedMaterial[bin[1]][bin[0]].trackVariance(slabReference);
    }
  } else {
    // Touched bins are not known: Run over all bins
    for (auto& matVec : m_accumulatedMaterial) {
      for (auto& mat : matVec) {
        mat.trackVariance(slabReference);
      }
    }
  }
}

// Void average for vacuum assignment
void Acts::AccumulatedSurfaceMaterial::trackAverage(const Vector3& gp,
                                                    bool emptyHit) {
  if (m_axes.empty()) {
    m_accumulatedMaterial[0][0].trackAverage(emptyHit);
    return;
  }

  std::array<std::size_t, 3> bTriple =
      Acts::ProtoAxisHelpers::binTripleFromProtoAxes(m_axes, gp);
  std::vector<std::array<std::size_t, 3>> trackBins = {bTriple};
  trackAverage(trackBins, emptyHit);
}

// Average the information accumulated during one event
void Acts::AccumulatedSurfaceMaterial::trackAverage(
    const std::vector<std::array<std::size_t, 3>>& trackBins, bool emptyHit) {
  // the homogeneous material case
  // the homogeneous material case
  if (m_axes.size() == 0) {
    m_accumulatedMaterial[0][0].trackAverage(emptyHit);
    return;
  }

  // The touched bins are known, so you can access them directly
  if (!trackBins.empty()) {
    for (auto bin : trackBins) {
      m_accumulatedMaterial[bin[1]][bin[0]].trackAverage(emptyHit);
    }
  } else {
    // Touched bins are not known: Run over all bins
    for (auto& matVec : m_accumulatedMaterial) {
      for (auto& mat : matVec) {
        mat.trackAverage(emptyHit);
      }
    }
  }
}

/// Total average creates SurfaceMaterial
std::unique_ptr<const Acts::ISurfaceMaterial>
Acts::AccumulatedSurfaceMaterial::totalAverage() {
  std::cout << "AccumulatedSurfaceMaterial::totalAverage called with "
            << m_axes.size() << " axes." << std::endl;
  if (ProtoAxisHelpers::totalBinsFromProtoAxes(m_axes) <= 1) {
    // Return HomogeneousSurfaceMaterial
    return std::make_unique<HomogeneousSurfaceMaterial>(
        m_accumulatedMaterial[0][0].totalAverage().first, m_splitFactor);
  }

  // number of bins per axis from DirectedProtoAxis
  std::size_t bins0 = m_axes[0].getAxis().getNBins();
  std::size_t bins1 = m_axes[1].getAxis().getNBins();

  // build the material-property matrix and fill from accumulated data
  MaterialSlabMatrix mpMatrix(
      bins1, MaterialSlabVector(bins0, MaterialSlab::Nothing()));
  for (std::size_t ib1 = 0; ib1 < bins1; ++ib1) {
    for (std::size_t ib0 = 0; ib0 < bins0; ++ib0) {
      mpMatrix[ib1][ib0] = m_accumulatedMaterial[ib1][ib0].totalAverage().first;
    }
  }
  return std::make_unique<const BinnedSurfaceMaterial>(
      m_axes, std::move(mpMatrix), m_splitFactor);
}
