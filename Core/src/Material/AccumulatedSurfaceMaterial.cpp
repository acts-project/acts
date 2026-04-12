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

#include <stdexcept>
#include <utility>

// Default Constructor - for homogeneous material
Acts::AccumulatedSurfaceMaterial::AccumulatedSurfaceMaterial(double splitFactor)
    : m_splitFactor(splitFactor) {
  AccumulatedVector accMat = {{AccumulatedMaterialSlab()}};
  m_accumulatedMaterial = {{accMat}};
}

// Binned Material constructor with split factor
Acts::AccumulatedSurfaceMaterial::AccumulatedSurfaceMaterial(
    std::vector<DirectedProtoAxis> directedProtoAxes,
    Transform3 globalToLocalTransform, double splitFactor)
    : m_directedProtoAxes(std::move(directedProtoAxes)),
      m_globalToLocalTransform(std::move(globalToLocalTransform)),
      m_splitFactor(splitFactor) {
  if (m_directedProtoAxes.size() > 2u) {
    throw std::invalid_argument(
        "AccumulatedSurfaceMaterial supports only 0D/1D/2D binning.");
  }
  std::size_t bins0 = m_directedProtoAxes.empty()
                          ? 1u
                          : m_directedProtoAxes[0u].getAxis().getNBins();
  std::size_t bins1 = m_directedProtoAxes.size() > 1u
                          ? m_directedProtoAxes[1u].getAxis().getNBins()
                          : 1u;
  AccumulatedVector accVec(bins0, AccumulatedMaterialSlab());
  m_accumulatedMaterial = AccumulatedMatrix(bins1, accVec);
}

// Assign a material properties object
std::array<std::size_t, 3> Acts::AccumulatedSurfaceMaterial::accumulate(
    const Vector2& lp, const MaterialSlab& mp, double pathCorrection) {
  if (m_directedProtoAxes.empty()) {
    m_accumulatedMaterial[0][0].accumulate(mp, pathCorrection);
    return {0, 0, 0};
  }
  std::array<std::size_t, 3> bins = lookupBins(Vector3(lp[0], lp[1], 0.));
  m_accumulatedMaterial[bins[1]][bins[0]].accumulate(mp, pathCorrection);
  return bins;
}

// Assign a material properties object
std::array<std::size_t, 3> Acts::AccumulatedSurfaceMaterial::accumulate(
    const Vector3& gp, const MaterialSlab& mp, double pathCorrection) {
  if (m_directedProtoAxes.empty()) {
    m_accumulatedMaterial[0][0].accumulate(mp, pathCorrection);
    return {0, 0, 0};
  }
  std::array<std::size_t, 3> bins = lookupBins(m_globalToLocalTransform * gp);
  m_accumulatedMaterial[bins[1]][bins[0]].accumulate(mp, pathCorrection);
  return bins;
}

// Void average for vacuum assignment
void Acts::AccumulatedSurfaceMaterial::trackVariance(const Vector3& gp,
                                                     MaterialSlab slabReference,
                                                     bool emptyHit) {
  if (m_directedProtoAxes.empty()) {
    m_accumulatedMaterial[0][0].trackVariance(slabReference, emptyHit);
    return;
  }
  std::array<std::size_t, 3> bins = lookupBins(m_globalToLocalTransform * gp);
  std::vector<std::array<std::size_t, 3>> trackBins = {bins};
  trackVariance(trackBins, slabReference);
}

// Average the information accumulated during one event
void Acts::AccumulatedSurfaceMaterial::trackVariance(
    const std::vector<std::array<std::size_t, 3>>& trackBins,
    MaterialSlab slabReference, bool emptyHit) {
  // the homogeneous material case
  if (m_directedProtoAxes.empty()) {
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
  if (m_directedProtoAxes.empty()) {
    m_accumulatedMaterial[0][0].trackAverage(emptyHit);
    return;
  }

  std::array<std::size_t, 3> bins = lookupBins(m_globalToLocalTransform * gp);
  std::vector<std::array<std::size_t, 3>> trackBins = {bins};
  trackAverage(trackBins, emptyHit);
}

// Average the information accumulated during one event
void Acts::AccumulatedSurfaceMaterial::trackAverage(
    const std::vector<std::array<std::size_t, 3>>& trackBins, bool emptyHit) {
  // the homogeneous material case
  if (m_directedProtoAxes.empty()) {
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
  const std::size_t bins0 = m_accumulatedMaterial[0u].size();
  const std::size_t bins1 = m_accumulatedMaterial.size();
  if (bins0 * bins1 == 1u) {
    // Return HomogeneousSurfaceMaterial
    return std::make_unique<HomogeneousSurfaceMaterial>(
        m_accumulatedMaterial[0][0].totalAverage().first, m_splitFactor);
  }

  if (m_directedProtoAxes.size() == 1u) {
    MaterialSlabVector mpVector(bins0, MaterialSlab::Nothing());
    for (std::size_t ib0 = 0; ib0 < bins0; ++ib0) {
      mpVector[ib0] = m_accumulatedMaterial[0u][ib0].totalAverage().first;
    }
    return std::make_unique<const BinnedSurfaceMaterial>(
        m_directedProtoAxes[0u], std::move(mpVector), m_splitFactor);
  }

  // Create the properties matrix
  MaterialSlabMatrix mpMatrix(
      bins1, MaterialSlabVector(bins0, MaterialSlab::Nothing()));
  for (std::size_t ib1 = 0; ib1 < bins1; ++ib1) {
    for (std::size_t ib0 = 0; ib0 < bins0; ++ib0) {
      mpMatrix[ib1][ib0] = m_accumulatedMaterial[ib1][ib0].totalAverage().first;
    }
  }
  return std::make_unique<const BinnedSurfaceMaterial>(
      std::array<DirectedProtoAxis, 2u>{m_directedProtoAxes[0u],
                                        m_directedProtoAxes[1u]},
      std::move(mpMatrix), m_splitFactor);
}
