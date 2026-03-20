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
  /// initializing the three vectors needed for Detailed Material description:
  m_totalThickness = std::vector<std::vector<double>>(
	1, std::vector<double>(1, 0.0));
  m_elementZ = std::vector<std::vector<std::vector<unsigned int>>>(
	1, std::vector<std::vector<unsigned int>>(1));
  m_weightedFrac = std::vector<std::vector<std::vector<float>>>(
	1, std::vector<std::vector<float>>(1));
}

// Binned Material constructor with split factor
Acts::AccumulatedSurfaceMaterial::AccumulatedSurfaceMaterial(
    const BinUtility& binUtility, double splitFactor)
    : m_binUtility(binUtility), m_splitFactor(splitFactor) {
  std::size_t bins0 = m_binUtility.bins(0);
  std::size_t bins1 = m_binUtility.bins(1);
  AccumulatedVector accVec(bins0, AccumulatedMaterialSlab());
  m_accumulatedMaterial = AccumulatedMatrix(bins1, accVec);
  m_totalThickness = std::vector<std::vector<double>>(
	bins1, std::vector<double>(bins0, 0.0));
  m_elementZ = std::vector<std::vector<std::vector<unsigned int>>>(
	bins1, std::vector<std::vector<unsigned int>>(bins0));
  m_weightedFrac = std::vector<std::vector<std::vector<float>>>(
	bins1, std::vector<std::vector<float>>(bins0));
}

// Assign a material properties object
std::array<std::size_t, 3> Acts::AccumulatedSurfaceMaterial::accumulate(
    const Vector2& lp, const MaterialSlab& mp, double pathCorrection,
    const std::vector<unsigned int>& elementZ,
    const std::vector<float>& elementFrac) {
  if (m_binUtility.dimensions() == 0) {
    m_accumulatedMaterial[0][0].accumulate(mp, pathCorrection);
    accumulateElementData(0, 0, mp, pathCorrection, elementZ, elementFrac);
    return {0, 0, 0};
  }
  std::size_t bin0 = m_binUtility.bin(lp, 0);
  std::size_t bin1 = m_binUtility.bin(lp, 1);
  m_accumulatedMaterial[bin1][bin0].accumulate(mp, pathCorrection);
  accumulateElementData(bin0, bin1, mp, pathCorrection, elementZ, elementFrac);
  return {bin0, bin1, 0};
}


// Assign a material properties object
std::array<std::size_t, 3> Acts::AccumulatedSurfaceMaterial::accumulate(
    const Vector3& gp, const MaterialSlab& mp, double pathCorrection,
    const std::vector<unsigned int>& elementZ,
    const std::vector<float>& elementFrac) {
  if (m_binUtility.dimensions() == 0) {
    m_accumulatedMaterial[0][0].accumulate(mp, pathCorrection);
    accumulateElementData(0, 0, mp, pathCorrection, elementZ, elementFrac);
    return {0, 0, 0};
  }
  std::array<std::size_t, 3> bTriple = m_binUtility.binTriple(gp);
  m_accumulatedMaterial[bTriple[1]][bTriple[0]].accumulate(mp, pathCorrection);
  accumulateElementData(bTriple[0],bTriple[1], mp, pathCorrection, elementZ, elementFrac);
  return bTriple;
}

// Helper- new element list and fraction accumulation calculation
void Acts::AccumulatedSurfaceMaterial::accumulateElementData(
    std::size_t bin0, size_t bin1, 
    const MaterialSlab& mp, double pathCorrection,
    const std::vector<unsigned int>& elementZ,
    const std::vector<float>& elementFrac) {
  if (elementZ.empty() || elementZ.size() != elementFrac.size()) {
    return;
  }

  double stepThickness = 
	  static_cast<double>(mp.thickness()) / pathCorrection;

  if (stepThickness <= 0.) {
	  return;
  }

  m_totalThickness[bin1][bin0] += stepThickness;
  for (std::size_t i = 0; i <elementZ.size(); ++i) {
    unsigned int z = elementZ[i];
    float frac = elementFrac[i];
    /// check to see if this element is already in the accumulated list
    auto it = std::find(m_elementZ[bin1][bin0].begin(),
			  m_elementZ[bin1][bin0].end(), z);
    if (it != m_elementZ[bin1][bin0].end()) {
      std::size_t idx =
	  std::distance(m_elementZ[bin1][bin0].begin(), it);
      m_weightedFrac[bin1][bin0][idx] += 
	      frac * static_cast<float>(stepThickness);
    } else{
      m_elementZ[bin1][bin0].push_back(z);
      m_weightedFrac[bin1][bin0].push_back(
		      frac * static_cast<float>(stepThickness));
    }
  }
}



// Void average for vacuum assignment
void Acts::AccumulatedSurfaceMaterial::trackVariance(const Vector3& gp,
                                                     MaterialSlab slabReference,
                                                     bool emptyHit) {
  if (m_binUtility.dimensions() == 0) {
    m_accumulatedMaterial[0][0].trackVariance(slabReference, emptyHit);
    return;
  }
  std::array<std::size_t, 3> bTriple = m_binUtility.binTriple(gp);
  std::vector<std::array<std::size_t, 3>> trackBins = {bTriple};
  trackVariance(trackBins, slabReference);
}

// Average the information accumulated during one event
void Acts::AccumulatedSurfaceMaterial::trackVariance(
    const std::vector<std::array<std::size_t, 3>>& trackBins,
    MaterialSlab slabReference, bool emptyHit) {
  // the homogeneous material case
  if (m_binUtility.dimensions() == 0) {
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
  if (m_binUtility.dimensions() == 0) {
    m_accumulatedMaterial[0][0].trackAverage(emptyHit);
    return;
  }

  std::array<std::size_t, 3> bTriple = m_binUtility.binTriple(gp);
  std::vector<std::array<std::size_t, 3>> trackBins = {bTriple};
  trackAverage(trackBins, emptyHit);
}

// Average the information accumulated during one event
void Acts::AccumulatedSurfaceMaterial::trackAverage(
    const std::vector<std::array<std::size_t, 3>>& trackBins, bool emptyHit) {
  // the homogeneous material case
  if (m_binUtility.dimensions() == 0) {
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
  if (m_binUtility.bins() == 1) {
    // Return HomogeneousSurfaceMaterial
    return std::make_unique<HomogeneousSurfaceMaterial>(
        m_accumulatedMaterial[0][0].totalAverage().first, m_splitFactor);
  }
  // Create the properties matrix
  MaterialSlabMatrix mpMatrix(
      m_binUtility.bins(1),
      MaterialSlabVector(m_binUtility.bins(0), MaterialSlab::Nothing()));
  // Loop over and fill
  for (std::size_t ib1 = 0; ib1 < m_binUtility.bins(1); ++ib1) {
    for (std::size_t ib0 = 0; ib0 < m_binUtility.bins(0); ++ib0) {
      mpMatrix[ib1][ib0] = m_accumulatedMaterial[ib1][ib0].totalAverage().first;
    }
  }
  // Build the element composition grids
  BinnedSurfaceMaterial::ElementZMatrix finalElementZ(
	m_binUtility.bins(1),
	std::vector<std::vector<unsigned int>>(m_binUtility.bins(0)));
  BinnedSurfaceMaterial::ElementFracMatrix finalElementFrac(
	m_binUtility.bins(1),
	std::vector<std::vector<float>>(m_binUtility.bins(0)));
  for (std::size_t ib1 = 0; ib1 < m_binUtility.bins(1); ++ib1){
    for (std::size_t ib0 = 0; ib0 < m_binUtility.bins(0); ++ib0){
      finalElementZ[ib1][ib0] = m_elementZ[ib1][ib0];
      double totalThickness = m_totalThickness[ib1][ib0];
      if (totalThickness > 0. && !m_weightedFrac[ib1][ib0].empty()) {
        finalElementFrac[ib1][ib0].resize(m_weightedFrac[ib1][ib0].size());
	for (std::size_t i = 0; i < m_weightedFrac[ib1][ib0].size(); ++i){
	  finalElementFrac[ib1][ib0][i] = 
		m_weightedFrac[ib1][ib0][i] / static_cast<float>(totalThickness);
	}
      }
    }
  }


  // Now return the BinnedSurfaceMaterial
  return std::make_unique<const BinnedSurfaceMaterial>(
      m_binUtility, std::move(mpMatrix), m_splitFactor,
      MappingType::Default,
      std::move(finalElementZ), std::move(finalElementFrac));
}
