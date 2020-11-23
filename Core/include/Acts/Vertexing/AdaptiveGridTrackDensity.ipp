// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include "Acts/Vertexing/VertexingError.hpp"

#include <algorithm>

template <int trkGridSize>
Acts::Result<float>
Acts::AdaptiveGridTrackDensity<trkGridSize>::getMaxZPosition(const std::vector<float>& mainGridDensity,
    const std::vector<int>& mainGridZValues) const {
  if (mainGridDensity.empty() || mainGridZValues.empty()) {
    return VertexingError::EmptyInput;
  }

  int zbin = -1;
  if (!m_cfg.useHighestSumZPosition) {
    
    int zGridPos = std::distance(mainGridDensity.begin(), std::max_element(mainGridDensity.begin(),mainGridDensity.end()));

    zbin = mainGridZValues[zGridPos];

  } else {
    // Get z position with highest density sum
    // of surrounding bins
    //zbin = getHighestSumZPosition(mainGrid);
  }

  // Derive corresponding z value
  int sign = (zbin > 0) ? +1 : -1;
  return (zbin + sign * 0.5f) * m_cfg.binSize;
}

template <int trkGridSize>
Acts::Result<std::pair<float, float>>
Acts::AdaptiveGridTrackDensity<trkGridSize>::
    getMaxZPositionAndWidth(const std::vector<float>& mainGridDensity,
    const std::vector<int>& mainGridZValues) const {
  // Get z maximum value
  auto maxZRes = getMaxZPosition(mainGridDensity, mainGridZValues);
  if (not maxZRes.ok()) {
    return maxZRes.error();
  }
  float maxZ = *maxZRes;

  // Get seed width estimate
  auto widthRes = estimateSeedWidth(mainGridDensity, mainGridZValues,maxZ);
  if (not widthRes.ok()) {
    return widthRes.error();
  }
  float width = *widthRes;
  std::pair<float, float> returnPair{maxZ, width};
  return returnPair;
}


template <int trkGridSize>
void
Acts::AdaptiveGridTrackDensity<trkGridSize>::addTrack(
    const Acts::BoundTrackParameters& trk,
    std::vector<float>& mainGridDensity,
    std::vector<int>& mainGridZValues) const {
  SymMatrix2D cov = trk.covariance()->block<2, 2>(0, 0);
  float d0 = trk.parameters()[0];
  float z0 = trk.parameters()[1];

  // Calculate offset in d direction to central bin at z-axis
  int dOffset = std::floor(d0 / m_cfg.binSize - 0.5) + 1;
  // Calculate bin in z
  int zBin = int(z0 / m_cfg.binSize);

  std::cout << "z0/zbin: " << z0 << "," << zBin << std::endl;

  // Calculate the positions of the bin centers
  float binCtrD = dOffset * m_cfg.binSize;

  int sign = (z0 > 0) ? +1 : -1;
  float binCtrZ = (zBin + sign * 0.5f) * m_cfg.binSize;

  std::cout << "binCtrZ: " << binCtrZ << std::endl;

  // Calculate the distance between IP values and their
  // corresponding bin centers
  float distCtrD = d0 - binCtrD;
  float distCtrZ = z0 - binCtrZ;

  // Check if current track does affect grid density
  // in central bins at z-axis
  if ((std::abs(dOffset) > trkGridSize - 1) / 2.) {
    return;
  }

  // Create the track grid
  ActsVectorF<trkGridSize> trackGrid =
      createTrackGrid(dOffset, cov, distCtrD, distCtrZ);

  std::vector<int> zBinValues;

  int startEnd = int(trkGridSize - 1)/2;

  std::cout << "current z values: " << std::endl;
  for(int i = 0; i < trkGridSize; i++){
    std::cout << "zbins: " << int(zBin + (i-startEnd)) << std::endl;
    zBinValues.push_back(int(zBin + (i-startEnd)));
  }

  std::cout << "mainGridZValues before: " <<std::endl;

  for(auto& bla : mainGridZValues){
    std::cout << bla << std::endl;
  }

  for(int i = 0; i < trkGridSize; i++){
    int z = zBinValues[i];
    
    std::cout << "testing z value: " << z << std::endl;

    auto findIter = std::find(mainGridZValues.begin(), mainGridZValues.end(), z);
    //bool exists = std::binary_search(mainGridZValues.begin(), mainGridZValues.end(), z);
    // TODO: wanted to check if element exists with findIter != mainGridZValues.end()
    // but if element is actually last element in list, caused problems maybe? Check if
    // results are always the same as if we're using bineary search here
    if(findIter != mainGridZValues.end()){
    // Z bin already exists
      std::cout << "z exists..." << std::endl;
      mainGridDensity[std::distance(mainGridZValues.begin(), findIter)] += trackGrid[i];
    }
    else{
      std::cout << "z did not exist..." << std::endl;
    // Create new z bin
      auto it = std::upper_bound( mainGridZValues.begin(), mainGridZValues.end(), z );
      mainGridDensity.insert(mainGridDensity.begin() + std::distance(mainGridZValues.begin(), it), trackGrid[i]);
      mainGridZValues.insert(it, z);
    } 
  }
}

// 
template <int trkGridSize>
Acts::ActsVectorF<trkGridSize>
Acts::AdaptiveGridTrackDensity<trkGridSize>::createTrackGrid(
    int offset, const Acts::SymMatrix2D& cov, float distCtrD,
    float distCtrZ) const {
  ActsVectorF<trkGridSize> trackGrid(ActsVectorF<trkGridSize>::Zero());

  int i = (trkGridSize - 1) / 2 + offset;
  float d = (i - static_cast<float>(trkGridSize) / 2 + 0.5f) * m_cfg.binSize;

  // Loop over columns
  for (int j = 0; j < trkGridSize; j++) {
    float z = (j - static_cast<float>(trkGridSize) / 2 + 0.5f) * m_cfg.binSize;
    trackGrid(j) = normal2D(d + distCtrD, z + distCtrZ, cov);
  }
  return trackGrid;
}

template <int trkGridSize>
Acts::Result<float>
Acts::AdaptiveGridTrackDensity<trkGridSize>::estimateSeedWidth(
    const std::vector<float>& mainGridDensity,
    const std::vector<int>& mainGridZValues, float maxZ) const {
  if (mainGridDensity.empty() || mainGridZValues.empty()) {
    return VertexingError::EmptyInput;
  }
  // Get z bin of max density z value
  int sign = (maxZ > 0) ? +1 : -1;
  int zMaxGridBin  = int(maxZ / m_cfg.binSize - sign * 0.5f);

  std::cout << "seed width estimation... max z bin: " << zMaxGridBin << std::endl; 

  // Find location in mainGridZValues
  auto findIter = std::find(mainGridZValues.begin(), mainGridZValues.end(), zMaxGridBin);
  int zBin = std::distance(mainGridZValues.begin(), findIter);

  const float maxValue = mainGridDensity[zBin];
  float gridValue = mainGridDensity[zBin];

  // Find right half-maximum bin
  int rhmBin = zBin;
  while (gridValue > maxValue / 2) {
    // Check if we are still operating on continous z values
    std::cout << "R:cont z values?: " << zMaxGridBin + (rhmBin - zBin) << "," << mainGridZValues[rhmBin] << std::endl;
    if( (zMaxGridBin + (rhmBin - zBin)) != mainGridZValues[rhmBin])
    {
      break;
    }
    rhmBin += 1;
    gridValue = mainGridDensity[rhmBin];
  }

  // Use linear approximation to find better z value for FWHM between bins
  float deltaZ1 = (maxValue / 2 - mainGridDensity[rhmBin - 1]) *
                  (m_cfg.binSize / (mainGridDensity[rhmBin - 1] - mainGridDensity[rhmBin]));
  // Find left half-maximum bin
  int lhmBin = zBin;
  gridValue = mainGridDensity[zBin];
  while (gridValue > maxValue / 2) {
    std::cout << "L:cont z values?: " << zMaxGridBin + (rhmBin - zBin) << "," << mainGridZValues[rhmBin] << std::endl;
    // Check if we are still operating on continous z values
    if( (zMaxGridBin + (rhmBin - zBin)) != mainGridZValues[rhmBin])
    {
      break;
    }
    lhmBin -= 1;
    gridValue = mainGridDensity[lhmBin];
  }

  // Use linear approximation to find better z value for FWHM between bins
  float deltaZ2 = (maxValue / 2 - mainGridDensity[lhmBin + 1]) *
                  (m_cfg.binSize / (mainGridDensity[rhmBin + 1] - mainGridDensity[rhmBin]));

  // Approximate FWHM
  float fwhm =
      rhmBin * m_cfg.binSize - deltaZ1 - lhmBin * m_cfg.binSize - deltaZ2;

  // FWHM = 2.355 * sigma
  float width = fwhm / 2.355f;

  return std::isnormal(width) ? width : 0.0f;
}


template <int trkGridSize>
float Acts::AdaptiveGridTrackDensity<trkGridSize>::normal2D(
    float d, float z, const Acts::SymMatrix2D& cov) const {
  float det = cov.determinant();
  float coef = 1 / (2 * M_PI * std::sqrt(det));
  float expo =
      -1 / (2 * det) *
      (cov(1, 1) * d * d - d * z * (cov(0, 1) + cov(1, 0)) + cov(0, 0) * z * z);
  return coef * std::exp(expo);
}

