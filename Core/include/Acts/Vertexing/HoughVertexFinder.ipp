// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

template <typename spacepoint_t>
Acts::HoughVertexFinder<spacepoint_t>::HoughVertexFinder(
    Config cfg, std::unique_ptr<const Logger> lgr)
    : m_cfg(std::move(cfg)), m_logger(std::move(lgr)) {
  if (m_cfg.absEtaFractions.size() != m_cfg.absEtaRanges.size()) {
    ACTS_ERROR("size of the absEtaFractions is "
               << m_cfg.absEtaFractions.size()
               << " but size of the absEtaRanges vector is "
               << m_cfg.absEtaRanges.size() << "; these two have to be equal.");
  }

  if (m_cfg.rangeIterZ.size() != m_cfg.nBinsZIterZ.size()) {
    ACTS_ERROR("size of the rangeIterZ is "
               << m_cfg.rangeIterZ.size()
               << " but size of the nBinsZIterZ vector is "
               << m_cfg.nBinsZIterZ.size() << "; these two have to be equal.");
  }

  if (m_cfg.rangeIterZ.size() != m_cfg.nBinsCotThetaIterZ.size()) {
    ACTS_ERROR("size of the rangeIterZ is "
               << m_cfg.rangeIterZ.size()
               << " but size of the nBinsCotThetaIterZ vector is "
               << m_cfg.nBinsCotThetaIterZ.size()
               << "; these two have to be equal.");
  }
}

template <typename spacepoint_t>
Acts::Result<Acts::Vertex> Acts::HoughVertexFinder<spacepoint_t>::find(
    const std::vector<spacepoint_t>& spacepoints) const {
  if (spacepoints.empty()) {
    return Acts::Result<Acts::Vertex>::failure(std::error_code());
  }

  float absEtaRange = m_cfg.maxAbsEta;
  float totalFrac = 0.;
  for (unsigned int r = 0; r < m_cfg.absEtaRanges.size(); ++r) {
    float addToTotalFrac = m_cfg.absEtaFractions.at(r);
    if ((totalFrac + addToTotalFrac) * spacepoints.size() > m_cfg.targetSPs) {
      float needOnly = (m_cfg.targetSPs - totalFrac * spacepoints.size()) /
                       (addToTotalFrac * spacepoints.size());
      absEtaRange = (r != 0u ? m_cfg.absEtaRanges.at(r - 1) : 0.) +
                    (m_cfg.absEtaRanges.at(r) -
                     (r != 0u ? m_cfg.absEtaRanges.at(r - 1) : 0.)) *
                        needOnly;
      totalFrac += needOnly * addToTotalFrac;
      break;
    }

    totalFrac += addToTotalFrac;
  }
  if (absEtaRange > m_cfg.maxAbsEta) {
    absEtaRange = m_cfg.maxAbsEta;
  }
  if (absEtaRange < m_cfg.minAbsEta) {
    absEtaRange = m_cfg.minAbsEta;
  }

  const float maxCotTheta = std::sinh(absEtaRange);
  const float minCotTheta = -1. * maxCotTheta;

  float binsNumDecrease = 1.;
  if (spacepoints.size() * totalFrac < m_cfg.targetSPs) {
    binsNumDecrease =
        std::pow(m_cfg.binsCotThetaDecrease,
                 std::log(m_cfg.targetSPs / (spacepoints.size() * totalFrac)));
  }

  Acts::Vertex vtx{m_cfg.defVtxPosition};
  for (std::uint32_t iter = 0; iter < m_cfg.rangeIterZ.size(); ++iter) {
    auto vtxNew = findHoughVertex(
        spacepoints, vtx, m_cfg.rangeIterZ.at(iter), m_cfg.nBinsZIterZ.at(iter),
        minCotTheta, maxCotTheta,
        static_cast<std::uint32_t>(m_cfg.nBinsCotThetaIterZ.at(iter) /
                                   binsNumDecrease));

    if (!vtxNew.ok()) {
      // vertex not found
      return Acts::Result<Acts::Vertex>::failure(std::error_code());
    }

    vtx = vtxNew.value();
  }

  return Acts::Result<Acts::Vertex>::success(vtx);
}

template <typename spacepoint_t>
Acts::Result<Acts::Vertex>
Acts::HoughVertexFinder<spacepoint_t>::findHoughVertex(
    const std::vector<spacepoint_t>& spacepoints, const Acts::Vertex& vtxOld,
    float rangeZ, std::uint32_t numZBins, float minCotTheta, float maxCotTheta,
    std::uint32_t numCotThetaBins) const {
  const float zBinSize = 2. * rangeZ / numZBins;
  const float invCotThetaBinSize =
      numCotThetaBins / (maxCotTheta - minCotTheta);
  const float minZ = vtxOld.position().z() - rangeZ,
              maxZ = vtxOld.position().z() + rangeZ;
  const float vtxOldX = vtxOld.position().x(), vtxOldY = vtxOld.position().y();

  HoughHist houghHist(HoughAxis(minZ, maxZ, numZBins),
                      HoughAxis(minCotTheta, maxCotTheta, numCotThetaBins));

  std::vector<std::uint32_t> houghZProjection(numZBins, 0);

  std::vector<float> vtxZPositions;
  for (std::uint32_t zBin = 0; zBin < numZBins; zBin++) {
    vtxZPositions.push_back(
        Acts::HoughTransformUtils::binCenter(minZ, maxZ, numZBins, zBin));
  }

  for (const auto& sp : spacepoints) {
    if (sp.z() > maxZ) {
      if ((sp.z() - maxZ) / sp.r() > maxCotTheta) {
        continue;
      }
    } else if (sp.z() < minZ) {
      if ((sp.z() - minZ) / sp.r() < minCotTheta) {
        continue;
      }
    }

    float sp_invr = 1. / std::sqrt((sp.x() - vtxOldX) * (sp.x() - vtxOldX) +
                                   (sp.y() - vtxOldY) * (sp.y() - vtxOldY));

    std::uint32_t zFrom = static_cast<std::uint32_t>(
        std::max(((sp.z() - maxCotTheta / sp_invr) - minZ) / zBinSize + 1, 0.));
    std::uint32_t zTo = static_cast<std::uint32_t>(std::min(
        ((sp.z() - minCotTheta / sp_invr) - minZ) / zBinSize, 1. * numZBins));

    for (std::uint32_t zBin = zFrom; zBin < zTo; zBin++) {
      float cotTheta = (sp.z() - vtxZPositions[zBin]) * sp_invr;

      std::uint32_t cotThetaBin = static_cast<std::uint32_t>(
          (cotTheta - minCotTheta) * invCotThetaBinSize);

      std::uint32_t cotThetaFrom =
          std::max<std::uint32_t>(cotThetaBin - m_cfg.fillNeighbours, 0);
      std::uint32_t cotThetaTo =
          std::min(cotThetaBin + m_cfg.fillNeighbours + 1, numCotThetaBins);

      for (std::uint32_t cotBin = cotThetaFrom; cotBin < cotThetaTo; ++cotBin) {
        ++houghHist.atLocalBins({zBin, cotBin});
      }
    }
  }

  for (std::uint32_t zBin = 0; zBin < numZBins; zBin++) {
    for (std::uint32_t cotBin = 0; cotBin < numCotThetaBins; ++cotBin) {
      auto rhs = houghHist.atLocalBins({zBin, cotBin});
      houghZProjection[zBin] += rhs * (rhs >= m_cfg.minHits);
    }
  }

  auto vtxNewZ = findHoughPeak(houghZProjection, vtxZPositions, numZBins);
  if (vtxNewZ.ok()) {
    Acts::Vertex newVertex(Vector3{vtxOldX, vtxOldY, vtxNewZ.value()});
    return Acts::Result<Acts::Vertex>::success(newVertex);
  }

  return Acts::Result<Acts::Vertex>::failure(std::error_code());
}

template <typename spacepoint_t>
Acts::Result<float> Acts::HoughVertexFinder<spacepoint_t>::findHoughPeak(
    const std::vector<std::uint32_t>& houghZProjection,
    const std::vector<float>& vtxZPositions, std::uint32_t numZBins) const {
  std::uint32_t maxZBin = 0;
  for (std::uint32_t zBin = 0; zBin < numZBins; zBin++) {
    if (houghZProjection[zBin] > houghZProjection.at(maxZBin)) {
      maxZBin = zBin;
    }
  }

  float avg =
      std::accumulate(houghZProjection.begin(), houghZProjection.end(), 0.) /
      houghZProjection.size();
  float sumEntries = 0;
  float meanZPeak = 0.;

  for (std::uint32_t zBin =
           std::max<std::uint32_t>(maxZBin - m_cfg.peakWidth, 0);
       zBin <= std::min(numZBins - 1, maxZBin + m_cfg.peakWidth); ++zBin) {
    sumEntries += std::max(houghZProjection.at(zBin) - avg, 0.f);
    meanZPeak +=
        vtxZPositions[zBin] * std::max(houghZProjection.at(zBin) - avg, 0.f);
  }

  if (sumEntries != 0.) {
    meanZPeak /= sumEntries;
    return Acts::Result<float>::success(meanZPeak);
  }

  // vertex not found; Hough image empty
  return Acts::Result<float>::failure(std::error_code());
}
