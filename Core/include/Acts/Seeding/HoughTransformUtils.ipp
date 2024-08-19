// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <algorithm>

template <class identifier_t>
template <class PointType>
void Acts::HoughTransformUtils::HoughPlane<identifier_t>::fill(
    const PointType& measurement, const HoughAxisRanges& axisRanges,
    LineParametrisation<PointType> linePar,
    LineParametrisation<PointType> widthPar, const identifier_t& identifier,
    unsigned layer, YieldType weight) {
  // loop over all bins in the first coordinate to populate the line
  for (std::size_t xBin = 0; xBin < m_cfg.nBinsX; xBin++) {
    // get the x-coordinate for the given bin
    auto x = binCenter(axisRanges.xMin, axisRanges.xMax, m_cfg.nBinsX, xBin);
    // now evaluate the line equation provided by the user
    CoordType y = linePar(x, measurement);
    CoordType dy = widthPar(x, measurement);
    // translate the y-coordinate range to a bin range
    int yBinDown =
        binIndex(axisRanges.yMin, axisRanges.yMax, m_cfg.nBinsY, y - dy);
    int yBinUp =
        binIndex(axisRanges.yMin, axisRanges.yMax, m_cfg.nBinsY, y + dy);
    // now we can loop over the bin range to fill the corresponding cells
    for (int yBin = yBinDown; yBin <= yBinUp; ++yBin) {
      // skip 'out of bounds' cases
      if (yBin >= static_cast<int>(m_cfg.nBinsY) || yBin < 0) {
        continue;
      }
      fillBin(xBin, yBin, identifier, layer, weight);
    }
  }
}

template <class identifier_t>
void Acts::HoughTransformUtils::HoughCell<identifier_t>::fill(
    const identifier_t& identifier, unsigned layer, YieldType weight) {
  // add the hit to the list of hits in the cell
  if (m_hits.insert(identifier).second) {
    // if new, increment the weighted hit counter
    m_nHits += weight;
  }
  // and to the same for layer counts
  if (m_layers.insert(layer).second) {
    m_nLayers += weight;
  }
}
template <class identifier_t>
void Acts::HoughTransformUtils::HoughCell<identifier_t>::reset() {
  // avoid expensive clear calls on empty cells
  if (m_nLayers > 0) {
    m_layers.clear();
  }
  if (m_nHits > 0) {
    m_hits.clear();
  }
  m_nHits = 0;
  m_nLayers = 0;
}

template <class identifier_t>
Acts::HoughTransformUtils::HoughPlane<identifier_t>::HoughPlane(
    const HoughPlaneConfig& cfg)
    : m_cfg(cfg) {
  // instantiate our histogram.
  m_houghHist = HoughHist(m_cfg.nBinsX, m_cfg.nBinsY);
}
template <class identifier_t>
void Acts::HoughTransformUtils::HoughPlane<identifier_t>::fillBin(
    std::size_t binX, std::size_t binY, const identifier_t& identifier,
    unsigned layer, double w) {
  // mark that this bin was filled with non trivial content
  m_touchedBins.insert({binX, binY});
  // add content to the cell
  m_houghHist(binX, binY).fill(identifier, layer, w);
  // and update our cached maxima
  YieldType nLayers = m_houghHist(binX, binY).nLayers();
  YieldType nHits = m_houghHist(binX, binY).nHits();
  if (nLayers > m_maxLayers) {
    m_maxLayers = nLayers;
    m_maxLocLayers = {binX, binY};
  }
  if (nHits > m_maxHits) {
    m_maxHits = nHits;
    m_maxLocHits = {binX, binY};
  }
}

template <class identifier_t>
void Acts::HoughTransformUtils::HoughPlane<identifier_t>::reset() {
  // reset all bins that were previously filled
  // avoid calling this on empty cells to save time
  for (auto bin : m_touchedBins) {
    m_houghHist(bin).reset();
  }
  // don't forget to reset our cached maxima
  m_maxHits = 0.;
  m_maxLayers = 0.;
  // and reset the list of nontrivial bins
  m_touchedBins.clear();
}

template <class identifier_t>
Acts::HoughTransformUtils::PeakFinders::LayerGuidedCombinatoric<identifier_t>::
    LayerGuidedCombinatoric(const LayerGuidedCombinatoricConfig& cfg)
    : m_cfg(cfg) {}

template <class identifier_t>
std::vector<typename Acts::HoughTransformUtils::PeakFinders::
                LayerGuidedCombinatoric<identifier_t>::Maximum>
Acts::HoughTransformUtils::PeakFinders::LayerGuidedCombinatoric<
    identifier_t>::findPeaks(const HoughPlane<identifier_t>& plane) const {
  // book the vector for the maxima
  std::vector<PeakFinders::LayerGuidedCombinatoric<identifier_t>::Maximum>
      maxima;
  // loop over the non empty bins
  for (auto [x, y] : plane.getNonEmptyBins()) {
    // and look for the ones that represent a maximum
    if (passThreshold(plane, x, y)) {
      // write out a maximum
      Maximum max;
      max.hitIdentifiers = plane.hitIds(x, y);
      maxima.push_back(max);
    }
  }
  return maxima;
}

template <class identifier_t>
bool Acts::HoughTransformUtils::PeakFinders::LayerGuidedCombinatoric<
    identifier_t>::passThreshold(const HoughPlane<identifier_t>& plane,
                                 std::size_t xBin, std::size_t yBin) const {
  // Check if we have sufficient layers for a maximum
  if (plane.nLayers(xBin, yBin) < m_cfg.threshold) {
    return false;
  }
  // if no local maximum is required, this is enough to classify as a max
  if (m_cfg.localMaxWindowSize == 0) {
    return true;
  }
  // otherwise, check for local maxima by looking at the surrounding cells

  /// loop over neighbours
  for (int dXbin = -m_cfg.localMaxWindowSize; dXbin <= m_cfg.localMaxWindowSize;
       dXbin++) {
    for (int dYbin = -m_cfg.localMaxWindowSize;
         dYbin <= m_cfg.localMaxWindowSize; dYbin++) {
      // exclude the candidate itself
      if (dYbin == 0 && dXbin == 0) {
        continue;
      }
      // out of bounds -> no need to test this bin
      if (static_cast<int>(xBin) + dXbin < 0 ||
          static_cast<int>(yBin) + dYbin < 0) {
        continue;
      }
      if (xBin + dXbin >= plane.nBinsX() || yBin + dYbin >= plane.nBinsY()) {
        continue;
      }
      // if a neighbour with more layers exist, this is not a minimum
      if (plane.nLayers(xBin + dXbin, yBin + dYbin) >
          plane.nLayers(xBin, yBin)) {
        return false;
      }
      // if the neighbour has fewer hits, we can move to the next neighbour
      if (plane.nLayers(xBin + dXbin, yBin + dYbin) <
          plane.nLayers(xBin, yBin)) {
        continue;
      }

      // we can still be in a plateau. In this case, resolve by counting hits

      // if the other bin with the same hit count has seen more clusters (
      // double hits in one layer), keep that one and reject the current
      if (plane.nHits(xBin + dXbin, yBin + dYbin) > plane.nHits(xBin, yBin)) {
        return false;
      }
      // and if we have completely identical hit and layer count, prefer bins in
      // the bottom (left) direction
      if (plane.nHits(xBin + dXbin, yBin + dYbin) == plane.nHits(xBin, yBin) &&
          (dYbin < 0 || (dYbin == 0 && dXbin < 0))) {
        return false;
      }
    }
  }
  return true;
}
template <class identifier_t>
Acts::HoughTransformUtils::PeakFinders::IslandsAroundMax<
    identifier_t>::IslandsAroundMax(const IslandsAroundMaxConfig& cfg)
    : m_cfg(cfg) {}

template <class identifier_t>
std::vector<typename Acts::HoughTransformUtils::PeakFinders::IslandsAroundMax<
    identifier_t>::Maximum>
Acts::HoughTransformUtils::PeakFinders::IslandsAroundMax<
    identifier_t>::findPeaks(const HoughPlane<identifier_t>& plane,
                             const HoughAxisRanges& ranges) {
  // check the global maximum hit count in the plane
  YieldType max = plane.maxHits();
  // and obtain the fraction of the max that is our cutoff for island formation
  YieldType min = std::max(m_cfg.threshold, m_cfg.fractionCutoff * max);
  // book a list for the candidates and the maxima
  std::vector<std::pair<std::size_t, std::size_t>> candidates;
  std::vector<Maximum> maxima;
  // keep track of the yields in each non empty cell
  std::map<std::pair<std::size_t, std::size_t>, YieldType> yieldMap;

  // now loop over all non empty bins
  for (auto [x, y] : plane.getNonEmptyBins()) {
    // we only consider cells above threshold from now on.
    // each is a potential candidate to seed or join an island
    // We also add each to our yield map
    if (plane.nHits(x, y) > min) {
      candidates.push_back({x, y});
      yieldMap[{x, y}] = plane.nHits(x, y);
    }
  }
  // sort the candidate cells descending in content
  std::sort(candidates.begin(), candidates.end(),
            [&plane](const std::pair<std::size_t, std::size_t>& bin1,
                     const std::pair<std::size_t, std::size_t>& bin2) {
              return (plane.nHits(bin1.first, bin1.second) >
                      plane.nHits(bin2.first, bin2.second));
            });

  // now we build islands from the candidate cells, starting with the most
  // populated one
  std::vector<std::pair<std::size_t, std::size_t>> toExplore;
  std::vector<std::pair<std::size_t, std::size_t>> solution;

  // loop over candidate cells
  for (auto& cand : candidates) {
    // check the content again (if used in a previous island, this will now
    // report empty)
    if (yieldMap[cand] < min) {
      continue;
    }
    // translate to parameter space for overlap veto
    CoordType xCand =
        binCenter(ranges.xMin, ranges.xMax, plane.nBinsX(), cand.first);
    CoordType yCand =
        binCenter(ranges.yMin, ranges.yMax, plane.nBinsY(), cand.second);
    // check if we are too close to a previously found maximum
    bool goodSpacing = true;
    for (auto& found : maxima) {
      if (std::abs(xCand - found.x) < m_cfg.minSpacingBetweenPeaks.first ||
          std::abs(yCand - found.y) < m_cfg.minSpacingBetweenPeaks.second) {
        goodSpacing = false;
        break;
      }
    }
    if (!goodSpacing) {
      continue;
    }
    // now we can extend the candidate into an island
    toExplore.clear();
    solution.clear();
    // initially, pass the candidate into the list of "neighbours" to explore
    // around an empty island
    toExplore.push_back(cand);
    // and incrementally add neighbours, filling the solution vector
    while (!toExplore.empty()) {
      extendMaximum(solution, toExplore, min, yieldMap);
    }
    // nothing found? Next candidate!
    if (solution.empty()) {
      continue;
    }
    // We found an island
    Maximum maximum;
    CoordType max_x = 0;
    CoordType max_y = 0;
    CoordType pos_den = 0;
    CoordType ymax = -std::numeric_limits<CoordType>::max();
    CoordType ymin = std::numeric_limits<CoordType>::max();
    CoordType xmax = -std::numeric_limits<CoordType>::max();
    CoordType xmin = std::numeric_limits<CoordType>::max();
    // loop over cells in the island and get the weighted mean position.
    // Also collect all hit identifiers in the island and the maximum
    // extent (within the count threshold!) of the island
    for (auto& [xBin, yBin] : solution) {
      const auto& hidIds = plane.hitIds(xBin, yBin);
      maximum.hitIdentifiers.insert(hidIds.cbegin(), hidIds.cend());
      CoordType xHit =
          binCenter(ranges.xMin, ranges.xMax, plane.nBinsX(), xBin);
      CoordType yHit =
          binCenter(ranges.yMin, ranges.yMax, plane.nBinsY(), yBin);
      YieldType nHits = plane.nHits(xBin, yBin);
      max_x += xHit * nHits;
      max_y += yHit * nHits;
      pos_den += nHits;
      if (xHit > xmax) {
        xmax = xHit;
      }
      if (yHit > ymax) {
        ymax = yHit;
      }
      if (xHit < xmin) {
        xmin = xHit;
      }
      if (yHit < ymin) {
        ymin = yHit;
      }
    }
    // calculate mean position
    maximum.x = max_x / pos_den;
    maximum.y = max_y / pos_den;
    // calculate width as symmetrised band
    maximum.wx = 0.5 * (xmax - xmin);
    maximum.wy = 0.5 * (ymax - ymin);
    maxima.push_back(maximum);
  }
  return maxima;
}

template <class identifier_t>
void Acts::HoughTransformUtils::PeakFinders::IslandsAroundMax<identifier_t>::
    extendMaximum(
        std::vector<std::pair<std::size_t, std::size_t>>& inMaximum,
        std::vector<std::pair<std::size_t, std::size_t>>& toExplore,
        YieldType threshold,
        std::map<std::pair<std::size_t, std::size_t>, YieldType>& yieldMap) {
  // in this call, we explore the last element of the toExplore list.
  // Fetch it and pop it from the vector.
  auto nextCand = toExplore.back();
  toExplore.pop_back();
  // check if we are above threshold. Don't add this cell to the island if not
  if (yieldMap[nextCand] < threshold) {
    return;
  }
  // This candidate is above threshold and should go on the island!

  // add it to the cell list for the island
  inMaximum.push_back(nextCand);
  // and "veto" the hit for further use via the yield map
  yieldMap[nextCand] = -1.0f;

  // now we have to collect the non empty neighbours of this cell and check them
  // as well
  for (auto step : m_stepDirections) {
    auto newCand = std::make_pair(nextCand.first + step.first,
                                  nextCand.second + step.second);
    // no need for bounds-checking, as the map structure will default to "empty"
    // yields for OOB

    // if the cell is above threshold, add it to our list of neighbours to
    // explore
    if (yieldMap[newCand] > threshold) {
      toExplore.push_back(newCand);
    }
  }
}
