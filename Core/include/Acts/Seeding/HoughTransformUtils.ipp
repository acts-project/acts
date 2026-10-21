// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/HoughTransformUtils.hpp"

#include <algorithm>
#include <tuple>

template <class identifier_t>
template <class PointType>
void Acts::HoughTransformUtils::HoughPlane<identifier_t>::fill(
    const PointType& measurement, const HoughAxisRanges& axisRanges,
    const LineParametrisation<PointType>& linePar,
    const LineParametrisation<PointType>& widthPar,
    const identifier_t& identifier, unsigned layer, YieldType weight) {
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
std::span<const identifier_t, std::dynamic_extent>
Acts::HoughTransformUtils::HoughCell<identifier_t>::getHits() const {
  std::span<const identifier_t, std::dynamic_extent> hits(m_hits.begin(),
                                                          m_iHit);
  return hits;
}

template <class identifier_t>
std::span<const unsigned, std::dynamic_extent>
Acts::HoughTransformUtils::HoughCell<identifier_t>::getLayers() const {
  std::span<const unsigned, std::dynamic_extent> layers(m_layers.begin(),
                                                        m_iLayer);

  return layers;
}

template <class identifier_t>
void Acts::HoughTransformUtils::HoughCell<identifier_t>::fill(
    const identifier_t& identifier, unsigned layer, YieldType weight) {
  // add the hit to the list of hits in the cell

  if (m_iHit != 0 && m_hits[m_iHit - 1] == identifier) {
    return;
  }

  if (m_iHit == m_hits.size()) {
    m_hits.resize(m_hits.size() + m_assignBatch);
  }

  m_hits[m_iHit] = identifier;
  m_iHit += 1;
  m_nHits += weight;

  // Check for duplicate layers
  if (m_iLayer != 0 && m_layers[m_iLayer - 1] == layer) {
    return;
  }

  if (m_iLayer == m_layers.size()) {
    m_layers.resize(m_layers.size() + m_assignBatch);
  }

  m_layers[m_iLayer] = layer;
  m_iLayer += 1;
  m_nLayers += weight;
}

template <class identifier_t>
void Acts::HoughTransformUtils::HoughCell<identifier_t>::reset() {
  m_nHits = 0;
  m_nLayers = 0;

  m_iLayer = 0;
  m_iHit = 0;
}

template <class identifier_t>
Acts::HoughTransformUtils::HoughPlane<identifier_t>::HoughPlane(
    const HoughPlaneConfig& cfg)
    : m_cfg(cfg),
      m_houghHist(Axis(0, m_cfg.nBinsX, m_cfg.nBinsX),
                  Axis(0, m_cfg.nBinsY, m_cfg.nBinsY)) {
  // instantiate our histogram.
  // m_houghHist = HoughHist(Axis(0, m_cfg.nBinsX, m_cfg.nBinsX), Axis(0,
  // m_cfg.nBinsY, m_cfg.nBinsY));
}
template <class identifier_t>
void Acts::HoughTransformUtils::HoughPlane<identifier_t>::fillBin(
    std::size_t binX, std::size_t binY, const identifier_t& identifier,
    unsigned layer, double w) {
  // mark that this bin was filled with non trivial content

  m_touchedBins.insert(globalBin({binX, binY}));

  // add content to the cell
  m_houghHist.atLocalBins({binX, binY}).fill(identifier, layer, w);
  // and update our cached maxima
  YieldType layers = nLayers(binX, binY);
  YieldType hits = nHits(binX, binY);
  if (layers > m_maxLayers) {
    m_maxLayers = layers;
    m_maxLocLayers = {binX, binY};
  }
  if (hits > m_maxHits) {
    m_maxHits = hits;
    m_maxLocHits = {binX, binY};
  }
}

template <class identifier_t>
void Acts::HoughTransformUtils::HoughPlane<identifier_t>::reset() {
  // reset all bins that were previously filled
  // avoid calling this on empty cells to save time
  for (auto bin : getNonEmptyBins()) {
    m_houghHist.at(bin).reset();
  }
  // don't forget to reset our cached maxima
  m_maxHits = 0.;
  m_maxLayers = 0.;
  // and reset the list of nontrivial bins
  m_iBin = 0;
  m_touchedBins.clear();
}
template <class identifier_t>
void Acts::HoughTransformUtils::HoughPlane<identifier_t>::checkIndices(
    std::size_t xBin, std::size_t yBin) const {
  if (xBin >= nBinsX()) {
    throw std::out_of_range("When accessing HoughPlane, X index " +
                            std::to_string(xBin) +
                            " is >= " + std::to_string(nBinsX()));
  }
  if (yBin >= nBinsY()) {
    throw std::out_of_range("When accessing HoughPlane, Y index " +
                            std::to_string(yBin) +
                            " is >= " + std::to_string(nBinsY()));
  }
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
  for (auto nbin : plane.getNonEmptyBins()) {
    // and look for the ones that represent a maximum

    std::array<std::size_t, 2> xy = plane.axisBins(nbin);

    if (passThreshold(plane, xy[0], xy[1])) {
      // write out a maximum
      Maximum max;
      auto hitIds = plane.hitIds(xy[0], xy[1]);
      max.hitIdentifiers.insert(std::make_move_iterator(hitIds.begin()),
                                std::make_move_iterator(hitIds.end()));
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
  const std::unordered_set<std::size_t>& nonEmptyBins{plane.getNonEmptyBins()};
  std::vector<std::size_t> candidates;
  candidates.reserve(nonEmptyBins.size());
  std::vector<Maximum> maxima;
  maxima.reserve(nonEmptyBins.size());
  // keep track of the yields in each non empty cell
  std::unordered_map<std::size_t, YieldType> yieldMap;

  // now loop over all non empty bins
  for (const std::size_t nbin : nonEmptyBins) {
    //  we only consider cells above threshold from now on.
    //  each is a potential candidate to seed or join an island
    //  We also add each to our yield map

    if (plane.nHits(nbin) > min) {
      candidates.push_back(nbin);
      yieldMap[nbin] = plane.nHits(nbin);
    }
  }
  // sort the candidate cells descending in content
  std::ranges::sort(candidates, std::greater{}, [&yieldMap](std::size_t c) {
    return std::make_tuple(yieldMap[c], c);
  });

  // now we build islands from the candidate cells, starting with the most
  // populated one
  std::vector<std::size_t> toExplore;
  std::vector<std::array<std::size_t, 2>> solution;

  // loop over candidate cells
  for (auto& cand : candidates) {
    // check the content again (if used in a previous island, this will now
    // report empty)
    if (yieldMap[cand] < min) {
      continue;
    }
    // CALL AXIS BINS HERE
    std::array<std::size_t, 2> xy = plane.axisBins(cand);
    // translate to parameter space for overlap veto
    CoordType xCand =
        binCenter(ranges.xMin, ranges.xMax, plane.nBinsX(), xy[0]);
    CoordType yCand =
        binCenter(ranges.yMin, ranges.yMax, plane.nBinsY(), xy[1]);
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
      extendMaximum(plane, solution, toExplore, min, yieldMap);
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
    std::vector<identifier_t> allHits;
    for (auto& [xBin, yBin] : solution) {
      auto hidIds = plane.hitIds(xBin, yBin);
      allHits.insert(allHits.end(), std::make_move_iterator(hidIds.begin()),
                     std::make_move_iterator(hidIds.end()));
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
    maximum.hitIdentifiers.insert(std::make_move_iterator(allHits.begin()),
                                  std::make_move_iterator(allHits.end()));
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
        const Acts::HoughTransformUtils::HoughPlane<identifier_t>& houghPlane,
        std::vector<std::array<std::size_t, 2>>& inMaximum,
        std::vector<std::size_t>& toExplore, YieldType threshold,
        std::unordered_map<std::size_t, YieldType>& yieldMap) {
  // in this call, we explore the last element of the toExplore list.
  // Fetch it and pop it from the vector.

  const std::size_t candidate = toExplore.back();
  YieldType& yield{yieldMap[candidate]};
  toExplore.pop_back();
  // check if we are above threshold. Don't add this cell to the island if not
  if (yield < threshold) {
    return;
  }
  // This candidate is above threshold and should go on the island!
  auto nextCand = houghPlane.axisBins(candidate);
  // add it to the cell list for the island
  inMaximum.push_back(nextCand);
  // and "veto" the hit for further use via the yield map
  yieldMap[candidate] = -1.0f;

  // now we have to collect the non empty neighbours of this cell and check them
  // as well
  for (auto step : m_stepDirections) {
    std::array<unsigned long, 2> newCandxy = {nextCand[0] + step.first,
                                              nextCand[1] + step.second};
    // check if we are moving out of the bounds
    if (newCandxy[0] >= houghPlane.nBinsX() ||
        newCandxy[1] >= houghPlane.nBinsY()) {
      continue;
    }
    std::size_t newCand = houghPlane.globalBin({newCandxy[0], newCandxy[1]});

    // if the cell is above threshold, add it to our list of neighbours to
    // explore
    if (yieldMap[newCand] > threshold) {
      toExplore.push_back(newCand);
    }
  }
}

namespace {

template <typename identifier_t>
bool passWindow(
    const typename Acts::HoughTransformUtils::HoughPlane<identifier_t>::Index
        index,
    const Acts::HoughTransformUtils::HoughPlane<identifier_t>& plane,
    const Acts::HoughTransformUtils::PeakFinders::SlidingWindowConfig& config) {
  using YieldType = Acts::HoughTransformUtils::YieldType;
  auto [xmax, ymax] = index;
  const YieldType max = plane.nHits(xmax, ymax);
  // window loopdefininf
  // this loop needs to be smarter to take care of wrapping
  for (std::size_t x = xmax - config.xWindowSize;
       x <= xmax + config.xWindowSize; ++x) {
    for (std::size_t y = ymax - config.yWindowSize;
         y <= ymax + config.yWindowSize; ++y) {
      const float xdist = static_cast<float>(x) - static_cast<float>(xmax);
      const float ydist = static_cast<float>(y) - static_cast<float>(ymax);
      const bool above = ydist + 0.1 > xdist;  // upper right corner
      const YieldType numOfHits = plane.nHits(x, y);

      if (above && numOfHits > max) {
        return false;
      }

      if (!above && numOfHits >= max) {
        return false;
      }
    }
  }
  return true;
}

template <typename identifier_t>
typename Acts::HoughTransformUtils::HoughPlane<identifier_t>::Index
slidingWindowRecenter(
    const typename Acts::HoughTransformUtils::HoughPlane<identifier_t>::Index
        index,
    const Acts::HoughTransformUtils::HoughPlane<identifier_t>& plane,
    const Acts::HoughTransformUtils::PeakFinders::SlidingWindowConfig& config) {
  using YieldType = Acts::HoughTransformUtils::YieldType;
  YieldType xw = 0;
  YieldType yw = 0;
  YieldType tot = 0;
  auto [xcenter, ycenter] = index;
  YieldType maxValue = plane.nHits(xcenter, ycenter);
  for (std::size_t x = xcenter - config.xRecenterSize;
       x <= xcenter + config.xRecenterSize; ++x) {
    for (std::size_t y = ycenter - config.yRecenterSize;
         y <= ycenter + config.yRecenterSize; ++y) {
      const YieldType numOfHits = plane.nHits(x, y);
      if (numOfHits >= maxValue) {
        xw += x * numOfHits;
        yw += y * numOfHits;
        tot += numOfHits;
      }
    }
  }
  return {static_cast<std::size_t>(xw / tot),
          static_cast<std::size_t>(yw / tot)};
}

}  // namespace

template <typename identifier_t>
std::vector<typename Acts::HoughTransformUtils::HoughPlane<identifier_t>::Index>
Acts::HoughTransformUtils::PeakFinders::slidingWindowPeaks(
    const Acts::HoughTransformUtils::HoughPlane<identifier_t>& plane,
    const Acts::HoughTransformUtils::PeakFinders::SlidingWindowConfig& config) {
  using IndexType = Acts::HoughTransformUtils::HoughPlane<identifier_t>::Index;
  std::vector<IndexType> output;

  for (std::size_t x = config.xWindowSize;
       x < plane.nBinsX() - config.xWindowSize; ++x) {
    for (std::size_t y = config.yWindowSize;
         y < plane.nBinsY() - config.yWindowSize; ++y) {
      if (plane.nHits(x, y) >= config.threshold) {
        const IndexType index({x, y});
        if (passWindow(index, plane, config)) {
          output.push_back(config.recenter
                               ? slidingWindowRecenter(index, plane, config)
                               : index);
        }
      }
    }
  }
  return output;
}

template <typename identifier_t, typename pixel_value_t>
std::vector<pixel_value_t>
Acts::HoughTransformUtils::PeakFinders::hitsCountImage(
    const Acts::HoughTransformUtils::HoughPlane<identifier_t>& plane,
    typename Acts::HoughTransformUtils::HoughPlane<identifier_t>::Index index,
    std::size_t xSize, std::size_t ySize,
    const std::function<pixel_value_t(const HoughPlane<identifier_t>&, int,
                                      int)>& summaryFunction) {
  std::vector<pixel_value_t> output;
  output.reserve((xSize) * (ySize));

  auto [xcenter, ycenter] = index;

  auto isInside = [&plane](int x, int y) -> bool {
    return 0 <= x && x < static_cast<int>(plane.nBinsX()) && 0 <= y &&
           y < static_cast<int>(plane.nBinsY());
  };

  for (std::size_t x = 0; x < xSize; ++x) {
    for (std::size_t y = 0; y < ySize; ++y) {
      int xPlaneCoord = x + xcenter - xSize / 2;
      int yPlaneCoord = y + ycenter - ySize / 2;
      output.push_back(isInside(xPlaneCoord, yPlaneCoord)
                           ? summaryFunction(plane, xPlaneCoord, yPlaneCoord)
                           : pixel_value_t{});
    }
  }
  return output;
}
