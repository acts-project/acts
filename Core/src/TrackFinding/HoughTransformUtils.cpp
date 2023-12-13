#include "Acts/Seeding/HoughTransformUtils.hpp"
namespace Acts {
namespace HoughTransformUtils {

void HoughCell::fill(idType identifier, unsigned layer, yieldType weight) {
  if (m_hits.insert(identifier).second) {
    m_nHits += weight;
  }
  if (m_layers.insert(layer).second) {
    m_nLayers += weight;
  }
}

void HoughCell::reset() {
  m_layers.clear();
  m_hits.clear();
  m_nHits = 0;
  m_nLayers = 0;
}

HoughPlane::HoughPlane(const HoughPlane::Config& cfg) : m_cfg(cfg) {
  m_houghHist = HoughHist(m_cfg.nBinsX, m_cfg.nBinsY);
}
void HoughPlane::fillBin(int binX, int binY, idType identifier, unsigned layer,
                         double w) {
  m_houghHist(binX, binY).fill(identifier, layer, w);
  yieldType nLayers = m_houghHist(binX, binY).nLayers();
  yieldType nHits = m_houghHist(binX, binY).nHits();
  if (nLayers > m_maxLayers) {
    m_maxLayers = nLayers;
    m_maxLocLayers = {binX, binY};
  }
  if (nHits > m_maxHits) {
    m_maxHits = nHits;
    m_maxLocHits = {binX, binY};
  }
}

void HoughPlane::reset() {
  m_houghHist.reset(HoughCell{});
}
// NYI
void HoughPlane::drawHoughHist(std::string const&) {}

HoughPeakFinder_LayerGuidedCombinatoric::
    HoughPeakFinder_LayerGuidedCombinatoric(
        const HoughPeakFinder_LayerGuidedCombinatoric::Config& cfg)
    : m_cfg(cfg) {}

std::vector<HoughPeakFinder_LayerGuidedCombinatoric::Maximum>
HoughPeakFinder_LayerGuidedCombinatoric::findPeaks(
    const HoughPlane& plane) const {
  std::vector<HoughPeakFinder_LayerGuidedCombinatoric::Maximum> maxima;

  for (int y = 0; y < plane.nBinsY(); y++) {
    for (int x = 0; x < plane.nBinsX(); x++) {
      if (passThreshold(plane, x, y)) {
        Maximum max;
        max.hitIdentifiers = plane.hitIds(x, y);
        maxima.push_back(max);
      }
    }
  }
  return maxima;
}

bool HoughPeakFinder_LayerGuidedCombinatoric::passThreshold(
    const HoughPlane& plane, int xBin, int yBin) const {
  // Pass window threshold
  if (plane.nLayers(xBin, yBin) < m_cfg.threshold)
    return false;

  if (m_cfg.localMaxWindowSize == 0)
    return true;

  // Pass local-maximum check, if used
  for (int j = -m_cfg.localMaxWindowSize; j <= m_cfg.localMaxWindowSize; j++) {
    for (int i = -m_cfg.localMaxWindowSize; i <= m_cfg.localMaxWindowSize;
         i++) {
      if (i == 0 && j == 0) {
        continue;
      }
      // out of bounds -> no need to test this bin
      if (xBin + j >= plane.nBinsX() || yBin + i >= plane.nBinsY())
        continue;
      // we are not in a minimum
      if (plane.nLayers(xBin + j, yBin + i) > plane.nLayers(xBin, yBin)) {
        return false;
      }
      if (plane.nLayers(xBin + j, yBin + i) < plane.nLayers(xBin, yBin)) {
        continue;  // this neighbour has fewer hits -> safe to ignore
      }

      // we can still be in a plateau. In this case, resolve by counting hits

      // if the other bin with the same hit count has seen more clusters (==
      // double hits in one layer), keep that one and reject the current
      if (plane.nHits(xBin + j, yBin + i) > plane.nHits(xBin, yBin)) {
        return false;
      }
      // and if we have completely identical points, bias the minimum towards
      // the bottom (left)
      if (plane.nHits(xBin + j, yBin + i) == plane.nHits(xBin, yBin) &&
          (i < 0 || (i == 0 && j < 0))) {
        return false;
      }
    }
  }
  return true;
}
std::vector<std::vector<int>>
HoughPeakFinder_LayerGuidedCombinatoric::getComboIndices(
    std::vector<size_t>& sizes) const {
  size_t nCombs = 1;
  std::vector<size_t> nCombs_prior(sizes.size());
  std::vector<int> temp(sizes.size(), 0);

  for (size_t i = 0; i < sizes.size(); i++) {
    if (sizes[i] > 0) {
      nCombs_prior[i] = nCombs;
      nCombs *= sizes[i];
    } else {
      temp[i] = -1;
    }
  }

  std::vector<std::vector<int>> combos(nCombs, temp);

  for (size_t icomb = 0; icomb < nCombs; icomb++) {
    size_t index = icomb;
    for (size_t isize = sizes.size() - 1; isize < sizes.size(); isize--) {
      if (sizes[isize] == 0) {
        continue;
      }
      combos[icomb][isize] = static_cast<int>(index / nCombs_prior[isize]);
      index = index % nCombs_prior[isize];
    }
  }

  return combos;
}

HoughPeakFinder_IslandsAroundMax::HoughPeakFinder_IslandsAroundMax(
    const HoughPeakFinder_IslandsAroundMax::Config& cfg)
    : m_cfg(cfg) {}

std::vector<HoughPeakFinder_IslandsAroundMax::Maximum>
HoughPeakFinder_IslandsAroundMax::findPeaks(const HoughPlane& plane, const houghAxisRanges & ranges) {
  yieldType max = plane.maxHits();
  yieldType min = m_cfg.fractionCutoff * max;
  std::vector<std::pair<int, int>> candidates;
  std::vector<Maximum> maxima;
  vector2D<yieldType> yieldMap(plane.nBinsX(), plane.nBinsY());
  for (int x = 0; x < plane.nBinsX(); ++x) {
    for (int y = 0; y < plane.nBinsY(); ++y) {
      if (plane.nHits(x, y) > min) {
        candidates.push_back({x, y});
        yieldMap(x, y) = plane.nHits(x, y);
      }
    }
  }
  std::sort(candidates.begin(), candidates.end(),
            [&plane](const std::pair<int, int>& bin1,
                     const std::pair<int, int>& bin2) {
              return (plane.nHits(bin1.first, bin1.second) >
                      plane.nHits(bin2.first, bin2.second));
  });
  std::vector<std::pair<int, int>> toExplore;
  std::vector<std::pair<int, int>> solution;
  for (auto& cand : candidates) {
    if (yieldMap(cand.first, cand.second) < min)
      continue;
    toExplore.clear();
    solution.clear();
    toExplore.push_back(cand);
    while (!toExplore.empty()) {
      extendMaximum(solution, toExplore, min, yieldMap);
    }
    if (solution.empty())
      continue;
    Maximum maximum;
    coordType max_x = 0;
    coordType max_y = 0;
    coordType pos_den = 0;
    coordType ymax = 0;
    coordType ymin = std::numeric_limits<coordType>::max();
    coordType xmax = 0;
    coordType xmin = std::numeric_limits<coordType>::max();
    for (auto& s : solution) {
      maximum.hitIdentifiers.insert(plane.hitIds(s.first, s.second).cbegin(),
                                plane.hitIds(s.first, s.second).cend());
      coordType xHit = plane.binCenter(ranges.xMin, ranges.xMax, plane.nBinsX(), s.first);
      coordType yHit = plane.binCenter(ranges.yMin, ranges.yMax, plane.nBinsY(), s.second);
      max_x += xHit * plane.nHits(s.first, s.second);
      max_y += yHit * plane.nHits(s.first, s.second);
      pos_den += plane.nHits(s.first, s.second);
      if (xHit > xmax)
        xmax = xHit;
      if (yHit > ymax)
        ymax = yHit;
      if (xHit < xmin)
        xmin = xHit;
      if (yHit < ymin)
        ymin = yHit;
    }
    maximum.x = max_x / pos_den;
    maximum.y = max_y / pos_den;
    maximum.wx = 0.5 * (xmax - xmin);
    maximum.wy = 0.5 * (ymax - ymin);
    maxima.push_back(maximum);
  }
  return maxima;
}

void HoughPeakFinder_IslandsAroundMax::extendMaximum(
    std::vector<std::pair<int, int>>& inMaximum,
    std::vector<std::pair<int, int>>& toExplore, yieldType threshold,
    vector2D<yieldType>& yieldMap) {
  /// check the next candidate
  auto nextCand = toExplore.back();
  toExplore.pop_back();
  if (yieldMap(nextCand.first, nextCand.second) < threshold)
    return;
  // try to extend our candidate using neighbouring bins
  inMaximum.push_back(nextCand);
  // this hit is now listed as "used" - reset it on the yield map
  yieldMap(nextCand.first, nextCand.second) = -1.0f;
  for (auto& step : m_stepDirections) {
    if (yieldMap(nextCand.first + step.first, nextCand.second + step.second) >
        threshold) {
      toExplore.push_back(
          {nextCand.first + step.first, nextCand.second + step.second});
    }
  }
}

}  // namespace HoughTransformUtils
}  // namespace Acts
