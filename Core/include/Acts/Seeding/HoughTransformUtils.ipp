
template <class identifier_t> template <class PointType>
void Acts::HoughTransformUtils::HoughPlane<identifier_t>::fill(
    const PointType& measurement, const houghAxisRanges& axisRanges,
    lineParametrisation<PointType> linePar,
    lineParametrisation<PointType> widthPar, const identifier_t & identifier, unsigned layer,
    yieldType weight) {
  for (int xBin = 0; xBin < m_cfg.nBinsX; xBin++) {
    auto x = binCenter(axisRanges.xMin, axisRanges.xMax, m_cfg.nBinsX, xBin);
    // now evaluate the line equation provided by the user
    coordType y = linePar(x, measurement);
    coordType dy = widthPar(x, measurement);
    auto yBinDown =
        binIndex(axisRanges.yMin, axisRanges.yMax, m_cfg.nBinsY, y - dy);
    auto yBinUp =
        binIndex(axisRanges.yMin, axisRanges.yMax, m_cfg.nBinsY, y + dy);
    for (int yBin = yBinDown; yBin <= yBinUp; ++yBin) {
      if (yBin >= m_cfg.nBinsY || yBin < 0)
        continue;  // out of bounds
      fillBin(xBin, yBin, identifier, layer, weight);
    }
  }
}


template <class identifier_t> void Acts::HoughTransformUtils::HoughCell<identifier_t>::fill(const identifier_t & identifier, unsigned layer, yieldType weight) {
  if (m_hits.insert(identifier).second) {
    m_nHits += weight;
  }
  if (m_layers.insert(layer).second) {
    m_nLayers += weight;
  }
}
template <class identifier_t> 
void Acts::HoughTransformUtils::HoughCell<identifier_t>::reset() {
  if (m_nLayers > 0) m_layers.clear();
  if (m_nHits > 0) m_hits.clear();
  m_hits.clear();
  m_nHits = 0;
  m_nLayers = 0;
}

template <class identifier_t> Acts::HoughTransformUtils::HoughPlane<identifier_t>::HoughPlane(const HoughPlaneConfig& cfg) : m_cfg(cfg) {
  m_houghHist = HoughHist(m_cfg.nBinsX, m_cfg.nBinsY);
}
template <class identifier_t> void Acts::HoughTransformUtils::HoughPlane<identifier_t>::fillBin(int binX, int binY, const identifier_t & identifier, unsigned layer,
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

template <class identifier_t> void Acts::HoughTransformUtils::HoughPlane<identifier_t>::reset() {
  for (size_t i =0; i < m_cfg.nBinsX; ++i){
    for (size_t j =0; j < m_cfg.nBinsY;++j){
      m_houghHist(i,j).reset();
    }
  }
  m_maxHits = 0.;
  m_maxLayers = 0.;
}
// NYI
template <class identifier_t> void Acts::HoughTransformUtils::HoughPlane<identifier_t>::drawHoughHist(std::string const&) {}

template <class identifier_t> Acts::HoughTransformUtils::HoughPeakFinder_LayerGuidedCombinatoric<identifier_t>::
    HoughPeakFinder_LayerGuidedCombinatoric(
        const HoughPeakFinder_LayerGuidedCombinatoricConfig& cfg)
    : m_cfg(cfg) {}

template <class identifier_t> std::vector<typename Acts::HoughTransformUtils::HoughPeakFinder_LayerGuidedCombinatoric<identifier_t>::Maximum>
Acts::HoughTransformUtils::HoughPeakFinder_LayerGuidedCombinatoric<identifier_t>::findPeaks(
    const HoughPlane<identifier_t>& plane) const {
  std::vector<HoughPeakFinder_LayerGuidedCombinatoric<identifier_t>::Maximum> maxima;

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

template <class identifier_t> bool Acts::HoughTransformUtils::HoughPeakFinder_LayerGuidedCombinatoric<identifier_t>::passThreshold(
    const HoughPlane<identifier_t>& plane, int xBin, int yBin) const {
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
template <class identifier_t> std::vector<std::vector<int>>
Acts::HoughTransformUtils::HoughPeakFinder_LayerGuidedCombinatoric<identifier_t>::getComboIndices(
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

template <class identifier_t> Acts::HoughTransformUtils::HoughPeakFinder_IslandsAroundMax<identifier_t>::HoughPeakFinder_IslandsAroundMax(
    const HoughPeakFinder_IslandsAroundMaxConfig& cfg)
    : m_cfg(cfg) {}

template <class identifier_t> std::vector<typename Acts::HoughTransformUtils::HoughPeakFinder_IslandsAroundMax<identifier_t>::Maximum>
Acts::HoughTransformUtils::HoughPeakFinder_IslandsAroundMax<identifier_t>::findPeaks(const HoughPlane<identifier_t>& plane, const houghAxisRanges & ranges) {
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
      std::cout <<" adding a SP with x,y bins "<<s.first<<" "<<s.second<<std::endl;
      maximum.hitIdentifiers.insert(plane.hitIds(s.first, s.second).cbegin(),
                                plane.hitIds(s.first, s.second).cend());
      coordType xHit = binCenter(ranges.xMin, ranges.xMax, plane.nBinsX(), s.first);
      coordType yHit = binCenter(ranges.yMin, ranges.yMax, plane.nBinsY(), s.second);
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

template <class identifier_t> void Acts::HoughTransformUtils::HoughPeakFinder_IslandsAroundMax<identifier_t>::extendMaximum(
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
    int xnew = nextCand.first + step.first; 
    int ynew = nextCand.second + step.second;
    if ((xnew < 0 || xnew >=  (int)yieldMap.size(0)) 
        || (ynew < 0 || ynew >= (int) yieldMap.size(1))) continue; 
    if (yieldMap(xnew, ynew) >
        threshold) {
      toExplore.push_back(
          {xnew, ynew});
    }
  }
}
