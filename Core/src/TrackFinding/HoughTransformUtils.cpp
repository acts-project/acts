#include "Acts/Seeding/HoughTransformUtils.hpp"
namespace Acts {
namespace HoughTransformUtils {

HoughPlane::HoughPlane(HoughPlaneConfig& cfg) : m_cfg(cfg) {
  for (int k = 0; k < m_cfg.nBinsX; ++k) {
    double low = unquant(m_cfg.xMin, m_cfg.xMax, m_cfg.nBinsX, k);
    double up = m_cfg.xMax;
    if (k < m_cfg.nBinsX - 1) {
      up = unquant(m_cfg.xMin, m_cfg.xMax, m_cfg.nBinsX, k + 1);
    }
    m_binCentersX.push_back(0.5 * (up + low));
  }
  for (int k = 0; k < m_cfg.nBinsY; ++k) {
    double low = unquant(m_cfg.yMin, m_cfg.yMax, m_cfg.nBinsY, k);
    double up = m_cfg.xMax;
    if (k < m_cfg.nBinsY - 1) {
      up = unquant(m_cfg.yMin, m_cfg.yMax, m_cfg.nBinsY, k + 1);
    }
    m_binCentersY.push_back(0.5 * (up + low));
  }
  m_layerHists.clear();
  m_layerHists.reserve(m_cfg.nLayers); 
  /// vector2D not compatible with resize()
  while(m_layerHists.size() < m_cfg.nLayers){
    m_layerHists.push_back(HoughHist(m_cfg.nBinsX, m_cfg.nBinsY));
  }
  m_totalHoughHist = HoughHist(m_cfg.nBinsX, m_cfg.nBinsY);
}
void HoughPlane::mergeLayers() {
  if (m_mergedLayers)
    return;
  for (int y = 0; y < m_cfg.nBinsY; y++) {
    for (int x = 0; x < m_cfg.nBinsX; x++) {
      m_totalHoughHist(x, y).first = 0;
      m_totalHoughHist(x, y).second.clear();
      for (auto& layer : m_layerHists) {
        if (layer(x,y).first > 0) {
          m_totalHoughHist(x,y).first += 1;
          m_totalHoughHist(x,y).second.insert(layer(x,y).second.begin(),
                                               layer(x,y).second.end());
        }
      }
    }
    m_mergedLayers = true;
  }
}
std::vector<HoughMaximum> HoughPlane::getMaxima() {
  if (!m_mergedLayers)
    mergeLayers();
  std::vector<HoughMaximum> maxima;
  for (int yBin = 0; yBin < m_cfg.nBinsY; yBin++) {
    for (int xBin = 0; xBin < m_cfg.nBinsX; xBin++) {
      if (passThreshold(m_totalHoughHist, xBin, yBin)) {
        HoughMaximum maximum;
        maximum.x = m_binCentersX.at(xBin);
        maximum.y = m_binCentersY.at(yBin);
        // find the size of our maximum
        const static std::array<std::pair<int, int>, 8> stepDirections{
            std::make_pair(-1, -1), std::make_pair(0, -1),
            std::make_pair(1, -1),  std::make_pair(-1, 0),
            std::make_pair(1, 0),   std::make_pair(-1, 1),
            std::make_pair(0, 1),   std::make_pair(1, 1)};
        std::set<std::pair<int, int>> explored;
        std::vector<std::pair<int, int>> toSearch{std::make_pair(xBin, yBin)};
        int xLeft = xBin;
        int xRight = xBin;
        int yBottom = yBin;
        int yTop = yBin;
        while (!toSearch.empty()) {
          std::pair<int, int> nextPair = toSearch.back();
          toSearch.pop_back();
          if (!explored.insert(nextPair).second)
            continue;
          if (m_totalHoughHist(nextPair.first, nextPair.second).first <
              m_cfg.fracForWidth * m_totalHoughHist(xBin, yBin).first)
            continue;
          if (nextPair.first < xLeft)
            xLeft = nextPair.first;
          if (nextPair.first > xRight)
            xRight = nextPair.first;
          if (nextPair.second < yBottom)
            yBottom = nextPair.second;
          if (nextPair.second > yTop)
            yTop = nextPair.second;
          for (int neighCell = 0; neighCell < 8; ++neighCell) {
            std::pair<int, int> newPair = std::make_pair(
                nextPair.first + stepDirections[neighCell].first,
                nextPair.second + stepDirections[neighCell].second);
            if (newPair.first < 0 || newPair.first >= m_cfg.nBinsX ||
                newPair.second < 0 || newPair.second >= m_cfg.nBinsY)
              continue;
            toSearch.push_back(newPair);
          }
        }
        maximum.wx = std::max(std::abs(m_binCentersX[xLeft] - maximum.x),
                              std::abs(m_binCentersX[xRight] - maximum.x));
        maximum.wy = std::max(std::abs(m_binCentersY[yTop] - maximum.y),
                              std::abs(m_binCentersY[yBottom] - maximum.y));
        maximum.nHits = m_totalHoughHist(xBin, yBin).second.size();
        maxima.push_back(maximum);
      }
    }
  }
  return maxima;
}

void HoughPlane::reset() {
  m_layerHists.clear();
}
bool HoughPlane::passThreshold(HoughHist const& houghHist, unsigned x,
                               unsigned y) const {
  // Pass window threshold
  if (houghHist(x,y).first < m_cfg.threshold)
    return false;

  // Pass local-maximum check, if used
  if (m_cfg.localMaxWindowSize != 0) {
    for (int j = -m_cfg.localMaxWindowSize; j <= m_cfg.localMaxWindowSize;
         j++) {
      for (int i = -m_cfg.localMaxWindowSize; i <= m_cfg.localMaxWindowSize;
           i++) {
        if (i == 0 && j == 0) {
          continue;
        }
        if (x + j < houghHist.size(0) && y + i < houghHist.size(1)) {
          if (houghHist(x + j, y + i).first > houghHist(x,y).first) {
            return false;
          }
          if (houghHist(x + j, y + i).first == houghHist(x,y).first) {
            if (houghHist(x + j, y + i).second.size() >
                houghHist(x,y).second.size()) {
              return false;
            }
            if (houghHist(x + j,  y + i).second.size() ==
                    houghHist(x,y).second.size() &&
                j <= 0 && i <= 0) {
              return false;  // favor bottom-left (low phi, low neg q/pt)
            }
          }
        }
      }
    }
  }

  return true;
}
void HoughPlane::drawHoughHist(HoughHist const&, std::string const&) {}

std::vector<std::vector<int>> HoughPlane::getComboIndices(
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
}  // namespace HoughTransformUtils
}  // namespace Acts
