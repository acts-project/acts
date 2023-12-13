// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// @file HoughTransformSeeder.hpp
// @author Max Goblirsch, based on work by Riley Xu and Jahred Adelman
// @brief Implements a (hopefully) generic set of tools to
/// implement a hough transform.

#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <array>
#include <optional>
#include <set>
#include <unordered_set>

#include "HoughVectors.hpp"

#pragma once

namespace Acts {
namespace HoughTransformUtils {

/// this type is responsible for encoding the coordinates of our hough hist
using coordType = double;
// this type is used as an identifier proxy
using idType = long long int;
using yieldType = float;

/// @brief this function represents a mapping of a coordinate point in detector space to a line in
/// hough space. Given the value of the first hough coordinate, it will return
/// the corresponding second coordinate according to the line parametrisation.
/// @tparam PointType: The class representing a point in detector space (can differ between implementations)
template <class PointType>
using lineParametrisation =
    std::function<coordType(coordType, const PointType&)>;
struct peakFinderConfig {
  yieldType threshold = 0;
  double fracForWidth = 0.5;
  yieldType deltaToGlobal =
      10000;  // max yield diff to global max for secondary maxima
  int localMaxWindowSize = 2;  // Only create candidates from a local maximum
};
struct houghAxisRanges {
  coordType xMin = 0.0f;
  coordType xMax = 0.0f;
  coordType yMin = 0.0f;
  coordType yMax = 0.0f;
};

class HoughCell {
 public:
  HoughCell() = default;
  yieldType nLayers() const { return m_nLayers; }
  yieldType nHits() const { return m_nHits; }
  const std::unordered_set<unsigned>& layers() const { return m_layers; }
  const std::unordered_set<idType>& hits() const { return m_hits; }
  void fill(idType identifier, unsigned layer, yieldType weight = 1.);
  void reset();

 private:
  yieldType m_nLayers = 0;
  yieldType m_nHits = 0;
  std::unordered_set<unsigned> m_layers = {};
  std::unordered_set<idType> m_hits = {};
};

using HoughHist = vector2D<HoughCell>;
class HoughPlane {
 public:
  struct Config {
    int nBinsX = 0;  // number of bins in the first coordinate
    int nBinsY = 0;  // number of bins in the second coordinate
  };
  /// @brief instantiate the (empty) hough plane
  /// @param cfg: configuration
  HoughPlane(const HoughPlane::Config& cfg);
  /// @brief add one measurement to the hough plane
  /// @tparam PointType: Type of the measurement
  /// @param measurement: The measurement to add
  /// @param linePar: The function y(x) parametrising the hough space line for a given measurement
  /// @param widthPar: The function dy(x) parametrising the width of the y(x) curve
  ///                   for a given measurement
  /// @param idGetter: A function that maps our measurement to a (unique) identifier.
  /// @param layerGetter: A function that maps our measurement to a layer index
  template <class PointType>
  void fill(const PointType& measurement, const houghAxisRanges& axisRanges,
            lineParametrisation<PointType> linePar,
            lineParametrisation<PointType> widthPar, idType identifier,
            unsigned layer = 0, yieldType weight = 1.0f);
  /// resets the contents of the grid, starting from scratch
  void reset();
  const std::unordered_set<unsigned>& layers(int x, int y) const {
    return m_houghHist(x, y).layers();
  }
  yieldType nLayers(int x, int y) const { return m_houghHist(x, y).nLayers(); }
  const std::unordered_set<Acts::HoughTransformUtils::idType>& hitIds(
      int x, int y) const {
    return m_houghHist(x, y).hits();
  }
  yieldType nHits(int x, int y) const { return m_houghHist(x, y).nHits(); }
  void fillBin(int binX, int binY, idType identifier, unsigned layer,
               double w = 1.0f);

  static inline int binIndex(double min, double max, unsigned nSteps,
                             double val) {
    return static_cast<int>((val - min) / (max - min) * nSteps);
  }
  // Returns the lower bound of the bin specified by step
  static inline double lowerBinEdge(double min, double max, unsigned nSteps,
                                    int step) {
    return min + (max - min) * step / nSteps;
  }
  // Returns the lower bound of the bin specified by step
  static inline double binCenter(double min, double max, unsigned nSteps,
                                    int step) {
    return min + (max - min) * 0.5 * (2 * step + 1) / nSteps;
  }
  int nBinsX() const { return m_cfg.nBinsX; }
  int nBinsY() const { return m_cfg.nBinsY; }

  yieldType maxHits() const { return m_maxHits; }
  std::pair<int, int> locMaxHits() const { return m_maxLocHits; }
  yieldType maxLayers() const { return m_maxLayers; }
  std::pair<int, int> locMaxLayers() const { return m_maxLocLayers; }

 private:
  void drawHoughHist(std::string const& name);  // for making pretty plots

  yieldType m_maxHits = 0.0f;
  yieldType m_maxLayers = 0.0f;
  std::pair<int, int> m_maxLocHits = {0, 0};
  std::pair<int, int> m_maxLocLayers = {0, 0};
  Config m_cfg;
  HoughHist m_houghHist;
};

/// example peak finders

class HoughPeakFinder_LayerGuidedCombinatoric {
 public:
  struct Config {
    yieldType threshold = 3.0f;  // min number of layers to obtain a maximum
    int localMaxWindowSize = 0;  // Only create candidates from a local maximum
  };
  struct Maximum {
    std::unordered_set<idType> hitIdentifiers =
        {};  // identifiers of contributing hits
  };
  HoughPeakFinder_LayerGuidedCombinatoric(
      const HoughPeakFinder_LayerGuidedCombinatoric::Config& cfg);
  std::vector<Maximum> findPeaks(const HoughPlane& plane) const;

 private:
  bool passThreshold(const HoughPlane& plane, int xBin,
                     int yBin) const;  // did we pass extensions?
  std::vector<std::vector<int>> getComboIndices(
      std::vector<size_t>& sizes) const;

  HoughPeakFinder_LayerGuidedCombinatoric::Config m_cfg;
};

class HoughPeakFinder_IslandsAroundMax {
 public:
  struct Config {
    yieldType threshold = 3.0f;    // min number of layers to obtain a maximum
    yieldType fractionCutoff = 0;  // Fraction of the global max to cut at
  };
  struct Maximum {
    coordType x = 0;   // x value of the maximum
    coordType y = 0;   // y value of the maximum
    coordType wx = 0;  // x width of the maximum
    coordType wy = 0;  // y width of the maximum
    std::unordered_set<idType> hitIdentifiers =
        {};  // identifiers of contributing hits
  };
  HoughPeakFinder_IslandsAroundMax(
      const HoughPeakFinder_IslandsAroundMax::Config& cfg);
  std::vector<Maximum> findPeaks(const HoughPlane& plane, const houghAxisRanges & ranges);

 private:
  void extendMaximum(
    std::vector<std::pair<int, int>>& inMaximum, std::vector<std::pair<int, int>>& toExplore, yieldType threshold,
      vector2D<yieldType>& yieldMap);
  HoughPeakFinder_IslandsAroundMax::Config m_cfg;
  const std::array<std::pair<int, int>, 8> m_stepDirections {
      std::make_pair(-1, -1), std::make_pair(0, -1), std::make_pair(1, -1),
      std::make_pair(-1, 0),  std::make_pair(1, 0),  std::make_pair(-1, 1),
      std::make_pair(0, 1),   std::make_pair(1, 1)
    };
};
}  // namespace HoughTransformUtils
};  // namespace Acts

#include "HoughTransformUtils.ipp"
