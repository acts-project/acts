// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// @file HoughTransformSeeder.hpp
// @author Max Goblirsch, based on work by Riley Xu and Jahred Adelman
// @brief Implements a set of tools to
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

/// this type is responsible for encoding the parameters of our hough space
using coordType = double;

// this type is used to encode hit counts. 
// Floating point to allow hit weights to be applied
using yieldType = float;

/// @brief this function represents a mapping of a coordinate point in detector space to a line in
/// hough space. Given the value of the first hough coordinate, it should return 
/// the corresponding second coordinate according to the line parametrisation.
/// Should be implemented by the user. 
/// @tparam PointType: The class representing a point in detector space (can differ between implementations)
template <class PointType>
using lineParametrisation =
    std::function<coordType(coordType, const PointType&)>;

/// @brief struct to define the ranges of the hough histogram. 
/// Used to move between parameter and bin index coordinates. 
/// Disconnected from the hough plane binning to be able to re-use
/// a plane with a given binning for several parameter ranges
struct houghAxisRanges {
  coordType xMin = 0.0f;  // minimum value of the first coordinate
  coordType xMax = 0.0f;  // maximim value of the first coordinate 
  coordType yMin = 0.0f;  // minimum value of the second coordinate
  coordType yMax = 0.0f;  // maximum value of the second coordinate
};

/// convenience functions to link bin indices to axis coordinates

/// @brief find the bin index corresponding to a certain abscissa 
/// of the coordinate axis, based on the axis limits and binning. 
/// @param min: Start of axis range 
/// @param max: End of axis range 
/// @param nSteps: Number of bins in axis
/// @param val: value to find the corresponding bin for 
/// @return the bin number. 
/// @note: No special logic to prevent over-/underflow, checking these is 
/// left to the caller 
int binIndex(double min, double max, unsigned nSteps,
                            double val) {
  return static_cast<int>((val - min) / (max - min) * nSteps);
}
// Returns the lower bound of the bin specified by step
double lowerBinEdge(double min, double max, unsigned nSteps,
                                  int step) {
  return min + (max - min) * step / nSteps;
}
// Returns the lower bound of the bin specified by step
double binCenter(double min, double max, unsigned nSteps,
                                  int step) {
  return min + (max - min) * 0.5 * (2 * step + 1) / nSteps;
}

/// @brief data class for the information to store for each 
/// cell of the hough histogram. 
/// @tparam identifier_t: Type of the identifier to associate to the hits 
///                       Should be sortable. Used to uniquely identify each 
///                       hit and to eventually return the list of hits per cell
template <class identifier_t> class HoughCell {
 public:
  /// @brief construct the cell as empty 
  HoughCell() = default;
  /// @brief add an entry to this cell
  /// @param identifier: Identifier of the hit (used to distinguish hits from another) 
  /// @param layer: Layer of the hit (used when counting layers)
  /// @param weight: Optional weight to assign to the hit 
  void fill(const identifier_t & identifier, unsigned layer, yieldType weight = 1.);
  /// @brief access the number of layers with hits compatible with this cell 
  yieldType nLayers() const { return m_nLayers; }
  /// @brief access the number of unique hits compatible with this cell 
  yieldType nHits() const { return m_nHits; }
  /// @brief access the set of layers compatible with this cell
  const std::unordered_set<unsigned>& layers() const { return m_layers; }
  /// @brief access the set of unique hits compatible with this cell
  const std::unordered_set<identifier_t>& hits() const { return m_hits; }
  /// @brief reset this cell, removing any existing content. 
  void reset();

 private:
  /// data members 

  yieldType m_nLayers = 0; // (weighted) number of layers with hits on this cell
  yieldType m_nHits = 0;   // (weighted) number of unique hits on this cell
  std::unordered_set<unsigned> m_layers = {};  // set of layers with hits on this cell
  std::unordered_set<identifier_t> m_hits = {};      // set of unique hits on this cell
};

/// @brief Configuration - number of bins in each axis. 
/// The Hough plane is agnostic of how the bins map to 
/// coordinates, allowing to re-use a plane for several 
/// (sub) detectors of different dimensions if the bin number 
/// remains applicable
struct HoughPlaneConfig {
  int nBinsX = 0;  // number of bins in the first coordinate
  int nBinsY = 0;  // number of bins in the second coordinate
};

/// @brief Representation of the hough plane - the histogram used
/// for the hough transform with methods to fill and evaluate 
/// the histogram. Templated to a class used as identifier for the hits 
template <class identifier_t> class HoughPlane {
 public:
  /// @brief hough histogram representation as a 2D-indexable vector of hough cells 
  using HoughHist = vector2D<HoughCell<identifier_t>>;

  /// @brief instantiate the (empty) hough plane
  /// @param cfg: configuration
  HoughPlane(const HoughPlaneConfig& cfg);

  /// fill and reset methods to modify the grid content 

  /// @brief add one measurement to the hough plane
  /// @tparam PointType: Type of the objects to use when adding measurements (e.g. experiment EDM object)
  /// @param measurement: The measurement to add
  /// @param axisRanges: Ranges of the hough axes, used to map the experimental EDM coordinates 
  ///                    to bin numbers
  /// @param linePar: The function y(x) parametrising the hough space line for a given measurement
  /// @param widthPar: The function dy(x) parametrising the width of the y(x) curve
  ///                   for a given measurement
  /// @param idGetter: A function that maps our measurement to a (unique) identifier.
  /// @param layerGetter: A function that maps our measurement to a layer index
  template <class PointType>
  void fill(const PointType& measurement, const houghAxisRanges& axisRanges,
            lineParametrisation<PointType> linePar,
            lineParametrisation<PointType> widthPar, const identifier_t & identifier,
            unsigned layer = 0, yieldType weight = 1.0f);
  /// @brief resets the contents of the grid. Can be used to avoid reallocating the histogram
  /// when switching regions / (sub)detectors
  void reset();

  //// user-facing accessors 

  /// @brief get the layers with hits in one cell of the histogram 
  /// @param x: bin index in the first coordinate
  /// @param y: bin index in the second coordinate
  /// @return the set of layer indices that have hits for this cell
  const std::unordered_set<unsigned>& layers(int x, int y) const {
    return m_houghHist(x, y).layers();
  }

  /// @brief get the (weighted) number of layers  with hits in one cell of the histogram 
  /// @param x: bin index in the first coordinate
  /// @param y: bin index in the second coordinate
  /// @return the (weighed) number of layers that have hits for this cell
  yieldType nLayers(int x, int y) const { return m_houghHist(x, y).nLayers(); }

  /// @brief get the identifiers of all hits in one cell of the histogram 
  /// @param x: bin index in the first coordinate
  /// @param y: bin index in the second coordinate
  /// @return the set of identifiers of the hits for this cell
  const std::unordered_set<identifier_t>& hitIds (
      int x, int y) const {
    return m_houghHist(x, y).hits();
  }
  /// @brief get the (weighted) number of hits in one cell of the histogram 
  /// @param x: bin index in the first coordinate
  /// @param y: bin index in the second coordinate
  /// @return the (weighted) number of hits for this cell
  yieldType nHits(int x, int y) const { return m_houghHist(x, y).nHits(); }

  /// @brief get the number of bins on the first coordinate 
  int nBinsX() const { return m_cfg.nBinsX; }
  /// @brief get the number of bins on the second coordinate 
  int nBinsY() const { return m_cfg.nBinsY; }

  /// @brief get the maximum number of (weighted) hits seen in a single 
  /// cell across the entire histrogram. 
  yieldType maxHits() const { return m_maxHits; }

  /// @brief get the bin indices of the cell containing the largest number
  /// of (weighted) hits across the entire histogram 
  std::pair<int, int> locMaxHits() const { return m_maxLocHits; }

  /// @brief get the maximum number of (weighted) layers with hits  seen 
  /// in a single cell across the entire histrogram. 
  yieldType maxLayers() const { return m_maxLayers; }

  /// @brief get the bin indices of the cell containing the largest number
  /// of (weighted) layers with hits across the entire histogram 
  std::pair<int, int> locMaxLayers() const { return m_maxLocLayers; }


 private:
  void fillBin(int binX, int binY, const identifier_t & identifier, unsigned layer,
               double w = 1.0f);
  void drawHoughHist(std::string const& name);  // for making pretty plots

  yieldType m_maxHits = 0.0f;
  yieldType m_maxLayers = 0.0f;
  std::pair<int, int> m_maxLocHits = {0, 0};
  std::pair<int, int> m_maxLocLayers = {0, 0};
  HoughPlaneConfig m_cfg;
  HoughHist m_houghHist;
};

/// example peak finders

struct HoughPeakFinder_LayerGuidedCombinatoricConfig {
  yieldType threshold = 3.0f;  // min number of layers to obtain a maximum
  int localMaxWindowSize = 0;  // Only create candidates from a local maximum
};
template <class identifier_t> class HoughPeakFinder_LayerGuidedCombinatoric {
 public:
  struct Maximum {
    std::unordered_set<identifier_t> hitIdentifiers =
        {};  // identifiers of contributing hits
  };
  HoughPeakFinder_LayerGuidedCombinatoric(
      const HoughPeakFinder_LayerGuidedCombinatoricConfig& cfg);
  std::vector<Maximum> findPeaks(const HoughPlane<identifier_t>& plane) const;

 private:
  bool passThreshold(const HoughPlane<identifier_t>& plane, int xBin,
                     int yBin) const;  // did we pass extensions?
  std::vector<std::vector<int>> getComboIndices(
      std::vector<size_t>& sizes) const;

  HoughPeakFinder_LayerGuidedCombinatoricConfig m_cfg;
};
struct HoughPeakFinder_IslandsAroundMaxConfig {
  yieldType threshold = 3.0f;    // min number of layers to obtain a maximum
  yieldType fractionCutoff = 0;  // Fraction of the global max to cut at
};
template <class identifier_t> class HoughPeakFinder_IslandsAroundMax {
 public:
  struct Maximum {
    coordType x = 0;   // x value of the maximum
    coordType y = 0;   // y value of the maximum
    coordType wx = 0;  // x width of the maximum
    coordType wy = 0;  // y width of the maximum
    std::unordered_set<identifier_t> hitIdentifiers =
        {};  // identifiers of contributing hits
  };
  HoughPeakFinder_IslandsAroundMax(
      const HoughPeakFinder_IslandsAroundMaxConfig& cfg);
  std::vector<Maximum> findPeaks(const HoughPlane<identifier_t>& plane, const houghAxisRanges & ranges);

 private:
  void extendMaximum(
    std::vector<std::pair<int, int>>& inMaximum, std::vector<std::pair<int, int>>& toExplore, yieldType threshold,
      vector2D<yieldType>& yieldMap);
  HoughPeakFinder_IslandsAroundMaxConfig m_cfg;
  const std::array<std::pair<int, int>, 8> m_stepDirections {
      std::make_pair(-1, -1), std::make_pair(0, -1), std::make_pair(1, -1),
      std::make_pair(-1, 0),  std::make_pair(1, 0),  std::make_pair(-1, 1),
      std::make_pair(0, 1),   std::make_pair(1, 1)
    };
};
}  // namespace HoughTransformUtils
};  // namespace Acts

#include "HoughTransformUtils.ipp"
