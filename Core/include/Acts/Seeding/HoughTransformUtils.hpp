// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Grid.hpp"

#include <array>
#include <span>
#include <stdexcept>
#include <unordered_set>

namespace Acts::HoughTransformUtils {

/// this type is responsible for encoding the parameters of our hough space
using CoordType = double;

/// Type alias for hit count/weight values in Hough space
// this type is used to encode hit counts.
// Floating point to allow hit weights to be applied
using YieldType = float;

/// @brief this function represents a mapping of a coordinate point in detector space to a line in
/// hough space. Given the value of the first hough coordinate, it shall return
/// the corresponding second coordinate according to the line parametrisation.
/// Should be implemented by the user.
/// @tparam PointType: The class representing a point in detector space (can differ between implementations)
template <class PointType>
using LineParametrisation =
    std::function<CoordType(CoordType, const PointType&)>;

/// @brief struct to define the ranges of the hough histogram.
/// Used to move between parameter and bin index coordinates.
/// Disconnected from the hough plane binning to be able to reuse
/// a plane with a given binning for several parameter ranges
struct HoughAxisRanges {
  /// Minimum value of the first hough coordinate
  CoordType xMin = 0.0f;  // minimum value of the first coordinate
  /// Maximum value of the first hough coordinate
  CoordType xMax = 0.0f;  // maximum value of the first coordinate
  /// Minimum value of the second hough coordinate
  CoordType yMin = 0.0f;  // minimum value of the second coordinate
  /// Maximum value of the second hough coordinate
  CoordType yMax = 0.0f;  // maximum value of the second coordinate
};

/// convenience functions to link bin indices to axis coordinates

/// @brief find the bin index corresponding to a certain abscissa
/// of the coordinate axis, based on the axis limits and binning.
/// @param min: Start of axis range
/// @param max: End of axis range
/// @param nSteps: Number of bins in axis
/// @param val: value to find the corresponding bin for
/// @return the bin number.
/// No special logic to prevent over-/underflow, checking these is
/// left to the caller
inline int binIndex(double min, double max, unsigned nSteps, double val) {
  return static_cast<int>((val - min) / (max - min) * nSteps);
}
// Returns the lower bound of the bin specified by step
/// @param min: Start of axis range
/// @param max: End of axis range
/// @param nSteps: Number of bins in axis
/// @param binIndex: The index of the bin
/// @return the parameter value at the lower bin edge.
/// No special logic to prevent over-/underflow, checking these is
/// left to the caller
inline double lowerBinEdge(double min, double max, unsigned nSteps,
                           std::size_t binIndex) {
  return min + (max - min) * binIndex / nSteps;
}
// Returns the lower bound of the bin specified by step
/// @param min: Start of axis range
/// @param max: End of axis range
/// @param nSteps: Number of bins in axis
/// @param binIndex: The index of the bin
/// @return the parameter value at the bin center.
/// No special logic to prevent over-/underflow, checking these is
/// left to the caller
inline double binCenter(double min, double max, unsigned nSteps,
                        std::size_t binIndex) {
  return min + (max - min) * 0.5 * (2 * binIndex + 1) / nSteps;
}

/// @brief data class for the information to store for each
/// cell of the hough histogram.
/// @tparam identifier_t: Type of the identifier to associate to the hits
///                       Should be sortable. Used to uniquely identify each
///                       hit and to eventually return the list of hits per cell
template <class identifier_t>
class HoughCell {
 public:
  /// @brief construct the cell as empty
  HoughCell() = default;
  /// @brief add an entry to this cell
  /// @param identifier: Identifier of the hit (used to distinguish hits from another)
  /// @param layer: Layer of the hit (used when counting layers)
  /// @param weight: Optional weight to assign to the hit
  void fill(const identifier_t& identifier, unsigned int layer,
            YieldType weight = 1.);
  /// @brief access the number of layers with hits compatible with this cell
  /// @return The (weighted) number of layers with hits in this cell
  YieldType nLayers() const { return m_nLayers; }
  /// @brief access the number of unique hits compatible with this cell
  /// @return The (weighted) number of unique hits in this cell
  YieldType nHits() const { return m_nHits; }
  /// @brief access the span of layers compatible with this cell
  /// @return Span containing the layer indices with hits in this cell
  std::span<const unsigned, std::dynamic_extent> getLayers() const;
  /// Access the span of unique hits compatible with this cell.
  /// @return Span containing the identifiers of hits in this cell
  std::span<const identifier_t, std::dynamic_extent> getHits() const;

  /// @brief reset this cell, removing any existing content.
  void reset();

 private:
  /// (weighted) number of layers with hits on this cell
  YieldType m_nLayers{0};
  /// (weighted) number of unique hits on this cell
  YieldType m_nHits{0};

  /// index for the hits -- keeps track of vector's size
  std::size_t m_iHit{0};
  /// index for the layers -- keeps track of vector's size
  std::size_t m_iLayer{0};

  /// a batch to resize the vector of the hits or the layers
  std::size_t m_assignBatch{20};

  /// vector of layers with hits on this cell
  std::vector<unsigned> m_layers{std::vector<unsigned>(m_assignBatch)};

  /// vector of hits on this cell
  std::vector<identifier_t> m_hits{std::vector<identifier_t>(m_assignBatch)};
};

/// @brief Configuration - number of bins in each axis.
/// The Hough plane is agnostic of how the bins map to
/// coordinates, allowing to reuse a plane for several
/// (sub) detectors of different dimensions if the bin number
/// remains applicable
struct HoughPlaneConfig {
  /// Number of bins in the first hough coordinate
  std::size_t nBinsX = 0;  // number of bins in the first coordinate
  /// Number of bins in the second hough coordinate
  std::size_t nBinsY = 0;  // number of bins in the second coordinate
};

/// @brief Representation of the hough plane - the histogram used
/// for the hough transform with methods to fill and evaluate
/// the histogram. Templated to a class used as identifier for the hits
template <class identifier_t>
class HoughPlane {
 public:
  /// @brief hough histogram representation as a 2D-indexable vector of hough cells
  using Axis =
      Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Bound>;
  /// Type alias for Hough histogram grid
  using HoughHist = Grid<HoughCell<identifier_t>, Axis, Axis>;
  /// Type alias for histogram index type
  using Index = typename HoughHist::index_t;

  /// @brief instantiate the (empty) hough plane
  /// @param cfg: configuration
  explicit HoughPlane(const HoughPlaneConfig& cfg);

  /// fill and reset methods to modify the grid content

  /// @brief add one measurement to the hough plane
  /// @tparam PointType: Type of the objects to use when adding measurements (e.g. experiment EDM object)
  /// @param measurement: The measurement to add
  /// @param axisRanges: Ranges of the hough axes, used to map the bin numbers to parameter values
  /// @param linePar: The function y(x) parametrising the hough space line for a given measurement
  /// @param widthPar: The function dy(x) parametrising the width of the y(x) curve
  ///                   for a given measurement
  /// @param identifier: The unique identifier for the given hit
  /// @param layer: A layer index for this hit
  /// @param weight: An optional weight to assign to this hit
  template <class PointType>
  void fill(const PointType& measurement, const HoughAxisRanges& axisRanges,
            const LineParametrisation<PointType>& linePar,
            const LineParametrisation<PointType>& widthPar,
            const identifier_t& identifier, unsigned layer = 0,
            YieldType weight = 1.0f);
  /// @brief resets the contents of the grid. Can be used to avoid reallocating the histogram
  /// when switching regions / (sub)detectors
  void reset();

  //// user-facing accessors

  /// @brief get the layers with hits in one cell of the histogram
  /// @param xBin: bin index in the first coordinate
  /// @param yBin: bin index in the second coordinate
  /// @return the layer indices that have hits for this cell
  /// @throws out of range if indices are not within plane limits
  std::span<const unsigned, std::dynamic_extent> layers(
      std::size_t xBin, std::size_t yBin) const {
    checkIndices(xBin, yBin);
    return m_houghHist.atLocalBins({xBin, yBin}).getLayers();
  }

  /// @brief get the (weighted) number of layers  with hits in one cell of the histogram
  /// @param xBin: bin index in the first coordinate
  /// @param yBin: bin index in the second coordinate
  /// @return the (weighed) number of layers that have hits for this cell
  /// @throws out of range if indices are not within plane limits
  YieldType nLayers(std::size_t xBin, std::size_t yBin) const {
    checkIndices(xBin, yBin);
    return m_houghHist.atLocalBins({xBin, yBin}).nLayers();
  }

  /// @brief get the identifiers of all hits in one cell of the histogram
  /// @param xBin: bin index in the first coordinate
  /// @param yBin: bin index in the second coordinate
  /// @return the list of identifiers of the hits for this cell
  /// Can include duplicates if a hit was filled more than once
  /// @throws out of range if indices are not within plane limits
  std::span<const identifier_t, std::dynamic_extent> hitIds(
      std::size_t xBin, std::size_t yBin) const {
    checkIndices(xBin, yBin);
    return m_houghHist.atLocalBins({xBin, yBin}).getHits();
  }
  /// @brief get the identifiers of all hits in one cell of the histogram
  /// @param xBin: bin index in the first coordinate
  /// @param yBin: bin index in the second coordinate
  /// @return the list of identifiers of the hits for this cell
  /// Guaranteed to not duplicate identifiers
  /// @throws out of range if indices are not within plane limits
  std::unordered_set<const identifier_t> uniqueHitIds(std::size_t xBin,
                                                      std::size_t yBin) const {
    checkIndices(xBin, yBin);
    const auto hits_span = m_houghHist.atLocalBins({xBin, yBin}).getHits();
    return std::unordered_set<identifier_t>(hits_span.begin(), hits_span.end());
  }
  /// @brief access the (weighted) number of hits in one cell of the histogram from bin's coordinates
  /// @param xBin: bin index in the first coordinate
  /// @param yBin: bin index in the second coordinate
  /// @return the (weighted) number of hits for this cell
  /// @throws out of range if indices are not within plane limits
  YieldType nHits(std::size_t xBin, std::size_t yBin) const {
    checkIndices(xBin, yBin);
    return m_houghHist.atLocalBins({xBin, yBin}).nHits();
  }

  /// @brief access the (weighted) number of hits in one cell of the histogram from globalBin index
  /// @param globalBin: global bin index
  /// @return the (weighted) number of hits for this cell
  YieldType nHits(std::size_t globalBin) const {
    return m_houghHist.at(globalBin).nHits();
  }

  /// @brief get the number of bins on the first coordinate
  /// @return Number of bins in the X direction
  std::size_t nBinsX() const { return m_cfg.nBinsX; }
  /// @brief get the number of bins on the second coordinate
  /// @return Number of bins in the Y direction
  std::size_t nBinsY() const { return m_cfg.nBinsY; }

  /// @brief get the maximum number of (weighted) hits seen in a single
  /// cell across the entire histrogram.
  /// @return Maximum number of hits found in any single cell
  YieldType maxHits() const { return m_maxHits; }

  /// @brief get the list of cells with non-zero content.
  /// Useful for peak-finders in sparse data
  /// to avoid looping over all cells
  /// @return Reference to set of global bin indices with non-zero content
  const std::unordered_set<std::size_t>& getNonEmptyBins() const {
    return m_touchedBins;
  }

  /// @brief get the coordinates of the bin given the global bin index
  /// @param globalBin Global bin index to convert to coordinates
  /// @return Local bin coordinates (x,y) corresponding to global bin index
  Index axisBins(std::size_t globalBin) const {
    return m_houghHist.localBinsFromGlobalBin(globalBin);
  }

  /// @brief get the globalBin index given the coordinates of the bin
  /// @param indexBin Bin coordinates to convert to global index
  /// @return Global bin index corresponding to local bin coordinates
  std::size_t globalBin(Index indexBin) const {
    return m_houghHist.globalBinFromLocalBins(indexBin);
  }

  /// @brief get the bin indices of the cell containing the largest number
  /// of (weighted) hits across the entire histogram
  /// @return Pair of (x,y) bin indices where maximum hits are found
  std::pair<std::size_t, std::size_t> locMaxHits() const {
    return m_maxLocHits;
  }

  /// @brief get the maximum number of (weighted) layers with hits  seen
  /// in a single cell across the entire histrogram.
  /// @return Maximum number of layers found in any single cell
  YieldType maxLayers() const { return m_maxLayers; }

  /// @brief get the bin indices of the cell containing the largest number
  /// of (weighted) layers with hits across the entire histogram
  /// @return Pair of (x,y) bin indices where maximum layers are found
  std::pair<std::size_t, std::size_t> locMaxLayers() const {
    return m_maxLocLayers;
  }

  /// @brief Helper method to fill a bin of the hough histogram.
  /// Updates the internal helper data structures (maximum tracker etc).
  /// @param binX: bin number along x
  /// @param binY: bin number along y
  /// @param identifier: hit identifier
  /// @param layer: layer index
  /// @param w: optional hit weight
  void fillBin(std::size_t binX, std::size_t binY,
               const identifier_t& identifier, unsigned layer, double w = 1.0f);

 private:
  YieldType m_maxHits = 0.0f;    // track the maximum number of hits seen
  YieldType m_maxLayers = 0.0f;  // track the maximum number of layers seen

  /// track the location of the maximum in hits
  std::pair<std::size_t, std::size_t> m_maxLocHits = {0, 0};
  /// track the location of the maximum in layers
  std::pair<std::size_t, std::size_t> m_maxLocLayers = {0, 0};

  std::size_t m_assignBatch{20};

  /// track the bins with non-trivial content
  std::unordered_set<std::size_t> m_touchedBins{};

  std::size_t m_iBin = 0;

  HoughPlaneConfig m_cfg;  // the configuration object
  HoughHist m_houghHist;   // the histogram data object

  /// @brief check if indices are are valid
  void checkIndices(std::size_t x, std::size_t y) const;
};

/// example peak finders.
namespace PeakFinders {
/// configuration for the LayerGuidedCombinatoric peak finder
struct LayerGuidedCombinatoricConfig {
  /// Minimum number of layers required to form a peak
  YieldType threshold = 3.0f;  // min number of layers to obtain a maximum
  /// Size of window for local maximum detection
  int localMaxWindowSize = 0;  // Only create candidates from a local maximum
};

/// @brief Peak finder inspired by ATLAS ITk event filter developments.
/// Builds peaks based on layer counts and allows for subsequent resolution
/// of the combinatorics by building multiple candidates per peak if needed.
/// In flat regions, peak positions are moved towards smaller values of the
/// second and first coordinate.
/// @tparam identifier_t: The identifier type to use. Should match the one used in the hough plane.
template <class identifier_t>
class LayerGuidedCombinatoric {
 public:
  /// @brief data class representing the found maxima.
  /// Here, we just have a list of cluster identifiers
  struct Maximum {
    /// Set of hit identifiers contributing to this peak
    std::unordered_set<identifier_t> hitIdentifiers =
        {};  // identifiers of contributing hits
  };
  /// @brief constructor
  /// @param cfg: Configuration object
  explicit LayerGuidedCombinatoric(const LayerGuidedCombinatoricConfig& cfg);

  /// @brief main peak finder method.
  /// @param plane: Filled hough plane to search
  /// @return vector of found maxima
  std::vector<Maximum> findPeaks(const HoughPlane<identifier_t>& plane) const;

 private:
  /// @brief check if a given bin is a local maximum.
  /// @param plane: The filled hough plane
  /// @param xBin: x bin index
  /// @param yBin: y bin index
  /// @return true if a maximum, false otherwise
  bool passThreshold(const HoughPlane<identifier_t>& plane, std::size_t xBin,
                     std::size_t yBin) const;  // did we pass extensions?

  LayerGuidedCombinatoricConfig m_cfg;  // configuration data object
};
/// @brief Configuration for the IslandsAroundMax peak finder
struct IslandsAroundMaxConfig {
  /// Minimum number of weighted hits required for a peak
  YieldType threshold =
      3.0f;  // min number of weigted hits required in a maximum
  /// Fraction of global maximum below which peaks are ignored
  YieldType fractionCutoff =
      0;  // Fraction of the global maximum at which to cut off maxima
  /// Minimum spacing between peaks in parameter space
  std::pair<CoordType, CoordType> minSpacingBetweenPeaks = {
      0.0f, 0.0f};  // minimum distance of a new peak from existing peaks in
                    // parameter space
};
/// @brief Peak finder inspired by ATLAS muon reconstruction.
/// Looks for regions above a given fraction of the global maximum
/// hit count and connects them into islands comprising adjacent bins
/// above the threshold. Peak positions are averaged across cells in the island,
/// weighted by hit counts
/// @tparam identifier_t: The identifier type to use. Should match the one used in the hough plane.
template <class identifier_t>
class IslandsAroundMax {
 public:
  /// @brief data struct representing a local maximum.
  /// Comes with a position estimate and a list of hits within the island
  struct Maximum {
    /// X coordinate of the peak maximum
    CoordType x = 0;
    /// Y coordinate of the peak maximum
    CoordType y = 0;
    /// Width of the peak in X direction
    CoordType wx = 0;
    /// Width of the peak in Y direction
    CoordType wy = 0;
    /// Set of hit identifiers contributing to this peak
    std::unordered_set<identifier_t> hitIdentifiers =
        {};  // identifiers of contributing hits
  };
  /// @brief constructor.
  /// @param cfg: configuration object
  explicit IslandsAroundMax(const IslandsAroundMaxConfig& cfg);

  /// @brief main peak finder method.
  /// @param plane: The filled hough plane to search
  /// @param ranges: The axis ranges used for mapping between parameter space and bins.
  /// @return List of the found maxima
  std::vector<Maximum> findPeaks(const HoughPlane<identifier_t>& plane,
                                 const HoughAxisRanges& ranges);

 private:
  /// @brief method to incrementally grow an island by adding adjacent cells
  /// Performs a breadth-first search for neighbours above threshold and adds
  /// them to candidate. Stops when no suitable neighbours are left.
  /// @param houghPlane: The current hough Plane we are looking for maxima
  /// @param inMaximum: List of cells found in the island. Incrementally populated by calls to the method
  /// @param toExplore: List of the global Bin indices of neighbour cell candidates left to explore. Method will not do anything once this is empty
  /// @param threshold: the threshold to apply to check if a cell should be added to an island
  /// @param yieldMap: A map of the hit content of above-threshold cells. Used cells will be set to empty content to avoid reuse by subsequent calls
  void extendMaximum(const HoughPlane<identifier_t>& houghPlane,
                     std::vector<std::array<std::size_t, 2>>& inMaximum,
                     std::vector<std::size_t>& toExplore, YieldType threshold,
                     std::unordered_map<std::size_t, YieldType>& yieldMap);

  IslandsAroundMaxConfig m_cfg;  // configuration data object

  /// @brief array of steps to consider when exploring neighbouring cells.
  const std::array<std::pair<int, int>, 8> m_stepDirections{
      std::make_pair(-1, -1), std::make_pair(0, -1), std::make_pair(1, -1),
      std::make_pair(-1, 0),  std::make_pair(1, 0),  std::make_pair(-1, 1),
      std::make_pair(0, 1),   std::make_pair(1, 1)};
};

/// @brief Peak finder using sliding window algorithm.
/// First it finds peaks by scanning all space for cells with number of hits
/// above threshold. Then applies sliding window (SW) logic to eliminate peaks
/// when maxima are adjacent leaving only one of them in a window. This SW
/// implementation requires that none on the upper right corner are above peak
/// and none in bottom left corner is below or equal to the peak. It can be
/// illustrated as follows for window size of 1:
///
///
///  <= <= <=
///   <  O <=
///   <  <  <
///
/// Then the algorithm collects maxima in a window (possibly of different size)
/// and calculates peak position using weighted average.

struct SlidingWindowConfig {
  /// peak threshold, cell content is compared with it using >= operator
  std::size_t threshold = 3;
  /// size of the window in x direction for sliding window
  std::size_t xWindowSize = 2;
  /// size of the window in y direction for sliding window
  std::size_t yWindowSize = 2;
  /// perform re-centering
  bool recenter = true;
  /// size of the window in x direction for recentering, this should be
  /// typically <= window size
  std::size_t xRecenterSize = 3;
  /// size of the window in y direction for recentering, this should be
  /// typically <= window size
  std::size_t yRecenterSize = 3;
};

/// @brief Obtain peaks list in Hough space using Sliding Window algorithm
/// @tparam identifier_t Hough plane content
/// @param plane Hough plane to work on
/// @param config algorithm configuration
/// @return list of indices (pairs of numbers)
template <typename identifier_t>
std::vector<typename HoughPlane<identifier_t>::Index> slidingWindowPeaks(
    const HoughPlane<identifier_t>& plane, const SlidingWindowConfig& config);

/// @brief Obtain an image around the peak
/// @tparam identifier_t Hough plane content
/// @param plane Hough plane to work on
/// @param index peak center
/// @param xSize number of cells around the peak in x direction
/// @param ySize number of cells around the peak in y direction
/// @param summaryFunction function constructing pixel content, default is just number of hits in cell
///        Other implementations may take into account layers
/// @return the vector with count of hits starting from lower left to upper right corner of rectangular window
/// The "image" is always of the same size, if it would happen to be outside of
/// Hough plane the content is padded with zeros
template <typename identifier_t, typename pixel_value_t = unsigned char>
std::vector<pixel_value_t> hitsCountImage(
    const HoughPlane<identifier_t>& plane,
    typename HoughPlane<identifier_t>::Index index, std::size_t xSize,
    std::size_t ySize,
    const std::function<pixel_value_t(const HoughPlane<identifier_t>&, int,
                                      int)>& summaryFunction =
        [](const HoughPlane<identifier_t>& plane, int x, int y) {
          return static_cast<pixel_value_t>(plane.nHits(x, y));
        });

}  // namespace PeakFinders
}  // namespace Acts::HoughTransformUtils

#include "HoughTransformUtils.ipp"
