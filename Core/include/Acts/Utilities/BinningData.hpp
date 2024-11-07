// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace Acts {

/// @class BinningData
///
///   This class holds all the data necessary for the bin calculation
///
///   phi has a very particular behaviour:
///   - there's the change around +/- PI
///
///   - it can be multiplicative or additive
///   multiplicative : each major bin has the same sub structure
///                    i.e. first binnning
///
/// structure is Equidistant
///   additive : sub structure replaces one bin (and one bin only)
///
///
class BinningData {
 public:
  AxisType axisType{};  ///< binning type: Equidistant, Variable
  AxisBoundaryType
      axisBoundaryType{};  ///< axis boundary type: (Open,) Bound, Closed
  AxisDirection
      axisDirection{};  ///< binning value: AxisX, AxisY, AxisZ, AxisR ...
  float min{};          ///< minimum value
  float max{};          ///< maximum value
  float step{};         ///< binning step
  bool zdim{};          ///< zero dimensional binning : direct access

  /// sub structure: describe some sub binning
  std::unique_ptr<const BinningData> subBinningData;
  /// sub structure: additive or multipicative
  bool subBinningAdditive{};

  /// Constructor for 0D binning
  ///
  /// @param aDir is the binning value: AxisX, AxisY, etc.
  /// @param bMin is the minimum value
  /// @param bMax is the maximum value
  BinningData(AxisDirection aDir, float bMin, float bMax)
      : axisType(AxisType::Equidistant),
        axisBoundaryType(AxisBoundaryType::Bound),
        axisDirection(aDir),
        min(bMin),
        max(bMax),
        step((bMax - bMin)),
        zdim(true),
        subBinningData(nullptr),
        m_bins(1),
        m_boundaries({{min, max}}),
        m_totalBins(1),
        m_totalBoundaries(std::vector<float>()),
        m_functionPtr(&searchEquidistantWithBoundary) {}

  /// Constructor for Equidistant binning
  /// and axisBoundaryTypeal sub structure can be
  /// multiplicative or additive
  ///
  /// @param abType is the binning axisBoundaryType : open, AxisBoundaryType::Closed
  /// @param aDir is the binning value: AxisX, AxisY, etc.
  /// @param bBins is number of equidistant bins
  /// @param bMin is the minimum value
  /// @param bMax is the maximum value
  /// @param sBinData is (axisBoundaryTypeal) sub structure
  /// @param sBinAdditive is the prescription for the sub structure
  BinningData(AxisBoundaryType abType, AxisDirection aDir, std::size_t bBins,
              float bMin, float bMax,
              std::unique_ptr<const BinningData> sBinData = nullptr,
              bool sBinAdditive = false)
      : axisType(AxisType::Equidistant),
        axisBoundaryType(abType),
        axisDirection(aDir),
        min(bMin),
        max(bMax),
        step((bMax - bMin) / bBins),
        zdim(bBins == 1 ? true : false),
        subBinningData(std::move(sBinData)),
        subBinningAdditive(sBinAdditive),
        m_bins(bBins),
        m_boundaries(std::vector<float>()),
        m_totalBins(bBins),
        m_totalBoundaries(std::vector<float>()) {
    // set to Equidistant search
    m_functionPtr = &searchEquidistantWithBoundary;
    // fill the boundary vector for fast access to center & boundaries
    m_boundaries.reserve(m_bins + 1);
    for (std::size_t ib = 0; ib < m_bins + 1; ++ib) {
      m_boundaries.push_back(min + ib * step);
    }
    // the binning data has sub structure - multiplicative or additive
    checkSubStructure();
  }

  /// Constructor for non-Equidistant binning
  ///
  /// @param abType is the binning axisBoundaryType : (Open), Bound, Closed
  /// @param aDir is the binning value : AxisX, AxisY, etc.
  /// @param bBoundaries are the bin boundaries
  /// @param sBinData is (axisBoundaryTypeal) sub structure
  BinningData(AxisBoundaryType abType, AxisDirection aDir,
              const std::vector<float>& bBoundaries,
              std::unique_ptr<const BinningData> sBinData = nullptr)
      : axisType(AxisType::Variable),
        axisBoundaryType(abType),
        axisDirection(aDir),
        zdim(bBoundaries.size() == 2 ? true : false),
        subBinningData(std::move(sBinData)),
        subBinningAdditive(true),
        m_bins(bBoundaries.size() - 1),
        m_boundaries(bBoundaries),
        m_totalBins(bBoundaries.size() - 1),
        m_totalBoundaries(bBoundaries) {
    // assert a no-size case
    throw_assert(m_boundaries.size() > 1, "Must have more than one boundary");
    min = m_boundaries[0];
    max = m_boundaries[m_boundaries.size() - 1];
    // set to Equidistant search
    m_functionPtr = &searchInVectorWithBoundary;
    // the binning data has sub structure - multiplicative
    checkSubStructure();
  }

  /// Copy constructor
  ///
  /// @param bdata is the source object
  BinningData(const BinningData& bdata)
      : axisType(bdata.axisType),
        axisBoundaryType(bdata.axisBoundaryType),
        axisDirection(bdata.axisDirection),
        min(bdata.min),
        max(bdata.max),
        step(bdata.step),
        zdim(bdata.zdim),
        subBinningData(nullptr),
        subBinningAdditive(bdata.subBinningAdditive),
        m_bins(bdata.m_bins),
        m_boundaries(bdata.m_boundaries),
        m_totalBins(bdata.m_totalBins),
        m_totalBoundaries(bdata.m_totalBoundaries) {
    // get the binning data
    subBinningData =
        bdata.subBinningData
            ? std::make_unique<const BinningData>(*bdata.subBinningData)
            : nullptr;
    // set the correct function pointer
    if (axisType == AxisType::Equidistant) {
      m_functionPtr = &searchEquidistantWithBoundary;
    } else {
      m_functionPtr = &searchInVectorWithBoundary;
    }
  }

  /// Assignment operator
  ///
  /// @param bdata is the source object
  BinningData& operator=(const BinningData& bdata) {
    if (this != &bdata) {
      axisType = bdata.axisType;
      axisBoundaryType = bdata.axisBoundaryType;
      axisDirection = bdata.axisDirection;
      min = bdata.min;
      max = bdata.max;
      step = bdata.step;
      zdim = bdata.zdim;
      subBinningAdditive = bdata.subBinningAdditive;
      subBinningData =
          bdata.subBinningData
              ? std::make_unique<const BinningData>(*bdata.subBinningData)
              : nullptr;
      m_bins = bdata.m_bins;
      m_boundaries = bdata.m_boundaries;
      m_totalBins = bdata.m_totalBins;
      m_totalBoundaries = bdata.m_totalBoundaries;
      // set the correct function pointer
      if (axisType == AxisType::Equidistant) {
        m_functionPtr = &searchEquidistantWithBoundary;
      } else {
        m_functionPtr = &searchInVectorWithBoundary;
      }
    }
    return (*this);
  }

  BinningData() = default;
  ~BinningData() = default;

  /// Equality operator
  ///
  /// @param bData is the binning data to be checked against
  ///
  /// @return a boolean indicating if they are the same
  bool operator==(const BinningData& bData) const {
    return (axisType == bData.axisType &&
            axisBoundaryType == bData.axisBoundaryType &&
            axisDirection == bData.axisDirection && min == bData.min &&
            max == bData.max && step == bData.step && zdim == bData.zdim &&
            ((subBinningData == nullptr && bData.subBinningData == nullptr) ||
             (subBinningData != nullptr && bData.subBinningData != nullptr &&
              (*subBinningData == *bData.subBinningData))) &&
            subBinningAdditive == bData.subBinningAdditive);
  }

  /// Return the number of bins - including sub bins
  std::size_t bins() const { return m_totalBins; }

  /// Return the boundaries  - including sub boundaries
  /// @return vector of floats indicating the boundary values
  const std::vector<float>& boundaries() const {
    if (subBinningData) {
      return m_totalBoundaries;
    }
    return m_boundaries;
  }

  /// Take the right float value
  ///
  /// @param lposition assumes the correct local position expression
  ///
  /// @return float value according to the binning setup
  float value(const Vector2& lposition) const {
    // ordered after occurrence
    if (axisDirection == AxisDirection::AxisR ||
        axisDirection == AxisDirection::AxisRPhi ||
        axisDirection == AxisDirection::AxisX ||
        axisDirection == AxisDirection::AxisTheta) {
      return lposition[0];
    }

    return lposition[1];
  }

  /// Take the right float value
  ///
  /// @param position is the global position
  ///
  /// @return float value according to the binning setup
  float value(const Vector3& position) const {
    using VectorHelpers::eta;
    using VectorHelpers::perp;
    using VectorHelpers::phi;
    // ordered after occurrence
    if (axisDirection == AxisDirection::AxisR ||
        axisDirection == AxisDirection::AxisTheta) {
      return (perp(position));
    }
    if (axisDirection == AxisDirection::AxisRPhi) {
      return (perp(position) * phi(position));
    }
    if (axisDirection == AxisDirection::AxisEta) {
      return (eta(position));
    }
    if (toUnderlying(axisDirection) < 3) {
      return static_cast<float>(position[toUnderlying(axisDirection)]);
    }
    // phi gauging
    return phi(position);
  }

  /// Get the center value of a bin
  ///
  /// @param bin is the bin for which the center value is requested
  ///
  /// @return float value according to the bin center
  float center(std::size_t bin) const {
    const std::vector<float>& bvals = boundaries();
    // take the center between bin boundaries
    float value =
        bin < (bvals.size() - 1) ? 0.5 * (bvals[bin] + bvals[bin + 1]) : 0.;
    return value;
  }

  /// Get the width of a bin
  ///
  /// @param bin is the bin for which the width is requested
  ///
  /// @return float value of width
  float width(std::size_t bin) const {
    const std::vector<float>& bvals = boundaries();
    // take the center between bin boundaries
    float value = bin < (bvals.size() - 1) ? bvals[bin + 1] - bvals[bin] : 0.;
    return value;
  }

  /// Check if bin is inside from Vector3
  ///
  /// @param position is the search position in global coordinated
  ///
  /// @return boolean if this is inside() method is true
  bool inside(const Vector3& position) const {
    // closed one is always inside
    if (axisBoundaryType == AxisBoundaryType::Closed) {
      return true;
    }
    // all other axisBoundaryTypes
    // @todo remove hard-coded tolerance parameters
    float val = value(position);
    return (val > min - 0.001 && val < max + 0.001);
  }

  /// Check if bin is inside from Vector2
  ///
  /// @param lposition is the search position in global coordinated
  ///
  /// @return boolean if this is inside() method is true
  bool inside(const Vector2& lposition) const {
    // closed one is always inside
    if (axisBoundaryType == AxisBoundaryType::Closed) {
      return true;
    }
    // all other axisBoundaryTypes
    // @todo remove hard-coded tolerance parameters
    float val = value(lposition);
    return (val > min - 0.001 && val < max + 0.001);
  }

  /// Generic search from a 2D position
  /// -- corresponds to local coordinate schema
  /// @param lposition is the search position in local coordinated
  ///
  /// @return bin according tot this
  std::size_t searchLocal(const Vector2& lposition) const {
    if (zdim) {
      return 0;
    }
    return search(value(lposition));
  }

  /// Generic search from a 3D position
  /// -- corresponds to global coordinate schema
  /// @param position is the search position in global coordinated
  ///
  /// @return bin according tot this
  std::size_t searchGlobal(const Vector3& position) const {
    if (zdim) {
      return 0;
    }
    return search(value(position));
  }

  /// Generic search - forwards to correct function pointer
  ///
  /// @param value is the searchvalue as float
  ///
  /// @return bin according tot this
  std::size_t search(float value) const {
    if (zdim) {
      return 0;
    }
    assert(m_functionPtr != nullptr);
    return (!subBinningData) ? (*m_functionPtr)(value, *this)
                             : searchWithSubStructure(value);
  }

  ///  Generic search with sub structure
  /// - forwards to correct function pointer
  ///
  /// @param value is the searchvalue as float
  ///
  /// @return bin according tot this
  std::size_t searchWithSubStructure(float value) const {
    // find the masterbin with the correct function pointer
    std::size_t masterbin = (*m_functionPtr)(value, *this);
    // additive sub binning -
    if (subBinningAdditive) {
      // no gauging done, for additive sub structure
      return masterbin + subBinningData->search(value);
    }
    // gauge the value to the subBinData
    float gvalue =
        value - masterbin * (subBinningData->max - subBinningData->min);
    // now go / additive or multiplicative
    std::size_t subbin = subBinningData->search(gvalue);
    // now return
    return masterbin * subBinningData->bins() + subbin;
  }

  /// Layer next direction is needed
  ///
  /// @param position is the start search position
  /// @param dir is the direction
  /// @todo check if this can be changed
  ///
  /// @return integer that indicates which direction to move
  int nextDirection(const Vector3& position, const Vector3& dir) const {
    if (zdim) {
      return 0;
    }
    float val = value(position);
    Vector3 probe = position + dir.normalized();
    float nextval = value(probe);
    return (nextval > val) ? 1 : -1;
  }

  /// access to the center value
  /// this uses the bin boundary vector, it also works with sub structure
  ///
  /// @param bin is the bin for which the value is requested, if bin > nbins
  /// it is set to max
  ///
  /// @return the center value of the bin is given
  float centerValue(std::size_t bin) const {
    if (zdim) {
      return 0.5 * (min + max);
    }
    float bmin = m_boundaries[bin];
    float bmax = bin < m_boundaries.size() ? m_boundaries[bin + 1] : max;
    return 0.5 * (bmin + bmax);
  }

 private:
  std::size_t m_bins{};             ///< number of bins
  std::vector<float> m_boundaries;  ///< vector of holding the bin boundaries
  std::size_t m_totalBins{};        ///< including potential substructure
  std::vector<float> m_totalBoundaries;  ///< including potential substructure

  std::size_t (*m_functionPtr)(float,
                               const BinningData&){};  /// function pointer

  /// helper method to set the sub structure
  void checkSubStructure() {
    // sub structure is only checked when sBinData is defined
    if (subBinningData) {
      m_totalBoundaries.clear();
      // (A) additive sub structure
      if (subBinningAdditive) {
        // one bin is replaced by the sub bins
        m_totalBins = m_bins + subBinningData->bins() - 1;
        // the tricky one - exchange one bin by many others
        m_totalBoundaries.reserve(m_totalBins + 1);
        // get the sub bin boundaries
        const std::vector<float>& subBinBoundaries =
            subBinningData->boundaries();
        float sBinMin = subBinBoundaries[0];
        // get the min value of the sub bin boundaries
        std::vector<float>::const_iterator maDir = m_boundaries.begin();
        for (; maDir != m_boundaries.end(); ++maDir) {
          // should define numerically stable
          if (std::abs((*maDir) - sBinMin) < 10e-10) {
            // copy the sub bin boundaries into the vector
            m_totalBoundaries.insert(m_totalBoundaries.begin(),
                                     subBinBoundaries.begin(),
                                     subBinBoundaries.end());
            ++maDir;
          } else {
            m_totalBoundaries.push_back(*maDir);
          }
        }
      } else {  // (B) multiplicative sub structure
        // every bin is just replaced by the sub binning structure
        m_totalBins = m_bins * subBinningData->bins();
        m_totalBoundaries.reserve(m_totalBins + 1);
        // get the sub bin boundaries if there are any
        const std::vector<float>& subBinBoundaries =
            subBinningData->boundaries();
        // create the boundary vector
        m_totalBoundaries.push_back(min);
        for (std::size_t ib = 0; ib < m_bins; ++ib) {
          float offset = ib * step;
          for (std::size_t isb = 1; isb < subBinBoundaries.size(); ++isb) {
            m_totalBoundaries.push_back(offset + subBinBoundaries[isb]);
          }
        }
      }
      // sort the total boundary vector
      std::ranges::sort(m_totalBoundaries);
    }
  }

  // Equidistant search
  // - fastest method
  static std::size_t searchEquidistantWithBoundary(float value,
                                                   const BinningData& bData) {
    // vanilla

    int bin = static_cast<int>((value - bData.min) / bData.step);
    // special treatment of the 0 bin for closed
    if (bData.axisBoundaryType == AxisBoundaryType::Closed) {
      if (value < bData.min) {
        return (bData.m_bins - 1);
      }
      if (value > bData.max) {
        return 0;
      }
    }
    // if outside boundary : return boundary for open, opposite bin for closed
    bin = bin < 0 ? ((bData.axisBoundaryType == AxisBoundaryType::Bound)
                         ? 0
                         : (bData.m_bins - 1))
                  : bin;
    return static_cast<std::size_t>(
        (bin <= static_cast<int>(bData.m_bins - 1))
            ? bin
            : ((bData.axisBoundaryType == AxisBoundaryType::Bound)
                   ? (bData.m_bins - 1)
                   : 0));
  }

  // Search in arbitrary boundary
  static std::size_t searchInVectorWithBoundary(float value,
                                                const BinningData& bData) {
    // lower boundary
    if (value <= bData.m_boundaries[0]) {
      return (bData.axisBoundaryType == AxisBoundaryType::Closed)
                 ? (bData.m_bins - 1)
                 : 0;
    }
    // higher boundary
    if (value >= bData.max) {
      return (bData.axisBoundaryType == AxisBoundaryType::Closed)
                 ? 0
                 : (bData.m_bins - 1);
    }

    auto lb = std::lower_bound(bData.m_boundaries.begin(),
                               bData.m_boundaries.end(), value);
    return static_cast<std::size_t>(
        std::distance(bData.m_boundaries.begin(), lb) - 1);
  }

 public:
  /// String screen output method
  /// @param indent the current indentation
  /// @return a string containing the screen information
  std::string toString(const std::string& indent = "") const {
    std::stringstream sl;
    sl << indent << "BinningData object:" << '\n';
    sl << indent << "  - type       : " << axisTypeToString(axisType) << '\n';
    sl << indent << "  - axisBoundaryType     : "
       << axisBoundaryTypeToString(axisBoundaryType) << '\n';
    sl << indent << "  - value      : " << axisDirectionToString(axisDirection)
       << '\n';
    sl << indent << "  - bins       : " << bins() << '\n';
    sl << indent << "  - min/max    : " << min << " / " << max << '\n';
    if (axisType == AxisType::Equidistant) {
      sl << indent << "  - step       : " << step << '\n';
    }
    sl << indent << "  - boundaries : | ";
    for (const auto& b : boundaries()) {
      sl << b << " | ";
    }
    sl << '\n';
    return sl.str();
  }
};
}  // namespace Acts
