// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/NeighborHoodIndices.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace Acts {

/// @brief calculate bin indices for an equidistant binning
///
/// This class provides some basic functionality for calculating bin indices
/// for a given equidistant binning.
template <AxisBoundaryType bdt>
class Axis<AxisType::Equidistant, bdt> : public IAxis {
 public:
  /// Static type identifier for this equidistant axis specialization
  static constexpr AxisType type = AxisType::Equidistant;

  /// Divide the range \f$[\text{xmin},\text{xmax})\f$ into \f$\text{nBins}\f$
  /// equidistant bins.
  ///
  /// @param xmin lower boundary of axis range
  /// @param xmax upper boundary of axis range
  /// @param nBins number of bins to divide the axis range into
  /// @param direction optional direction of the axis
  Axis(double xmin, double xmax, std::size_t nBins,
       std::optional<AxisDirection> direction = std::nullopt)
      : IAxis(direction),
        m_min(xmin),
        m_max(xmax),
        m_width((xmax - xmin) / nBins),
        m_bins(nBins) {
    if (m_min >= m_max) {
      std::string msg = "Axis: Invalid axis range'";
      msg += "', min edge (" + std::to_string(m_min) + ") ";
      msg += " needs to be smaller than max edge (";
      msg += std::to_string(m_max) + ").";
      throw std::invalid_argument(msg);
    }
    if (m_bins < 1u) {
      throw std::invalid_argument(
          "Axis: Invalid binning, at least one bin is needed.");
    }
  }

  /// Divide the range \f$[\text{xmin},\text{xmax})\f$ into \f$\text{nBins}\f$
  /// equidistant bins.
  ///
  /// @param typeTag boundary type tag
  /// @param xmin lower boundary of axis range
  /// @param xmax upper boundary of axis range
  /// @param nBins number of bins to divide the axis range into
  /// @param direction optional direction of the axis
  Axis(AxisBoundaryTypeTag<bdt> typeTag, double xmin, double xmax,
       std::size_t nBins, std::optional<AxisDirection> direction = std::nullopt)
      : Axis(xmin, xmax, nBins, direction) {
    static_cast<void>(typeTag);
  }

  /// returns whether the axis is equidistant
  /// @return bool is equidistant
  bool isEquidistant() const final { return true; }

  /// returns whether the axis is variable
  /// @return bool is variable
  bool isVariable() const final { return false; }

  /// returns the type of the axis
  /// @return @c AxisType of this axis
  AxisType getType() const final { return type; }

  /// returns the boundary type set in the template param
  /// @return @c AxisBoundaryType of this axis
  AxisBoundaryType getBoundaryType() const final { return bdt; }

  /// Get #size bins which neighbor the one given. Generic overload with
  /// symmetric size.
  /// @param idx requested bin index
  /// @param size how many neighboring bins (up/down)
  /// @return Set of neighboring bin indices (global)
  NeighborHoodIndices neighborHoodIndices(std::size_t idx,
                                          std::size_t size = 1) const {
    return neighborHoodIndices(idx,
                               std::make_pair(-static_cast<int>(size), size));
  }

  /// Get #size bins which neighbor the one given. This is the version for Open.
  /// @param idx requested bin index
  /// @param sizes how many neighboring bins (up/down)
  /// @return Set of neighboring bin indices (global)
  /// @note Open varies given bin and allows 0 and NBins+1 (underflow, overflow) as neighbors
  NeighborHoodIndices neighborHoodIndices(std::size_t idx,
                                          std::pair<int, int> sizes = {-1,
                                                                       1}) const
    requires(bdt == AxisBoundaryType::Open)
  {
    constexpr int min = 0;
    const int max = getNBins() + 1;
    const int itmin = std::clamp(static_cast<int>(idx + sizes.first), min, max);
    const int itmax =
        std::clamp(static_cast<int>(idx + sizes.second), min, max);
    return NeighborHoodIndices(itmin, itmax + 1);
  }

  /// Get #size bins which neighbor the one given. This is the version for
  /// Bound.
  /// @param idx requested bin index
  /// @param sizes how many neighboring bins (up/down)
  /// @return Set of neighboring bin indices (global)
  /// @note Bound varies given bin and allows 1 and NBins (regular bins) as neighbors
  NeighborHoodIndices neighborHoodIndices(std::size_t idx,
                                          std::pair<int, int> sizes = {-1,
                                                                       1}) const
    requires(bdt == AxisBoundaryType::Bound)
  {
    if (idx <= 0 || idx >= (getNBins() + 1)) {
      return NeighborHoodIndices();
    }
    constexpr int min = 1;
    const int max = getNBins();
    const int itmin = std::clamp(static_cast<int>(idx) + sizes.first, min, max);
    const int itmax =
        std::clamp(static_cast<int>(idx) + sizes.second, min, max);
    return NeighborHoodIndices(itmin, itmax + 1);
  }

  /// Get #size bins which neighbor the one given. This is the version for
  /// Closed (i.e. Wrapping).
  /// @param idx requested bin index
  /// @param sizes how many neighboring bins (up/down)
  /// @return Set of neighboring bin indices (global)
  /// @note Closed varies given bin and allows bins on the opposite
  ///       side of the axis as neighbors. (excludes underflow / overflow)
  NeighborHoodIndices neighborHoodIndices(std::size_t idx,
                                          std::pair<int, int> sizes = {-1,
                                                                       1}) const
    requires(bdt == AxisBoundaryType::Closed)
  {
    // Handle invalid indices
    if (idx <= 0 || idx >= (getNBins() + 1)) {
      return NeighborHoodIndices();
    }

    // Handle corner case where user requests more neighbours than the number
    // of bins on the axis. All bins are returned in this case.

    const int max = getNBins();
    sizes.first = std::clamp(sizes.first, -max, max);
    sizes.second = std::clamp(sizes.second, -max, max);
    if (std::abs(sizes.first - sizes.second) >= max) {
      sizes.first = 1 - idx;
      sizes.second = max - idx;
    }

    // If the entire index range is not covered, we must wrap the range of
    // targeted neighbor indices into the range of valid bin indices. This may
    // split the range of neighbor indices in two parts:
    //
    // Before wraparound - [        XXXXX]XXX
    // After wraparound  - [ XXXX   XXXX ]
    //
    const int itmin = idx + sizes.first;
    const int itmax = idx + sizes.second;
    const std::size_t itfirst = wrapBin(itmin);
    const std::size_t itlast = wrapBin(itmax);
    if (itfirst <= itlast) {
      return NeighborHoodIndices(itfirst, itlast + 1);
    } else {
      return NeighborHoodIndices(itfirst, max + 1, 1, itlast + 1);
    }
  }

  /// Converts bin index into a valid one for this axis.
  /// @note Open: bin index is clamped to [0, nBins+1]
  /// @param bin The bin to wrap
  /// @return valid bin index
  std::size_t wrapBin(int bin) const
    requires(bdt == AxisBoundaryType::Open)
  {
    return std::max(std::min(bin, static_cast<int>(getNBins()) + 1), 0);
  }

  /// Converts bin index into a valid one for this axis.
  /// @note Bound: bin index is clamped to [1, nBins]
  /// @param bin The bin to wrap
  /// @return valid bin index
  std::size_t wrapBin(int bin) const
    requires(bdt == AxisBoundaryType::Bound)
  {
    return std::max(std::min(bin, static_cast<int>(getNBins())), 1);
  }

  /// Converts bin index into a valid one for this axis.
  /// @note Closed: bin index wraps around to other side
  /// @param bin The bin to wrap
  /// @return valid bin index
  std::size_t wrapBin(int bin) const
    requires(bdt == AxisBoundaryType::Closed)
  {
    const int w = getNBins();
    return 1 + (w + ((bin - 1) % w)) % w;
    // return int(bin<1)*w - int(bin>w)*w + bin;
  }

  /// get corresponding bin index for given coordinate
  /// @param x input coordinate
  /// @return index of bin containing the given value
  /// @note Bin intervals are defined with closed lower bounds and open upper
  ///       bounds, that is \f$l <= x < u\f$ if the value @c x lies within a
  ///       bin with lower bound @c l and upper bound @c u.
  /// @note Bin indices start at @c 1. The underflow bin has the index @c 0
  ///       while the index <tt>nBins + 1</tt> indicates the overflow bin.
  std::size_t getBin(double x) const final {
    return wrapBin(
        static_cast<int>(std::floor((x - getMin()) / getBinWidth()) + 1));
  }

  /// get bin width
  /// @return constant width for all bins
  double getBinWidth(std::size_t /*bin*/) const final { return m_width; }

  /// get bin width
  /// @return constant width for all bins
  double getBinWidth() const { return getBinWidth(0); }

  /// get lower bound of bin
  /// @param bin index of bin
  /// @return lower bin boundary
  ///
  /// @pre @c bin must be a valid bin index (excluding the underflow bin),
  ///      i.e. \f$1 \le \text{bin} \le \text{nBins} + 1\f$
  ///
  /// @note Bin intervals have a closed lower bound, i.e. the lower boundary
  ///       belongs to the bin with the given bin index.
  double getBinLowerBound(std::size_t bin) const final {
    return getMin() + (bin - 1) * getBinWidth();
  }

  /// get upper bound of bin
  /// @param bin index of bin
  /// @return upper bin boundary
  /// @pre @c bin must be a valid bin index (excluding the overflow bin),
  ///      i.e. \f$0 \le \text{bin} \le \text{nBins}\f$
  /// @note Bin intervals have an open upper bound, i.e. the upper boundary
  ///       does @b not belong to the bin with the given bin index.
  double getBinUpperBound(std::size_t bin) const final {
    return getMin() + bin * getBinWidth();
  }

  /// get bin center
  /// @param bin index of bin
  /// @return bin center position
  /// @pre @c bin must be a valid bin index (excluding under-/overflow bins),
  ///      i.e. \f$1 \le \text{bin} \le \text{nBins}\f$
  double getBinCenter(std::size_t bin) const final {
    return getMin() + (bin - 0.5) * getBinWidth();
  }

  /// get maximum of binning range
  /// @return maximum of binning range
  double getMax() const final { return m_max; }

  /// get minimum of binning range
  /// @return minimum of binning range
  double getMin() const final { return m_min; }

  /// get total number of bins
  /// @return total number of bins (excluding under-/overflow bins)
  std::size_t getNBins() const final { return m_bins; }

  /// check whether value is inside axis limits
  /// @param x The value to check
  /// @return @c true if \f$\text{xmin} \le x < \text{xmax}\f$, otherwise
  ///         @c false
  /// @post If @c true is returned, the bin containing the given value is a
  ///       valid bin, i.e. it is neither the underflow nor the overflow bin.
  bool isInside(double x) const final { return (m_min <= x) && (x < m_max); }

  /// Return a vector of bin edges
  /// @return Vector which contains the bin edges
  std::vector<double> getBinEdges() const final {
    std::vector<double> binEdges;
    for (std::size_t i = 1; i <= m_bins; i++) {
      binEdges.push_back(getBinLowerBound(i));
    }
    binEdges.push_back(getBinUpperBound(m_bins));
    return binEdges;
  }

  friend std::ostream& operator<<(std::ostream& os, const Axis& axis) {
    os << "Axis<Equidistant, " << bdt << ">(";
    os << axis.m_min << ", ";
    os << axis.m_max << ", ";
    os << axis.m_bins << ", ";
    if (axis.getDirection().has_value()) {
      os << *axis.getDirection();
    } else {
      os << "Undefined";
    }
    os << ")";
    return os;
  }

 protected:
  void toStream(std::ostream& os) const final { os << *this; }

 private:
  /// minimum of binning range
  double m_min{};
  /// maximum of binning range
  double m_max{};
  /// constant bin width
  double m_width{};
  /// number of bins (excluding under-/overflow bins)
  std::size_t m_bins{};
};

/// calculate bin indices for a variable binning
///
/// This class provides some basic functionality for calculating bin indices
/// for a given binning with variable bin sizes.
template <AxisBoundaryType bdt>
class Axis<AxisType::Variable, bdt> : public IAxis {
 public:
  /// Static type identifier for this variable-width axis specialization
  static constexpr AxisType type = AxisType::Variable;

  /// Create a binning structure with @c nBins variable-sized bins from the
  /// given bin boundaries. @c nBins is given by the number of bin edges
  /// reduced by one.
  /// @param binEdges vector of bin edges
  /// @param direction optional direction of the axis
  /// @pre @c binEdges must be strictly sorted in ascending order.
  /// @pre @c binEdges must contain at least two entries.
  explicit Axis(std::vector<double> binEdges,
                std::optional<AxisDirection> direction = std::nullopt)
      : IAxis(direction), m_binEdges(std::move(binEdges)) {
    if (m_binEdges.size() < 2) {
      throw std::invalid_argument(
          "Axis: Invalid binning, at least two bin edges are needed.");
    }
    if (!std::ranges::is_sorted(m_binEdges)) {
      throw std::invalid_argument(
          "Axis: Invalid binning, bin edges are not sorted.");
    }
  }

  /// Create a binning structure with @c nBins variable-sized bins from the
  /// given bin boundaries. @c nBins is given by the number of bin edges
  /// reduced by one.
  /// @param typeTag boundary type tag
  /// @param binEdges vector of bin edges
  /// @param direction optional direction of the axis
  /// @pre @c binEdges must be strictly sorted in ascending order.
  /// @pre @c binEdges must contain at least two entries.
  Axis(AxisBoundaryTypeTag<bdt> typeTag, std::vector<double> binEdges,
       std::optional<AxisDirection> direction = std::nullopt)
      : Axis(std::move(binEdges), direction) {
    static_cast<void>(typeTag);
  }

  /// returns whether the axis is equidistante
  /// @return bool is equidistant
  bool isEquidistant() const final { return false; }

  /// returns whether the axis is variable
  /// @return bool is variable
  bool isVariable() const final { return true; }

  /// returns the type of the axis
  /// @return @c AxisType of this axis
  AxisType getType() const final { return type; }

  /// returns the boundary type set in the template param
  /// @return @c AxisBoundaryType of this axis
  AxisBoundaryType getBoundaryType() const final { return bdt; }

  /// Get #size bins which neighbor the one given. Generic overload with
  /// symmetric size.
  /// @param idx requested bin index
  /// @param size how many neighboring bins
  /// @return Set of neighboring bin indices (global)
  NeighborHoodIndices neighborHoodIndices(std::size_t idx,
                                          std::size_t size = 1) const {
    return neighborHoodIndices(idx,
                               std::make_pair(-static_cast<int>(size), size));
  }

  /// Get #size bins which neighbor the one given. This is the version for Open.
  /// @param idx requested bin index
  /// @param sizes how many neighboring bins (up/down)
  /// @return Set of neighboring bin indices (global)
  /// @note Open varies given bin and allows 0 and NBins+1 (underflow, overflow) as neighbors
  NeighborHoodIndices neighborHoodIndices(std::size_t idx,
                                          std::pair<int, int> sizes = {-1,
                                                                       1}) const
    requires(bdt == AxisBoundaryType::Open)
  {
    constexpr int min = 0;
    const int max = getNBins() + 1;
    const int itmin = std::max(min, static_cast<int>(idx) + sizes.first);
    const int itmax = std::min(max, static_cast<int>(idx) + sizes.second);
    return NeighborHoodIndices(itmin, itmax + 1);
  }

  /// Get #size bins which neighbor the one given. This is the version for
  /// Bound.
  /// @param idx requested bin index
  /// @param sizes how many neighboring bins (up/down)
  /// @return Set of neighboring bin indices (global)
  /// @note Bound varies given bin and allows 1 and NBins (regular bins) as neighbors
  NeighborHoodIndices neighborHoodIndices(std::size_t idx,
                                          std::pair<int, int> sizes = {-1,
                                                                       1}) const
    requires(bdt == AxisBoundaryType::Bound)
  {
    if (idx <= 0 || idx >= (getNBins() + 1)) {
      return NeighborHoodIndices();
    }
    constexpr int min = 1;
    const int max = getNBins();
    const int itmin = std::max(min, static_cast<int>(idx) + sizes.first);
    const int itmax = std::min(max, static_cast<int>(idx) + sizes.second);
    return NeighborHoodIndices(itmin, itmax + 1);
  }

  /// Get #size bins which neighbor the one given. This is the version for
  /// Closed.
  /// @param idx requested bin index
  /// @param sizes how many neighboring bins (up/down)
  /// @return Set of neighboring bin indices (global)
  /// @note Closed varies given bin and allows bins on the opposite
  ///       side of the axis as neighbors. (excludes underflow / overflow)
  NeighborHoodIndices neighborHoodIndices(std::size_t idx,
                                          std::pair<int, int> sizes = {-1,
                                                                       1}) const
    requires(bdt == AxisBoundaryType::Closed)
  {
    // Handle invalid indices
    if (idx <= 0 || idx >= (getNBins() + 1)) {
      return NeighborHoodIndices();
    }

    // Handle corner case where user requests more neighbours than the number
    // of bins on the axis. All bins are returned in this case

    const int max = getNBins();
    sizes.first = std::clamp(sizes.first, -max, max);
    sizes.second = std::clamp(sizes.second, -max, max);
    if (std::abs(sizes.first - sizes.second) >= max) {
      sizes.first = 1 - idx;
      sizes.second = max - idx;
    }

    // If the entire index range is not covered, we must wrap the range of
    // targeted neighbor indices into the range of valid bin indices. This may
    // split the range of neighbor indices in two parts:
    //
    // Before wraparound - [        XXXXX]XXX
    // After wraparound  - [ XXXX   XXXX ]
    //
    const int itmin = idx + sizes.first;
    const int itmax = idx + sizes.second;
    const std::size_t itfirst = wrapBin(itmin);
    const std::size_t itlast = wrapBin(itmax);
    if (itfirst <= itlast) {
      return NeighborHoodIndices(itfirst, itlast + 1);
    } else {
      return NeighborHoodIndices(itfirst, max + 1, 1, itlast + 1);
    }
  }

  /// Converts bin index into a valid one for this axis.
  /// @note Open: bin index is clamped to [0, nBins+1]
  /// @param bin The bin to wrap
  /// @return valid bin index
  std::size_t wrapBin(int bin) const
    requires(bdt == AxisBoundaryType::Open)
  {
    return std::max(std::min(bin, static_cast<int>(getNBins()) + 1), 0);
  }

  /// Converts bin index into a valid one for this axis.
  /// @note Bound: bin index is clamped to [1, nBins]
  /// @param bin The bin to wrap
  /// @return valid bin index
  std::size_t wrapBin(int bin) const
    requires(bdt == AxisBoundaryType::Bound)
  {
    return std::max(std::min(bin, static_cast<int>(getNBins())), 1);
  }

  /// Converts bin index into a valid one for this axis.
  /// @note Closed: bin index wraps around to other side
  /// @param bin The bin to wrap
  /// @return valid bin index
  std::size_t wrapBin(int bin) const
    requires(bdt == AxisBoundaryType::Closed)
  {
    const int w = getNBins();
    return 1 + (w + ((bin - 1) % w)) % w;
    // return int(bin<1)*w - int(bin>w)*w + bin;
  }

  /// get corresponding bin index for given coordinate
  /// @param x input coordinate
  /// @return index of bin containing the given value
  /// @note Bin intervals are defined with closed lower bounds and open upper
  ///       bounds, that is \f$l <= x < u\f$ if the value @c x lies within a
  ///       bin with lower bound @c l and upper bound @c u.
  /// @note Bin indices start at @c 1. The underflow bin has the index @c 0
  ///       while the index <tt>nBins + 1</tt> indicates the overflow bin .
  std::size_t getBin(double x) const final {
    const auto it = std::ranges::upper_bound(m_binEdges, x);
    return wrapBin(
        static_cast<int>(std::ranges::distance(m_binEdges.begin(), it)));
  }

  /// get bin width
  /// @param bin index of bin
  /// @return width of given bin
  /// @pre @c bin must be a valid bin index (excluding under-/overflow bins),
  ///      i.e. \f$1 \le \text{bin} \le \text{nBins}\f$
  double getBinWidth(std::size_t bin) const final {
    return m_binEdges.at(bin) - m_binEdges.at(bin - 1);
  }

  /// get lower bound of bin
  /// @param bin index of bin
  /// @return lower bin boundary
  /// @pre @c bin must be a valid bin index (excluding the underflow bin),
  ///      i.e. \f$1 \le \text{bin} \le \text{nBins} + 1\f$
  /// @note Bin intervals have a closed lower bound, i.e. the lower boundary
  ///       belongs to the bin with the given bin index.
  double getBinLowerBound(std::size_t bin) const final {
    return m_binEdges.at(bin - 1);
  }

  /// get upper bound of bin
  /// @param bin index of bin
  /// @return upper bin boundary
  /// @pre @c bin must be a valid bin index (excluding the overflow bin),
  ///      i.e. \f$0 \le \text{bin} \le \text{nBins}\f$
  /// @note Bin intervals have an open upper bound, i.e. the upper boundary
  ///       does @b not belong to the bin with the given bin index.
  double getBinUpperBound(std::size_t bin) const final {
    return m_binEdges.at(bin);
  }

  /// get bin center
  /// @param bin index of bin
  /// @return bin center position
  /// @pre @c bin must be a valid bin index (excluding under-/overflow bins),
  ///      i.e. \f$1 \le \text{bin} \le \text{nBins}\f$
  double getBinCenter(std::size_t bin) const final {
    return 0.5 * (getBinLowerBound(bin) + getBinUpperBound(bin));
  }

  /// get maximum of binning range
  /// @return maximum of binning range
  double getMax() const final { return m_binEdges.back(); }

  /// get minimum of binning range
  /// @return minimum of binning range
  double getMin() const final { return m_binEdges.front(); }

  /// get total number of bins
  /// @return total number of bins (excluding under-/overflow bins)
  std::size_t getNBins() const final { return m_binEdges.size() - 1; }

  /// check whether value is inside axis limits
  /// @param x The value to check
  /// @return @c true if \f$\text{xmin} \le x < \text{xmax}\f$, otherwise @c false
  /// @post If @c true is returned, the bin containing the given value is a
  ///       valid bin, i.e. it is neither the underflow nor the overflow bin.
  bool isInside(double x) const final {
    return (m_binEdges.front() <= x) && (x < m_binEdges.back());
  }

  /// Return a vector of bin edges
  /// @return Vector which contains the bin edges
  std::vector<double> getBinEdges() const final { return m_binEdges; }

  friend std::ostream& operator<<(std::ostream& os, const Axis& axis) {
    os << "Axis<Variable, " << bdt << ">({";
    os << axis.m_binEdges.front();
    for (std::size_t i = 1; i < axis.m_binEdges.size(); ++i) {
      os << ", " << axis.m_binEdges.at(i);
    }
    os << "}, ";
    if (axis.getDirection().has_value()) {
      os << *axis.getDirection();
    } else {
      os << "Undefined";
    }
    os << ")";
    return os;
  }

 protected:
  void toStream(std::ostream& os) const final { os << *this; }

 private:
  /// vector of bin edges (sorted in ascending order)
  std::vector<double> m_binEdges;
};

}  // namespace Acts
