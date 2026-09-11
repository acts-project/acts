// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <functional>
#include <vector>

namespace Acts::Experimental {

/// @brief helper class for adaptive traversal of HT space
class HoughAccumulatorSection {
 public:
  /// @brief Defines the fate of section during traversal
  enum class Decision {
    Discard,  //< the section is not to be explored further
    Accept,   //< the section should be accepted as solution without further
              //< exploration
    Drill,    //< the section should be expred further by splitting according to
            //< binning definition (split into 4 or 2 left-right or top-bottom)
    DrillAndExpand,  //< the section should be source of subsections as in the
                     //< case of drill & but the sections will be made larger
    //< size increase is configured in opt by relative factors @see expandX, @see expandY
  };

  HoughAccumulatorSection() = default;

  /// @brief Construct the section
  /// @param xw - width in x direction
  /// @param yw - widths in y direction
  /// @param xBegin - location of left side of the section
  /// @param yBegin - location of bottom side of the section
  /// @param div - division level (will be incremented during each division)
  /// @param indices - indices of measurements that generate lines passing this section
  /// @param history - storage for arbitrary information related to the section
  HoughAccumulatorSection(float xw, float yw, float xBegin, float yBegin,
                          int div = 0,
                          const std::vector<std::uint32_t> &indices = {},
                          const std::vector<float> &history = {});

  /// @brief decision designated for this section in next traversal step
  /// @return currently set decision
  Decision decision() const { return m_decision; }

  /// @brief updates decision designated for this step
  /// If this is not set specifically the default is to split (called Drill)
  /// @param d - new decision for the section
  /// @return value of the decision set
  Decision updateDecision(Decision d) { return m_decision = d; }

  /// @brief keep indices and update parameters of the box
  /// This method is useful when changing direction of the search
  /// @param xw - new x width
  /// @param yw - new x width
  /// @param xBegin - new left side
  /// @param yBegin - new bottom side
  void updateDimensions(float xw, float yw, float xBegin, float yBegin);

  /// @brief keep indices and update parameters of the box by scalling
  /// @param xs - scale in x direction, if bigger than 1 the size increases
  /// @param ys - scale in y direction
  /// The box is recentred
  void expand(float xs, float ys);

  /// @brief indices of measurements that result in lines in this section
  /// @return reference to the vector of indices
  inline const std::vector<std::uint32_t> &indices() const { return m_indices; }

  /// @brief mutable version of @see indices
  /// @return mutable access to indices
  inline std::vector<std::uint32_t> &indices() { return m_indices; }

  /// @brief number of lines in the section
  /// @return number of lines
  std::size_t count() const { return m_indices.size(); }

  /// @brief create bottom section that is result of splitting this one into 2
  /// @param copyIndices - copies indices from the parent if true
  /// +------+
  /// |      |
  /// +------+
  /// |     <|-- this part
  /// +------+
  /// @return - new section
  HoughAccumulatorSection bottom(bool copyIndices = false) const;

  /// @brief create top section that is result of splitting this one into 2
  /// @param copyIndices - copies indices from the parent if true
  /// +------+
  /// |      <|-- this part
  /// +------+
  /// |      |
  /// +------+
  /// @return - new section
  HoughAccumulatorSection top(bool copyIndices = false) const;

  /// @brief create left section that is result of splitting this one into 2
  /// @param copyIndices - copies indices from the parent if true
  /// +---+---+
  /// |  <|---|-- this part
  /// +---+---+
  /// @return - new section
  HoughAccumulatorSection left(bool copyIndices = false) const;

  /// @brief create right section that is result of splitting this one into 2
  /// @param copyIndices - copies indices from the parent if true
  /// +---+---+
  /// |   |  <|-- this part
  /// +---+---+
  /// @return - new section
  HoughAccumulatorSection right(bool copyIndices = false) const;

  /// @brief create section that is result of splitting this one into 4
  /// @arg copyIndices - copies indices from the parent
  /// create section that is bottom left corner of this this one
  /// by default the section is divided into 4 quadrants,
  /// if parameters are provided the quadrants size can be adjusted
  /// +---+---+
  /// |   |   |
  /// +---+---+
  /// |   |  <|-- this part
  /// +---+---+
  /// @return new accumulator section (with new dimensions and location)
  /// @param copyIndices - copies indices from the parent if true
  HoughAccumulatorSection bottomRight(bool copyIndices = false) const;

  /// @brief create section that is result of splitting this one into 4
  /// @param copyIndices - copies indices from the parent
  /// @see bottomRight
  /// @return new accumulator section (with new dimensions and location)
  HoughAccumulatorSection bottomLeft(bool copyIndices = false) const;

  /// @brief create section that is result of splitting this one into 4
  /// @param copyIndices - copies indices from the parent
  /// @see bottomRight
  /// @return new accumulator section (with new dimensions and location)
  HoughAccumulatorSection topLeft(bool copyIndices = false) const;

  /// @brief create section that is result of splitting this one into 4
  /// @param copyIndices - copies indices from the parent
  /// @see bottomRight
  /// @return new accumulator section (with new dimensions and location)
  HoughAccumulatorSection topRight(bool copyIndices = false) const;

  /// @brief true if the line defined by given parameters passes the section
  /// @param function is callable used to check crossing at the edges
  /// @return true if the line passes the section
  template <typename F>
  bool isLineInside(F &&function) const &
    requires std::invocable<F, float>;

  /// @brief check if the lines cross inside the section
  /// @param line1 - functional form of line 1
  /// @param line2 - functional form of line 2
  /// @warning note that this function is assuming that these are lines and the derivative is positive.
  /// @return true if the two lines cross in the section
  template <typename F>
  bool isCrossingInside(F &&line1, F &&line2) const &
    requires std::invocable<F, float>;

  /// @brief size accessor
  /// @return size in x direction
  float xSize() const { return m_xSize; }
  /// @brief size accessor
  /// @return size in x direction
  float ySize() const { return m_ySize; }
  /// @brief location accessor
  /// @return left side position
  float xBegin() const { return m_xBegin; }
  /// @brief location accessor
  /// @return bottom side position
  float yBegin() const { return m_yBegin; }
  /// @brief number of divisions that lead to this section
  /// @return integer > 0
  std::uint32_t divisionLevel() const { return m_divisionLevel; }

  /// store additional (arbitrary) info in indexed array
  /// @param index - identifier
  /// @param value - value to store
  void setHistory(std::size_t index, float value) {
    m_history.resize(std::max(index + 1, m_history.size()));
    m_history.at(index) = value;
  }
  /// @brief retrieve history info
  /// @param index - item index
  /// @return value stored by @see setHistory
  float history(std::uint32_t index) const { return m_history.at(index); }

 private:
  Decision m_decision = Decision::Drill;
  float m_xSize = 0;
  float m_ySize = 0;
  float m_xBegin = 0;
  float m_yBegin = 0;
  /// number of times the starting section was already divided
  std::uint32_t m_divisionLevel = 0;
  /// indices of measurements contributing to this section
  std::vector<std::uint32_t> m_indices;
  /// additional record where an arbitrary information can be stored
  std::vector<float> m_history;
};

template <typename F>
inline bool HoughAccumulatorSection::isLineInside(F &&function) const &
  requires std::invocable<F, float>
{
  const float yB = function(m_xBegin);
  const float yE = function(m_xBegin + m_xSize);
  return (yE > yB) ? yB < m_yBegin + m_ySize && yE > m_yBegin
                   : yB > m_yBegin && yE < m_yBegin + m_ySize;
}

template <typename F>
inline bool HoughAccumulatorSection::isCrossingInside(F &&line1,
                                                      F &&line2) const &
  requires std::invocable<F, float>
{
  // this microalgorithm idea is illustrated below
  // section left section right
  // example with crossing
  //                                       |            +2
  // line 1 crossing left section edge     +1          _|
  // left edge mid point                   |_           |
  //                                       |            +1
  // line 2crossing left section           +2           |
  //
  // example with no crossing
  //                                       |            +1
  // line 1 crossing left section edge     +1          _|
  // left edge mid point                   |_           |
  //                                       |            +2
  // line 2crossing left section           +2           |
  // The above covers most of the cases.
  // Additional precautions are made when both lines cross
  // left & right (x) bounds outside of vertical (y) bounds.

  const float xL = xBegin();
  const float xR = xBegin() + xSize();

  // Evaluate both lines at section boundaries
  const float y1L = line1(xL);
  const float y1R = line1(xR);
  const float y2L = line2(xL);
  const float y2R = line2(xR);

  // --- Step 1: Quick rejection based on ordering ---
  // If one line is consistently above the other → no crossing
  const float dL = y1L - y2L;
  const float dR = y1R - y2R;

  if (dL * dR > 0) {
    return false;
  }

  // --- Step 2: Check if either line fully spans inside vertically ---
  const auto checkVerticalOverlap = [this](float y) {
    return yBegin() < y && y < yBegin() + ySize();
  };

  if (checkVerticalOverlap(y1L) && checkVerticalOverlap(y1R)) {
    return true;
  }

  if (checkVerticalOverlap(y2L) && checkVerticalOverlap(y2R)) {
    return true;
  }

  // --- Step 3: Approximate crossing position ---
  // Linear interpolation assuming near-linear behavior
  const float abs_dL = std::abs(dL);
  const float abs_dR = std::abs(dR);

  // Avoid division by zero (parallel & overlapping edge case)
  if (abs_dL + abs_dR == 0.0f) {
    return false;
  }

  const float t = abs_dL / (abs_dL + abs_dR);  // fraction from left
  const float xCross = std::lerp(xL, xR, t);

  // --- Step 4: Check if crossing lies inside vertical bounds ---
  const float yCross1 = line1(xCross);
  const float yCross2 = line2(xCross);
  const float yCross = 0.5f * (yCross1 + yCross2);  // intersection approx

  return checkVerticalOverlap(yCross);
}

/// @brief structure that keeps HT space traversal options
/// @tparam Measurement - measurement type which are used in generating lines in HT space
template <typename Measurement>
struct HoughExplorationOptions {
  /// minimum bin size in x direction, beyond that
  /// value the sections are not split
  float xMinBinSize = 1.0f;

  /// minimum bin size in y direction, beyond that
  /// value the sections are not split
  float yMinBinSize = 1.0f;

  /// expand in x direction (default by 10%) if Expand
  /// Decision is made
  float expandX = 1.1f;

  /// expand in y direction (default by 10%) if Expand
  /// Decision is made
  float expandY = 1.1f;

  /// @brief functional, that given measurement and argument in HT space x return y
  using LineFunctor = std::function<float(const Measurement &, float)>;

  /// functional that, given measurement and
  /// "x" coordinate of Hough space return "y" coordinate
  LineFunctor lineFunctor;
  /// @brief functional type, that given section and measurements decides the evolution of section
  using DecisionFunctor = std::function<HoughAccumulatorSection::Decision(
      const HoughAccumulatorSection &, const std::vector<Measurement> &)>;

  /// function deciding if the accumulator section should be, discarded,
  /// split further (and how), or is a solution
  DecisionFunctor decisionFunctor;
};

/// @brief function that walks through the section splitting them (depth-first) and
/// @tparam Measurement - type of measurements
/// @tparam Functor - for translation from measurements to Hough space
/// @param sectionsStack - the stacks to consider
/// @param measurements - measurements to which indices are kept in HoughAccumulatorSection
/// @param opt - exploration directives
/// @param results - sectoins that satisfied acceptance criteria
template <typename Measurement>
void exploreHoughParametersSpace(
    std::vector<HoughAccumulatorSection> &sectionsStack,
    const std::vector<Measurement> &measurements,
    const HoughExplorationOptions<Measurement> &opt,
    std::vector<HoughAccumulatorSection> &results) {
  while (!sectionsStack.empty()) {
    HoughAccumulatorSection thisSection = std::move(sectionsStack.back());
    sectionsStack.pop_back();

    std::vector<HoughAccumulatorSection> newSections;
    const bool splitX = thisSection.xSize() > opt.xMinBinSize;
    const bool splitY = thisSection.ySize() > opt.yMinBinSize;
    if (splitX && splitY) {
      // Split into 4 sections
      newSections.push_back(thisSection.bottomLeft());
      newSections.push_back(thisSection.topLeft());
      newSections.push_back(thisSection.bottomRight());
      newSections.push_back(thisSection.topRight());
    } else if (!splitX && splitY) {
      // Split into 2 sections horizontally
      newSections.push_back(thisSection.bottom());
      newSections.push_back(thisSection.top());
    } else if (splitX && !splitY) {
      // Split into 2 sections vertically
      newSections.push_back(thisSection.left());
      newSections.push_back(thisSection.right());
    }

    if (thisSection.decision() ==
        HoughAccumulatorSection::Decision::DrillAndExpand) {
      for (HoughAccumulatorSection &s : newSections) {
        s.expand(opt.expandX, opt.expandY);
      }
    }

    for (const std::uint32_t idx : thisSection.indices()) {
      const auto &m = measurements[idx];
      const auto line = [&](float x) { return opt.lineFunctor(m, x); };

      for (HoughAccumulatorSection &s : newSections) {
        if (s.isLineInside(line)) {
          s.indices().push_back(idx);
        }
      }
    }
    for (HoughAccumulatorSection &s : newSections) {
      s.updateDecision(opt.decisionFunctor(s, measurements));
      if (s.decision() == HoughAccumulatorSection::Decision::Accept) {
        results.push_back(std::move(s));
      } else if (s.decision() == HoughAccumulatorSection::Decision::Drill ||
                 s.decision() ==
                     HoughAccumulatorSection::Decision::DrillAndExpand) {
        sectionsStack.push_back(std::move(s));
      }
    }
  }
}

/// @brief Helper for fast check if there is enough line crossings within HoughAccumulatorSection
/// @tparam Measurement
/// @tparam Functor
/// @param section - section to rest against
/// @param measurements - vector of measurements
/// @param lineFunctor - the function defining line
/// @param threshold - number of crossings required for positive decision (it is provided so early exit can be implemented)
/// @return true if the check passes
template <typename Measurement, typename Functor>
bool passIntersectionsCheck(const HoughAccumulatorSection &section,
                            const std::vector<Measurement> &measurements,
                            const Functor &lineFunctor,
                            const std::uint32_t threshold) {
  const std::size_t count = section.count();
  const float xLeft = section.xBegin();
  const float xRight = xLeft + section.xSize();
  const auto &indices = section.indices();

  // Small Buffer Optimization
  constexpr std::size_t kMaxStackLines = 64;

  if (count <= kMaxStackLines) {
    std::array<float, kMaxStackLines> yLeft{};
    std::array<float, kMaxStackLines> yRight{};

    for (std::size_t i = 0; i < count; ++i) {
      const auto &m = measurements[indices[i]];
      yLeft[i] = lineFunctor(m, xLeft);
      yRight[i] = lineFunctor(m, xRight);
    }

    std::uint32_t inside = 0;
    for (std::uint32_t i = 0; i < count; ++i) {
      for (std::uint32_t j = i + 1; j < count; ++j) {
        if ((yLeft[i] - yLeft[j]) * (yRight[i] - yRight[j]) < 0.0f) {
          inside++;
          if (inside >= threshold) {
            return true;  // Early exit
          }
        }
      }
    }
    return false;
  } else {
    std::vector<float> yLeft(count);
    std::vector<float> yRight(count);
    for (std::uint32_t i = 0; i < count; ++i) {
      const auto &m = measurements[indices[i]];
      yLeft[i] = lineFunctor(m, xLeft);
      yRight[i] = lineFunctor(m, xRight);
    }
    std::uint32_t inside = 0;
    for (std::uint32_t i = 0; i < count; ++i) {
      for (std::uint32_t j = i + 1; j < count; ++j) {
        if ((yLeft[i] - yLeft[j]) * (yRight[i] - yRight[j]) < 0.0f) {
          inside++;
          if (inside >= threshold) {
            return true;  // Early exit
          }
        }
      }
    }
    return inside >= threshold;
  }
}

}  // namespace Acts::Experimental
