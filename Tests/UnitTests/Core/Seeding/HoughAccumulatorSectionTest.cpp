// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Seeding/HoughAccumulatorSection.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <functional>
#include <map>
#include <vector>

using namespace Acts;
using namespace Acts::Experimental;

namespace ActsTests {

auto logger = getDefaultLogger("UnitTests", Logging::VERBOSE);

/// @brief Structure representing test parameters for a line in the accumulator space: y
/// = slope * x + intercept
struct LineParameters {
  float slope;
  float intercept;
};

// Structure for statistics to enable BOOST_CHECK on results
struct Stats {
  double area{};
  std::int32_t nSections{};
  std::int32_t nLines{};
  std::int32_t discardedByThresholdCut{};
  std::int32_t discardedByCrossingCut{};
};

template <typename measurement_t = LineParameters>
struct TestExplorationOptions : public HoughExplorationOptions<LineParameters> {
  std::uint32_t threshold =
      4;  // number of lines passing section for it to be still considered
  std::uint32_t noiseThreshold = 12;  // number of lines passing section at the
                                      // final split to consider it noise
  std::uint32_t min_division_level = 2;
};

// SUITE 1: HoughAccumulatorSection Test
BOOST_AUTO_TEST_SUITE(HoughAccumulatorSectionSuite)

using namespace ActsTests;

BOOST_AUTO_TEST_CASE(construct) {
  HoughAccumulatorSection s(10., 100., -5., -50.0);
  BOOST_CHECK_EQUAL(s.xSize(), 10.);
  BOOST_CHECK_EQUAL(s.ySize(), 100.);
  BOOST_CHECK_EQUAL(s.xBegin(), -5.);
  BOOST_CHECK_EQUAL(s.yBegin(), -50.);
}

BOOST_AUTO_TEST_CASE(split_vertical) {
  HoughAccumulatorSection s(10., 100., -5., -50.0);
  HoughAccumulatorSection t = s.top();
  HoughAccumulatorSection b = s.bottom();
  BOOST_CHECK_EQUAL(s.xSize(), t.xSize());
  BOOST_CHECK_EQUAL(s.xSize(), b.xSize());
  BOOST_CHECK_EQUAL(0.5f * s.ySize(), t.ySize());
  BOOST_CHECK_EQUAL(0.5f * s.ySize(), b.ySize());
  BOOST_CHECK_EQUAL(s.xBegin(), t.xBegin());
  BOOST_CHECK_EQUAL(s.xBegin(), b.xBegin());
  BOOST_CHECK_EQUAL(t.yBegin(), 0.0);
  BOOST_CHECK_EQUAL(b.yBegin(), -50.0);
}

BOOST_AUTO_TEST_CASE(split_horizontal) {
  HoughAccumulatorSection s(10., 100., -5., -50.0);
  HoughAccumulatorSection l = s.left();
  HoughAccumulatorSection r = s.right();
  BOOST_CHECK_EQUAL(s.ySize(), l.ySize());
  BOOST_CHECK_EQUAL(s.ySize(), r.ySize());
  BOOST_CHECK_EQUAL(0.5f * s.xSize(), l.xSize());
  BOOST_CHECK_EQUAL(0.5f * s.xSize(), r.xSize());
  BOOST_CHECK_EQUAL(l.xBegin(), -5.0);
  BOOST_CHECK_EQUAL(r.xBegin(), 0.0);
}

BOOST_AUTO_TEST_CASE(split_4) {
  HoughAccumulatorSection s(10., 100., 5., 50.0);
  HoughAccumulatorSection tl = s.topLeft();
  HoughAccumulatorSection tr = s.topRight();
  HoughAccumulatorSection bl = s.bottomLeft();
  HoughAccumulatorSection br = s.bottomRight();

  BOOST_CHECK_EQUAL(tl.xBegin(), 5.0);
  BOOST_CHECK_EQUAL(tl.yBegin(), 100.0);
  BOOST_CHECK_EQUAL(tl.xSize(), 5.0);
  BOOST_CHECK_EQUAL(tl.ySize(), 50.0);

  BOOST_CHECK_EQUAL(tr.xBegin(), 10.0);
  BOOST_CHECK_EQUAL(tr.yBegin(), 100.0);
  BOOST_CHECK_EQUAL(tr.xSize(), 5.0);
  BOOST_CHECK_EQUAL(tr.ySize(), 50.0);

  BOOST_CHECK_EQUAL(bl.xBegin(), 5.0);
  BOOST_CHECK_EQUAL(bl.yBegin(), 50.0);
  BOOST_CHECK_EQUAL(bl.xSize(), 5.0);
  BOOST_CHECK_EQUAL(bl.ySize(), 50.0);

  BOOST_CHECK_EQUAL(br.xBegin(), 10.0);
  BOOST_CHECK_EQUAL(br.yBegin(), 50.0);
  BOOST_CHECK_EQUAL(br.xSize(), 5.0);
  BOOST_CHECK_EQUAL(br.ySize(), 50.0);
}

BOOST_AUTO_TEST_CASE(is_line_inside) {
  HoughAccumulatorSection s(10., 2., -5., -1.0);
  BOOST_CHECK(!s.isLineInside([](float x) { return 0.5f * x + 5.f; }));
  BOOST_CHECK(s.isLineInside([](float x) { return 0.5f * x + 3.f; }));
  BOOST_CHECK(s.isLineInside([](float x) { return 0.1f * x; }));
  BOOST_CHECK(s.isLineInside([](float x) { return 1.0f * x - 5.f; }));
  BOOST_CHECK(!s.isLineInside([](float x) { return 0.5f * x - 5.f; }));
}

BOOST_AUTO_TEST_CASE(is_crossing_inside) {
  HoughAccumulatorSection s(10., 6., -5., -6.);
  std::function<float(float)> l1 = [](float x) { return 1.0f * x + 3.f; };
  std::function<float(float)> l2 = [](float x) { return 0.5f * x + 1.f; };
  std::function<float(float)> l3 = [](float x) { return 0.5f * x - 1.f; };
  std::function<float(float)> l4 = [](float x) { return 1.0f * x - 3.f; };

  BOOST_CHECK(s.isCrossingInside(l1, l2));
  BOOST_CHECK(!s.isCrossingInside(l2, l3));
  BOOST_CHECK(!s.isCrossingInside(l1, l3));
  BOOST_CHECK(!s.isCrossingInside(l4, l3));
}

BOOST_AUTO_TEST_CASE(history_handing) {
  HoughAccumulatorSection s;
  s.setHistory(2, 0.6f);  // 2 is index
  s.setHistory(0, 1.0f);  // 0 is index

  BOOST_CHECK_EQUAL(s.history(0), 1.0f);
  BOOST_CHECK_EQUAL(s.history(2), 0.6f);
}

BOOST_AUTO_TEST_SUITE_END()

// SUITE 2: exploreParametersSpace Test
BOOST_AUTO_TEST_SUITE(ExploreParametersSpaceSuite)

BOOST_AUTO_TEST_CASE(test_extra_split_x_min_div_0) {
  // DESCRIPTION:
  // This test verifies the asymmetric binary split logic (!splitX && splitY /
  // splitX && !splitY) We define 3 lines with very low slopes that cross in the
  // center of the X-axis but remain entirely in the lower half of the Y-axis.
  // The grid has a yMinBinSize (10.0) and a smaller xMinBinSize (6.0)

  std::vector<LineParameters> measurements = {
      {0.05f, 1.0f},   // y = 0.05x + 1
      {0.08f, 1.01f},  // y = 0.08x + 1.01
      {-0.10f, 3.0f}   // y = -0.1x + 3
  };

  TestExplorationOptions<LineParameters> opt;
  opt.xMinBinSize = 6.0f;
  opt.yMinBinSize = 10.0f;
  opt.noiseThreshold = 3;
  opt.threshold = 2;
  opt.min_division_level = 0;
  // Linear function definition: y = slope * x + intercept
  opt.lineFunctor = [](const LineParameters &p, float arg) {
    return p.slope * arg + p.intercept;
  };

  HoughAccumulatorSection section(20.0f, 20.0f, 0.0f, 0.0f, 0);
  section.indices() = {0, 1, 2};

  std::map<int, Stats> sStat;

  opt.decisionFunctor = [&sStat, &opt](const HoughAccumulatorSection &sec,
                                       const std::vector<LineParameters> &mes) {
    using enum HoughAccumulatorSection::Decision;

    if (sec.count() < opt.threshold) {
      sStat[sec.divisionLevel()].discardedByThresholdCut += 1;
      return Discard;
    }
    if (sec.count() < 3 * opt.threshold &&
        !passIntersectionsCheck(sec, mes, opt.lineFunctor,
                                opt.threshold * (opt.threshold - 1))) {
      sStat[sec.divisionLevel()].discardedByCrossingCut += 1;
      return Discard;
    }

    sStat[sec.divisionLevel()].area += sec.xSize() * sec.ySize();
    sStat[sec.divisionLevel()].nSections += 1;
    sStat[sec.divisionLevel()].nLines += sec.count();

    if (sec.divisionLevel() <= opt.min_division_level) {
      return Drill;
    }
    if (sec.count() <= opt.noiseThreshold && sec.xSize() <= opt.xMinBinSize &&
        sec.ySize() <= opt.yMinBinSize) {
      return Accept;
    }

    return Drill;
  };

  std::vector<HoughAccumulatorSection> sectionsStack;
  sectionsStack.push_back(std::move(section));
  std::vector<HoughAccumulatorSection> results;

  exploreHoughParametersSpace(sectionsStack, measurements, opt, results);

  BOOST_CHECK(sectionsStack.empty());

  BOOST_REQUIRE_EQUAL(results.size(), 1);
  BOOST_CHECK_EQUAL(results[0].xSize(), 5.0);
  BOOST_CHECK_EQUAL(results[0].ySize(), 10.0);

  // ==========================================
  // TREE EXECUTION TRACE
  // ==========================================

  // LEVEL 0:
  // this section is not tested with decisionFunctor and so the statistics is
  // not collected

  // LEVEL 1: Split into 4 (10x10). The lines are in the bottom half.
  // Section(10, 10, 0, 10) [tL] - Discarded by threshold cut (0 lines).
  // Section(10, 10, 10, 10) [tR] - Discarded by threshold cut (0 lines).
  // Section(10, 10, 0, 0) [bL] - Discarded by crossing cut (has 3 lines, but no
  // intersections inside). Section(10, 10, 10, 0) [bR] - Contains crossings ->
  // Drill.
  BOOST_CHECK_EQUAL(sStat[1].discardedByThresholdCut, 2);
  BOOST_CHECK_EQUAL(sStat[1].discardedByCrossingCut, 1);
  BOOST_CHECK_EQUAL(sStat[1].nSections, 1);
  BOOST_CHECK_EQUAL(sStat[1].area, 100.0);

  // LEVEL 2: bR section splits. ySize(10) is NO longer > yMinBinSize(10).
  // splitX=true, splitY=false. Splits vertically into 2 (5x10).
  // Section(5, 10, 15, 0) [right] - Discarded by crossing cut (intersections
  // are around x=13). Section(5, 10, 10, 0) [left] - Contains crossings.
  // xSize(5) <= xMinBinSize(6). -> Accept
  BOOST_CHECK_EQUAL(sStat[2].nSections, 1);
  BOOST_CHECK_EQUAL(sStat[2].area, 50.0);
  BOOST_CHECK_EQUAL(sStat[2].discardedByThresholdCut, 0);
  BOOST_CHECK_EQUAL(sStat[2].discardedByCrossingCut, 1);
}

BOOST_AUTO_TEST_CASE(test_with_min_div_lvl_is_1) {
  // DESCRIPTION:
  // This test verifies the standard QuadTree exploration (splitX && splitY).
  // We inject 5 direct line parameter sets (y = slope*x + intercept) into a
  // 10x20 area. min_division_level is set to 1, forcing an initial blind split
  // into 4 sections, after which the Early Culling and Small Buffer
  // Optimization (SBO) logic takes over.

  // Basic config
  std::vector<LineParameters> measurements = {
      {5.0, 3.0},   // Line 0
      {4.0, -4.0},  // Lina 1
      {1.0, 2.0},   // Line 2
      {0.5, 7.0},   // Line 3
      {0.2, 1.01}   // Line 4
  };

  TestExplorationOptions<LineParameters> opt;
  opt.xMinBinSize = 2.51;
  opt.yMinBinSize = 2.51;
  opt.noiseThreshold = 3;
  opt.threshold = 2;
  opt.min_division_level = 1;
  // Linear function definition: y = a*x + b, where LineParameters stores (a, b)
  // parameters for test purposes.
  opt.lineFunctor = [](const LineParameters &p, float arg) {
    return p.slope * arg + p.intercept;
  };

  // Create initial section and populate it with indices mapping to our 3
  // measurements
  HoughAccumulatorSection section(10.0, 20.0, -5.0, -10.0);
  section.indices() = {0, 1, 2, 3, 4};

  // Statistics map for verification in tests
  std::map<int, Stats> sStat;

  // Optimized decision functor
  opt.decisionFunctor = [&sStat, &opt](const HoughAccumulatorSection &sec,
                                       const std::vector<LineParameters> &mes) {
    using enum HoughAccumulatorSection::Decision;
    if (sec.count() < opt.threshold) {
      sStat[sec.divisionLevel()].discardedByThresholdCut += 1;
      return Discard;
    }
    if (sec.count() < 3 * opt.threshold &&
        !passIntersectionsCheck(sec, mes, opt.lineFunctor,
                                opt.threshold * (opt.threshold - 1))) {
      sStat[sec.divisionLevel()].discardedByCrossingCut += 1;
      return Discard;
    }
    // Stats to save info about not discarded accumulator space
    sStat[sec.divisionLevel()].area += sec.xSize() * sec.ySize();
    sStat[sec.divisionLevel()].nSections += 1;
    sStat[sec.divisionLevel()].nLines += sec.count();

    if (sec.divisionLevel() <= opt.min_division_level) {
      return Drill;
    }
    if (sec.count() <= opt.noiseThreshold && sec.xSize() <= opt.xMinBinSize &&
        sec.ySize() <= opt.yMinBinSize) {
      std::cerr << "Accepted section: \n";
      return Accept;
    }

    return Drill;
  };

  std::vector<HoughAccumulatorSection> sectionsStack;
  sectionsStack.push_back(std::move(section));  // pushing root node

  std::vector<HoughAccumulatorSection> results;

  // Core engine
  exploreHoughParametersSpace(sectionsStack, measurements, opt, results);

  // Initial Assertions to verify the engine explored the space
  BOOST_CHECK(!sStat.empty());  // Verifies the tree was drilled

  // Empty section stack after algorithm
  BOOST_CHECK(sectionsStack.empty());

  // Results vector tests
  BOOST_REQUIRE_EQUAL(results.size(), 2);
  BOOST_CHECK_EQUAL(results[0].xSize(), 2.5);
  BOOST_CHECK_EQUAL(results[0].count(), 3);
  BOOST_CHECK_EQUAL(results[1].ySize(), 2.5);
  BOOST_CHECK_EQUAL(results[1].count(), 3);

  // ==========================================
  // TREE EXECUTION TRACE
  // ==========================================

  // LEVEL 0:
  // this section is not tested with decisionFunctor and so the statistics is
  // not collected

  // LEVEL 1: Splits into 4 (5.0 x 10.0).
  // 1 section discarded by threshold cut.
  // 1 section discarded by crossing cut.
  // 2 sections survive -> Drill.
  BOOST_CHECK_EQUAL(sStat[1].nSections, 2);
  BOOST_CHECK_EQUAL(sStat[1].area, 100.0);
  BOOST_CHECK_EQUAL(sStat[1].nLines, 9);
  BOOST_CHECK_EQUAL(sStat[1].discardedByThresholdCut, 1);
  BOOST_CHECK_EQUAL(sStat[1].discardedByCrossingCut, 1);

  // LEVEL 2: 2 sections split into 4 each (Total 8 generated sections of 2.5
  // x 5.0). 2 sections discarded by threshold cut. 4 sections discarded by
  // crossing cut. 2 sections survive -> Drill.
  BOOST_CHECK_EQUAL(sStat[2].nSections, 2);
  BOOST_CHECK_EQUAL(sStat[2].area, 25.0);
  BOOST_CHECK_EQUAL(sStat[2].nLines, 7);
  BOOST_CHECK_EQUAL(sStat[2].discardedByThresholdCut, 2);
  BOOST_CHECK_EQUAL(sStat[2].discardedByCrossingCut, 4);

  // LEVEL 3: 2 sections split. xSize(2.5) is NO longer > xMinBinSize(2.5).
  // splitX=false, splitY=true. Splits horizontally into 2 each (Total 4
  // sections of 2.5 x 2.5). 1 section discarded by threshold cut. 1 section
  // discarded by crossing cut. 2 sections survive. Size is (2.5 x 2.5) <=
  // minBins. -> Accept
  BOOST_CHECK_EQUAL(sStat[3].nSections, 2);
  BOOST_CHECK_EQUAL(sStat[3].area, 12.5);
  BOOST_CHECK_EQUAL(sStat[3].nLines, 6);
  BOOST_CHECK_EQUAL(sStat[3].discardedByThresholdCut, 1);
  BOOST_CHECK_EQUAL(sStat[3].discardedByCrossingCut, 1);
}

BOOST_AUTO_TEST_CASE(test_drill_and_expand_logic) {
  // DESCRIPTION:
  // This test verifies how new ExploreParametersSpace algorithm works with
  // DrillAndExpand option in the decision functor.
  // Bottom-right section is the found solution thanks to expanding sections

  std::vector<LineParameters> measurements = {
      {0.0f, -3.5f},  // y = -3.5
      {2.0f, -4.0f},  // y = 2x - 4
      {0.0f, 0.5f}    //  y = 0.5
  };

  TestExplorationOptions<LineParameters> opt;
  opt.xMinBinSize = 6.0f;
  opt.yMinBinSize = 6.0f;
  opt.noiseThreshold = 3;
  opt.threshold = 2;
  opt.min_division_level = 0;
  opt.expandX = 1.5;
  opt.expandY = 1.5;
  opt.lineFunctor = [](const LineParameters &p, float arg) {
    return p.slope * arg + p.intercept;
  };

  HoughAccumulatorSection section(8.0f, 8.0f, -4.0f, -4.0f, 0);
  section.updateDecision(HoughAccumulatorSection::Decision::DrillAndExpand);
  section.indices() = {0, 1, 2};

  std::map<int, Stats> sStat;

  opt.decisionFunctor = [&sStat, &opt](const HoughAccumulatorSection &sec,
                                       const std::vector<LineParameters> &mes) {
    using enum HoughAccumulatorSection::Decision;
    if (sec.count() < opt.threshold) {
      sStat[sec.divisionLevel()].discardedByThresholdCut += 1;
      return Discard;
    }
    if (sec.count() < 3 * opt.threshold &&
        !passIntersectionsCheck(sec, mes, opt.lineFunctor,
                                opt.threshold * (opt.threshold - 1))) {
      sStat[sec.divisionLevel()].discardedByCrossingCut += 1;
      return Discard;
    }

    sStat[sec.divisionLevel()].area += sec.xSize() * sec.ySize();
    sStat[sec.divisionLevel()].nSections += 1;
    sStat[sec.divisionLevel()].nLines += sec.count();

    if (sec.divisionLevel() <= opt.min_division_level) {
      return DrillAndExpand;
    }
    if (sec.count() <= opt.noiseThreshold && sec.xSize() <= opt.xMinBinSize &&
        sec.ySize() <= opt.yMinBinSize) {
      return Accept;
    }
    return DrillAndExpand;
  };

  std::vector<HoughAccumulatorSection> sectionsStack;
  sectionsStack.push_back(std::move(section));
  std::vector<HoughAccumulatorSection> results;

  exploreHoughParametersSpace(sectionsStack, measurements, opt, results);

  BOOST_CHECK(sectionsStack.empty());

  BOOST_REQUIRE_EQUAL(results.size(), 1);
  BOOST_CHECK_EQUAL(results[0].xSize(), 6.0);
  BOOST_CHECK_EQUAL(results[0].ySize(), 6.0);

  // ==========================================
  // TREE EXECUTION TRACE
  // ==========================================

  // LEVEL 0:
  // this section is not tested with decisionFunctor and so the statistics is
  // not collected

  // LEVEL 1: Splits into 4 sections, but each expands by 1.5 in both directions
  // (4.0 -> 6.0). Section(6.0, 6.0, -5.0, 1.0) [tL] - Discarded by threshold
  // cut (only y=0.5 reaches it). Section(6.0, 6.0, -1.0, 1.0) [tR] - Discarded
  // by crossing cut (has 2 lines, but no crossings inside). Section(6.0, 6.0,
  // -5.0, -5.0) [bL] - Discarded by crossing cut (has enough lines but misses
  // intersection). Section(6.0, 6.0, -1.0, -5.0) [bR] - Intersection captured
  // inside thanks to expansion -> Accept
  BOOST_CHECK_EQUAL(sStat[1].discardedByThresholdCut, 1);
  BOOST_CHECK_EQUAL(sStat[1].discardedByCrossingCut, 2);
  BOOST_CHECK_EQUAL(sStat[1].nSections, 1);
  BOOST_CHECK_EQUAL(sStat[1].area, 36.0);
  BOOST_CHECK_EQUAL(sStat[1].nLines, 3);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
