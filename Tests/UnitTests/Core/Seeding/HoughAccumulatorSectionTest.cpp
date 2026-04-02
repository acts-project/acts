// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>
#include "Acts/Seeding2/HoughAccumulatorSection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <vector>
#include <functional>
#include <map>

using namespace Acts;

namespace ActsTests {
    
auto logger = getDefaultLogger("UnitTests", Logging::VERBOSE);

struct Point2D {
    float x;
    float y;
};

// Structure for statistics to enable BOOST_CHECK on results
struct Stats {
    double area{};
    int nSections{};
    int nLines{};
    int discardedByThresholdCut{};
    int discardedByCrossingCut{};
};

template <typename measurement_t = Point2D>
struct TestExplorationOptions {
    float xMinBinSize = 1.0f;
    float yMinBinSize = 1.0f;
    float expandX = 1.1f;
    float expandY = 1.1f;
    unsigned threshold = 4;  // number of lines passing section for it to be still considered
    unsigned noiseThreshold = 12;  // number of lines passing section at the final split to consider it noise
    unsigned min_division_level = 2;

    enum class Decision { Discard, Accept, Drill, DrillAndExpand };

    using DecisionFunctor = std::function<Decision(
        const HoughAccumulatorSection &section,
        const std::vector<measurement_t> &measurements)>;
        
    DecisionFunctor decisionFunctor;
};

template <typename SectionType, typename MeasurementType, typename Functor>
bool passIntersectionsCheck(
    const SectionType &section,
    const std::vector<MeasurementType> &measurements,
    Functor lineFunctor, 
    const unsigned threshold) {
    
    const unsigned count = section.count();
    const float xLeft = section.xBegin();
    const float xRight = xLeft + section.xSize();
    const auto& indices = section.indices();

    // Small Buffer Optimization
    constexpr unsigned kMaxStackLines = 64; 
    
    if (count <= kMaxStackLines) {
        float yLeft[kMaxStackLines];
        float yRight[kMaxStackLines];        
        
        for (unsigned i = 0; i < count; ++i) {
            const auto &m = measurements[indices[i]];
            yLeft[i] = lineFunctor(m, xLeft);
            yRight[i] = lineFunctor(m, xRight);
        }
        
        unsigned inside = 0;
        for (unsigned i = 0; i < count; ++i) {
            for (unsigned j = i + 1; j < count; ++j) {
                if ((yLeft[i] - yLeft[j]) * (yRight[i] - yRight[j]) < 0.0f) {
                    inside++;
                    if (inside >= threshold) return true; // Early exit
                }
            }
        }
        return false;
    } 
    else {
        std::vector<float> yLeft(count);
        std::vector<float> yRight(count);
        for (unsigned i = 0; i < count; ++i) {
            const auto &m = measurements[indices[i]];
            yLeft[i] = lineFunctor(m, xLeft);
            yRight[i] = lineFunctor(m, xRight);
        }
        unsigned inside = 0;
        for (unsigned i = 0; i < count; ++i) {
            for (unsigned j = i + 1; j < count; ++j) {
                if ((yLeft[i] - yLeft[j]) * (yRight[i] - yRight[j]) < 0.0f) {
                    inside++;
                    if (inside >= threshold) return true; // Early exit
                }
            }
        }
        return inside >= threshold;
    }
}

// SUITE 1: HoughAccumulatorSection Test
BOOST_AUTO_TEST_SUITE(HoughAccumulatorSectionSuite)

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
    BOOST_CHECK_EQUAL(0.5f*s.ySize(), t.ySize());
    BOOST_CHECK_EQUAL(0.5f*s.ySize(), b.ySize());
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
    BOOST_CHECK_EQUAL(0.5f*s.xSize(), l.xSize());
    BOOST_CHECK_EQUAL(0.5f*s.xSize(), r.xSize());
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
    BOOST_CHECK( ! s.isLineInside( [](float x){ return 0.5f*x+5.f; }) );
    BOOST_CHECK( s.isLineInside( [](float x){ return 0.5f*x+3.f; }) );
    BOOST_CHECK( s.isLineInside( [](float x){ return 0.1f*x; }) );
    BOOST_CHECK( s.isLineInside( [](float x){ return 1.0f*x-5.f; }) );
    BOOST_CHECK( ! s.isLineInside( [](float x){ return 0.5f*x-5.f; }) );
}

BOOST_AUTO_TEST_CASE(is_crossing_inside) {
    HoughAccumulatorSection s(10., 6., -5., -6.);
    std::function<float(float)> l1 = [](float x){ return 1.0f*x+3.f; };
    std::function<float(float)> l2 = [](float x){ return 0.5f*x+1.f; };
    std::function<float(float)> l3 = [](float x){ return 0.5f*x-1.f; };
    std::function<float(float)> l4 = [](float x){ return 1.0f*x-3.f; };

    BOOST_CHECK(s.isCrossingInside(l1, l2));
    BOOST_CHECK(! s.isCrossingInside(l2, l3));
    BOOST_CHECK(! s.isCrossingInside(l1, l3));
    BOOST_CHECK(! s.isCrossingInside(l4, l3));
}

BOOST_AUTO_TEST_SUITE_END()

// Core Engine
template <typename M, typename Options, typename Functor>
void exploreParametersSpace(std::vector<HoughAccumulatorSection> &sectionsStack,
                            const std::vector<M> &measurements,
                            const Options &opt,
                            Functor lineFunctor,
                            std::vector<HoughAccumulatorSection> &results) {
    
    using Decision = typename Options::Decision;
    
    while (!sectionsStack.empty()) {
        //ACTS_VERBOSE("Stack size " << sectionsStack.size());
        
        HoughAccumulatorSection thisSection = std::move(sectionsStack.back());
        sectionsStack.pop_back();

        Decision whatNext = opt.decisionFunctor(thisSection, measurements);

        if (whatNext == Decision::Discard) {
            continue;
        } else if (whatNext == Decision::Accept) {
            results.push_back(std::move(thisSection));
        } else {
            const float xMin = thisSection.xBegin();
            const float xMax = xMin + thisSection.xSize();
            const float yMin = thisSection.yBegin();
            const float yMax = yMin + thisSection.ySize();
            const float yMid = yMin + thisSection.ySize() * 0.5;

            bool splitX = thisSection.xSize() > opt.xMinBinSize;
            bool splitY = thisSection.ySize() > opt.yMinBinSize;

            if (splitX && splitY) {
                // Split into 4 sections
                HoughAccumulatorSection bL = thisSection.bottomLeft();
                HoughAccumulatorSection tL = thisSection.topLeft();
                HoughAccumulatorSection bR = thisSection.bottomRight();
                HoughAccumulatorSection tR = thisSection.topRight();

                if (whatNext == Decision::DrillAndExpand) {
                    bL.expand(opt.expandX, opt.expandY);
                    tL.expand(opt.expandX, opt.expandY);
                    bR.expand(opt.expandX, opt.expandY);
                    tR.expand(opt.expandX, opt.expandY);
                }
                
                // checking if lines are crossing 4 new sections 
                for (unsigned idx : thisSection.indices()) {
                    const auto &m = measurements[idx];
                    float yL = lineFunctor(m, xMin);
                    float yR = lineFunctor(m, xMax);
                    float yM = (yL + yR) * 0.5;

                    float minL = std::min(yL, yM);
                    float maxL = std::max(yL, yM);
                    float minR = std::min(yM, yR);
                    float maxR = std::max(yM, yR);

                    if (maxL > yMin && minL < yMid) bL.indices().push_back(idx);
                    if (maxL > yMid && minL < yMax) tL.indices().push_back(idx);
                    if (maxR > yMin && minR < yMid) bR.indices().push_back(idx);
                    if (maxR > yMid && minR < yMax) tR.indices().push_back(idx);
                }
                sectionsStack.push_back(std::move(bL));
                sectionsStack.push_back(std::move(tL));
                sectionsStack.push_back(std::move(bR));
                sectionsStack.push_back(std::move(tR));

            } else if (!splitX && splitY) {
                // Split into 2 horizontal sections
                HoughAccumulatorSection b = thisSection.bottom();
                HoughAccumulatorSection t = thisSection.top();

                if (whatNext == Decision::DrillAndExpand) {
                    b.expand(opt.expandX, opt.expandY);
                    t.expand(opt.expandX, opt.expandY);
                }

                for (unsigned idx : thisSection.indices()) {
                    const auto &m = measurements[idx];
                    float yL = lineFunctor(m, xMin);
                    float yR = lineFunctor(m, xMax);
                    float minAll = std::min(yL, yR);
                    float maxAll = std::max(yL, yR);

                    if (maxAll > yMin && minAll < yMid) b.indices().push_back(idx);
                    if (maxAll > yMid && minAll < yMax) t.indices().push_back(idx);
                }
                sectionsStack.push_back(std::move(b));
                sectionsStack.push_back(std::move(t));

            } else if (splitX && !splitY) {
                // Split into 2 vertical sections
                HoughAccumulatorSection l = thisSection.left();
                HoughAccumulatorSection r = thisSection.right();

                if (whatNext == Decision::DrillAndExpand) {
                    l.expand(opt.expandX, opt.expandY);
                    r.expand(opt.expandX, opt.expandY);
                }

                for (unsigned idx : thisSection.indices()) {
                    const auto &m = measurements[idx];
                    float yL = lineFunctor(m, xMin);
                    float yR = lineFunctor(m, xMax);
                    float yM = (yL + yR) * 0.5f;

                    float minL = std::min(yL, yM);
                    float maxL = std::max(yL, yM);
                    float minR = std::min(yM, yR);
                    float maxR = std::max(yM, yR);

                    if (maxL > yMin && minL < yMax) l.indices().push_back(idx);
                    if (maxR > yMin && minR < yMax) r.indices().push_back(idx);
                }
                sectionsStack.push_back(std::move(l));
                sectionsStack.push_back(std::move(r));
            }
        }
    }
}


// SUITE 2: exploreParametersSpace Test
BOOST_AUTO_TEST_SUITE(ExploreParametersSpaceSuite)

BOOST_AUTO_TEST_CASE(test_extra_split_x_min_div_0) {
    
    // Linear function definition: y = a*x + b
    auto lineFunctor = [](const Point2D& p, float arg) { return p.x * arg + p.y; };

    std::vector<Point2D> measurements = {
        { 0.05, 1.0},  // y = 0.05x + 1 
        { 0.08, 1.01},  // y = 0.08x + 1.01
        {-0.1, 3.0}   // y = -0.1x + 3
    };

    TestExplorationOptions<Point2D> opt;
    opt.xMinBinSize = 6.0f;   
    opt.yMinBinSize = 10.0f; 
    opt.noiseThreshold = 3;
    opt.threshold = 2;
    opt.min_division_level = 0; 
    
    HoughAccumulatorSection section(20.0f, 20.0f, 0.0f, 0.0f, 0);
    section.indices() = {0, 1, 2}; 

    std::map<int, Stats> sStat;

    opt.decisionFunctor = [&sStat, &opt, lineFunctor](
                                const HoughAccumulatorSection &sec,
                                const std::vector<Point2D> &mes) {
        
        using enum TestExplorationOptions<Point2D>::Decision;

        if (sec.count() < opt.threshold) {
            sStat[sec.divisionLevel()].discardedByThresholdCut += 1;
            return Discard;
        }
        if (sec.count() < 3 * opt.threshold &&
            !passIntersectionsCheck(sec, mes, lineFunctor,
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
        if (sec.count() <= opt.noiseThreshold &&
            sec.xSize() <= opt.xMinBinSize && 
            sec.ySize() <= opt.yMinBinSize) {
            return Accept;
        }

        return Drill;
    };

    std::vector<HoughAccumulatorSection> sectionsStack;
    sectionsStack.push_back(std::move(section));
    std::vector<HoughAccumulatorSection> results;

    exploreParametersSpace(sectionsStack, measurements, opt, lineFunctor, results);
    
    BOOST_CHECK(sectionsStack.empty()); 

    // result dimensions check
    BOOST_CHECK(results.size() == 1);
    BOOST_CHECK_EQUAL(results[0].xSize(), 5.0);
    BOOST_CHECK_EQUAL(results[0].ySize(), 10.0);

    // stats check:
    // level 0:
    BOOST_CHECK_EQUAL(sStat[0].nSections, 1);
    BOOST_CHECK_EQUAL(sStat[0].area, 400.0);
    
    // level 1: Split into 4 (10x10). The lines are in the bottom half. 
    // Top-left and Top-right get 0 lines -> early discarded!.
    BOOST_CHECK_EQUAL(sStat[1].discardedByThresholdCut, 2);
    BOOST_CHECK_EQUAL(sStat[1].discardedByCrossingCut, 1);
    BOOST_CHECK_EQUAL(sStat[1].nSections, 1); 
    BOOST_CHECK_EQUAL(sStat[1].area, 100.0);

    // level 2:
    BOOST_CHECK_EQUAL(sStat[2].nSections, 1);
    BOOST_CHECK_EQUAL(sStat[2].area, 50.0); // one 5x10 box
    BOOST_CHECK_EQUAL(sStat[2].discardedByThresholdCut, 0);
    BOOST_CHECK_EQUAL(sStat[2].discardedByCrossingCut, 1);
}


BOOST_AUTO_TEST_CASE(test_with_min_div_lvl_is_1) {
    
    // Linear function definition: y = a*x + b, where Point2D stores (a, b) parameters for test purposes.
    auto lineFunctor = [](const Point2D& p, float arg) { return p.x * arg + p.y; };

    // Basic config
    std::vector<Point2D> measurements = {
        {5.0, 3.0}, 
        {4.0, -4.0},
        {1.0, 2.0},
        {0.5, 7.0}, 
        {0.2, 1.0}
    };

    TestExplorationOptions<Point2D> opt;
    opt.xMinBinSize = 2.5;
    opt.yMinBinSize = 2.5;
    opt.noiseThreshold = 3;
    opt.threshold = 2;
    opt.min_division_level = 1;  
    
    // Create initial section and populate it with indices mapping to our 3 measurements
    HoughAccumulatorSection section(10.0, 20.0, -5.0, -10.0);
    section.indices() = {0, 1, 2, 3, 4}; 

    // Statistics map for verification in tests
    std::map<int, Stats> sStat;

    // Optimized decision functor
    opt.decisionFunctor = [&sStat, &opt, lineFunctor](
                                const HoughAccumulatorSection &sec,
                                const std::vector<Point2D> &mes) {
        
        using enum TestExplorationOptions<Point2D>::Decision;

        if (sec.count() < opt.threshold) {
            sStat[sec.divisionLevel()].discardedByThresholdCut += 1;
            return Discard;
        }
        if (sec.count() < 3 * opt.threshold &&
            !passIntersectionsCheck(sec, mes, lineFunctor,
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
        if (sec.count() <= opt.noiseThreshold &&
            sec.xSize() <= opt.xMinBinSize && 
            sec.ySize() <= opt.yMinBinSize) {
            return Accept;
        }

        return Drill;
    };

    std::vector<HoughAccumulatorSection> sectionsStack;
    sectionsStack.push_back(std::move(section)); // pushing root node
    
    std::vector<HoughAccumulatorSection> results;

    // Core engine
    exploreParametersSpace(sectionsStack, measurements, opt, lineFunctor, results);

    // Initial Assertions to verify the engine explored the space
    BOOST_CHECK(sStat.size() > 0); // Verifies the tree was drilled
    
    // Empty section stack after algorithm
    BOOST_CHECK(sectionsStack.empty()); 
    
    // Results vector tests
    BOOST_CHECK_EQUAL(results.size(), 2);
    BOOST_CHECK_EQUAL(results[0].xSize(), 2.5);
    BOOST_CHECK_EQUAL(results[0].count(), 3);
    BOOST_CHECK_EQUAL(results[1].ySize(), 2.5);
    BOOST_CHECK_EQUAL(results[1].count(), 3);
    
    // nSections tests
    BOOST_CHECK_EQUAL(sStat[0].nSections, 1);
    BOOST_CHECK_EQUAL(sStat[1].nSections, 2);
    BOOST_CHECK_EQUAL(sStat[2].nSections, 2);
    BOOST_CHECK_EQUAL(sStat[3].nSections, 2);

    // area tests
    BOOST_CHECK_EQUAL(sStat[0].area, 200.0);
    BOOST_CHECK_EQUAL(sStat[1].area, 100.0);
    BOOST_CHECK_EQUAL(sStat[2].area, 25.0);
    BOOST_CHECK_EQUAL(sStat[3].area, 12.5);

    // nLines tests
    BOOST_CHECK_EQUAL(sStat[0].nLines, 5);
    BOOST_CHECK_EQUAL(sStat[1].nLines, 9);
    BOOST_CHECK_EQUAL(sStat[2].nLines, 7);
    BOOST_CHECK_EQUAL(sStat[3].nLines, 6);

    // discardedByThresholdCut & discardedByCrossingCut tests
    BOOST_CHECK_EQUAL(sStat[0].discardedByThresholdCut, 0);
    BOOST_CHECK_EQUAL(sStat[1].discardedByThresholdCut, 1);
    BOOST_CHECK_EQUAL(sStat[2].discardedByThresholdCut, 2);
    BOOST_CHECK_EQUAL(sStat[3].discardedByThresholdCut, 1);
    BOOST_CHECK_EQUAL(sStat[0].discardedByCrossingCut, 0);
    BOOST_CHECK_EQUAL(sStat[1].discardedByCrossingCut, 1);
    BOOST_CHECK_EQUAL(sStat[2].discardedByCrossingCut, 4);
    BOOST_CHECK_EQUAL(sStat[3].discardedByCrossingCut, 1);
}

BOOST_AUTO_TEST_SUITE_END()

} // namespace ActsTests