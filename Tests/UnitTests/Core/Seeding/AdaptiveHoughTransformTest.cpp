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
using namespace Acts;

namespace ActsTests {
    
auto logger = getDefaultLogger("UnitTests", Logging::VERBOSE);

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
    BOOST_CHECK_EQUAL(0.5*s.ySize(), t.ySize());
    BOOST_CHECK_EQUAL(0.5*s.ySize(), b.ySize());
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

    BOOST_CHECK_EQUAL(0.5*s.xSize(), l.xSize());
    BOOST_CHECK_EQUAL(0.5*s.xSize(), r.xSize());

    BOOST_CHECK_EQUAL(l.xBegin(), -5.0);
    BOOST_CHECK_EQUAL(r.xBegin(), 0.0);
}

BOOST_AUTO_TEST_CASE(split_4) {
    // it is asymmetric and shifted
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
    // rather asymmetric section
    HoughAccumulatorSection s(10., 2., -5., -1.0);

    // line above
    BOOST_CHECK( ! s.isLineInside( [](float x){ return 0.5*x+5.; }) );

    // line above on the right
    BOOST_CHECK( s.isLineInside( [](float x){ return 0.5*x+3.; }) );

    // line that is totally inside
    BOOST_CHECK( s.isLineInside( [](float x){ return 0.1*x; }) );

    // line below on the right
    BOOST_CHECK( s.isLineInside( [](float x){ return 1.0*x-5.; }) );

    // line below
    BOOST_CHECK( ! s.isLineInside( [](float x){ return 0.5*x-5.; }) );
}

BOOST_AUTO_TEST_CASE(is_crossing_inside) {
    // rather asymmetric section again
    HoughAccumulatorSection s(10., 6., -5., -6.);
    // test lines (best to draw it)
    std::function<float(float)> l1 = [](float x){ return 1.0*x+3.; };
    std::function<float(float)> l2 = [](float x){ return 0.5*x+1.; };
    std::function<float(float)> l3 = [](float x){ return 0.5*x-1.; };
    std::function<float(float)> l4 = [](float x){ return 1.0*x-3.; };

    BOOST_CHECK(s.isCrossingInside(l1, l2));
    BOOST_CHECK(! s.isCrossingInside(l2, l3));
    BOOST_CHECK(! s.isCrossingInside(l1, l3));
    BOOST_CHECK(! s.isCrossingInside(l4, l3));
}

BOOST_AUTO_TEST_SUITE_END()

}