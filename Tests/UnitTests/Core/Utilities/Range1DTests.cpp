// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/RangeXD.hpp"

#include <limits>
#include <utility>
#include <vector>

namespace {
std::vector<int> v = {-100, -90, -80, -70, -60, -50, -40, -30, -20, 10, 0,
                      10,   20,  30,  40,  50,  60,  70,  80,  90,  100};
}

namespace Acts::Test {

BOOST_AUTO_TEST_SUITE(Utilities)

BOOST_AUTO_TEST_SUITE(Range1D)

BOOST_AUTO_TEST_CASE(infinite_range_int) {
  Acts::Range1D<int> r;

  BOOST_CHECK_EQUAL(r.min(), std::numeric_limits<int>::lowest());
  BOOST_CHECK_EQUAL(r.max(), std::numeric_limits<int>::max());
}

BOOST_AUTO_TEST_CASE(infinite_range_long double) {
  Acts::Range1D<long double> r;

  BOOST_CHECK_EQUAL(r.min(), std::numeric_limits<long double>::lowest());
  BOOST_CHECK_EQUAL(r.max(), std::numeric_limits<long double>::max());
}

BOOST_AUTO_TEST_CASE(constructor_range_int) {
  Acts::Range1D<int> r(-11, 2);

  BOOST_CHECK_EQUAL(r.min(), -11);
  BOOST_CHECK_EQUAL(r.max(), 2);
}

BOOST_AUTO_TEST_CASE(constructor_range_long double) {
  Acts::Range1D<long double> r(-11.0, 2.0);

  BOOST_CHECK_EQUAL(r.min(), -11.0);
  BOOST_CHECK_EQUAL(r.max(), 2.0);
}

BOOST_AUTO_TEST_CASE(constructor_copy_long double) {
  Acts::Range1D<long double> q(-11.0, 2.0);
  Acts::Range1D<long double> r(q);

  BOOST_CHECK_EQUAL(r.min(), -11.0);
  BOOST_CHECK_EQUAL(r.max(), 2.0);
}

BOOST_AUTO_TEST_CASE(constructor_tuple_long double) {
  Acts::Range1D<long double> r(std::pair<long double, long double>(-11.0, 2.0));

  BOOST_CHECK_EQUAL(r.min(), -11.0);
  BOOST_CHECK_EQUAL(r.max(), 2.0);
}

BOOST_AUTO_TEST_CASE(assign_long double) {
  Acts::Range1D<long double> q(-11.0, 2.0);
  Acts::Range1D<long double> r = q;

  BOOST_CHECK_EQUAL(r.min(), -11.0);
  BOOST_CHECK_EQUAL(r.max(), 2.0);
}

BOOST_AUTO_TEST_CASE(shrinkMin_long double) {
  Acts::Range1D<long double> r;

  r.shrinkMin(-11.0);

  BOOST_CHECK_EQUAL(r.min(), -11.0);
  BOOST_CHECK_EQUAL(r.max(), std::numeric_limits<long double>::max());
}

BOOST_AUTO_TEST_CASE(shrinkMax_long double) {
  Acts::Range1D<long double> r;

  r.shrinkMax(2.0);

  BOOST_CHECK_EQUAL(r.min(), std::numeric_limits<long double>::lowest());
  BOOST_CHECK_EQUAL(r.max(), 2.0);
}

BOOST_AUTO_TEST_CASE(shrink_both_long double) {
  Acts::Range1D<long double> r;

  r.shrinkMin(-11.0);
  r.shrinkMax(2.0);

  BOOST_CHECK_EQUAL(r.min(), -11.0);
  BOOST_CHECK_EQUAL(r.max(), 2.0);
}

BOOST_AUTO_TEST_CASE(shrink_long double) {
  Acts::Range1D<long double> r;

  r.shrink(-11.0, 2.0);

  BOOST_CHECK_EQUAL(r.min(), -11.0);
  BOOST_CHECK_EQUAL(r.max(), 2.0);
}

BOOST_AUTO_TEST_CASE(shrink_twice_long double) {
  Acts::Range1D<long double> r;

  r.shrink(-100.0, 20.0);
  r.shrink(-11.0, 2.0);

  BOOST_CHECK_EQUAL(r.min(), -11.0);
  BOOST_CHECK_EQUAL(r.max(), 2.0);
}

BOOST_AUTO_TEST_CASE(shrink_noop_long double) {
  Acts::Range1D<long double> r;

  r.shrink(-11.0, 2.0);
  r.shrink(-100.0, 20.0);

  BOOST_CHECK_EQUAL(r.min(), -11.0);
  BOOST_CHECK_EQUAL(r.max(), 2.0);
}

BOOST_AUTO_TEST_CASE(shrink_noop_min_long double) {
  Acts::Range1D<long double> r;

  r.shrink(-11.0, 2.0);
  r.shrink(-100.0, 1.0);

  BOOST_CHECK_EQUAL(r.min(), -11.0);
  BOOST_CHECK_EQUAL(r.max(), 1.0);
}

BOOST_AUTO_TEST_CASE(shrink_noop_max_long double) {
  Acts::Range1D<long double> r;

  r.shrink(-11.0, 2.0);
  r.shrink(-10.0, 20.0);

  BOOST_CHECK_EQUAL(r.min(), -10.0);
  BOOST_CHECK_EQUAL(r.max(), 2.0);
}

BOOST_AUTO_TEST_CASE(shrinkMin_noop_long double) {
  Acts::Range1D<long double> r;

  r.shrinkMin(-11.0);
  r.shrinkMin(-100.0);

  BOOST_CHECK_EQUAL(r.min(), -11.0);
  BOOST_CHECK_EQUAL(r.max(), std::numeric_limits<long double>::max());
}

BOOST_AUTO_TEST_CASE(shrinkMax_noop_long double) {
  Acts::Range1D<long double> r;

  r.shrinkMax(2.0);
  r.shrinkMax(4.0);

  BOOST_CHECK_EQUAL(r.min(), std::numeric_limits<long double>::lowest());
  BOOST_CHECK_EQUAL(r.max(), 2.0);
}

BOOST_AUTO_TEST_CASE(setMin_long double) {
  Acts::Range1D<long double> r;

  r.setMin(2.0);

  BOOST_CHECK_EQUAL(r.min(), 2.0);
  BOOST_CHECK_EQUAL(r.max(), std::numeric_limits<long double>::max());
}

BOOST_AUTO_TEST_CASE(setMin_twice_long double) {
  Acts::Range1D<long double> r;

  r.setMin(-2.0);
  r.setMin(-4.0);

  BOOST_CHECK_EQUAL(r.min(), -4.0);
  BOOST_CHECK_EQUAL(r.max(), std::numeric_limits<long double>::max());
}

BOOST_AUTO_TEST_CASE(setMax_long double) {
  Acts::Range1D<long double> r;

  r.setMax(2.0);

  BOOST_CHECK_EQUAL(r.min(), std::numeric_limits<long double>::lowest());
  BOOST_CHECK_EQUAL(r.max(), 2.0);
}

BOOST_AUTO_TEST_CASE(setMax_twice_long double) {
  Acts::Range1D<long double> r;

  r.setMax(2.0);
  r.setMax(4.0);

  BOOST_CHECK_EQUAL(r.min(), std::numeric_limits<long double>::lowest());
  BOOST_CHECK_EQUAL(r.max(), 4.0);
}

BOOST_AUTO_TEST_CASE(expandMin_long double) {
  Acts::Range1D<long double> r(0.0, 0.0);

  r.expandMin(-11.0);

  BOOST_CHECK_EQUAL(r.min(), -11.0);
  BOOST_CHECK_EQUAL(r.max(), 0.0);
}

BOOST_AUTO_TEST_CASE(expandMax_long double) {
  Acts::Range1D<long double> r(0.0, 0.0);

  r.expandMax(2.0);

  BOOST_CHECK_EQUAL(r.min(), 0.0);
  BOOST_CHECK_EQUAL(r.max(), 2.0);
}

BOOST_AUTO_TEST_CASE(expand_both_long double) {
  Acts::Range1D<long double> r(0.0, 0.0);

  r.expandMin(-11.0);
  r.expandMax(2.0);

  BOOST_CHECK_EQUAL(r.min(), -11.0);
  BOOST_CHECK_EQUAL(r.max(), 2.0);
}

BOOST_AUTO_TEST_CASE(expand_long double) {
  Acts::Range1D<long double> r(0.0, 0.0);

  r.expand(-11.0, 2.0);

  BOOST_CHECK_EQUAL(r.min(), -11.0);
  BOOST_CHECK_EQUAL(r.max(), 2.0);
}

BOOST_AUTO_TEST_CASE(expand_twice_long double) {
  Acts::Range1D<long double> r(0.0, 0.0);

  r.expand(-11.0, 2.0);
  r.expand(-100.0, 20.0);

  BOOST_CHECK_EQUAL(r.min(), -100.0);
  BOOST_CHECK_EQUAL(r.max(), 20.0);
}

BOOST_AUTO_TEST_CASE(expand_noop_long double) {
  Acts::Range1D<long double> r(0.0, 0.0);

  r.expand(-100.0, 20.0);
  r.expand(-11.0, 2.0);

  BOOST_CHECK_EQUAL(r.min(), -100.0);
  BOOST_CHECK_EQUAL(r.max(), 20.0);
}

BOOST_AUTO_TEST_CASE(expand_noop_min_long double) {
  Acts::Range1D<long double> r(0.0, 0.0);

  r.expand(-100.0, 1.0);
  r.expand(-11.0, 2.0);

  BOOST_CHECK_EQUAL(r.min(), -100.0);
  BOOST_CHECK_EQUAL(r.max(), 2.0);
}

BOOST_AUTO_TEST_CASE(expand_noop_max_long double) {
  Acts::Range1D<long double> r(0.0, 0.0);

  r.expand(-10.0, 20.0);
  r.expand(-11.0, 2.0);

  BOOST_CHECK_EQUAL(r.min(), -11.0);
  BOOST_CHECK_EQUAL(r.max(), 20.0);
}

BOOST_AUTO_TEST_CASE(expandMin_noop_long double) {
  Acts::Range1D<long double> r(0.0, 0.0);

  r.expandMin(-100.0);
  r.expandMin(-11.0);

  BOOST_CHECK_EQUAL(r.min(), -100.0);
  BOOST_CHECK_EQUAL(r.max(), 0.0);
}

BOOST_AUTO_TEST_CASE(expandMax_noop_long double) {
  Acts::Range1D<long double> r(0.0, 0.0);

  r.expandMax(4.0);
  r.expandMax(2.0);

  BOOST_CHECK_EQUAL(r.min(), 0.0);
  BOOST_CHECK_EQUAL(r.max(), 4.0);
}

BOOST_AUTO_TEST_CASE(size_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);

  BOOST_CHECK_EQUAL(r.size(), 35.0);
}

BOOST_AUTO_TEST_CASE(size_zero_long double) {
  Acts::Range1D<long double> r(-10.0, -10.0);

  BOOST_CHECK_EQUAL(r.size(), 0.0);
}

BOOST_AUTO_TEST_CASE(size_zero2_long double) {
  Acts::Range1D<long double> r(-10.0, -50.0);

  BOOST_CHECK_EQUAL(r.size(), 0.0);
}

BOOST_AUTO_TEST_CASE(degenerate_false_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);

  BOOST_CHECK(!r.degenerate());
}

BOOST_AUTO_TEST_CASE(degenerate_true_long double) {
  Acts::Range1D<long double> r(-10.0, -25.0);

  BOOST_CHECK(r.degenerate());
}

BOOST_AUTO_TEST_CASE(contains_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);

  BOOST_CHECK(!r.contains(-11.0));
  BOOST_CHECK(!r.contains(-30.0));
  BOOST_CHECK(r.contains(-5.0));
  BOOST_CHECK(r.contains(0.0));
  BOOST_CHECK(r.contains(10.0));
}

BOOST_AUTO_TEST_CASE(contains_degenerate_long double) {
  Acts::Range1D<long double> r(-10.0, -25.0);

  BOOST_CHECK(!r.contains(-11.0));
  BOOST_CHECK(!r.contains(-30.0));
  BOOST_CHECK(!r.contains(-5.0));
  BOOST_CHECK(!r.contains(0.0));
  BOOST_CHECK(!r.contains(10.0));
}

BOOST_AUTO_TEST_CASE(intersect_true1_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);
  Acts::Range1D<long double> q(5.0, 50.0);

  BOOST_CHECK((r && q));
}

BOOST_AUTO_TEST_CASE(intersect_true2_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);
  Acts::Range1D<long double> q(-100.0, 50.0);

  BOOST_CHECK((r && q));
}

BOOST_AUTO_TEST_CASE(intersect_true3_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);
  Acts::Range1D<long double> q(-5.0, 5.0);

  BOOST_CHECK((r && q));
}

BOOST_AUTO_TEST_CASE(intersect_false1_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);
  Acts::Range1D<long double> q(-50.0, -15.0);

  BOOST_CHECK(!(r && q));
}

BOOST_AUTO_TEST_CASE(intersect_false2_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);
  Acts::Range1D<long double> q(50.0, 55.0);

  BOOST_CHECK(!(r && q));
}

BOOST_AUTO_TEST_CASE(equals_true_int) {
  Acts::Range1D<int> r(-5, 5);
  Acts::Range1D<int> q(-5, 5);

  BOOST_CHECK((r == q));
}

BOOST_AUTO_TEST_CASE(equals_false_int) {
  Acts::Range1D<int> r(-5, 5);
  Acts::Range1D<int> q(-6, 4);

  BOOST_CHECK(!(r == q));
}

BOOST_AUTO_TEST_CASE(subset1_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);
  Acts::Range1D<long double> q(5.0, 50.0);

  BOOST_CHECK(!(r <= q));
  BOOST_CHECK(!(q <= r));
}

BOOST_AUTO_TEST_CASE(subset2_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);
  Acts::Range1D<long double> q(-100.0, 50.0);

  BOOST_CHECK((r <= q));
  BOOST_CHECK(!(q <= r));
}

BOOST_AUTO_TEST_CASE(subset3_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);
  Acts::Range1D<long double> q(-5.0, 5.0);

  BOOST_CHECK(!(r <= q));
  BOOST_CHECK((q <= r));
}

BOOST_AUTO_TEST_CASE(subset4_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);
  Acts::Range1D<long double> q(-50.0, -15.0);

  BOOST_CHECK(!(r <= q));
  BOOST_CHECK(!(q <= r));
}

BOOST_AUTO_TEST_CASE(subset5_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);
  Acts::Range1D<long double> q(50.0, 55.0);

  BOOST_CHECK(!(r <= q));
  BOOST_CHECK(!(q <= r));
}

BOOST_AUTO_TEST_CASE(superset1_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);
  Acts::Range1D<long double> q(5.0, 50.0);

  BOOST_CHECK(!(r >= q));
  BOOST_CHECK(!(q >= r));
}

BOOST_AUTO_TEST_CASE(superset2_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);
  Acts::Range1D<long double> q(-100.0, 50.0);

  BOOST_CHECK(!(r >= q));
  BOOST_CHECK((q >= r));
}

BOOST_AUTO_TEST_CASE(superset3_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);
  Acts::Range1D<long double> q(-5.0, 5.0);

  BOOST_CHECK((r >= q));
  BOOST_CHECK(!(q >= r));
}

BOOST_AUTO_TEST_CASE(superset4_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);
  Acts::Range1D<long double> q(-50.0, -15.0);

  BOOST_CHECK(!(r >= q));
  BOOST_CHECK(!(q >= r));
}

BOOST_AUTO_TEST_CASE(superset5_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);
  Acts::Range1D<long double> q(50.0, 55.0);

  BOOST_CHECK(!(r >= q));
  BOOST_CHECK(!(q >= r));
}

BOOST_AUTO_TEST_CASE(intersection1_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);
  Acts::Range1D<long double> q(5.0, 50.0);
  Acts::Range1D<long double> i = r & q;

  BOOST_CHECK_EQUAL(i.min(), 5.0);
  BOOST_CHECK_EQUAL(i.max(), 25.0);
}

BOOST_AUTO_TEST_CASE(intersection2_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);
  Acts::Range1D<long double> q(-100.0, 50.0);
  Acts::Range1D<long double> i = r & q;

  BOOST_CHECK_EQUAL(i.min(), -10.0);
  BOOST_CHECK_EQUAL(i.max(), 25.0);
}

BOOST_AUTO_TEST_CASE(intersection3_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);
  Acts::Range1D<long double> q(-5.0, 5.0);
  Acts::Range1D<long double> i = r & q;

  BOOST_CHECK_EQUAL(i.min(), -5.0);
  BOOST_CHECK_EQUAL(i.max(), 5.0);
}

BOOST_AUTO_TEST_CASE(intersection4_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);
  Acts::Range1D<long double> q(-50.0, -15.0);
  Acts::Range1D<long double> i = r & q;

  BOOST_CHECK(i.degenerate());
}

BOOST_AUTO_TEST_CASE(intersection5_long double) {
  Acts::Range1D<long double> r(-10.0, 25.0);
  Acts::Range1D<long double> q(50.0, 55.0);
  Acts::Range1D<long double> i = r & q;

  BOOST_CHECK(i.degenerate());
}

BOOST_AUTO_TEST_CASE(intersects_edge_int_positive) {
  Acts::Range1D<int> r(-10, 11);
  Acts::Range1D<int> q(10, 20);

  BOOST_CHECK((r && q));
}

BOOST_AUTO_TEST_CASE(intersects_edge_int_negative) {
  Acts::Range1D<int> r(-10, 10);
  Acts::Range1D<int> q(10, 20);

  BOOST_CHECK(!(r && q));
}

BOOST_AUTO_TEST_CASE(single_value_not_degenerate_int_positive) {
  Acts::Range1D<int> r(10, 10);

  BOOST_CHECK(r.degenerate());
}

BOOST_AUTO_TEST_CASE(single_value_not_degenerate_int_negative) {
  Acts::Range1D<int> r(10, 11);

  BOOST_CHECK(!r.degenerate());
}

BOOST_AUTO_TEST_CASE(intersects_implies_non_degenerate_intersection) {
  for (int i1 : v) {
    for (int i2 : v) {
      for (int j1 : v) {
        for (int j2 : v) {
          Acts::Range1D<int> r(i1, i2);
          Acts::Range1D<int> q(j1, j2);

          if (r.degenerate() || q.degenerate()) {
            continue;
          }

          if (r && q) {
            BOOST_CHECK(!((r & q).degenerate()));
          } else {
            BOOST_CHECK(((r & q).degenerate()));
          }
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(degeneracy_implies_size_zero_intersection) {
  for (int i1 : v) {
    for (int i2 : v) {
      Acts::Range1D<int> r(i1, i2);

      if (r.degenerate()) {
        BOOST_CHECK_EQUAL(r.size(), 0);
      } else {
        BOOST_CHECK_GE(r.size(), 0);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(range_inclusive_left) {
  Acts::Range1D<int> r(10, 20);

  BOOST_CHECK(!r.contains(9));
  BOOST_CHECK(r.contains(10));
  BOOST_CHECK(r.contains(11));
}

BOOST_AUTO_TEST_CASE(range_exclusive_right) {
  Acts::Range1D<int> r(10, 20);

  BOOST_CHECK(r.contains(19));
  BOOST_CHECK(!r.contains(20));
  BOOST_CHECK(!r.contains(21));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
