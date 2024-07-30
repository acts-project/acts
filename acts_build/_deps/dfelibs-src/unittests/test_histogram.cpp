// SPDX-License-Identifier: MIT

/// \file
/// \brief Unit tests for dfe::Histogram

#include <boost/test/unit_test.hpp>
#include <iostream>

#include "dfe/dfe_histogram.hpp"

// storage tests

using dfe::histogram_impl::NArray;

BOOST_AUTO_TEST_CASE(histogram_narray2f_init) {
  NArray<float, 2> s({10, 9}, 0.0f);

  BOOST_TEST(s.size().size() == 2);
  BOOST_TEST(s.size()[0] == 10);
  BOOST_TEST(s.size()[1] == 9);

  for (size_t i = 0; i < 10; ++i) {
    for (size_t j = 0; j < 9; ++j) {
      auto x = s[{i, j}];
      BOOST_TEST(x == 0.0f);
    }
  }
}

BOOST_AUTO_TEST_CASE(histogram_narray2f_at) {
  NArray<float, 2> s({10, 9}, 12.0f);

  BOOST_TEST(s.at({0, 0}) == 12.0f);
  BOOST_TEST(s.at({0, 8}) == 12.0f);
  BOOST_TEST(s.at({9, 0}) == 12.0f);
  BOOST_TEST(s.at({9, 8}) == 12.0f);
  BOOST_CHECK_THROW(s.at({0, 9}), std::out_of_range);
  BOOST_CHECK_THROW(s.at({10, 0}), std::out_of_range);
  BOOST_CHECK_THROW(s.at({10, 9}), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(histogram_narray3f_init) {
  NArray<float, 3> s({10, 9, 8}, 23.0f);

  BOOST_TEST(s.size().size() == 3);
  BOOST_TEST(s.size()[0] == 10);
  BOOST_TEST(s.size()[1] == 9);
  BOOST_TEST(s.size()[2] == 8);

  for (size_t i = 0; i < 10; ++i) {
    for (size_t j = 0; j < 9; ++j) {
      for (size_t k = 0; k < 8; ++k) {
        auto x = s[{i, j}];
        BOOST_TEST(x == 23.0f);
      }
    }
  }
}

// histogram tests

BOOST_AUTO_TEST_CASE(histogram_uniform1) {
  using H1 = dfe::Histogram<double, dfe::UniformAxis<double>>;

  H1 h({0.0, 1.0, 8});
  BOOST_CHECK_THROW(h.fill(-0.12), std::out_of_range);
  // upper boundary is alway exclusive
  BOOST_CHECK_THROW(h.fill(1.0), std::out_of_range);
  // lower boundary is inclusive
  BOOST_CHECK_NO_THROW(h.fill(0.0)); // bin 0 lower edge
  BOOST_CHECK_NO_THROW(h.fill(0.125)); // bin 1 lower edge
  BOOST_CHECK_NO_THROW(h.fill(0.2)); // bin 1 inside
  BOOST_CHECK_NO_THROW(h.fill(0.25)); // bin 2 lower edge
  BOOST_CHECK_NO_THROW(h.fill(0.99));
  BOOST_TEST(h.value({0}) == 1);
  BOOST_TEST(h.value({1}) == 2);
  BOOST_TEST(h.value({2}) == 1);
  BOOST_TEST(h.value({3}) == 0);
  BOOST_TEST(h.value({4}) == 0);
  BOOST_TEST(h.value({5}) == 0);
  BOOST_TEST(h.value({6}) == 0);
  BOOST_TEST(h.value({7}) == 1);
}

BOOST_AUTO_TEST_CASE(histogram_uniform2f) {
  using H2 =
    dfe::Histogram<double, dfe::UniformAxis<double>, dfe::UniformAxis<double>>;

  H2 h({0.0, 1.0, 4}, {-10.0, 10, 4});
  BOOST_CHECK_THROW(h.fill(-1.2, 0.5), std::out_of_range);
  BOOST_CHECK_THROW(h.fill(23.0, 0.5), std::out_of_range);
  BOOST_CHECK_THROW(h.fill(0.5, -16.2), std::out_of_range);
  BOOST_CHECK_THROW(h.fill(0.5, 24.1), std::out_of_range);
  BOOST_CHECK_THROW(h.fill(-1.2, -16.2), std::out_of_range);
  BOOST_CHECK_THROW(h.fill(-1.2, 24.1), std::out_of_range);
  BOOST_CHECK_THROW(h.fill(17.2, -16.2), std::out_of_range);
  BOOST_CHECK_THROW(h.fill(17.2, 24.1), std::out_of_range);
  BOOST_CHECK_NO_THROW(h.fill(0.2, 1, 0.5));
  BOOST_CHECK_NO_THROW(h.fill(0.3, 2));
  BOOST_CHECK_NO_THROW(h.fill(0.5, -5.4));
  BOOST_TEST(h.value({0, 0}) == 0);
  BOOST_TEST(h.value({0, 1}) == 0);
  BOOST_TEST(h.value({0, 2}) == 0.5);
  BOOST_TEST(h.value({0, 3}) == 0);
  BOOST_TEST(h.value({1, 0}) == 0);
  BOOST_TEST(h.value({1, 1}) == 0);
  BOOST_TEST(h.value({1, 2}) == 1);
  BOOST_TEST(h.value({1, 3}) == 0);
  BOOST_TEST(h.value({2, 0}) == 1);
  BOOST_TEST(h.value({2, 1}) == 0);
  BOOST_TEST(h.value({2, 2}) == 0);
  BOOST_TEST(h.value({2, 3}) == 0);
  BOOST_TEST(h.value({3, 0}) == 0);
  BOOST_TEST(h.value({3, 1}) == 0);
  BOOST_TEST(h.value({3, 2}) == 0);
  BOOST_TEST(h.value({3, 3}) == 0);
}

BOOST_AUTO_TEST_CASE(histogram_overflow1) {
  using H1 = dfe::Histogram<double, dfe::OverflowAxis<double>>;

  H1 h({0.0, 1.0, 8});
  BOOST_CHECK_NO_THROW(h.fill(std::numeric_limits<double>::lowest()));
  BOOST_CHECK_NO_THROW(h.fill(-0.12));
  BOOST_CHECK_NO_THROW(h.fill(std::numeric_limits<double>::max()));
  BOOST_CHECK_NO_THROW(h.fill(std::numeric_limits<double>::infinity()));
  // upper boundary is alway exclusive
  BOOST_CHECK_NO_THROW(h.fill(1.0));
  BOOST_CHECK_NO_THROW(h.fill(1.2));
  // lower boundary is inclusive
  BOOST_CHECK_NO_THROW(h.fill(0.0)); // bin 0 lower edge
  BOOST_CHECK_NO_THROW(h.fill(0.125)); // bin 1 lower edge
  BOOST_CHECK_NO_THROW(h.fill(0.2)); // bin 1 inside
  BOOST_CHECK_NO_THROW(h.fill(0.25)); // bin 2 lower edge
  BOOST_CHECK_NO_THROW(h.fill(0.99));
  BOOST_TEST(h.value({0}) == 2); // underflow
  BOOST_TEST(h.value({1}) == 1);
  BOOST_TEST(h.value({2}) == 2);
  BOOST_TEST(h.value({3}) == 1);
  BOOST_TEST(h.value({4}) == 0);
  BOOST_TEST(h.value({5}) == 0);
  BOOST_TEST(h.value({6}) == 0);
  BOOST_TEST(h.value({7}) == 0);
  BOOST_TEST(h.value({8}) == 1);
  BOOST_TEST(h.value({9}) == 4); // overflow
}

BOOST_AUTO_TEST_CASE(histogram_variable1d) {
  using H1 = dfe::Histogram<double, dfe::VariableAxis<double>>;

  // not enough edges
  BOOST_CHECK_THROW(H1({1.0}), std::invalid_argument);
  // edges not sorted
  BOOST_CHECK_THROW(H1({1.0, -2.0, 3.0}), std::invalid_argument);
  // edges are not unique
  BOOST_CHECK_THROW(H1({1.0, 2.0, 3.0, 3.0}), std::invalid_argument);

  H1 h({1.0, 10.0, 100.0, 1000.0});
  BOOST_TEST(h.size() == H1::Index{3});
  // invalid ranges
  BOOST_CHECK_THROW(h.fill(0.0), std::out_of_range);
  BOOST_CHECK_THROW(h.fill(1000.0), std::out_of_range); // on upper limit
  BOOST_CHECK_THROW(h.fill(99999999.0), std::out_of_range);
  BOOST_CHECK_THROW(
    h.fill(std::numeric_limits<double>::lowest()), std::out_of_range);
  BOOST_CHECK_THROW(
    h.fill(std::numeric_limits<double>::max()), std::out_of_range);
  BOOST_CHECK_THROW(
    h.fill(std::numeric_limits<double>::infinity()), std::out_of_range);
  BOOST_CHECK_THROW(
    h.fill(std::numeric_limits<double>::quiet_NaN()), std::out_of_range);
  // example data
  BOOST_CHECK_NO_THROW(h.fill(1));
  BOOST_CHECK_NO_THROW(h.fill(5));
  BOOST_CHECK_NO_THROW(h.fill(7));
  BOOST_CHECK_NO_THROW(h.fill(10));
  BOOST_CHECK_NO_THROW(h.fill(100));
  BOOST_CHECK_NO_THROW(h.fill(125));
  BOOST_CHECK_NO_THROW(h.fill(130));
  BOOST_TEST(h.value({0}) == 3);
  BOOST_TEST(h.value({1}) == 1);
  BOOST_TEST(h.value({2}) == 3);
}
