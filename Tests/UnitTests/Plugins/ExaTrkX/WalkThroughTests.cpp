// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Plugins/ExaTrkX/BoostWalkthrough.hpp"

#include <algorithm>

auto testLogger =
    Acts::getDefaultLogger("test_walkthrough", Acts::Logging::VERBOSE);

// Convenience for test graph construction
template <typename T>
std::vector<T> &operator+=(std::vector<T> &first,
                           const std::vector<T> &second) {
  first.insert(first.end(), second.begin(), second.end());
  return first;
}

using Vi = std::vector<int>;
using Vf = std::vector<float>;

BOOST_AUTO_TEST_CASE(test_single_chain) {
  Vi sources = {0, 1, 2, 3, 4, 5};
  Vi targets = {1, 2, 3, 4, 5, 6};
  Vf scores(6, 0.9f);

  Acts::Walkthrough walkthrough({}, testLogger->clone());

  auto candidates = walkthrough(sources, targets, scores, 7);

  BOOST_CHECK(candidates.size() == 1);
  BOOST_CHECK(candidates.at(0) == (Vi{0, 1, 2, 3, 4, 5, 6}));
}

BOOST_AUTO_TEST_CASE(test_single_source_two_branches) {
  Vi sources = {0, 1, 2, 3, 4, 5};
  Vi targets = {1, 2, 3, 4, 5, 6};

  sources += Vi{2, 7, 8};
  targets += Vi{7, 8, 9};

  Vf scores(sources.size(), 0.9f);

  Acts::Walkthrough walkthrough({}, testLogger->clone());

  auto candidates = walkthrough(sources, targets, scores, 10);

  std::ranges::sort(candidates, std::less{}, [](auto &c) { return c.size(); });
  BOOST_CHECK_EQUAL(candidates.size(), 2);
  BOOST_CHECK(candidates.at(0) == (Vi{7, 8, 9}));
  BOOST_CHECK(candidates.at(1) == (Vi{0, 1, 2, 3, 4, 5, 6}));
}

BOOST_AUTO_TEST_CASE(test_single_source_weak_link) {
  Vi sources = {0, 1, 2, 3, 4, 5};
  Vi targets = {1, 2, 3, 4, 5, 6};

  sources += Vi{2, 7, 8};
  targets += Vi{7, 8, 9};

  Vf scores(sources.size(), 0.9f);
  // add a weak link here, so the other branch should be taken
  scores[2] = 0.3;

  Acts::Walkthrough walkthrough({}, testLogger->clone());

  auto candidates = walkthrough(sources, targets, scores, 10);

  std::ranges::sort(candidates, std::less{}, [](auto &c) { return c.size(); });
  BOOST_CHECK_EQUAL(candidates.size(), 2);
  BOOST_CHECK(candidates.at(0) == (Vi{3, 4, 5, 6}));
  BOOST_CHECK(candidates.at(1) == (Vi{0, 1, 2, 7, 8, 9}));
}

BOOST_AUTO_TEST_CASE(test_crossing) {
  Vi sources = {0, 1, 2, 3, 4, 5};
  Vi targets = {1, 2, 3, 4, 5, 6};

  sources += Vi{7, 8, 9, 10, 3, 11, 12};
  targets += Vi{8, 9, 10, 3, 11, 12, 13};

  Vf scores(sources.size(), 0.9f);

  Acts::Walkthrough walkthrough({}, testLogger->clone());

  auto candidates = walkthrough(sources, targets, scores, 14);

  std::ranges::sort(candidates, std::less{}, [](auto &c) { return c.size(); });
  BOOST_CHECK_EQUAL(candidates.size(), 3);
  BOOST_CHECK(candidates.at(0) == (Vi{11, 12, 13}));
  BOOST_CHECK(candidates.at(1) == (Vi{7, 8, 9, 10}));
  BOOST_CHECK(candidates.at(2) == (Vi{0, 1, 2, 3, 4, 5, 6}));
}
