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
#include <filesystem>
#include <fstream>
#include <ranges>

auto testLogger =
    Acts::getDefaultLogger("test_walkthrough", Acts::Logging::VERBOSE);

// Convenience for test graph construction
template <typename T>
std::vector<T> &operator+=(std::vector<T> &first,
                           const std::vector<T> &second) {
  first.insert(first.end(), second.begin(), second.end());
  return first;
}

template <typename T>
std::span<T> toSpan(std::vector<T> &v) {
  return {v.begin(), v.end()};
}

using Vi = std::vector<int>;
using Vf = std::vector<float>;

BOOST_AUTO_TEST_CASE(test_single_chain) {
  Vi sources = {0, 1, 2, 3, 4, 5};
  Vi targets = {1, 2, 3, 4, 5, 6};
  Vf scores(6, 0.9f);

  Acts::WalkthroughAlgorithm walkthrough({}, testLogger->clone());

  auto candidates =
      walkthrough(toSpan(sources), toSpan(targets), toSpan(scores), 7);

  BOOST_REQUIRE_EQUAL(candidates.size(), 1);
  BOOST_CHECK(candidates.at(0) == (Vi{0, 1, 2, 3, 4, 5, 6}));
}

BOOST_AUTO_TEST_CASE(test_single_source_two_branches) {
  Vi sources = {0, 1, 2, 3, 4, 5};
  Vi targets = {1, 2, 3, 4, 5, 6};

  sources += Vi{2, 7, 8};
  targets += Vi{7, 8, 9};

  Vf scores(sources.size(), 0.9f);

  Acts::WalkthroughAlgorithm walkthrough({}, testLogger->clone());

  auto candidates =
      walkthrough(toSpan(sources), toSpan(targets), toSpan(scores), 10);

  std::ranges::sort(candidates, std::less{}, [](auto &c) { return c.size(); });
  BOOST_REQUIRE_EQUAL(candidates.size(), 2);
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

  Acts::WalkthroughAlgorithm walkthrough({}, testLogger->clone());

  auto candidates =
      walkthrough(toSpan(sources), toSpan(targets), toSpan(scores), 10);

  std::ranges::sort(candidates, std::less{}, [](auto &c) { return c.size(); });
  BOOST_REQUIRE_EQUAL(candidates.size(), 2);
  BOOST_CHECK(candidates.at(0) == (Vi{3, 4, 5, 6}));
  BOOST_CHECK(candidates.at(1) == (Vi{0, 1, 2, 7, 8, 9}));
}

BOOST_AUTO_TEST_CASE(test_crossing) {
  Vi sources = {0, 1, 2, 3, 4, 5};
  Vi targets = {1, 2, 3, 4, 5, 6};

  sources += Vi{7, 8, 9, 10, 3, 11, 12};
  targets += Vi{8, 9, 10, 3, 11, 12, 13};

  Vf scores(sources.size(), 0.9f);

  auto numNodes = std::max(*std::ranges::max_element(sources),
                           *std::ranges::max_element(targets)) +
                  1;

  Acts::WalkthroughAlgorithm walkthrough({}, testLogger->clone());

  auto candidates =
      walkthrough(toSpan(sources), toSpan(targets), toSpan(scores), numNodes);

  std::ranges::sort(candidates, std::less{}, [](auto &c) { return c.size(); });
  BOOST_REQUIRE_EQUAL(candidates.size(), 3);
  BOOST_CHECK(candidates.at(0) == (Vi{0, 1, 2}));
  BOOST_CHECK(candidates.at(1) == (Vi{11, 12, 13}));
  BOOST_CHECK(candidates.at(2) == (Vi{7, 8, 9, 10, 3, 4, 5, 6}));
}

boost::test_tools::assertion_result files_available(
    boost::unit_test::test_unit_id) {
  namespace fs = std::filesystem;
  return fs::exists("edge_sources.txt") && fs::exists("edge_targets.txt") &&
         fs::exists("edge_scores.txt");
}

BOOST_AUTO_TEST_CASE(test_from_files,
                     *boost::unit_test::precondition(files_available)) {
  auto read = [](auto &vector, const auto &filename) {
    std::ifstream file(filename);
    std::decay_t<decltype(vector.front())> v;
    while (file >> v) {
      vector.push_back(v);
    }
  };

  Vi sources, targets;
  Vf scores;
  read(sources, "edge_sources.txt");
  read(targets, "edge_targets.txt");
  read(scores, "edge_scores.txt");

  auto numNodes = std::max(*std::ranges::max_element(sources),
                           *std::ranges::max_element(targets)) +
                  1;

  Acts::WalkthroughAlgorithm::Config cfg;
  cfg.ccScoreCut = 0.01;
  Acts::WalkthroughAlgorithm walkthrough(cfg, testLogger->clone());

  auto candidates =
      walkthrough(toSpan(sources), toSpan(targets), toSpan(scores), numNodes);
}
