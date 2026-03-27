// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsPlugins/Json/JsonKindDispatcher.hpp"

#include <stdexcept>

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(JsonKindDispatcherDispatchesByKind) {
  Acts::JsonKindDispatcher<int, int> dispatcher{"kind", "test payload"};
  dispatcher.registerKind("A", [](const nlohmann::json& j, int increment) {
    return j.at("value").get<int>() + increment;
  });

  nlohmann::json encoded = {{"kind", "A"}, {"value", 41}};
  BOOST_CHECK_EQUAL(dispatcher(encoded, 1), 42);
  BOOST_CHECK(dispatcher.hasKind("A"));
  BOOST_CHECK_EQUAL(dispatcher.size(), 1u);
}

BOOST_AUTO_TEST_CASE(JsonKindDispatcherRejectsDuplicateKind) {
  Acts::JsonKindDispatcher<int> dispatcher{"kind", "test payload"};
  dispatcher.registerKind("A", [](const nlohmann::json&) { return 1; });

  BOOST_CHECK_THROW(
      dispatcher.registerKind("A", [](const nlohmann::json&) { return 2; }),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(JsonKindDispatcherRejectsMissingKind) {
  Acts::JsonKindDispatcher<int> dispatcher{"kind", "test payload"};
  dispatcher.registerKind("A", [](const nlohmann::json&) { return 1; });

  nlohmann::json encoded = {{"value", 1}};
  BOOST_CHECK_THROW(dispatcher(encoded), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(JsonKindDispatcherRejectsNonStringKind) {
  Acts::JsonKindDispatcher<int> dispatcher{"kind", "test payload"};
  dispatcher.registerKind("A", [](const nlohmann::json&) { return 1; });

  nlohmann::json encoded = {{"kind", 1}};
  BOOST_CHECK_THROW(dispatcher(encoded), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(JsonKindDispatcherRejectsUnknownKind) {
  Acts::JsonKindDispatcher<int> dispatcher{"kind", "test payload"};
  dispatcher.registerKind("A", [](const nlohmann::json&) { return 1; });

  nlohmann::json encoded = {{"kind", "B"}};
  BOOST_CHECK_THROW(dispatcher(encoded), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
