// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/EventData/TrackStateProxy.hpp"
#include "Acts/Utilities/PointerTraits.hpp"

#include <iostream>
#include <memory>

using namespace Acts;

namespace ActsTests {

template <typename T>
bool testPointer(const T /*ptr*/) {
  BOOST_TEST_MESSAGE("Object of " << typeid(T).name()
                                  << " does not pass the pointer concept ");
  return false;
}

template <PointerConcept T>
bool testPointer(const T /*ptr*/) {
  BOOST_TEST_MESSAGE("Object of " << typeid(T).name()
                                  << " passes the pointer concept ");
  return true;
}

BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

BOOST_AUTO_TEST_CASE(testConceptPass) {
  int* raw_ptr{nullptr};
  BOOST_CHECK(testPointer(raw_ptr));
  int raw_val{0};
  BOOST_CHECK(!testPointer(raw_val));

  BOOST_CHECK(testPointer(std::unique_ptr<int>{nullptr}));
  BOOST_CHECK(testPointer(std::unique_ptr<const int>{nullptr}));

  BOOST_CHECK(testPointer(std::shared_ptr<int>{nullptr}));
  BOOST_CHECK(testPointer(std::shared_ptr<const int>{nullptr}));

  BOOST_CHECK(testPointer(detail_lt::TransitiveConstPointer<int>{nullptr}));
  BOOST_CHECK(
      testPointer(detail_lt::TransitiveConstPointer<const int>{nullptr}));
  // Class with partial pointer-like behavior
  struct PartialPointerLike {
    int* ptr = nullptr;
    int* operator->() const { return ptr; }
    // Missing operator*() and operator bool()
  };
  BOOST_CHECK(!testPointer(PartialPointerLike{}));

  /** Ensure that the remove_pointer_t trait is doing what's supposed to do */
  static_assert(std::is_same_v<RemovePointer_t<std::unique_ptr<int>>, int>);
  static_assert(std::is_same_v<RemovePointer_t<std::shared_ptr<int>>, int>);
  static_assert(std::is_same_v<RemovePointer_t<int>, int>);
  static_assert(std::is_same_v<RemovePointer_t<int*>, int>);
}
BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
