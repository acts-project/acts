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

BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

BOOST_AUTO_TEST_CASE(testConceptPass) {
  static_assert(PointerConcept<int*>, "int* is a pointer concept");
  static_assert(!PointerConcept<int>, "int is not a pointer concept");

  static_assert(PointerConcept<std::unique_ptr<int>>,
                "std::unique_ptr<int> is a pointer concept");
  static_assert(PointerConcept<std::unique_ptr<const int>>,
                "std::unique_ptr<const int> is a pointer concept");

  static_assert(PointerConcept<std::shared_ptr<int>>,
                "std::shared_ptr<int> is a pointer concept");
  static_assert(PointerConcept<std::shared_ptr<const int>>,
                "std::shared_ptr<const int> is a pointer concept");

  static_assert(PointerConcept<detail_lt::TransitiveConstPointer<int>>,
                "detail_lt::TransitiveConstPointer<int> is a pointer concept");
  static_assert(
      PointerConcept<detail_lt::TransitiveConstPointer<const int>>,
      "detail_lt::TransitiveConstPointer<const int> is a pointer concept");

  // Class with partial pointer-like behavior
  struct PartialPointerLike {
    int* ptr = nullptr;
    int* operator->() const { return ptr; }
    // Missing operator*() and operator bool()
  };
  static_assert(!PointerConcept<PartialPointerLike>,
                "PartialPointerLike is not a pointer concept");

  /** Ensure that the remove_pointer_t trait is doing what's supposed to do */
  static_assert(std::is_same_v<RemovePointer_t<std::unique_ptr<int>>, int>);
  static_assert(std::is_same_v<RemovePointer_t<std::shared_ptr<int>>, int>);
  static_assert(std::is_same_v<RemovePointer_t<int>, int>);
  static_assert(std::is_same_v<RemovePointer_t<int*>, int>);
}
BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
