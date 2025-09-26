// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "Acts/Utilities/detail/TransformComparator.hpp"

BOOST_AUTO_TEST_SUITE(TransformComparatorTests)

BOOST_AUTO_TEST_CASE(Vector3Sorting) {
  ///
  Acts::detail::TransformComparator sorter{};
  Acts::detail::TransformComparator sorterWithTol{0.5, 0.5};

  for (unsigned dimA = 0; dimA < 3; ++dimA) {
    const Acts::Vector3 aPlus = Acts::Vector3::Unit(dimA);
    const Acts::Vector3 aMinus = -aPlus;
    BOOST_CHECK_EQUAL(sorter.compare<3>(aPlus, aPlus), 0);
    BOOST_CHECK_EQUAL(sorter.compare<3>(aPlus, aMinus), 1);
    BOOST_CHECK_EQUAL(sorter.compare<3>(aMinus, aPlus), -1);
    BOOST_CHECK_NE(sorter(aPlus, aMinus), sorter(aMinus, aPlus));

    const Acts::Vector3 halfPlus =
        (0.5 + std::numeric_limits<double>::epsilon()) * aPlus;
    const Acts::Vector3 halfMinus =
        (0.5 - std::numeric_limits<double>::epsilon()) * aPlus;

    BOOST_CHECK_EQUAL(sorterWithTol.compare<3>(aPlus, halfPlus), 0);
    BOOST_CHECK_NE(sorterWithTol.compare<3>(aPlus, halfMinus), 0);
    for (unsigned dimB = 0; dimB < 3; ++dimB) {
      const Acts::Vector3 bPlus = Acts::Vector3::Unit(dimB);
      const int expectSign = (dimA == dimB) ? 0 : dimA > dimB ? -1 : 1;
      BOOST_CHECK_EQUAL(sorter.compare<3>(aPlus, bPlus), expectSign);
    }
  }
}

BOOST_AUTO_TEST_CASE(Vector2Sorting) {
  Acts::detail::TransformComparator sorter{};
  Acts::detail::TransformComparator sorterWithTol{0.5, 0.5};

  for (unsigned dimA = 0; dimA < 2; ++dimA) {
    const Acts::Vector2 aPlus = Acts::Vector2::Unit(dimA);
    const Acts::Vector2 aMinus = -aPlus;
    BOOST_CHECK_EQUAL(sorter.compare<2>(aPlus, aPlus), 0);
    BOOST_CHECK_EQUAL(sorter.compare<2>(aPlus, aMinus), 1);
    BOOST_CHECK_EQUAL(sorter.compare<2>(aMinus, aPlus), -1);
    BOOST_CHECK_NE(sorter(aPlus, aMinus), sorter(aMinus, aPlus));

    const Acts::Vector2 halfPlus =
        (0.5 + std::numeric_limits<double>::epsilon()) * aPlus;
    const Acts::Vector2 halfMinus =
        (0.5 - std::numeric_limits<double>::epsilon()) * aPlus;

    BOOST_CHECK_EQUAL(sorterWithTol.compare<2>(aPlus, halfPlus), 0);
    BOOST_CHECK_NE(sorterWithTol.compare<2>(aPlus, halfMinus), 0);
    for (unsigned dimB = 0; dimB < 2; ++dimB) {
      const Acts::Vector2 bPlus = Acts::Vector2::Unit(dimB);
      const int expectSign = (dimA == dimB) ? 0 : dimA > dimB ? -1 : 1;
      BOOST_CHECK_EQUAL(sorter.compare<2>(aPlus, bPlus), expectSign);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
