// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <any>
#include <sstream>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(EventDataSuite)

BOOST_AUTO_TEST_CASE(TestSourceLinkCoverage) {
  using detail::Test::TestSourceLink;

  Vector2 stddev(0.01, 0.1);
  SquareMatrix2 cov = stddev.cwiseProduct(stddev).asDiagonal();
  TestSourceLink l1(eBoundLoc0, 0.1, cov(0, 0), GeometryIdentifier(0x999), 0);
  TestSourceLink l2(l1);

  BOOST_CHECK(l1 == l2);     // testing the ==
  BOOST_CHECK(!(l1 != l2));  // testing the !=
  std::ostringstream str;
  str << l1;
}

struct MySourceLink {
  GeometryIdentifier m_geometryId;

  GeometryIdentifier geometryId() const { return m_geometryId; }
};

BOOST_AUTO_TEST_CASE(Construct) {
  MySourceLink msl;
  msl.m_geometryId = GeometryIdentifier().withSensitive(42);
  {
    SourceLink sl{msl};
    BOOST_CHECK_EQUAL(sl.get<MySourceLink>().geometryId(), msl.geometryId());
    BOOST_CHECK_THROW(sl.get<int>(), std::bad_any_cast);
  }
}

BOOST_AUTO_TEST_CASE(GetPtr) {
  {
    // correct type returns non-null pointer
    MySourceLink msl;
    msl.m_geometryId = GeometryIdentifier().withSensitive(7);
    SourceLink sl{msl};

    auto* p = sl.getPtr<MySourceLink>();
    BOOST_REQUIRE_NE(p, static_cast<MySourceLink*>(nullptr));
    BOOST_CHECK_EQUAL(p->geometryId(), msl.geometryId());

    // wrong type returns nullptr
    BOOST_CHECK_EQUAL(sl.getPtr<int>(), static_cast<int*>(nullptr));
    BOOST_CHECK_EQUAL(sl.getPtr<double>(), static_cast<double*>(nullptr));
  }

  {
    // const overload
    int value = 42;
    const SourceLink sl{value};

    const int* p = sl.getPtr<int>();
    BOOST_REQUIRE_NE(p, static_cast<const int*>(nullptr));
    BOOST_CHECK_EQUAL(*p, 42);

    BOOST_CHECK_EQUAL(sl.getPtr<float>(), static_cast<const float*>(nullptr));
  }

  {
    // mutation through pointer
    int value = 10;
    SourceLink sl{value};
    int* p = sl.getPtr<int>();
    BOOST_REQUIRE_NE(p, static_cast<int*>(nullptr));
    *p = 99;
    BOOST_CHECK_EQUAL(sl.get<int>(), 99);
  }
}

BOOST_AUTO_TEST_CASE(Reassign) {
  int value = 5;
  SourceLink sl{value};

  BOOST_CHECK_EQUAL(sl.get<int>(), value);
  BOOST_CHECK_THROW(sl.get<double>(), std::bad_any_cast);

  double otherValue = 42.42;

  // this changes the stored type
  sl = SourceLink{otherValue};
  BOOST_CHECK_EQUAL(sl.get<double>(), otherValue);
  BOOST_CHECK_THROW(sl.get<int>(), std::bad_any_cast);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
