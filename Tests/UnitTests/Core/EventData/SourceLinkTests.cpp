// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <any>

using namespace Acts::UnitLiterals;

BOOST_AUTO_TEST_SUITE(EventDataSourceLink)

struct MySourceLink {
  Acts::GeometryIdentifier m_geometryId;

  Acts::GeometryIdentifier geometryId() const { return m_geometryId; }
};

BOOST_AUTO_TEST_CASE(Construct) {
  MySourceLink msl;
  msl.m_geometryId.setSensitive(42);
  {
    Acts::SourceLink sl{msl};
    BOOST_CHECK_EQUAL(sl.get<MySourceLink>().geometryId(), msl.geometryId());
    BOOST_CHECK_THROW(sl.get<int>(), std::bad_any_cast);
  }
}

BOOST_AUTO_TEST_SUITE_END()
