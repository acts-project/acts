// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>
#include <iterator>
#include <stdexcept>

#include "Acts/Plugins/Json/JsonGeometricalHierarchyConverter.hpp"

#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"

namespace {

using Acts::GeometryID;
using json = nlohmann::json;

struct Thing {
  GeometryID id;
  double value = 1.0;

  Thing() {}
  Thing(GeometryID i, double v) : id(i), value(v) {}

  constexpr auto geometryId() const { return id; }
};

json thingToJson(Thing thing) {
  json jThing;
  std::ostringstream sID;
  sID << thing.id;

  jThing["id"] = sID.str();
  jThing["value"] = thing.value;
  return jThing;
}

Thing initThing(GeometryID id) {
  Thing thing;
  thing.id = id;
  thing.value = 2.0;
  return thing;
}

}  // namespace

BOOST_AUTO_TEST_CASE(Initialise_With_Geometry) {
  Acts::GeometryContext tgContext = Acts::GeometryContext();
  Acts::Test::CylindricalTrackingGeometry cGeometry(tgContext);
  auto tGeometry = cGeometry();

  Acts::JsonGeometricalHierarchyConverter<Thing>::Config cfg;
  Acts::JsonGeometricalHierarchyConverter<Thing> converter(cfg);

  json map =
      converter.trackingGeometryToJson(*tGeometry, thingToJson, initThing);

  BOOST_CHECK_EQUAL(map["volumes"].size(), 3);
  BOOST_CHECK_EQUAL(map["volumes"]["3"]["layers"]["8"]["sensitive"].size(),
                    1092);
  BOOST_CHECK_EQUAL(map["volumes"]["2"]["boundaries"].size(), 3);
}
