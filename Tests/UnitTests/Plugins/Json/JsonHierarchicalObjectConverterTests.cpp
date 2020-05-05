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

#include "Acts/Plugins/Json/JsonHierarchicalObjectConverter.hpp"

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

Thing jsonToThing(json map) {
  Thing thing;
  std::stringstream sID;
  sID << map["id"];
  int v, b, l, a, s;
  std::string tstr;
  sID >> tstr >> v >> tstr >> b >> tstr >> l >> tstr >> a >> tstr >> s >> tstr;
  GeometryID geoID = GeometryID();
  geoID.setVolume(v).setBoundary(b).setLayer(l).setApproach(a).setSensitive(s);

  thing.id = geoID;
  thing.value = map["value"];
  return thing;
}

Thing initThing(GeometryID id) {
  Thing thing;
  thing.id = id;
  thing.value = 2.0;
  return thing;
}

using Container = Acts::HierarchicalGeometryContainer<Thing>;

// helper function to create geometry ids
GeometryID makeId(int volume, int layer = 0, int sensitive = 0) {
  return GeometryID().setVolume(volume).setLayer(layer).setSensitive(sensitive);
}

// helper function to create geometry ids
GeometryID makeId(int volume, int boundary, int layer, int approach,
                  int sensitive) {
  return GeometryID()
      .setVolume(volume)
      .setBoundary(boundary)
      .setLayer(layer)
      .setApproach(approach)
      .setSensitive(sensitive);
}

}  // namespace

BOOST_AUTO_TEST_CASE(Convert_HierarchicalObject) {
  Container c({
      {makeId(1), 1.0},
      {makeId(1, 2), 1.2},
      {makeId(1, 3), 1.5},
      {makeId(1, 3, 1), 2.0},
      {makeId(1, 3, 2), 2.2},
      {makeId(2), 3.0},
      {makeId(2, 1, 0, 0, 0), 10.0},
      {makeId(2, 1), 4.0},
      {makeId(2, 1, 2), 5.0},
      {makeId(2, 2), 15.0},
      {makeId(2, 0, 2, 1, 0), 13.0},
      {makeId(2, 0, 2, 3, 0), 14.0},
      {makeId(2, 0, 2, 0, 2), 15.0},
  });
  Acts::JsonHierarchicalObjectConverter<Thing>::Config cfg;
  Acts::JsonHierarchicalObjectConverter<Thing> converter(cfg);
  json map = converter.hierarchicalObjectToJson(c, thingToJson);

  BOOST_CHECK_EQUAL(jsonToThing(map[cfg.volkey]["1"][cfg.datakey]).value, 1.0);
  BOOST_CHECK_EQUAL(
      jsonToThing(map[cfg.volkey]["1"][cfg.laykey]["2"][cfg.repkey]).value,
      1.2);
  BOOST_CHECK_EQUAL(
      jsonToThing(map[cfg.volkey]["1"][cfg.laykey]["3"][cfg.repkey]).value,
      1.5);
  BOOST_CHECK_EQUAL(
      jsonToThing(map[cfg.volkey]["1"][cfg.laykey]["3"][cfg.senkey]["1"]).value,
      2.0);
  BOOST_CHECK_EQUAL(
      jsonToThing(map[cfg.volkey]["1"][cfg.laykey]["3"][cfg.senkey]["2"]).value,
      2.2);

  BOOST_CHECK_EQUAL(jsonToThing(map[cfg.volkey]["2"][cfg.datakey]).value, 3.0);
  BOOST_CHECK_EQUAL(
      jsonToThing(map[cfg.volkey]["2"][cfg.laykey]["1"][cfg.repkey]).value,
      4.0);
  BOOST_CHECK_EQUAL(
      jsonToThing(map[cfg.volkey]["2"][cfg.laykey]["1"][cfg.senkey]["2"]).value,
      5.0);
  BOOST_CHECK_EQUAL(jsonToThing(map[cfg.volkey]["2"][cfg.boukey]["1"]).value,
                    10.0);
  BOOST_CHECK_EQUAL(
      jsonToThing(map[cfg.volkey]["2"][cfg.laykey]["2"][cfg.appkey]["1"]).value,
      13.0);
  BOOST_CHECK_EQUAL(
      jsonToThing(map[cfg.volkey]["2"][cfg.laykey]["2"][cfg.appkey]["3"]).value,
      14.0);
  BOOST_CHECK_EQUAL(
      jsonToThing(map[cfg.volkey]["2"][cfg.laykey]["2"][cfg.senkey]["2"]).value,
      15.0);

  Container c2 = converter.jsonToHierarchicalContainer(map, jsonToThing);

  BOOST_CHECK_EQUAL(c2.find(makeId(1))->value, 1.0);
  BOOST_CHECK_EQUAL(c2.find(makeId(1, 2))->value, 1.2);
  BOOST_CHECK_EQUAL(c2.find(makeId(1, 3))->value, 1.5);
  BOOST_CHECK_EQUAL(c2.find(makeId(1, 3, 1))->value, 2.0);
  BOOST_CHECK_EQUAL(c2.find(makeId(1, 3, 2))->value, 2.2);

  BOOST_CHECK_EQUAL(c2.find(makeId(2))->value, 3.0);
  BOOST_CHECK_EQUAL(c2.find(makeId(2, 1))->value, 4.0);
  BOOST_CHECK_EQUAL(c2.find(makeId(2, 1, 2))->value, 5.0);
  BOOST_CHECK_EQUAL(c2.find(makeId(2, 1, 0, 0, 0))->value, 10.0);
  BOOST_CHECK_EQUAL(c2.find(makeId(2, 0, 2, 1, 0))->value, 13.0);
  BOOST_CHECK_EQUAL(c2.find(makeId(2, 0, 2, 3, 0))->value, 14.0);
}

BOOST_AUTO_TEST_CASE(Initialise_With_Geometry) {
  Acts::GeometryContext tgContext = Acts::GeometryContext();
  Acts::Test::CylindricalTrackingGeometry cGeometry(tgContext);
  auto tGeometry = cGeometry();

  Acts::JsonHierarchicalObjectConverter<Thing>::Config cfg;
  Acts::JsonHierarchicalObjectConverter<Thing> converter(cfg);

  json map =
      converter.trackingGeometryToJson(*tGeometry, thingToJson, initThing);

  BOOST_CHECK_EQUAL(map[cfg.volkey].size(), 3);
  BOOST_CHECK_EQUAL(map[cfg.volkey]["3"][cfg.laykey]["8"][cfg.senkey].size(),
                    1092);
  BOOST_CHECK_EQUAL(map[cfg.volkey]["2"][cfg.boukey].size(), 3);
}
