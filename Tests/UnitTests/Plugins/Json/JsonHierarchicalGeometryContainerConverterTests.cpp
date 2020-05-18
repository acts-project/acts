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

#include "Acts/Plugins/Json/JsonHierarchicalGeometryContainerConverter.hpp"

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

Thing jsonToThing(GeometryID id, json map) {
  Thing thing;
  thing.id = id;
  thing.value = map["value"];
  return thing;
}

double jsonTovalue(json map) {
  return map["value"];
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
  Acts::JsonHierarchicalGeometryContainerConverter<Thing> converter;
  converter.datakey = "Thing";
  json map = converter.hierarchicalObjectToJson(c, thingToJson);

  BOOST_CHECK_EQUAL(jsonTovalue(map["volumes"]["1"]["Thing"]), 1.0);
  BOOST_CHECK_EQUAL(
      jsonTovalue(map["volumes"]["1"]["layers"]["2"]["representing"]["Thing"]),
      1.2);
  BOOST_CHECK_EQUAL(
      jsonTovalue(map["volumes"]["1"]["layers"]["3"]["representing"]["Thing"]),
      1.5);
  BOOST_CHECK_EQUAL(
      jsonTovalue(
          map["volumes"]["1"]["layers"]["3"]["sensitive"]["1"]["Thing"]),
      2.0);
  BOOST_CHECK_EQUAL(
      jsonTovalue(
          map["volumes"]["1"]["layers"]["3"]["sensitive"]["2"]["Thing"]),
      2.2);

  BOOST_CHECK_EQUAL(jsonTovalue(map["volumes"]["2"]["Thing"]), 3.0);
  BOOST_CHECK_EQUAL(
      jsonTovalue(map["volumes"]["2"]["layers"]["1"]["representing"]["Thing"]),
      4.0);

  BOOST_CHECK_EQUAL(
      jsonTovalue(
          map["volumes"]["2"]["layers"]["1"]["sensitive"]["2"]["Thing"]),
      5.0);
  BOOST_CHECK_EQUAL(
      jsonTovalue(map["volumes"]["2"]["boundaries"]["1"]["Thing"]), 10.0);
  BOOST_CHECK_EQUAL(
      jsonTovalue(map["volumes"]["2"]["layers"]["2"]["approach"]["1"]["Thing"]),
      13.0);
  BOOST_CHECK_EQUAL(
      jsonTovalue(map["volumes"]["2"]["layers"]["2"]["approach"]["3"]["Thing"]),
      14.0);
  BOOST_CHECK_EQUAL(
      jsonTovalue(
          map["volumes"]["2"]["layers"]["2"]["sensitive"]["2"]["Thing"]),
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
