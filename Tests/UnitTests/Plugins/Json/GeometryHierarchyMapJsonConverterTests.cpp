// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>
#include <fstream>
#include <nlohmann/json.hpp>
#include <ostream>

#include "Acts/Plugins/Json/GeometryHierarchyMapJsonConverter.hpp"
#include "Acts/Tests/CommonHelpers/DataDirectory.hpp"

namespace {

using Acts::GeometryID;
using nlohmann::json;

// helper function to create geometry ids.
GeometryID makeId(int volume = 0, int layer = 0, int sensitive = 0) {
  return GeometryID().setVolume(volume).setLayer(layer).setSensitive(sensitive);
}

// example element for the container

struct Thing {
  double x = 1.0;
  int y = 23;

  friend constexpr bool operator==(const Thing& lhs, const Thing& rhs) {
    return (lhs.x == rhs.x) and (lhs.y == rhs.y);
  }
};

// custom Json encoder/decoder. naming is mandated by nlohman::json and thus can
// not match our naming guidelines.

void to_json(json& j, const Thing& t) {
  j = json{{"x", t.x}, {"y", t.y}};
}

void from_json(const json& j, Thing& t) {
  j.at("x").get_to(t.x);
  j.at("y").get_to(t.y);
}

std::ostream& operator<<(std::ostream& os, const Thing& t) {
  os << nlohmann::json(t);
  return os;
}

using Container = Acts::GeometryHierarchyMap<Thing>;
using Converter = Acts::GeometryHierarchyMapJsonConverter<Thing>;

}  // namespace

BOOST_TEST_DONT_PRINT_LOG_VALUE(json::iterator)
BOOST_TEST_DONT_PRINT_LOG_VALUE(Container::Iterator)

BOOST_AUTO_TEST_SUITE(GeometryHierarchyMapJsonConverter)

BOOST_AUTO_TEST_CASE(ToJson) {
  Container c = {
      {makeId(1), {2.0, -3}},
      {makeId(2, 3), {-4.5, 5}},
      {makeId(4, 5, 6), {7.25, -8}},
  };
  json j = Converter("thing").toJson(c);

  BOOST_TEST(j.is_object());
  // check header
  auto header = j.find("acts-geometry-hierarchy-map");
  BOOST_TEST(header != j.end());
  BOOST_TEST(header->is_object());
  BOOST_TEST(header->at("format-version").is_number_integer());
  BOOST_TEST(header->at("format-version").get<int>() == 0);
  BOOST_TEST(header->at("value-identifier").is_string());
  BOOST_TEST(header->at("value-identifier").get<std::string>() == "thing");
  // check entries
  auto entries = j.find("entries");
  BOOST_TEST(entries != j.end());
  BOOST_TEST(entries->is_array());
  BOOST_TEST(entries->size() == 3u);
}

BOOST_AUTO_TEST_CASE(FromJson) {
  json j = {
      {
          "acts-geometry-hierarchy-map",
          {
              {"format-version", 0},
              {"value-identifier", "thing"},
          },
      },
      {
          "entries",
          {
              {
                  {"volume", 2},
                  {"layer", 3},
                  {"value", {{"x", 4.0}, {"y", 4}}},
              },
              {
                  {"volume", 5},
                  {"layer", 6},
                  {"sensitive", 7},
                  {"value", {{"x", 3.0}, {"y", 3}}},
              },
          },
      },
  };
  Container c = Converter("thing").fromJson(j);

  BOOST_TEST(not c.empty());
  BOOST_TEST(c.size() == 2);
  {
    auto it = c.find(makeId(2, 3));
    BOOST_TEST(it != c.end());
    BOOST_TEST(it->x == 4.0);
    BOOST_TEST(it->y == 4);
  }
  {
    auto it = c.find(makeId(5, 6, 7));
    BOOST_TEST(it != c.end());
    BOOST_TEST(it->x == 3.0);
    BOOST_TEST(it->y == 3);
  }
}

BOOST_AUTO_TEST_CASE(FromJsonMissingHeader) {
  json j = {
      {"entries", {}},
  };
  BOOST_CHECK_THROW(Converter("an-identifier").fromJson(j),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(FromJsonInvalidFormatVersion) {
  json j = {
      {
          "acts-geometry-hierarchy-map",
          {
              {"format-version", -1},
              {"value-identifier", "an-identifier"},
          },
      },
      {"entries", {}},
  };
  BOOST_CHECK_THROW(Converter("an-identifier").fromJson(j),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(FromJsonInvalidValueIdentifier) {
  json j = {
      {
          "acts-geometry-hierarchy-map",
          {
              {"format-version", 0},
              {"value-identifier", "an-identifier"},
          },
      },
      {"entries", {}},
  };
  BOOST_CHECK_THROW(Converter("not-the-identifier").fromJson(j),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(FromJsonMissingEntries) {
  json j = {
      {
          "acts-geometry-hierarchy-map",
          {
              {"format-version", 0},
              {"value-identifier", "an-identifier"},
          },
      },
  };
  BOOST_CHECK_THROW(Converter("an-identifier").fromJson(j),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(Roundtrip) {
  Container c0 = {
      {makeId(1), {2.0, -3}},
      {makeId(2, 3), {-4.5, 5}},
      {makeId(4, 5, 6), {7.25, -8}},
  };
  auto j = Converter("the-identifier").toJson(c0);
  auto c1 = Converter("the-identifier").fromJson(j);

  BOOST_TEST(c0.size() == c1.size());
  for (auto i = std::min(c0.size(), c1.size()); 0 < i--;) {
    BOOST_TEST(c0.idAt(i) == c1.idAt(i));
    BOOST_TEST(c0.valueAt(i) == c1.valueAt(i));
  }
}

BOOST_AUTO_TEST_CASE(FromFile) {
  // read json data from file
  auto path = Acts::Test::getDataPath("geometry-hierarchy-map.json");
  auto file = std::ifstream(path, std::ifstream::in | std::ifstream::binary);
  BOOST_TEST(file.good());
  json j;
  file >> j;
  BOOST_TEST(file.good());
  // convert json to container
  Container c = Converter("thing").fromJson(j);
  // check container content
  BOOST_TEST(not c.empty());
  BOOST_TEST(c.size() == 4);
  BOOST_TEST(c.find(makeId()) != c.end());
  BOOST_TEST(c.find(makeId(1, 2)) != c.end());
  BOOST_TEST(c.find(makeId(3)) != c.end());
  BOOST_TEST(c.find(makeId(3, 4)) != c.end());
}

BOOST_AUTO_TEST_SUITE_END()
