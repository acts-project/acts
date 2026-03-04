// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"
#include "ActsPlugins/Json/ProtoAxisJsonConverter.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(EquidistantProtoAxisJsonConversion) {
  using enum AxisBoundaryType;
  using enum AxisType;

  // Bound, equidistant axis
  ProtoAxis epab(Bound, 0.0, 1.0, 10);

  nlohmann::json jProtoAxis = ProtoAxisJsonConverter::toJson(epab);

  BOOST_CHECK(jProtoAxis.contains("axis"));
  BOOST_CHECK(jProtoAxis.contains("autorange"));

  ProtoAxis epabRead = ProtoAxisJsonConverter::fromJson(jProtoAxis);

  BOOST_CHECK_EQUAL(epabRead.getAxis(), epab.getAxis());
  BOOST_CHECK_EQUAL(epabRead.isAutorange(), epab.isAutorange());
  BOOST_CHECK_EQUAL(epabRead.toString(), epab.toString());
}

BOOST_AUTO_TEST_CASE(AutorangeProtoAxisJsonConversion) {
  using enum AxisBoundaryType;
  using enum AxisType;

  // Bound, equidistant axis, autorange
  ProtoAxis epa(Bound, 10);

  nlohmann::json jProtoAxis = ProtoAxisJsonConverter::toJson(epa);

  BOOST_CHECK(jProtoAxis.contains("axis"));
  BOOST_CHECK(jProtoAxis.contains("autorange"));

  ProtoAxis epaRead = ProtoAxisJsonConverter::fromJson(jProtoAxis);

  BOOST_CHECK_EQUAL(epaRead.getAxis(), epa.getAxis());
  BOOST_CHECK_EQUAL(epaRead.isAutorange(), epa.isAutorange());
  BOOST_CHECK_EQUAL(epaRead.toString(), epa.toString());
}

BOOST_AUTO_TEST_CASE(VariableProtoAxisJsonConversion) {
  using enum AxisBoundaryType;
  using enum AxisType;

  // Bound, variable axis
  ProtoAxis vpab(Bound, {0.0, 1.0, 10});

  nlohmann::json jProtoAxis = ProtoAxisJsonConverter::toJson(vpab);
  BOOST_CHECK(jProtoAxis.contains("axis"));
  BOOST_CHECK(jProtoAxis.contains("autorange"));

  ProtoAxis vpabRead = ProtoAxisJsonConverter::fromJson(jProtoAxis);

  BOOST_CHECK_EQUAL(vpabRead.getAxis(), vpab.getAxis());
  BOOST_CHECK_EQUAL(vpabRead.isAutorange(), vpab.isAutorange());
  BOOST_CHECK_EQUAL(vpabRead.toString(), vpab.toString());
}

BOOST_AUTO_TEST_CASE(InvalidAndValidInputJson) {
  // valid eq axis input
  nlohmann::json jValidEqAxis = {{"bins", 10},
                                 {"boundary_type", "Bound"},
                                 {"range", std::array<double, 2>{0.0, 1.0}},
                                 {"type", "Equidistant"}};

  // Valid input first
  nlohmann::json jValidEq = {{"axis", jValidEqAxis}, {"autorange", false}};

  BOOST_CHECK_NO_THROW(ProtoAxisJsonConverter::fromJson(jValidEq));

  // Invalid input - zero bins
  nlohmann::json jInvalidEqAxis = jValidEqAxis;
  jInvalidEqAxis["bins"] = 0;

  nlohmann::json jInvalidEq = {{"axis", jInvalidEqAxis}, {"autorange", false}};

  BOOST_CHECK_THROW(ProtoAxisJsonConverter::fromJson(jInvalidEq),
                    std::invalid_argument);

  // Invalid input - auto range without bins
  jInvalidEq = {{"axis", jInvalidEqAxis}, {"autorange", true}};
  BOOST_CHECK_THROW(ProtoAxisJsonConverter::fromJson(jInvalidEq),
                    std::invalid_argument);

  // Invalid input - min >= max
  jInvalidEqAxis = jValidEqAxis;
  jInvalidEqAxis["range"] = std::array<double, 2>{1.0, 0.0};

  jInvalidEq = {{"axis", jInvalidEqAxis}, {"autorange", false}};

  BOOST_CHECK_THROW(ProtoAxisJsonConverter::fromJson(jInvalidEq),
                    std::invalid_argument);

  nlohmann::json jValidVarAxis = {
      {"boundary_type", "Bound"},
      {"boundaries", std::vector<double>{0.0, 0.25, 0.75, 1.0}},
      {"type", "Variable"}};

  // Valid input first
  nlohmann::json jValidVar = {{"axis", jValidVarAxis}, {"autorange", false}};
  BOOST_CHECK_NO_THROW(ProtoAxisJsonConverter::fromJson(jValidVar));

  // Invalid input - less than two edges
  nlohmann::json jInvalidVarAxis = jValidVarAxis;
  jInvalidVarAxis["boundaries"] = std::vector<double>{0.0};

  nlohmann::json jInvalidVar = {{"axis", jInvalidVarAxis},
                                {"autorange", false}};
  BOOST_CHECK_THROW(ProtoAxisJsonConverter::fromJson(jInvalidVar),
                    std::invalid_argument);

  // Invalid input - non-increasing edges
  jInvalidVarAxis = jValidVarAxis;
  jInvalidVarAxis["boundaries"] = std::vector<double>{0.0, 0.75, 0.25, 1.0};

  jInvalidVar = {{"axis", jInvalidVarAxis}, {"autorange", false}};

  BOOST_CHECK_THROW(ProtoAxisJsonConverter::fromJson(jInvalidVar),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
