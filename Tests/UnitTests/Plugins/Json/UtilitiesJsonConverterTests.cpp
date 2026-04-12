// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"
#include "ActsPlugins/Json/UtilitiesJsonConverter.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>
#include <fstream>
#include <initializer_list>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "EqualityHelpers.hpp"

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(BinUtilityRoundTripTests) {
  BinUtility reference(2, 0., 4., open, AxisDirection::AxisR);

  std::ofstream out;

  // Test in one dimension
  nlohmann::json joneDimOut;
  to_json(joneDimOut, reference);
  out.open("BinUtility_1D.json");
  out << joneDimOut.dump(2);
  out.close();

  auto in = std::ifstream("BinUtility_1D.json",
                          std::ifstream::in | std::ifstream::binary);
  BOOST_CHECK(in.good());
  nlohmann::json joneDimIn;
  in >> joneDimIn;
  in.close();

  BinUtility test;
  from_json(joneDimIn, test);

  BOOST_CHECK(isEqual(reference, test, 0.0001));

  // Increase to two dimensions
  reference += BinUtility(10., -std::numbers::pi, std::numbers::pi, closed,
                          AxisDirection::AxisPhi);
  nlohmann::json jtwoDimOut;
  to_json(jtwoDimOut, reference);
  out.open("BinUtility_2D.json");
  out << jtwoDimOut.dump(2);
  out.close();

  in = std::ifstream("BinUtility_2D.json",
                     std::ifstream::in | std::ifstream::binary);
  BOOST_CHECK(in.good());
  nlohmann::json jtwoDimIn;
  in >> jtwoDimIn;
  in.close();

  test = BinUtility();
  from_json(jtwoDimIn, test);

  BOOST_CHECK(isEqual(reference, test, 0.0001));

  // Increase to three dimensions
  std::vector<float> boundaries = {-4., -1.5, 0., 10.};
  reference += BinUtility(boundaries, open, AxisDirection::AxisZ);
  nlohmann::json jthreeDimOut;
  to_json(jthreeDimOut, reference);
  out.open("BinUtility_3D.json");
  out << jthreeDimOut.dump(2);
  out.close();

  in = std::ifstream("BinUtility_3D.json",
                     std::ifstream::in | std::ifstream::binary);
  BOOST_CHECK(in.good());
  nlohmann::json jthreeDimIn;
  in >> jthreeDimIn;
  in.close();

  test = BinUtility();
  from_json(jthreeDimIn, test);

  BOOST_CHECK(isEqual(reference, test, 0.0001));

  // One with transform
  Transform3 t;
  t = Eigen::AngleAxis(0.12334, Vector3(1., 2., 3).normalized());
  t.pretranslate(Vector3(1., 2., 3.));

  auto bData = reference.binningData()[0];

  reference = BinUtility(bData, t);

  nlohmann::json jtransformOut;
  to_json(jtransformOut, reference);
  out.open("BinUtility_Transform.json");
  out << jtransformOut.dump(2);
  out.close();

  in = std::ifstream("BinUtility_Transform.json",
                     std::ifstream::in | std::ifstream::binary);
  BOOST_CHECK(in.good());
  nlohmann::json jtransformIn;
  in >> jtransformIn;
  in.close();

  test = BinUtility();
  from_json(jtransformIn, test);

  BOOST_CHECK(isEqual(reference, test, 0.0001));
}

BOOST_AUTO_TEST_CASE(Range1DRoundTrip) {
  Range1D<double> r(-10., 100.);

  nlohmann::json jrange;
  jrange["range"] = r;

  Range1D<double> rIn = jrange["range"];

  CHECK_CLOSE_ABS(rIn.min(), -10., 10e-5);
  CHECK_CLOSE_ABS(rIn.max(), 100., 10e-5);
}

BOOST_AUTO_TEST_CASE(ProtoAxisRoundTripTests) {
  ProtoAxis protoEq(AxisBoundaryType::Bound, -5., 5., 10u);
  nlohmann::json jProtoEq = protoEq;
  ProtoAxis protoEqRead(AxisBoundaryType::Bound, 1u);
  from_json(jProtoEq, protoEqRead);
  BOOST_CHECK_EQUAL(protoEqRead.getAxis(), protoEq.getAxis());
  BOOST_CHECK_EQUAL(protoEqRead.isAutorange(), protoEq.isAutorange());

  ProtoAxis protoVar(AxisBoundaryType::Closed,
                     std::vector<double>{-2., -1., 3.});
  nlohmann::json jProtoVar = protoVar;
  ProtoAxis protoVarRead(AxisBoundaryType::Bound, 1u);
  from_json(jProtoVar, protoVarRead);
  BOOST_CHECK_EQUAL(protoVarRead.getAxis(), protoVar.getAxis());
  BOOST_CHECK_EQUAL(protoVarRead.isAutorange(), protoVar.isAutorange());
}

BOOST_AUTO_TEST_CASE(DirectedProtoAxisRoundTripTests) {
  DirectedProtoAxis dProtoEq(AxisDirection::AxisPhi, AxisBoundaryType::Closed,
                             -std::numbers::pi, std::numbers::pi, 8u);
  nlohmann::json jdProtoEq = dProtoEq;
  DirectedProtoAxis dProtoEqRead(AxisDirection::AxisX, AxisBoundaryType::Bound,
                                 1u);
  from_json(jdProtoEq, dProtoEqRead);
  BOOST_CHECK_EQUAL(dProtoEqRead.getAxisDirection(),
                    dProtoEq.getAxisDirection());
  BOOST_CHECK_EQUAL(dProtoEqRead.getAxis(), dProtoEq.getAxis());
  BOOST_CHECK_EQUAL(dProtoEqRead.isAutorange(), dProtoEq.isAutorange());

  DirectedProtoAxis dProtoVar(AxisDirection::AxisR, AxisBoundaryType::Bound,
                              std::vector<double>{0., 1., 4., 9.});
  nlohmann::json jdProtoVar = dProtoVar;
  DirectedProtoAxis dProtoVarRead(AxisDirection::AxisX, AxisBoundaryType::Bound,
                                  1u);
  from_json(jdProtoVar, dProtoVarRead);
  BOOST_CHECK_EQUAL(dProtoVarRead.getAxisDirection(),
                    dProtoVar.getAxisDirection());
  BOOST_CHECK_EQUAL(dProtoVarRead.getAxis(), dProtoVar.getAxis());
  BOOST_CHECK_EQUAL(dProtoVarRead.isAutorange(), dProtoVar.isAutorange());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
