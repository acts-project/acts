// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/Json/AlgebraJsonConverter.hpp"

#include <fstream>
#include <string>
#include <utility>

#include <nlohmann/json.hpp>

using namespace Acts;

BOOST_AUTO_TEST_SUITE(AlgebraJsonConverter)

BOOST_AUTO_TEST_CASE(TransformRoundTripTests) {
  Transform3 reference = Transform3::Identity();

  std::ofstream out;

  // Test the identity transform
  nlohmann::json identityOut;
  to_json(identityOut, reference);
  out.open("Transform3_Identity.json");
  out << identityOut.dump(2);
  out.close();

  auto in = std::ifstream("Transform3_Identity.json",
                          std::ifstream::in | std::ifstream::binary);
  BOOST_CHECK(in.good());
  nlohmann::json identityIn;
  in >> identityIn;
  in.close();

  Transform3 test;
  from_json(identityIn, test);

  BOOST_CHECK(test.isApprox(reference));

  // Test a pure translation transform
  reference.pretranslate(Vector3(1., 2., 3.));

  nlohmann::json translationOut;
  to_json(translationOut, reference);
  out.open("Transform3_Translation.json");
  out << translationOut.dump(2);
  out.close();

  in = std::ifstream("Transform3_Translation.json",
                     std::ifstream::in | std::ifstream::binary);
  BOOST_CHECK(in.good());
  nlohmann::json translationIn;
  in >> translationIn;
  in.close();

  test = Transform3::Identity();
  from_json(translationIn, test);

  BOOST_CHECK(test.isApprox(reference));

  // Test a full transform
  reference = Eigen::AngleAxis(0.12334, Vector3(1., 2., 3).normalized());
  reference.pretranslate(Vector3(1., 2., 3.));

  nlohmann::json fullOut;
  to_json(fullOut, reference);
  out.open("Transform3_Full.json");
  out << fullOut.dump(2);
  out.close();

  in = std::ifstream("Transform3_Full.json",
                     std::ifstream::in | std::ifstream::binary);
  BOOST_CHECK(in.good());
  nlohmann::json fullIn;
  in >> fullIn;
  in.close();

  test = Transform3::Identity();
  from_json(fullIn, test);

  BOOST_CHECK(test.isApprox(reference));
}

BOOST_AUTO_TEST_SUITE_END()