// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Plugins/Json/AlgebraJsonConverter.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

using namespace Acts;

BOOST_AUTO_TEST_SUITE(AlgebraJsonConversion)

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

BOOST_AUTO_TEST_CASE(TransformNullIdentity) {
  // An identity matrix
  Transform3 reference = Transform3::Identity();

  // Test the identity transform with nulled
  Transform3JsonConverter::Options nulledOption{false, false};
  nlohmann::json nulledOut =
      Transform3JsonConverter::toJson(reference, nulledOption);
  BOOST_CHECK_EQUAL(nulledOut["translation"], nullptr);
  BOOST_CHECK_EQUAL(nulledOut["rotation"], nullptr);

  // Test with writing the identity
  Transform3JsonConverter::Options writtenOption{true, false};
  nlohmann::json writtenOut =
      Transform3JsonConverter::toJson(reference, writtenOption);
  BOOST_CHECK_NE(writtenOut["translation"], nullptr);
  BOOST_CHECK_NE(writtenOut["rotation"], nullptr);
}

BOOST_AUTO_TEST_CASE(TransformTranspose) {
  // An identity matrix
  Transform3 reference = Transform3::Identity();
  reference.pretranslate(Vector3(1., 2., 3.));
  reference.rotate(Eigen::AngleAxis(0.12334, Vector3(1., 2., 3).normalized()));

  std::vector<double> referenceT = {1., 2., 3.};
  std::vector<double> referenceR = {0.992946,   -0.0975562, 0.0673888,
                                    0.0997267,  0.994574,   -0.0296247,
                                    -0.0641331, 0.0361362,  0.997287};

  // Test standard writing
  Transform3JsonConverter::Options standardOptions{true, false};
  nlohmann::json standardOut =
      Transform3JsonConverter::toJson(reference, standardOptions);
  // Check translation read back in
  BOOST_CHECK(standardOut["translation"].get<std::vector<double>>() ==
              referenceT);
  // Check rotation read back in - not transposed
  std::vector<double> readR =
      standardOut["rotation"].get<std::vector<double>>();
  for (auto [i, rr] : Acts::enumerate(referenceR)) {
    CHECK_CLOSE_ABS(readR[i], rr, 1e-5);
  }

  // Test transposed writing
  Transform3JsonConverter::Options transposeOptions{true, true};
  nlohmann::json transposeOut =
      Transform3JsonConverter::toJson(reference, transposeOptions);
  // Check translation read back in
  BOOST_CHECK(transposeOut["translation"].get<std::vector<double>>() ==
              referenceT);

  // Check rotation read back in - transposed
  std::vector<std::size_t> transposedIndices = {0, 3, 6, 1, 4, 7, 2, 5, 8};
  readR = transposeOut["rotation"].get<std::vector<double>>();
  for (auto [i, rr] : Acts::enumerate(referenceR)) {
    CHECK_CLOSE_ABS(readR[transposedIndices[i]], rr, 1e-5);
  }
}

BOOST_AUTO_TEST_CASE(IdentifiedTransformsJson) {
  std::unordered_map<GeometryIdentifier, Transform3> transformMap;
  std::ofstream out;

  Transform3 reference0 = Transform3::Identity();
  reference0.pretranslate(Vector3(1., 2., 3.));
  reference0.rotate(Eigen::AngleAxis(0.12334, Vector3(1., 2., 3).normalized()));

  std::vector<double> referenceT0 = {1., 2., 3.};
  std::vector<double> referenceR0 = {0.992946,   -0.0975562, 0.0673888,
                                     0.0997267,  0.994574,   -0.0296247,
                                     -0.0641331, 0.0361362,  0.997287};

  Transform3 reference1 = Transform3::Identity();

  auto geometryId0 =
      GeometryIdentifier().withVolume(1).withLayer(2).withSensitive(3);
  auto geometryId1 = GeometryIdentifier().withVolume(10).withBoundary(4);

  transformMap[geometryId0] = reference0;
  transformMap[geometryId1] = reference1;

  // Write the standard map out
  IdentifiedTransform3JsonConverter::Options jOptions{};

  nlohmann::json itsOut =
      IdentifiedTransform3JsonConverter::toJson(transformMap, jOptions);

  out.open("IdentifiedTransforms.json");
  out << itsOut.dump(4);
  out.close();

  jOptions.compressIdentifier = true;
  nlohmann::json itsOutCompressed =
      IdentifiedTransform3JsonConverter::toJson(transformMap, jOptions);

  out.open("IdentifiedTransformsCompressed.json");
  out << itsOutCompressed.dump(4);
  out.close();

  auto itsIn = IdentifiedTransform3JsonConverter::fromJson(itsOut);

  BOOST_CHECK_EQUAL(itsIn.size(), transformMap.size());
  BOOST_CHECK(itsIn[geometryId0].isApprox(transformMap[geometryId0]));
  BOOST_CHECK(itsIn[geometryId1].isApprox(transformMap[geometryId1]));

  auto itsInCompressed =
      IdentifiedTransform3JsonConverter::fromJson(itsOutCompressed);

  BOOST_CHECK_EQUAL(itsInCompressed.size(), transformMap.size());
  BOOST_CHECK(itsInCompressed[geometryId0].isApprox(transformMap[geometryId0]));
  BOOST_CHECK(itsInCompressed[geometryId1].isApprox(transformMap[geometryId1]));
}

BOOST_AUTO_TEST_SUITE_END()
