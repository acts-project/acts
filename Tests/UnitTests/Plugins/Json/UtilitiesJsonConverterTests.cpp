// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Plugins/Json/UtilitiesJsonConverter.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningData.hpp"

#include <fstream>
#include <iostream>

using namespace Acts;

BOOST_AUTO_TEST_SUITE(UtilitiesJsonConverter)

namespace {
/// Check whether the BinningData objects are equal
///
/// @param ba The first BinningData object
/// @param bb The second BinningData object
/// @param tolerance a tolerance parameter
///
/// @return a boolean
bool isEqual(const BinningData& ba, const BinningData& bb, float tolerance) {
  bool equalBool = (ba.type == bb.type) and (ba.option == bb.option) and
                   (ba.binvalue == bb.binvalue) and (ba.zdim == bb.zdim) and
                   (ba.subBinningAdditive == bb.subBinningAdditive);

  bool equalRange = (std::abs(ba.min - bb.min) < tolerance) and
                    (std::abs(ba.max - bb.max) < tolerance) and
                    (std::abs(ba.step - bb.step) < tolerance);

  bool euqalStructure =
      (ba.subBinningData != nullptr)
          ? isEqual(*ba.subBinningData, *bb.subBinningData, tolerance)
          : (bb.subBinningData == nullptr);

  bool equalBoundaries = (ba.boundaries().size() == bb.boundaries().size());
  if (equalBoundaries) {
    for (size_t ib = 0; ib < ba.boundaries().size(); ++ib) {
      equalBoundaries =
          (std::abs(ba.boundaries()[ib] - bb.boundaries()[ib]) < tolerance);
      if (not equalBoundaries) {
        break;
      }
    }
  }
  return equalBool and equalRange and euqalStructure;
}

/// Check whether the BinUtility ojbects are equal
///
/// @param ba The first BinUtility object
/// @param bb the second BinUtility object
/// @param tolerance a tolerance parameter
///
/// @return a bollean if equal
bool isEqual(const BinUtility& ba, const BinUtility& bb, float tolerance) {
  bool equal = (ba.binningData().size() == bb.binningData().size());
  if (equal) {
    for (size_t ib = 0; ib < ba.binningData().size(); ++ib) {
      equal = isEqual(ba.binningData()[ib], bb.binningData()[ib], tolerance);
    }
  }
  return equal;
}

}  // namespace

BOOST_AUTO_TEST_CASE(BinUtilityRoundTripTests) {
  BinUtility reference(2, 0., 4., open, binR);

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
  reference += BinUtility(10., -M_PI, M_PI, closed, binPhi);
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
  reference += BinUtility(boundaries, open, binZ);
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

BOOST_AUTO_TEST_SUITE_END()
