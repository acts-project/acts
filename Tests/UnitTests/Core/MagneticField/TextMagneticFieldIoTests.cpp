// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/TextMagneticFieldIo.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "ActsTests/CommonHelpers/TemporaryDirectory.hpp"

#include <fstream>

namespace bdata = boost::unit_test::data;

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(MagneticFieldSuite)

BOOST_AUTO_TEST_CASE(InterpolatedBFieldMap_rz_from_text) {
  TemporaryDirectory tmp{};

  // This tests octant reading, different deliminators, different scales
  for (auto fieldScale : {1_T, 1_kGauss}) {
    for (auto [id, cdlm] : enumerate(std::vector<std::string>{";", ",", ""})) {
      for (auto fOctant : {true, false}) {
        std::string fieldName = "rz_magfield" + std::to_string(id) + "_" +
                                std::to_string(static_cast<int>(fOctant)) +
                                "fieldScale" + std::to_string(fieldScale) +
                                ".csv";

        std::string dlm = cdlm.empty() ? " " : cdlm;

        // Write out a simple csv file with r,z,Br,Bz
        std::ofstream csvFile(tmp.path() / fieldName);

        // header should be ignored if indicated by %
        csvFile << "% r,z,Br,Bz\n";
        // header should be ignored if indicated by #
        csvFile << "# this is a B-Field file in r/z\n";
        // Empty line should be ignored
        csvFile << '\n';

        // Innermost r-bin
        if (!fOctant) {
          csvFile << "0." << dlm << "-2." << dlm << "0." << dlm << "0.2\n";
          csvFile << "0." << dlm << "-1." << dlm << "0.1" << dlm << "1.8\n";
        }
        csvFile << "0." << dlm << "0." << dlm << "0." << dlm << "2.0\n";
        csvFile << "0." << dlm << "1." << dlm << "0.1" << dlm << "1.8\n";
        csvFile << "0." << dlm << "2." << dlm << "0." << dlm << "0.2\n";
        // Middle r bin
        if (!fOctant) {
          csvFile << "0.5" << dlm << "-2." << dlm << "0." << dlm << "0.1\n";
          csvFile << "0.5" << dlm << "-1." << dlm << "0.2" << dlm << "1.7\n";
        }
        csvFile << "0.5" << dlm << "0." << dlm << "0.1" << dlm << "1.9\n";
        csvFile << "0.5" << dlm << "1." << dlm << "0.2" << dlm << "1.7\n";
        csvFile << "0.5" << dlm << "2." << dlm << "0." << dlm << "0.1\n";
        // Outer r bin
        if (!fOctant) {
          csvFile << "1." << dlm << "-2." << dlm << "0." << dlm << "0.\n";
          csvFile << "1." << dlm << "-1." << dlm << "0.1" << dlm << "1.5\n";
        }
        csvFile << "1." << dlm << "0." << dlm << "0.1" << dlm << "1.7\n";
        csvFile << "1." << dlm << "1." << dlm << "0.2" << dlm << "1.5\n";
        csvFile << "1." << dlm << "2." << dlm << "0." << dlm << "0.\n";
        csvFile.close();

        auto rzField = makeMagneticFieldMapRzFromText(
            [](std::array<std::size_t, 2> binsRZ,
               std::array<std::size_t, 2> nBinsRZ) {
              return (binsRZ.at(0) * nBinsRZ.at(1) + binsRZ.at(1));
            },
            tmp.path() / fieldName, 1_mm, fieldScale, fOctant, dlm);

        // Check that the bfield is two dimensional
        auto nBins = rzField.getNBins();
        BOOST_CHECK_EQUAL(nBins.size(), 2u);

        // Check number of bins in r and z
        BOOST_CHECK_EQUAL(nBins.at(0), 3u);
        BOOST_CHECK_EQUAL(nBins.at(1), 5u);

        // Check that the bin edges are correct
        auto mins = rzField.getMin();
        auto maxs = rzField.getMax();
        BOOST_CHECK_EQUAL(mins.size(), 2u);
        BOOST_CHECK_EQUAL(maxs.size(), 2u);

        // Get the unchecked field in the middle
        auto centralField = rzField.getFieldUnchecked({0., .0, 0.0});
        BOOST_CHECK(centralField.isApprox(Vector3(0., 0., 2. * fieldScale)));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(InterpolatedBFieldMap_xyz_from_text) {
  TemporaryDirectory tmp{};
  // This tests octant reading, different deliminators, different scales
  for (auto fieldScale : {1_T, 1_kGauss}) {
    for (auto [id, cdlm] : enumerate(std::vector<std::string>{";", ",", ""})) {
      for (auto fOctant : {true, false}) {
        std::string fieldName = "xyz_magfield" + std::to_string(id) + "_" +
                                std::to_string(static_cast<int>(fOctant)) +
                                "fieldScale" + std::to_string(fieldScale) +
                                ".csv";

        std::string dlm = cdlm.empty() ? " " : cdlm;

        // Write out a simple csv file with x,y,z,Bx,By,Bz
        std::ofstream csvFile(fieldName);

        // header should be ignored if indicated by %
        csvFile << "% x,y,z,Bx,By,Bz\n";
        // header should be ignored if indicated by #
        csvFile << "# this is a B-Field file in x/y/z\n";
        // Empty line should be ignored
        csvFile << '\n';

        std::vector<double> xvals = {0., 0.5, 1.};
        std::vector<double> yvals = {0., 1., 2.};
        std::vector<double> zvals = {0., 1., 2., 3.};
        if (!fOctant) {
          xvals = {-1., -0.5, 0., 0.5, 1.};
          yvals = {-2., -1., 0., 1., 2.};
          zvals = {-3., -2., -1., 0., 1., 2., 3.};
        }
        for (auto xv : xvals) {
          for (auto yv : yvals) {
            for (auto zv : zvals) {
              csvFile << xv << dlm << yv << dlm << zv << dlm << 0. << dlm << 0.
                      << dlm << 2.1 - std::abs(zv) * 0.1 << '\n';
            }
          }
        }
        csvFile.close();

        auto xyzField = makeMagneticFieldMapXyzFromText(
            [](std::array<std::size_t, 3> binsXYZ,
               std::array<std::size_t, 3> nBinsXYZ) {
              return (binsXYZ.at(0) * nBinsXYZ.at(1) * nBinsXYZ.at(2) +
                      binsXYZ.at(1) * nBinsXYZ.at(2) + binsXYZ.at(2));
            },
            fieldName, 1_mm, fieldScale, true, dlm);

        // Check that the bfield is three-dimenstional
        auto nBins = xyzField.getNBins();
        BOOST_CHECK_EQUAL(nBins.size(), 3u);
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
