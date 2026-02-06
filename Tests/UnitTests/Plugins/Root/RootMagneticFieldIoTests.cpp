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
#include "Acts/Utilities/Enumerate.hpp"
#include "ActsPlugins/Root/RootMagneticFieldIo.hpp"
#include "ActsTests/CommonHelpers/TemporaryDirectory.hpp"

#include <fstream>

#include <TFile.h>
#include <TTree.h>

namespace bdata = boost::unit_test::data;

using namespace Acts;
using namespace ActsPlugins;
using namespace Acts::UnitLiterals;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(RootSuite)

BOOST_AUTO_TEST_CASE(InterpolatedBFieldMap_rz_from_root) {
  TemporaryDirectory tmp{};

  auto fileName = tmp.path() / "rz_magfield.root";
  std::string fieldName = "BFieldMap";
  // Create the file
  auto rzFile = TFile::Open(fileName.c_str(), "RECREATE");
  auto tree = new TTree(fieldName.c_str(), "BFieldMap");
  double r = 0.;
  double z = 0.;
  double Br = 0.;
  double Bz = 0.;
  tree->Branch("r", &r, "r/D");
  tree->Branch("z", &z, "z/D");
  tree->Branch("Br", &Br, "Br/D");
  tree->Branch("Bz", &Bz, "Bz/D");
  for (auto rv : {0., 0.1, 0.2, 0.3, 0.4, 0.5}) {
    for (auto zv : {-2., -1., 0., 1., 2.}) {
      r = rv;
      z = zv;
      Br = r * 0.1;
      Bz = z * 0.1 + 2.;
      tree->Fill();
    }
  }
  tree->Write();
  rzFile->Close();

  // Now read it back in
  auto rzField = makeMagneticFieldMapRzFromRoot(
      [](std::array<std::size_t, 2> binsRZ,
         std::array<std::size_t, 2> nBinsRZ) {
        return (binsRZ.at(0) * nBinsRZ.at(1) + binsRZ.at(1));
      },
      fileName, fieldName, 1_mm, 1_T, false);

  // Check that the bfield is two dimensional
  auto nBins = rzField.getNBins();
  BOOST_CHECK_EQUAL(nBins.size(), 2u);

  // Check number of bins in r and z
  BOOST_CHECK_EQUAL(nBins.at(0), 6u);
  BOOST_CHECK_EQUAL(nBins.at(1), 5u);

  // Check that the bin edges are correct
  auto mins = rzField.getMin();
  auto maxs = rzField.getMax();
  BOOST_CHECK_EQUAL(mins.size(), 2u);
  BOOST_CHECK_EQUAL(maxs.size(), 2u);
}

BOOST_AUTO_TEST_CASE(InterpolatedBFieldMap_xyz_from_root) {
  TemporaryDirectory tmp{};

  for (auto fOctant : {true, false}) {
    std::string octantSuffix = fOctant ? "octant" : "full";

    std::string baseFileName = "xyz_magfield" + octantSuffix + ".root";
    auto fileName = tmp.path() / baseFileName;
    std::string fieldName = "BFieldMap";

    // Create the file
    auto xyzFile = TFile::Open(fileName.c_str(), "RECREATE");
    auto tree = new TTree(fieldName.c_str(), "BFieldMap");

    double x = 0.;
    double y = 0.;
    double z = 0.;
    double Bx = 0.;
    double By = 0.;
    double Bz = 0.;
    tree->Branch("x", &x, "x/D");
    tree->Branch("y", &y, "y/D");
    tree->Branch("z", &z, "z/D");
    tree->Branch("Bx", &Bx, "Bx/D");
    tree->Branch("By", &By, "By/D");
    tree->Branch("Bz", &Bz, "Bz/D");
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
          x = xv;
          y = yv;
          z = zv;
          Bx = x * 0.1;
          By = y * 0.1;
          Bz = z * 0.1 + 2.;
          tree->Fill();
        }
      }
    }

    tree->Write();
    xyzFile->Close();

    // Now read it back in
    auto xyzField = makeMagneticFieldMapXyzFromRoot(
        [](std::array<std::size_t, 3> binsXYZ,
           std::array<std::size_t, 3> nBinsXYZ) {
          return (binsXYZ.at(0) * nBinsXYZ.at(1) * nBinsXYZ.at(2) +
                  binsXYZ.at(1) * nBinsXYZ.at(2) + binsXYZ.at(2));
        },
        fileName, fieldName, 1_mm, 1_T, fOctant);

    // Check that the bfield is three-dimenstional
    auto nBins = xyzField.getNBins();
    BOOST_CHECK_EQUAL(nBins.size(), 3u);

    // Check number of bins in x, y and z
    BOOST_CHECK_EQUAL(nBins.at(0), 5u);
    BOOST_CHECK_EQUAL(nBins.at(1), 5u);
    BOOST_CHECK_EQUAL(nBins.at(2), 7u);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
