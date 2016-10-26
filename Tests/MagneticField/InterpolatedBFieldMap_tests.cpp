// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file InterpolatedBFieldMap_tests.cpp

#define BOOST_TEST_MODULE Mapped magnetic field tests
#include <boost/test/included/unit_test.hpp>
#include <fstream>
#include <memory>
#include <string>

#include <ACTS/MagneticField/InterpolatedBFieldMap.hpp>
#include "ACTS/Utilities/Units.hpp"

namespace Acts {

namespace Test {

  const std::string Points1  = "0 0 0 1 0 0\n";
  const std::string Points2  = Points1 + "0 0 1 1 1 0\n";
  const std::string Points3  = Points2 + "0 1 0 1 1 1\n";
  const std::string Points4  = Points3 + "0 1 1 2 1 1\n";
  const std::string Points5  = Points4 + "1 0 0 2 2 1\n";
  const std::string Points6  = Points5 + "1 0 1 2 2 2\n";
  const std::string Points7  = Points6 + "1 1 0 3 2 2\n";
  const std::string Points8  = Points7 + "1 1 1 3 3 2\n";
  std::string       Points27 = "0 0 0 1 0 0\n\
0 0 1 1 1 0\n\
0 0 2 1 1 1\n\
0 1 0 2 1 1\n\
0 1 1 2 2 1\n\
0 1 2 2 2 2\n\
0 2 0 3 2 2\n\
0 2 1 3 3 2\n\
0 2 2 3 3 3\n\
1 0 0 4 3 3\n\
    \n\
1 0 1 4 4 3\n\
1 0 2 4 4 4\n\
1 1 0 5 4 4\n\
1 1 1 5 5 4\n\
1 1 2 5 5 5\n\
1 2 0 6 5 5\n\
1 2 1 6 6 5\n\
1 2 2 6 6 6\n\
2 0 0 7 6 6\n\
2 0 1 7 7 6\n\
2 0 2 7 7 7\n\
2 1 0 8 7 7\n\
2 1 1 8 8 7\n\
2 1 2 8 8 8\n\
2 2 0 9 8 8\n\
2 2 1 9 9 8\n\
2 2 2 9 9 9\n";

  /// @brief helper class for unit testing InterpolatedBFieldMap
  ///
  /// In order to unit test private methods, a friend helper class is needed.
  struct InterpolatedBFieldTester
  {
    /// @cond
    InterpolatedBFieldTester(const std::string& inFile)
    {
      InterpolatedBFieldMap::Config c;
      c.fieldMapFile = inFile;
      BField         = std::make_unique<InterpolatedBFieldMap>(std::move(c));
    }

    size_t
    NPointsX() const
    {
      return BField->m_NPointsX;
    }

    size_t
    NPointsY() const
    {
      return BField->m_NPointsY;
    }

    size_t
    NPointsZ() const
    {
      return BField->m_NPointsZ;
    }

    size_t
    NBFieldValues() const
    {
      return BField->m_BField.size();
    }

    void
    getBinNumbers(const Vector3D& pos,
                  size_t&         xBin,
                  size_t&         yBin,
                  size_t&         zBin) const
    {
      BField->getBinNumbers(pos, xBin, yBin, zBin);
    }

    size_t
    globalIndex(size_t xBin, size_t yBin, size_t zBin) const
    {
      return BField->globalIndex(xBin, yBin, zBin);
    }

    bool
    insideGrid(const Vector3D& pos) const
    {
      return BField->insideGrid(pos);
    }

    std::unique_ptr<InterpolatedBFieldMap> BField = nullptr;
    /// @endcond
  };

  /// @brief unit test for interpolated BField from map with a single point
  ///
  /// Tests the correct behaviour and consistency of
  /// -# InterpolatedBFieldMap::InterpolatedBFieldMap
  /// -# InterpolatedBFieldMap::initialize
  /// -# InterpolatedBFieldMap::insideGrid
  BOOST_AUTO_TEST_CASE(InterpolatedBFieldMap_1Point)
  {
    std::ofstream map_file("TmpField.txt", std::ios::out | std::ios::trunc);
    map_file << Points1;
    map_file.close();

    InterpolatedBFieldTester t("TmpField.txt");

    BOOST_TEST(t.NPointsX() == 1u);
    BOOST_TEST(t.NPointsY() == 1u);
    BOOST_TEST(t.NPointsZ() == 1u);

    BOOST_TEST(t.NBFieldValues() == 3u);

    BOOST_TEST(not t.insideGrid(Vector3D(0, 0, 0)));
    BOOST_TEST(not t.insideGrid(Vector3D(0, 0, 1)));
    BOOST_TEST(not t.insideGrid(Vector3D(0, 0.1, 0)));
    BOOST_TEST(not t.insideGrid(Vector3D(-1e-5, 0, 0)));
  }

  /// @brief unit test for interpolated BField from map with a single bin
  ///
  /// A single bin in 3D requires eight grid points.
  ///
  /// Tests the correct behaviour and consistency of
  /// -# InterpolatedBFieldMap::InterpolatedBFieldMap
  /// -# InterpolatedBFieldMap::getBinNumbers
  /// -# InterpolatedBFieldMap::getField
  /// -# InterpolatedBFieldMap::globalIndex
  /// -# InterpolatedBFieldMap::initialize
  /// -# InterpolatedBFieldMap::insideGrid
  BOOST_AUTO_TEST_CASE(InterpolatedBFieldMap_8Points)
  {
    std::ofstream map_file("TmpField.txt", std::ios::out | std::ios::trunc);
    map_file << Points8;
    map_file.close();

    InterpolatedBFieldTester t("TmpField.txt");
    size_t                   xBin = 99, yBin = 99, zBin = 99;

    BOOST_TEST(t.NPointsX() == 2u);
    BOOST_TEST(t.NPointsY() == 2u);
    BOOST_TEST(t.NPointsZ() == 2u);

    BOOST_TEST(t.NBFieldValues() == 24u);

    BOOST_TEST(t.insideGrid(Vector3D(0, 0, 0)));
    BOOST_TEST(t.insideGrid(Vector3D(0, 0.99, 0)));
    BOOST_TEST(t.insideGrid(Vector3D(0, 0, 0.99)));
    BOOST_TEST(t.insideGrid(Vector3D(0, 0.99, 0.99)));
    BOOST_TEST(t.insideGrid(Vector3D(0.99, 0, 0)));
    BOOST_TEST(t.insideGrid(Vector3D(0.99, 0.99, 0)));
    BOOST_TEST(t.insideGrid(Vector3D(0.99, 0, 0.99)));
    BOOST_TEST(not t.insideGrid(Vector3D(0, 0, 1)));
    BOOST_TEST(not t.insideGrid(Vector3D(0, 1, 0)));
    BOOST_TEST(not t.insideGrid(Vector3D(1, 0, 0)));
    BOOST_TEST(not t.insideGrid(Vector3D(1, 1, 1)));
    BOOST_TEST(not t.insideGrid(Vector3D(0, 0, -0.1)));
    BOOST_TEST(not t.insideGrid(Vector3D(0, -0.1, 0)));
    BOOST_TEST(not t.insideGrid(Vector3D(-0.1, 0, 0)));
    BOOST_TEST(not t.insideGrid(Vector3D(-0.1, 0.2, 0.3)));

    t.getBinNumbers(Vector3D(0, 0, 0), xBin, yBin, zBin);
    BOOST_TEST(xBin == 0u);
    BOOST_TEST(yBin == 0u);
    BOOST_TEST(zBin == 0u);

    t.getBinNumbers(Vector3D(0.5, 0.5, 0.5), xBin, yBin, zBin);
    BOOST_TEST(xBin == 0u);
    BOOST_TEST(yBin == 0u);
    BOOST_TEST(zBin == 0u);

    BOOST_TEST(t.globalIndex(0, 0, 0) == 0u);

    // clang-format off
    BOOST_TEST((t.BField->getField(Vector3D(0, 0, 0)) / units::_T).isApprox(Vector3D(1, 0, 0)));
    BOOST_TEST((t.BField->getField(Vector3D(0.2, 0, 0)) / units::_T).isApprox(Vector3D(1.2, 0.4, 0.2)));
    BOOST_TEST((t.BField->getField(Vector3D(0, 0.6, 0)) / units::_T).isApprox(Vector3D(1, 0.6, 0.6)));
    BOOST_TEST((t.BField->getField(Vector3D(0, 0, 0.7)) / units::_T).isApprox(Vector3D(1, 0.7, 0)));
    BOOST_TEST((t.BField->getField(Vector3D(0.5, 0.5, 0.5)) / units::_T).isApprox(Vector3D(1.875, 1.5, 1.125)));
    // clang-format on
  }

  /// @brief unit test for interpolated BField from map with 8 bins
  ///
  /// Eight bins in 3D requires 3x3x3 = 27 grid points.
  ///
  /// Tests the correct behaviour and consistency of
  /// -# InterpolatedBFieldMap::InterpolatedBFieldMap
  /// -# InterpolatedBFieldMap::getBinNumbers
  /// -# InterpolatedBFieldMap::getField
  /// -# InterpolatedBFieldMap::globalIndex
  /// -# InterpolatedBFieldMap::initialize
  /// -# InterpolatedBFieldMap::insideGrid
  BOOST_AUTO_TEST_CASE(InterpolatedBFieldMap_27Points)
  {
    std::ofstream map_file("TmpField.txt", std::ios::out | std::ios::trunc);
    map_file << Points27;
    map_file.close();

    InterpolatedBFieldTester t("TmpField.txt");
    size_t                   xBin = 99, yBin = 99, zBin = 99;

    BOOST_TEST(t.NPointsX() == 3u);
    BOOST_TEST(t.NPointsY() == 3u);
    BOOST_TEST(t.NPointsZ() == 3u);

    BOOST_TEST(t.NBFieldValues() == 81u);

    BOOST_TEST(t.insideGrid(Vector3D(0, 0, 0)));
    BOOST_TEST(t.insideGrid(Vector3D(0, 1.99, 0)));
    BOOST_TEST(t.insideGrid(Vector3D(0, 0, 1.99)));
    BOOST_TEST(t.insideGrid(Vector3D(0, 1.99, 1.99)));
    BOOST_TEST(t.insideGrid(Vector3D(1.99, 0, 0)));
    BOOST_TEST(t.insideGrid(Vector3D(1.99, 1.99, 0)));
    BOOST_TEST(t.insideGrid(Vector3D(1.99, 0, 1.99)));
    BOOST_TEST(not t.insideGrid(Vector3D(0, 0, 2)));
    BOOST_TEST(not t.insideGrid(Vector3D(0, 2, 0)));
    BOOST_TEST(not t.insideGrid(Vector3D(2, 0, 0)));
    BOOST_TEST(not t.insideGrid(Vector3D(2, 2, 2)));
    BOOST_TEST(not t.insideGrid(Vector3D(0, 0, -0.1)));
    BOOST_TEST(not t.insideGrid(Vector3D(0, -0.1, 0)));
    BOOST_TEST(not t.insideGrid(Vector3D(-0.1, 0, 0)));
    BOOST_TEST(not t.insideGrid(Vector3D(-0.1, 0.2, 0.3)));

    t.getBinNumbers(Vector3D(0, 0, 0), xBin, yBin, zBin);
    BOOST_TEST(xBin == 0u);
    BOOST_TEST(yBin == 0u);
    BOOST_TEST(zBin == 0u);

    t.getBinNumbers(Vector3D(0.2, 1, 0.5), xBin, yBin, zBin);
    BOOST_TEST(xBin == 0u);
    BOOST_TEST(yBin == 1u);
    BOOST_TEST(zBin == 0u);

    t.getBinNumbers(Vector3D(0.2, 0.6, 1), xBin, yBin, zBin);
    BOOST_TEST(xBin == 0u);
    BOOST_TEST(yBin == 0u);
    BOOST_TEST(zBin == 1u);

    t.getBinNumbers(Vector3D(0.5, 1.2, 1.5), xBin, yBin, zBin);
    BOOST_TEST(xBin == 0u);
    BOOST_TEST(yBin == 1u);
    BOOST_TEST(zBin == 1u);

    t.getBinNumbers(Vector3D(1.5, 0, 0), xBin, yBin, zBin);
    BOOST_TEST(xBin == 1u);
    BOOST_TEST(yBin == 0u);
    BOOST_TEST(zBin == 0u);

    t.getBinNumbers(Vector3D(1, 1.2, 0.3), xBin, yBin, zBin);
    BOOST_TEST(xBin == 1u);
    BOOST_TEST(yBin == 1u);
    BOOST_TEST(zBin == 0u);

    t.getBinNumbers(Vector3D(1, 0.2, 1.3), xBin, yBin, zBin);
    BOOST_TEST(xBin == 1u);
    BOOST_TEST(yBin == 0u);
    BOOST_TEST(zBin == 1u);

    t.getBinNumbers(Vector3D(1.5, 1, 1), xBin, yBin, zBin);
    BOOST_TEST(xBin == 1u);
    BOOST_TEST(yBin == 1u);
    BOOST_TEST(zBin == 1u);

    BOOST_TEST(t.globalIndex(0, 0, 0) == 0u);
    BOOST_TEST(t.globalIndex(0, 0, 1) == 3u);
    BOOST_TEST(t.globalIndex(0, 0, 2) == 6u);
    BOOST_TEST(t.globalIndex(0, 1, 0) == 9u);
    BOOST_TEST(t.globalIndex(0, 1, 1) == 12u);
    BOOST_TEST(t.globalIndex(0, 1, 2) == 15u);
    BOOST_TEST(t.globalIndex(0, 2, 0) == 18u);
    BOOST_TEST(t.globalIndex(0, 2, 1) == 21u);
    BOOST_TEST(t.globalIndex(0, 2, 2) == 24u);
    BOOST_TEST(t.globalIndex(1, 0, 0) == 27u);
    BOOST_TEST(t.globalIndex(1, 0, 1) == 30u);
    BOOST_TEST(t.globalIndex(1, 0, 2) == 33u);
    BOOST_TEST(t.globalIndex(1, 1, 0) == 36u);
    BOOST_TEST(t.globalIndex(1, 1, 1) == 39u);
    BOOST_TEST(t.globalIndex(1, 1, 2) == 42u);
    BOOST_TEST(t.globalIndex(1, 2, 0) == 45u);
    BOOST_TEST(t.globalIndex(1, 2, 1) == 48u);
    BOOST_TEST(t.globalIndex(1, 2, 2) == 51u);
    BOOST_TEST(t.globalIndex(2, 0, 0) == 54u);
    BOOST_TEST(t.globalIndex(2, 0, 1) == 57u);
    BOOST_TEST(t.globalIndex(2, 0, 2) == 60u);
    BOOST_TEST(t.globalIndex(2, 1, 0) == 63u);
    BOOST_TEST(t.globalIndex(2, 1, 1) == 66u);
    BOOST_TEST(t.globalIndex(2, 1, 2) == 69u);
    BOOST_TEST(t.globalIndex(2, 2, 0) == 72u);
    BOOST_TEST(t.globalIndex(2, 2, 1) == 75u);
    BOOST_TEST(t.globalIndex(2, 2, 2) == 78u);

    // clang-format off
    BOOST_TEST((t.BField->getField(Vector3D(0, 0, 0)) / units::_T).isApprox(Vector3D(1, 0, 0)));
    BOOST_TEST((t.BField->getField(Vector3D(0, 0, 1)) / units::_T).isApprox(Vector3D(1, 1, 0)));
    BOOST_TEST((t.BField->getField(Vector3D(0, 1, 0)) / units::_T).isApprox(Vector3D(2, 1, 1)));
    BOOST_TEST((t.BField->getField(Vector3D(0, 1, 1)) / units::_T).isApprox(Vector3D(2, 2, 1)));
    BOOST_TEST((t.BField->getField(Vector3D(1, 0, 0)) / units::_T).isApprox(Vector3D(4, 3, 3)));
    BOOST_TEST((t.BField->getField(Vector3D(1, 0, 1)) / units::_T).isApprox(Vector3D(4, 4, 3)));
    BOOST_TEST((t.BField->getField(Vector3D(1, 1, 0)) / units::_T).isApprox(Vector3D(5, 4, 4)));
    BOOST_TEST((t.BField->getField(Vector3D(1, 1, 1)) / units::_T).isApprox(Vector3D(5, 5, 4)));
    BOOST_TEST((t.BField->getField(Vector3D(1.2, 1, 1)) / units::_T).isApprox(Vector3D(5.6, 5.6, 4.6)));
    BOOST_TEST((t.BField->getField(Vector3D(0, 1.6, 0)) / units::_T).isApprox(Vector3D(2.6, 1.6, 1.6)));
    BOOST_TEST((t.BField->getField(Vector3D(1, 0, 0.7)) / units::_T).isApprox(Vector3D(4, 3.7, 3)));
    // clang-format on
  }
}  // namespace Test

}  // namespace Acts
