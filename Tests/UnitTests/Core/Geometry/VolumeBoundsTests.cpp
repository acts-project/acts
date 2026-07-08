// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/ConeVolumeBounds.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/DiamondVolumeBounds.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"

#include <iomanip>
#include <memory>
#include <numbers>
#include <sstream>
#include <vector>

using namespace Acts;

namespace ActsTests {

namespace {

struct StreamState {
  std::ios_base::fmtflags flags;
  std::streamsize precision;
  std::streamsize width;
  char fill;
};

StreamState setNonDefaultStreamState(std::ostringstream& stream) {
  stream << std::scientific << std::showpos << std::setfill('#')
         << std::setprecision(3);
  stream.width(17);
  return {stream.flags(), stream.precision(), stream.width(), stream.fill()};
}

void checkStreamStatePreserved(const VolumeBounds& bounds) {
  std::ostringstream stream;
  const auto state = setNonDefaultStreamState(stream);

  bounds.toStream(stream);

  BOOST_CHECK(stream.flags() == state.flags);
  BOOST_CHECK_EQUAL(stream.precision(), state.precision);
  BOOST_CHECK_EQUAL(stream.width(), state.width);
  BOOST_CHECK_EQUAL(stream.fill(), state.fill);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(GeometrySuite)

BOOST_AUTO_TEST_CASE(VolumeBoundsTest) {
  // Tests if the planes are correctly oriented
  // s_planeXY
  // s_planeYZ
  // s_planeZX

  Vector3 xaxis(1., 0., 0.);
  Vector3 yaxis(0., 1., 0.);
  Vector3 zaxis(0., 0., 1.);

  auto rotXY = s_planeXY.rotation();
  BOOST_CHECK(rotXY.col(0).isApprox(xaxis));
  BOOST_CHECK(rotXY.col(1).isApprox(yaxis));
  BOOST_CHECK(rotXY.col(2).isApprox(zaxis));

  auto rotYZ = s_planeYZ.rotation();
  BOOST_CHECK(rotYZ.col(0).isApprox(yaxis));
  BOOST_CHECK(rotYZ.col(1).isApprox(zaxis));
  BOOST_CHECK(rotYZ.col(2).isApprox(xaxis));

  auto rotZX = s_planeZX.rotation();
  BOOST_CHECK(rotZX.col(0).isApprox(zaxis));
  BOOST_CHECK(rotZX.col(1).isApprox(xaxis));
  BOOST_CHECK(rotZX.col(2).isApprox(yaxis));
}

BOOST_AUTO_TEST_CASE(VolumeBoundsToStreamPreservesStreamState) {
  std::vector<std::unique_ptr<VolumeBounds>> bounds;
  bounds.push_back(std::make_unique<ConeVolumeBounds>(0., 0., 0.45, 50., 50.,
                                                      0., std::numbers::pi));
  bounds.push_back(std::make_unique<CuboidVolumeBounds>(10., 20., 30.));
  bounds.push_back(std::make_unique<CylinderVolumeBounds>(10., 20., 30.));
  bounds.push_back(
      std::make_unique<DiamondVolumeBounds>(10., 12., 8., 5., 10., 2.));
  bounds.push_back(std::make_unique<TrapezoidVolumeBounds>(5., 10., 8., 4.));

  for (const auto& bound : bounds) {
    checkStreamStatePreserved(*bound);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
