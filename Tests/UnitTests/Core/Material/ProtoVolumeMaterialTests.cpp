// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/ProtoVolumeMaterial.hpp"
#include "ActsTests/CommonHelpers/PredefinedMaterials.hpp"

#include <array>
#include <sstream>
#include <utility>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(MaterialSuite)

/// Test the constructors
BOOST_AUTO_TEST_CASE(ProtoVolumeMaterial_construction_test) {
  const MultiAxisSpec3D binning(
      std::array{AxisSpec::DeferredEquidistant(10, AxisDirection::AxisR),
                 AxisSpec::DeferredEquidistant(20, AxisDirection::AxisPhi),
                 AxisSpec::DeferredVariable({0., 0.1, 1.}, std::nullopt,
                                            AxisDirection::AxisZ)});

  // Constructor from arguments
  ProtoVolumeMaterial vmp(binning);
  // Copy constructor
  ProtoVolumeMaterial vmpCopy(vmp);
  // Copy move constructor
  ProtoVolumeMaterial vmpCopyMoved(std::move(vmpCopy));
  BOOST_CHECK(vmpCopyMoved.binning() == binning);
  ProtoVolumeMaterial assigned;
  assigned = vmpCopyMoved;
  BOOST_CHECK(assigned.binning() == binning);
  assigned = ProtoVolumeMaterial(binning);
  BOOST_CHECK(assigned.binning() == binning);

  // The consumer supplies the missing geometry properties without changing
  // the bin counts or the normalized variable edges in the proxy.
  const MultiAxisSpec3D::Options options{
      {{.min = 1., .max = 11., .boundaryType = AxisBoundaryType::Bound},
       {.min = -3., .max = 3., .boundaryType = AxisBoundaryType::Closed},
       {.min = -10., .max = 10., .boundaryType = AxisBoundaryType::Bound}}};
  const auto axes = assigned.binning().buildMultiAxis(options);
  BOOST_CHECK_EQUAL(axes->getAxis(0).getNBins(), 10u);
  BOOST_CHECK_EQUAL(axes->getAxis(1).getNBins(), 20u);
  BOOST_CHECK_EQUAL(axes->getAxis(2).getNBins(), 2u);
  BOOST_CHECK_EQUAL(axes->getAxis(2).getBinUpperBound(1), -8.);
  BOOST_CHECK(assigned.binning().isDeferred());
}

BOOST_AUTO_TEST_CASE(ProtoVolumeMaterialIsVacuumAndDefaultsToHomogeneous) {
  const ProtoVolumeMaterial proxy;
  BOOST_CHECK_EQUAL(proxy.binning().size(), 3u);
  for (const auto& axis : proxy.binning().axisSpecs()) {
    BOOST_CHECK_EQUAL(axis.nBins(), 1u);
    BOOST_CHECK(axis.isDeferred());
  }
  IVolumeMaterial::Cache cache;
  const HomogeneousVolumeMaterial homogeneous(makeSilicon());
  homogeneous.material(Vector3::Zero(), cache);
  BOOST_CHECK(proxy.material(Vector3::Zero(), cache).isVacuum());
  BOOST_CHECK(proxy.material(Vector3(100., 200., -300.), cache).isVacuum());
  BOOST_REQUIRE(cache.cell.has_value());
  BOOST_CHECK(cache.cell->isInside(Vector3(-100., -200., 300.)));
  BOOST_CHECK(proxy.material(Vector3::Zero()).isVacuum());
  std::ostringstream stream;
  stream << proxy;
  BOOST_CHECK(stream.str().find("ProtoVolumeMaterial") != std::string::npos);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
