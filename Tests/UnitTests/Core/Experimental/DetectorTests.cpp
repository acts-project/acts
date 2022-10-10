// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Experimental/Detector.hpp"
#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/detail/NavigationStateUpdators.hpp"
#include "Acts/Experimental/detail/PortalGenerators.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <exception>
#include <memory>

using namespace Acts::Experimental;

/// Trial and error volume finder
///
/// @param gctx is the geometry context of this call
/// @param detector is the detector
/// @param position is the position
///
/// @return a detector volume pointer or null
const DetectorVolume* trialAndError(const Acts::GeometryContext& gctx,
                                    const Detector& detector,
                                    const Acts::Vector3& position) {
  for (const auto v : detector.volumes()) {
    if (v->inside(gctx, position)) {
      return v;
    }
  }
  return nullptr;
}

/// Unpack to shared - simply to test the getSharedPtr mechanism
///
/// @tparam referenced_type is the type of the referenced object
///
/// @param rt is the referenced object
///
/// @returns a shared pointer
template <typename referenced_type>
std::shared_ptr<referenced_type> unpackToShared(referenced_type& rt) {
  return rt.getSharedPtr();
}

Acts::GeometryContext tContext;

BOOST_AUTO_TEST_SUITE(Experimental)

BOOST_AUTO_TEST_CASE(DetectorConstruction) {
  Acts::ActsScalar r0 = 0.;
  Acts::ActsScalar r1 = 10.;
  Acts::ActsScalar r2 = 100.;
  Acts::ActsScalar r3 = 200.;
  Acts::ActsScalar zHalfL = 200.;

  Acts::Transform3 nominal = Acts::Transform3::Identity();

  // Create a bunch of volumes
  auto cyl0Bounds =
      std::make_unique<Acts::CylinderVolumeBounds>(r0, r1, zHalfL);

  auto cyl0BoundsCopy =
      std::make_unique<Acts::CylinderVolumeBounds>(r0, r1, zHalfL);

  auto cyl1Bounds =
      std::make_unique<Acts::CylinderVolumeBounds>(r1, r2, zHalfL);

  auto cyl2Bounds =
      std::make_unique<Acts::CylinderVolumeBounds>(r2, r3, zHalfL);

  auto portalGenerator = detail::defaultPortalGenerator();

  auto cyl0 = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Cyl0", nominal, std::move(cyl0Bounds),
      detail::defaultPortalProvider());

  auto cyl0nameDup = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Cyl0", nominal, std::move(cyl0BoundsCopy),
      detail::defaultPortalProvider());

  auto cyl1 = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Cyl1", nominal, std::move(cyl1Bounds),
      detail::defaultPortalProvider());

  auto cyl2 = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Cyl2", nominal, std::move(cyl2Bounds),
      detail::defaultPortalProvider());

  // A detector construction that should work
  DetectorVolumeFinder trialAndErrorFinder;
  trialAndErrorFinder.connect<&trialAndError>();

  std::vector<std::shared_ptr<DetectorVolume>> volumes012 = {cyl0, cyl1, cyl2};
  auto det012 = Detector::makeShared("Det012", volumes012, trialAndErrorFinder);

  // Check the basic return functions
  BOOST_CHECK(det012->name() == "Det012");
  BOOST_CHECK(det012->volumes().size() == 3u);
  BOOST_CHECK(det012->volumePtrs().size() == 3u);

  // Check the shared pointer mechanism
  BOOST_CHECK(det012 == unpackToShared<Detector>(*det012));
  BOOST_CHECK(det012 == unpackToShared<const Detector>(*det012));

  // Check the inside function with positions
  auto find0 = det012->findVolume(tContext, Acts::Vector3(5., 0., 0.));
  BOOST_CHECK(find0 == cyl0.get());

  auto find1 = det012->findVolume(tContext, Acts::Vector3(15., 0., 0.));
  BOOST_CHECK(find1 == cyl1.get());

  auto find2 = det012->findVolume(tContext, Acts::Vector3(150., 0., 0.));
  BOOST_CHECK(find2 == cyl2.get());

  auto findNull = det012->findVolume(tContext, Acts::Vector3(1500., 0., 0.));
  BOOST_CHECK(findNull == nullptr);

  /// Find by name
  find0 = det012->findVolume("Cyl0");
  BOOST_CHECK(find0 == cyl0.get());

  findNull = det012->findVolume("Null");
  BOOST_CHECK(findNull == nullptr);

  // Misconfigured - unkonnected finder
  DetectorVolumeFinder unconnected;
  BOOST_CHECK_THROW(
      Detector::makeShared("Det012_unconnected", volumes012, unconnected),
      std::invalid_argument);

  // Misconfigured - duplicate name
  std::vector<std::shared_ptr<DetectorVolume>> volumes002 = {cyl0, cyl0nameDup,
                                                             cyl2};
  BOOST_CHECK_THROW(Detector::makeShared("Det002_name_duplicate", volumes002,
                                         trialAndErrorFinder),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
