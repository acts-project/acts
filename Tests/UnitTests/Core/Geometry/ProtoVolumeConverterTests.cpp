// This file is part of the Acts project.
 //
 // Copyright (C) 2020 CERN for the benefit of the Acts project
 //
 // This Source Code Form is subject to the terms of the Mozilla Public
 // License, v. 2.0. If a copy of the MPL was not distributed with this
 // file, You can obtain one at http://mozilla.org/MPL/2.0/.

 #include <boost/test/data/test_case.hpp>
 #include <boost/test/unit_test.hpp>

 #include "Acts/Geometry/GeometryContext.hpp"
 #include "Acts/Geometry/ProtoDetector.hpp"
 #include "Acts/Geometry/ProtoVolumeConverter.hpp"
 #include "Acts/Geometry/VolumeBounds.hpp"
 #include "Acts/Geometry/detail/PortalGenerators.hpp"
 #include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
 #include "Acts/Utilities/BinningData.hpp"

 #include <cmath>
 #include <sstream>

 namespace Acts {

 // Create a test context
 GeometryContext tContext;

 using Values = CylinderVolumeBounds::BoundValues;

 namespace Experimental {

 namespace Test {

 BOOST_AUTO_TEST_SUITE(Geometry)

 BOOST_AUTO_TEST_CASE(CylinderVolumeShapeConversion) {
   ProtoVolume cylinder;
   cylinder.name = "cylinder";
   cylinder.extent.set(binR, 0., 200.);
   cylinder.extent.set(binZ, -300., 300.);

   const auto& [transform0, bounds0] =
       ConcentricCylinderConverter{cylinder}.create(tContext);

   CHECK_CLOSE_ABS(bounds0->get(Values::eMinR), 0., 1e-5);
   CHECK_CLOSE_ABS(bounds0->get(Values::eMaxR), 200., 1e-5);
   CHECK_CLOSE_ABS(bounds0->get(Values::eHalfLengthZ), 300., 1e-5);
   CHECK_CLOSE_ABS(bounds0->get(Values::eAveragePhi), 0., 1e-5);
   CHECK_CLOSE_ABS(bounds0->get(Values::eHalfPhiSector), M_PI, 1e-5);
   BOOST_CHECK(transform0.isApprox(Transform3::Identity()));

   ProtoVolume tube;
   tube.name = "tube";
   tube.extent.set(binR, 10., 200.);
   tube.extent.set(binZ, -400., 400.);

   auto [transform1, bounds1] =
       ConcentricCylinderConverter{tube}.create(tContext);

   CHECK_CLOSE_ABS(bounds1->get(Values::eMinR), 10., 1e-5);
   CHECK_CLOSE_ABS(bounds1->get(Values::eMaxR), 200., 1e-5);
   CHECK_CLOSE_ABS(bounds1->get(Values::eHalfLengthZ), 400., 1e-5);
   CHECK_CLOSE_ABS(bounds1->get(Values::eAveragePhi), 0., 1e-5);
   CHECK_CLOSE_ABS(bounds1->get(Values::eHalfPhiSector), M_PI, 1e-5);
   BOOST_CHECK(transform1.isApprox(Transform3::Identity()));

   ProtoVolume sectoralCylinder;
   sectoralCylinder.name = "sectoralCylinder";
   sectoralCylinder.extent.set(binR, 10., 200.);
   sectoralCylinder.extent.set(binZ, -400., 400.);
   sectoralCylinder.extent.set(binPhi, -0.75 * M_PI, -0.25 * M_PI);

   auto [transform2, bounds2] =
       ConcentricCylinderConverter{sectoralCylinder}.create(tContext);

   CHECK_CLOSE_ABS(bounds2->get(Values::eMinR), 10., 1e-5);
   CHECK_CLOSE_ABS(bounds2->get(Values::eMaxR), 200., 1e-5);
   CHECK_CLOSE_ABS(bounds2->get(Values::eHalfLengthZ), 400., 1e-5);
   ActsScalar refAvPhi = -0.5 * M_PI;
   CHECK_CLOSE_ABS(bounds2->get(Values::eAveragePhi), refAvPhi, 1e-5);
   CHECK_CLOSE_ABS(bounds2->get(Values::eHalfPhiSector), 0.25 * M_PI, 1e-5);
   BOOST_CHECK(transform2.isApprox(Transform3::Identity()));

   ProtoVolume shiftedCylinder;
   cylinder.name = "cylinder";
   cylinder.extent.set(binR, 0., 200.);
   cylinder.extent.set(binZ, -300., 300.);
 }

 BOOST_AUTO_TEST_CASE(EmptyInternalsTest) {
   ProtoVolume cylinder;
   cylinder.name = "cylinder";
   cylinder.extent.set(binR, 0., 200.);
   cylinder.extent.set(binZ, -300., 300.);

   auto [transform, bounds] =
       ConcentricCylinderConverter{cylinder}.create(tContext);

   auto [surfaces, volumes, candidateUpdator] =
       EmptyInternals{cylinder}.create(tContext);
   BOOST_CHECK(surfaces.empty());
   BOOST_CHECK(volumes.empty());
   BOOST_CHECK(candidateUpdator.connected());

   auto portalGenerator = detail::defaultPortalGenerator();

   // A cylinder - fully equipped
   auto cylinderVolume = DetectorVolumeFactory::construct(
       portalGenerator, tContext, "CylinderVolume", transform, std::move(bounds),
       std::move(candidateUpdator));
 }

 BOOST_AUTO_TEST_CASE(EmptyVolumeTest) {
   ProtoVolume cylinder;
   cylinder.name = "cylinder";
   cylinder.extent.set(binR, 0., 200.);
   cylinder.extent.set(binZ, -300., 300.);

   DetectorBlock dBlock;
   SingleBlockBuilder<>{cylinder}(dBlock, tContext);

   // One volume contained
   BOOST_CHECK(std::get<0>(dBlock).size() == 1u);
   // Three bounding surfaces define the shell/skin
   BOOST_CHECK(std::get<1>(dBlock).size() == 3u);
 }

 BOOST_AUTO_TEST_SUITE_END()

 }  // namespace Test

 }  // namespace Experimental
 }  // namespace Acts