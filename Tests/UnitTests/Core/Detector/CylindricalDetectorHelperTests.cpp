// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Detector/detail/CylindricalDetectorHelper.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CutoutCylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <map>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Acts {
namespace Experimental {
class Portal;
}  // namespace Experimental
}  // namespace Acts

using namespace Acts::Experimental;
using namespace Acts::Experimental::detail;
using namespace Acts::Experimental::detail::CylindricalDetectorHelper;

Acts::Logging::Level logLevel = Acts::Logging::VERBOSE;

Acts::GeometryContext tContext;
std::vector<std::shared_ptr<DetectorVolume>> eVolumes = {};

auto portalGenerator = defaultPortalGenerator();

BOOST_AUTO_TEST_SUITE(Experimental)

BOOST_AUTO_TEST_CASE(ConnectVolumeExceptions) {
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("Faulty setups", logLevel));

  auto cBounds0 = std::make_unique<Acts::CylinderVolumeBounds>(0., 100., 100);
  auto volume0 = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Volume0", Acts::Transform3::Identity(),
      std::move(cBounds0), tryAllPortals());

  auto cBounds1 =
      std::make_unique<Acts::CylinderVolumeBounds>(0., 100., 100, 0.2, 1.);
  auto volume1 = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Volume0", Acts::Transform3::Identity(),
      std::move(cBounds1), tryAllPortals());

  ACTS_INFO("*** Test: nullptr in the list of volumes");

  // Invalid arguments: nullptr
  std::vector<std::shared_ptr<DetectorVolume>> volumesWithNullptr = {
      volume0, nullptr, volume1};
  BOOST_CHECK_THROW(connectInR(tContext, volumesWithNullptr, {}, logLevel),
                    std::invalid_argument);

  ACTS_INFO("*** Test: non-cylinder in the list of volumes");

  auto cubeBounds = std::make_unique<Acts::CuboidVolumeBounds>(100., 100., 100);
  auto cube = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Cube", Acts::Transform3::Identity(),
      std::move(cubeBounds), tryAllPortals());

  // Invalid arguments: cube
  std::vector<std::shared_ptr<DetectorVolume>> volumesWithCube = {
      volume0, volume1, cube};
  BOOST_CHECK_THROW(connectInR(tContext, volumesWithCube, {}, logLevel),
                    std::invalid_argument);

  ACTS_INFO("*** Test: non-aligned volume in the list of volumes");
  Acts::Transform3 rotated = Acts::Transform3::Identity();
  Acts::AngleAxis3 rotX(0.1234, Acts::Vector3::UnitX());
  rotated *= rotX;

  auto cBounds2 = std::make_unique<Acts::CylinderVolumeBounds>(0., 100., 100);
  auto volume2 = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Volume2", rotated, std::move(cBounds2),
      tryAllPortals());

  // Invalid arguments: non-aligned
  std::vector<std::shared_ptr<DetectorVolume>> volumesWithNonaligned = {
      volume0, volume1, volume2};
  BOOST_CHECK_THROW(connectInR(tContext, volumesWithNonaligned, {}, logLevel),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(ConnectInR) {
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("Connect: R", logLevel));
  ACTS_INFO("*** Test: connect DetectorVolumes in R, create proto container");
  // Test with different opening angles
  std::vector<Acts::ActsScalar> testOpenings = {M_PI, 0.5 * M_PI};

  std::vector<Acts::ActsScalar> radii = {0., 10., 100., 200.};
  Acts::ActsScalar halfZ = 100.;

  // This should work for full cylinder and sector openings
  for (auto [io, opening] : Acts::enumerate(testOpenings)) {
    ACTS_INFO("    -> test  with phi openeing: " << opening);
    std::string opStr = "opening_" + std::to_string(io);
    std::vector<std::shared_ptr<DetectorVolume>> rVolumes = {};
    // Create the voluems
    for (auto [i, r] : Acts::enumerate(radii)) {
      if (i > 0) {
        auto cBounds = std::make_unique<Acts::CylinderVolumeBounds>(
            radii[i - 1u], r, halfZ, opening, 0.);
        rVolumes.push_back(DetectorVolumeFactory::construct(
            portalGenerator, tContext, "Cylinder_r" + std::to_string(i),
            Acts::Transform3::Identity(), std::move(cBounds), tryAllPortals()));
      }
    }

    auto protoContainer = connectInR(tContext, rVolumes, {}, logLevel);
    // Check the portal setup
    BOOST_CHECK(rVolumes[0u]->portalPtrs()[2u] ==
                rVolumes[1u]->portalPtrs()[3u]);
    BOOST_CHECK(rVolumes[1u]->portalPtrs()[2u] ==
                rVolumes[2u]->portalPtrs()[3u]);
    BOOST_CHECK(rVolumes[0u]->portalPtrs()[0u] ==
                rVolumes[1u]->portalPtrs()[0u]);
    BOOST_CHECK(rVolumes[1u]->portalPtrs()[0u] ==
                rVolumes[2u]->portalPtrs()[0u]);
    BOOST_CHECK(rVolumes[0u]->portalPtrs()[1u] ==
                rVolumes[1u]->portalPtrs()[1u]);
    BOOST_CHECK(rVolumes[1u]->portalPtrs()[1u] ==
                rVolumes[2u]->portalPtrs()[1u]);
    BOOST_CHECK(rVolumes[0u]->portalPtrs()[0u] == protoContainer[0u]);
    BOOST_CHECK(rVolumes[0u]->portalPtrs()[1u] == protoContainer[1u]);

    // A detector construction that should work
    auto detector =
        Detector::makeShared("DetectorInR", rVolumes, tryRootVolumes());

    // Make a rzphi grid
    const auto& volumes = detector->volumes();
    auto boundaries = rzphiBoundaries(tContext, volumes);
    const auto& rBoundaries = boundaries[0u];
    const auto& zBoundaries = boundaries[1u];

    // Check the radii
    std::vector<Acts::ActsScalar> zvalues = {-halfZ, halfZ};
    BOOST_CHECK(radii == rBoundaries);
    BOOST_CHECK(zvalues == zBoundaries);
  }

  // Invalid arguments
  ACTS_INFO("*** Test: faulty empty vector");
  BOOST_CHECK_THROW(connectInR(tContext, eVolumes, {}, logLevel),
                    std::invalid_argument);

  // Faulty setups, not matchint in R
  ACTS_INFO("*** Test: volumes are not matching in R");

  auto cBounds00 = std::make_unique<Acts::CylinderVolumeBounds>(0., 100., 100);
  auto volume00 = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Volume00", Acts::Transform3::Identity(),
      std::move(cBounds00), tryAllPortals());

  auto cBounds01 =
      std::make_unique<Acts::CylinderVolumeBounds>(101., 200., 100);
  auto volume01 = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Volume01", Acts::Transform3::Identity(),
      std::move(cBounds01), tryAllPortals());

  std::vector<std::shared_ptr<DetectorVolume>> volumesNotMatching = {volume00,
                                                                     volume01};
  BOOST_CHECK_THROW(connectInR(tContext, volumesNotMatching, {}, logLevel),
                    std::runtime_error);

  ACTS_INFO("*** Test: volume bounds are not aligned");
  Acts::Transform3 shifted = Acts::Transform3::Identity();
  shifted.pretranslate(Acts::Vector3(0., 0., 10.));

  auto cBounds10 = std::make_unique<Acts::CylinderVolumeBounds>(0., 100., 100);
  auto volume10 = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Volume10", shifted, std::move(cBounds10),
      tryAllPortals());

  auto cBounds11 = std::make_unique<Acts::CylinderVolumeBounds>(100., 200., 90);
  auto volume11 = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Volume11", shifted, std::move(cBounds11),
      tryAllPortals());

  std::vector<std::shared_ptr<DetectorVolume>> volumesNotAligned = {volume10,
                                                                    volume11};
  BOOST_CHECK_THROW(connectInR(tContext, volumesNotAligned, {}, logLevel),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ConnectInZ) {
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("Connect: Z", logLevel));
  ACTS_INFO("*** Test: connect DetectorVolumes in Z, create proto container");

  // @TODO: test with different transforms, this should work in, not used yet
  std::vector<Acts::Transform3> transforms = {Acts::Transform3::Identity()};
  std::vector<std::array<Acts::ActsScalar, 2>> radii = {{0., 100.},
                                                        {20., 120.}};
  std::vector<Acts::ActsScalar> zValues = {-100., -20, 10., 100., 200.};

  for (auto [it, t] : Acts::enumerate(transforms)) {
    ACTS_INFO("    -> test series with transfrom id " << it);

    std::string trfStr = "_transform_" + std::to_string(it);
    for (auto [ir, r] : Acts::enumerate(radii)) {
      ACTS_INFO("        -> test series with radii setup "
                << radii[ir][0u] << ", " << radii[ir][1u]);

      std::string radStr = "_radii_" + std::to_string(ir);
      std::vector<std::shared_ptr<DetectorVolume>> zVolumes = {};
      for (auto [i, z] : Acts::enumerate(zValues)) {
        if (i > 0) {
          auto cBounds = std::make_unique<Acts::CylinderVolumeBounds>(
              r[0], r[1], 0.5 * (z - zValues[i - 1u]));
          // z center
          Acts::ActsScalar zCenter = 0.5 * (z + zValues[i - 1u]);
          Acts::Transform3 ti = Acts::Transform3::Identity();
          ti.pretranslate(t.translation() +
                          zCenter * t.rotation().matrix().col(2));
          ti.prerotate(t.rotation());
          // create the volume
          zVolumes.push_back(DetectorVolumeFactory::construct(
              portalGenerator, tContext,
              "Cylinder_z" + std::to_string(i) + trfStr + radStr, ti,
              std::move(cBounds), tryAllPortals()));
        }
      }
      // Now call the connector
      auto protoContainer = connectInZ(tContext, zVolumes, {}, logLevel);

      // Check the portal setup.
      // Glued, remainders are outside skin
      BOOST_CHECK(zVolumes[0u]->portalPtrs()[1u] ==
                  zVolumes[1u]->portalPtrs()[0u]);
      BOOST_CHECK(zVolumes[1u]->portalPtrs()[1u] ==
                  zVolumes[2u]->portalPtrs()[0u]);
      BOOST_CHECK(zVolumes[2u]->portalPtrs()[1u] ==
                  zVolumes[3u]->portalPtrs()[0u]);
      BOOST_CHECK(protoContainer[0u] == zVolumes[0u]->portalPtrs()[0u]);
      BOOST_CHECK(protoContainer[1u] == zVolumes[3u]->portalPtrs()[1u]);

      // Covered with the same surface, shich is the outside skin
      std::vector<unsigned int> checkShared = {2u};
      if (radii[ir][0u] > 0.) {
        checkShared.push_back(3u);
      }

      for (const auto& ip : checkShared) {
        BOOST_CHECK(zVolumes[0u]->portalPtrs()[ip] ==
                    zVolumes[1u]->portalPtrs()[ip]);
        BOOST_CHECK(zVolumes[1u]->portalPtrs()[ip] ==
                    zVolumes[2u]->portalPtrs()[ip]);
        BOOST_CHECK(zVolumes[2u]->portalPtrs()[ip] ==
                    zVolumes[3u]->portalPtrs()[ip]);
        BOOST_CHECK(protoContainer[ip] == zVolumes[0u]->portalPtrs()[ip]);
      }

      auto detector =
          Detector::makeShared("DetectorInZ", zVolumes, tryRootVolumes());
    }
  }

  // Invalid arguments
  BOOST_CHECK_THROW(connectInZ(tContext, eVolumes, {}, logLevel),
                    std::invalid_argument);

  // Volumes have different radii - other bounds will be the same
  auto cBounds00 = std::make_unique<Acts::CylinderVolumeBounds>(0., 100., 100);
  auto volume00 = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Volume00",
      Acts::Transform3::Identity() * Acts::Translation3(0., 0., -100.),
      std::move(cBounds00), tryAllPortals());

  auto cBounds01 = std::make_unique<Acts::CylinderVolumeBounds>(0., 105., 100);
  auto volume01 = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Volume01",
      Acts::Transform3::Identity() * Acts::Translation3(0., 0., 100.),
      std::move(cBounds01), tryAllPortals());

  std::vector<std::shared_ptr<DetectorVolume>> volumesNonalignedBounds = {
      volume00, volume01};
  BOOST_CHECK_THROW(connectInZ(tContext, volumesNonalignedBounds, {}, logLevel),
                    std::runtime_error);

  // Volumes are not attached
  auto cBounds10 = std::make_unique<Acts::CylinderVolumeBounds>(0., 100., 100);
  auto volume10 = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Volume00",
      Acts::Transform3::Identity() * Acts::Translation3(0., 0., -105.),
      std::move(cBounds10), tryAllPortals());

  auto cBounds11 = std::make_unique<Acts::CylinderVolumeBounds>(0., 100., 100);
  auto volume11 = DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Volume01",
      Acts::Transform3::Identity() * Acts::Translation3(0., 0., 100.),
      std::move(cBounds11), tryAllPortals());

  std::vector<std::shared_ptr<DetectorVolume>> volumesNotAttached = {volume10,
                                                                     volume11};
  BOOST_CHECK_THROW(connectInZ(tContext, volumesNotAttached, {}, logLevel),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ConnectInPhi) {
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("Connect: Phi", logLevel));
  ACTS_INFO("*** Test: connect DetectorVolumes in Phi, create proto container");

  std::vector<Acts::Transform3> transforms = {Acts::Transform3::Identity()};
  unsigned int phiSectors = 5;
  Acts::ActsScalar phiHalfSector = M_PI / phiSectors;

  for (auto [it, t] : Acts::enumerate(transforms)) {
    ACTS_INFO("    -> test series with transfrom id " << it);

    std::vector<std::shared_ptr<DetectorVolume>> phiVolumes = {};
    for (unsigned int i = 0; i < phiSectors; ++i) {
      auto cBounds = std::make_unique<Acts::CylinderVolumeBounds>(
          10., 100., 100., phiHalfSector,
          -M_PI + (2u * i + 1u) * phiHalfSector);

      // create the volume
      phiVolumes.push_back(DetectorVolumeFactory::construct(
          portalGenerator, tContext, "Cylinder_phi" + std::to_string(i), t,
          std::move(cBounds), tryAllPortals()));
    }

    auto protoContainer = connectInPhi(tContext, phiVolumes, {}, logLevel);

    // All phiVolumes share : inner tube, outer cover, negative & positive disc
    std::vector<unsigned int> checkShared = {0u, 1u, 2u, 3u};
    for (auto [iv, v] : Acts::enumerate(phiVolumes)) {
      if (iv > 0u) {
        auto current = v;
        auto last = phiVolumes[iv - 1u];
        for (const auto& ch : checkShared) {
          BOOST_CHECK(current->portalPtrs()[ch] == last->portalPtrs()[ch]);
        }
      }
    }

    auto detector =
        Detector::makeShared("DetectorInPhi", phiVolumes, tryRootVolumes());
  }

  // Invalid arguments
  BOOST_CHECK_THROW(connectInPhi(tContext, eVolumes, {}, logLevel),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(WrapVolumeinRZ) {
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("Wrap: Z-R", logLevel));
  ACTS_INFO(
      "*** Test: wrap volume in Z-R with CutoutCylinderVolume, create proto "
      "container");

  // @TODO: test with different transforms, this should work in, not used yet
  std::vector<Acts::Transform3> transforms = {Acts::Transform3::Identity()};

  // Test with different inner radii
  std::vector<std::array<Acts::ActsScalar, 3u>> radii = {{0., 100., 500.},
                                                         {20., 120., 500.}};

  Acts::ActsScalar innerHalfZ = 150.;
  Acts::ActsScalar outerHalfZ = 175.;

  // Set up all the different tests
  for (auto [it, tf] : Acts::enumerate(transforms)) {
    ACTS_INFO("    Test series with transfrom id " << it);

    std::string trfStr = "_transform_" + std::to_string(it);
    for (auto [ir, r] : Acts::enumerate(radii)) {
      ACTS_INFO("    -> test series with radii setup " << radii[ir][0u] << ", "
                                                       << radii[ir][1u]);

      std::vector<std::shared_ptr<DetectorVolume>> volumes = {};

      std::string radStr = "_radii_" + std::to_string(ir);
      // Create the inner bounds
      auto iBounds = std::make_unique<Acts::CylinderVolumeBounds>(
          radii[ir][0u], radii[ir][1u], innerHalfZ);
      volumes.push_back(DetectorVolumeFactory::construct(
          portalGenerator, tContext, "InnerCylinder" + radStr + trfStr, tf,
          std::move(iBounds), tryAllPortals()));

      // Create the wrapping bounds
      auto wBounds = std::make_unique<Acts::CutoutCylinderVolumeBounds>(
          radii[ir][0u], radii[ir][1u], radii[ir][2u], outerHalfZ, innerHalfZ);

      volumes.push_back(DetectorVolumeFactory::construct(
          portalGenerator, tContext, "WrappingCylinder" + radStr + trfStr, tf,
          std::move(wBounds), tryAllPortals()));

      wrapInZR(tContext, volumes, logLevel);
    }
  }

  // Invalid arguments
  BOOST_CHECK_THROW(wrapInZR(tContext, eVolumes, logLevel),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(ProtoContainerZR) {
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("Container: Z-R", logLevel));
  ACTS_INFO("*** Test: create a container in Z-R.");

  auto transform = Acts::Transform3::Identity();

  std::vector<Acts::ActsScalar> innerMostRadii = {0., 2.};

  for (auto [ir, imr] : Acts::enumerate(innerMostRadii)) {
    ACTS_INFO("    -> test series innermost radius setup "
              << innerMostRadii[ir]);

    // A container in R
    std::vector<Acts::ActsScalar> radii = {25., 100., 200.};
    Acts::ActsScalar halfZ = 200;

    // An innermost Pipe
    auto bBounds =
        std::make_unique<Acts::CylinderVolumeBounds>(imr, radii[0u], halfZ);

    auto innerPipe = DetectorVolumeFactory::construct(
        portalGenerator, tContext, "InnerPipe", transform, std::move(bBounds),
        tryAllPortals());

    // Make a container representation out of it
    std::map<unsigned int, std::shared_ptr<Portal>> ipContainer;
    for (auto [ip, p] : Acts::enumerate(innerPipe->portalPtrs())) {
      ipContainer[ip] = p;
    }

    // Create the r - sorted volumes
    std::vector<std::shared_ptr<DetectorVolume>> rVolumes = {};
    // Create the voluems
    for (auto [i, r] : Acts::enumerate(radii)) {
      if (i > 0) {
        auto cBounds = std::make_unique<Acts::CylinderVolumeBounds>(
            radii[i - 1u], r, halfZ);
        rVolumes.push_back(DetectorVolumeFactory::construct(
            portalGenerator, tContext, "Cylinder_r" + std::to_string(i),
            transform, std::move(cBounds), tryAllPortals()));
      }
    }

    auto protoContainerInR = connectInR(tContext, rVolumes, {}, logLevel);

    std::vector<Acts::ActsScalar> zValues = {-200., -120, 10., 100., 200.};
    std::vector<std::shared_ptr<DetectorVolume>> zVolumes = {};
    for (auto [i, z] : Acts::enumerate(zValues)) {
      if (i > 0) {
        auto cBounds = std::make_unique<Acts::CylinderVolumeBounds>(
            200., 300., 0.5 * (z - zValues[i - 1u]));
        // z center
        Acts::ActsScalar zCenter = 0.5 * (z + zValues[i - 1u]);
        Acts::Transform3 ti = transform;
        ti.pretranslate(transform.translation() +
                        zCenter * transform.rotation().matrix().col(2));

        // create the volume
        zVolumes.push_back(DetectorVolumeFactory::construct(
            portalGenerator, tContext, "Cylinder_z" + std::to_string(i), ti,
            std::move(cBounds), tryAllPortals()));
      }
    }
    // Now call the connector
    auto protoContainerInZ = connectInZ(tContext, zVolumes, {}, logLevel);
    auto centralContainer = connectInR(
        tContext, {ipContainer, protoContainerInR, protoContainerInZ}, {},
        logLevel);

    // Let's make two endcaps
    // Nec
    auto necBounds =
        std::make_unique<Acts::CylinderVolumeBounds>(imr, 300., 50.);

    auto necTransform = Acts::Transform3::Identity();
    necTransform.pretranslate(Acts::Vector3(0., 0., -250));
    auto necVolume = DetectorVolumeFactory::construct(
        portalGenerator, tContext, "Nec", necTransform, std::move(necBounds),
        tryAllPortals());

    std::map<unsigned int, std::shared_ptr<Portal>> necContainer;
    for (auto [ip, p] : Acts::enumerate(necVolume->portalPtrs())) {
      necContainer[ip] = p;
    }

    // Pec container
    auto pecInnerBounds =
        std::make_unique<Acts::CylinderVolumeBounds>(imr, 175., 100.);

    auto pecOuterBounds =
        std::make_unique<Acts::CylinderVolumeBounds>(175., 300., 100.);

    auto pecTransform = Acts::Transform3::Identity();
    pecTransform.pretranslate(Acts::Vector3(0., 0., 300));
    auto pecInner = DetectorVolumeFactory::construct(
        portalGenerator, tContext, "PecInner", pecTransform,
        std::move(pecInnerBounds), tryAllPortals());
    auto pecOuter = DetectorVolumeFactory::construct(
        portalGenerator, tContext, "PecOuter", pecTransform,
        std::move(pecOuterBounds), tryAllPortals());

    std::vector<std::shared_ptr<DetectorVolume>> pecVolumes = {pecInner,
                                                               pecOuter};
    auto pecContainer = connectInR(tContext, pecVolumes, {}, logLevel);

    auto overallContainer = connectInZ(
        tContext, {necContainer, centralContainer, pecContainer}, {}, logLevel);

    //  Add them togeter
    std::vector<std::shared_ptr<DetectorVolume>> dVolumes;
    dVolumes.push_back(innerPipe);
    dVolumes.push_back(necVolume);
    dVolumes.insert(dVolumes.end(), rVolumes.begin(), rVolumes.end());
    dVolumes.insert(dVolumes.end(), zVolumes.begin(), zVolumes.end());
    dVolumes.push_back(pecInner);
    dVolumes.push_back(pecOuter);

    auto detector = Detector::makeShared("DetectorFromProtoContainer", dVolumes,
                                         tryRootVolumes());
  }  // test with different innermost radii
}

BOOST_AUTO_TEST_CASE(WrapContainernRZ) {
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("Container: Wrap", logLevel));
  ACTS_INFO("*** Test: create a container in Z-R by wrapping.");

  // Test with different inner radii
  std::vector<std::array<Acts::ActsScalar, 3u>> radii = {{0., 100., 500.},
                                                         {20., 120., 500.}};

  Acts::ActsScalar innerHalfZ = 150.;
  Acts::ActsScalar innerBarrelHalfZ = 75.;
  Acts::ActsScalar innerEndcapHalfZ = 0.5 * (innerHalfZ - innerBarrelHalfZ);
  Acts::ActsScalar outerHalfZ = 175.;

  Acts::Transform3 tf = Acts::Transform3::Identity();

  // Set up all the different tests
  for (auto [ir, r] : Acts::enumerate(radii)) {
    std::string radStr = "_radii_" + std::to_string(ir);
    ACTS_INFO("    -> test series innermost radius setup " << radii[ir][0u]);

    // Let's create the inner container first
    std::vector<std::shared_ptr<DetectorVolume>> iVolumes = {};

    auto iNecBounds = std::make_unique<Acts::CylinderVolumeBounds>(
        radii[ir][0u], radii[ir][1u], innerEndcapHalfZ);
    Acts::Transform3 ntf = tf;
    ntf.pretranslate(
        Acts::Vector3(0., 0., -innerBarrelHalfZ - innerEndcapHalfZ));
    iVolumes.push_back(DetectorVolumeFactory::construct(
        portalGenerator, tContext, "InnerNec" + radStr, ntf,
        std::move(iNecBounds), tryAllPortals()));

    auto iBarrelBounds = std::make_unique<Acts::CylinderVolumeBounds>(
        radii[ir][0u], radii[ir][1u], innerBarrelHalfZ);
    iVolumes.push_back(DetectorVolumeFactory::construct(
        portalGenerator, tContext, "InnerBarrel" + radStr, tf,
        std::move(iBarrelBounds), tryAllPortals()));

    auto iPecBounds = std::make_unique<Acts::CylinderVolumeBounds>(
        radii[ir][0u], radii[ir][1u], innerEndcapHalfZ);
    Acts::Transform3 ptf = tf;
    ptf.pretranslate(
        Acts::Vector3(0., 0., innerBarrelHalfZ + innerEndcapHalfZ));
    iVolumes.push_back(DetectorVolumeFactory::construct(
        portalGenerator, tContext, "InnerPec" + radStr, ptf,
        std::move(iPecBounds), tryAllPortals()));

    auto innerContainer = connectInZ(tContext, iVolumes, {}, logLevel);

    // Create the wrapping volume
    auto wBounds = std::make_unique<Acts::CutoutCylinderVolumeBounds>(
        radii[ir][0u], radii[ir][1u], radii[ir][2u], outerHalfZ, innerHalfZ);
    auto wVolume = DetectorVolumeFactory::construct(
        portalGenerator, tContext, "WrappingVolume" + radStr, tf,
        std::move(wBounds), tryAllPortals());

    std::vector<DetectorComponent::PortalContainer> containers;
    containers.push_back(innerContainer);

    DetectorComponent::PortalContainer outerContainer;
    for (auto [ip, p] : Acts::enumerate(wVolume->portalPtrs())) {
      outerContainer[ip] = p;
    }
    containers.push_back(outerContainer);

    auto detector = wrapInZR(tContext, containers, logLevel);
  }
}

BOOST_AUTO_TEST_SUITE_END()
