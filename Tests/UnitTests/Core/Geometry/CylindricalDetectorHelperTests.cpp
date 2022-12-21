// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/CylindricalDetectorHelper.hpp"
#include "Acts/Geometry/Detector.hpp"
#include "Acts/Geometry/DetectorVolume.hpp"
#include "Acts/Geometry/detail/DetectorVolumeFinders.hpp"
#include "Acts/Geometry/detail/NavigationStateUpdators.hpp"
#include "Acts/Geometry/detail/PortalGenerators.hpp"
#include "Acts/Geometry/detail/SurfaceCandidatesUpdators.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <exception>
#include <memory>

using namespace Acts::Experimental;

Acts::GeometryContext tContext;

BOOST_AUTO_TEST_SUITE(Experimental)

BOOST_AUTO_TEST_CASE(ConnectPerfectlyFitInR) {
  // Test with different transforms
  std::vector<Acts::Transform3> transforms = {Acts::Transform3::Identity()};

  // Test with different opening angles
  std::vector<Acts::ActsScalar> testOpenings = {M_PI};

  std::vector<Acts::ActsScalar> radii = {0., 10., 100., 200.};
  Acts::ActsScalar halfZ = 100.;

  auto portalGenerator = detail::defaultPortalGenerator();
  auto navigationStateUpdator = detail::allPortals();

  // This should work for generic transforms
  for (auto [it, tf] : Acts::enumerate(transforms)) {
    std::string trfStr = "transform_" + std::to_string(it);
    // This should work for full cylinder and sector openings
    for (auto [io, opening] : Acts::enumerate(testOpenings)) {
      std::string opStr = "_opening_" + std::to_string(it);

      std::vector<std::shared_ptr<DetectorVolume>> rVolumes = {};
      // Create the voluems
      for (auto [i, r] : Acts::enumerate(radii)) {
        if (i > 0) {
          auto cBounds = std::make_unique<Acts::CylinderVolumeBounds>(
              radii[i - 1u], r, halfZ, opening, 0.);
          rVolumes.push_back(DetectorVolumeFactory::construct(
              portalGenerator, tContext, "Cylinder_r" + std::to_string(i), tf,
              std::move(cBounds), detail::allPortals()));
        }
      }

      connectDetectorVolumesInR(tContext, rVolumes);

      // A detector construction that should work
      auto detector = Detector::makeShared("DetectorInR", rVolumes,
                                           detail::tryAllVolumes());

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
  }
}

BOOST_AUTO_TEST_CASE(PerfectFitAttachmentZ) {
  std::vector<Acts::Transform3> transforms = {Acts::Transform3::Identity()};
  std::vector<std::array<Acts::ActsScalar, 2>> radii = {{0., 100.},
                                                        {20., 120.}};
  std::vector<Acts::ActsScalar> zValues = {-100., -20, 10., 100., 200.};

  auto portalGenerator = detail::defaultPortalGenerator();
  auto navigationStateUpdator = detail::allPortals();

  for (auto [it, t] : Acts::enumerate(transforms)) {
    std::string trfStr = "transform_" + std::to_string(it);
    for (auto [ir, r] : Acts::enumerate(radii)) {
      std::string radStr = "radii_" + std::to_string(ir);
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
              portalGenerator, tContext, "Cylinder_z" + std::to_string(i), ti,
              std::move(cBounds), detail::allPortals()));
        }
      }
      // Now call the connector
      connectDetectorVolumesInZ(tContext, zVolumes);

      auto detector = Detector::makeShared("DetectorInZ", zVolumes,
                                           detail::tryAllVolumes());
    }
  }
}

BOOST_AUTO_TEST_CASE(PerfectFitPhiAttachment) {
  std::vector<Acts::Transform3> transforms = {Acts::Transform3::Identity()};
  unsigned int phiSectors = 5;
  Acts::ActsScalar phiHalfSector = M_PI / phiSectors;

  auto portalGenerator = detail::defaultPortalGenerator();
  auto navigationStateUpdator = detail::allPortals();

  for (const auto& t : transforms) {
    std::vector<std::shared_ptr<DetectorVolume>> phiVolumes = {};
    for (unsigned int i = 0; i < phiSectors; ++i) {
      auto cBounds = std::make_unique<Acts::CylinderVolumeBounds>(
          10., 100., 100., phiHalfSector, 0.);
      // Rotate in phi center
      const Acts::Vector3 colZ = t.rotation().matrix().col(2);
      Acts::Transform3 ti =
          Acts::AngleAxis3(-M_PI + 2 * i * phiHalfSector, colZ) * t;

      // create the volume
      phiVolumes.push_back(DetectorVolumeFactory::construct(
          portalGenerator, tContext, "Cylinder_phi" + std::to_string(i), ti,
          std::move(cBounds), detail::allPortals()));
    }

    connectDetectorVolumesInPhi(tContext, phiVolumes);
    auto detector = Detector::makeShared("DetectorInPhi", phiVolumes,
                                         detail::tryAllVolumes());
  }
}

BOOST_AUTO_TEST_CASE(ProtoContainerRZ) {

  auto transform = Acts::Transform3::Identity();

  std::vector<Acts::ActsScalar> innerMostRadii = {0., 2.};

  for (auto [ir, imr] : Acts::enumerate(innerMostRadii)) {
    auto portalGenerator = detail::defaultPortalGenerator();
    auto navigationStateUpdator = detail::allPortals();
    // A container in R
    std::vector<Acts::ActsScalar> radii = {25., 100., 200.};
    Acts::ActsScalar halfZ = 200;

    // An innermost Pipe
    auto bBounds =
        std::make_unique<Acts::CylinderVolumeBounds>(imr, radii[0u], halfZ);

    auto innerPipe = DetectorVolumeFactory::construct(
        portalGenerator, tContext, "InnerPipe", transform, std::move(bBounds),
        detail::allPortals());

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
            transform, std::move(cBounds), detail::allPortals()));
      }
    }

    auto protoContainerInR = connectDetectorVolumesInR(tContext, rVolumes);

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
            std::move(cBounds), detail::allPortals()));
      }
    }
    // Now call the connector
    auto protoContainerInZ = connectDetectorVolumesInZ(tContext, zVolumes);

    auto centralContainer = connectContainersInR(
        tContext, {ipContainer, protoContainerInR, protoContainerInZ});

    // Let's make two endcaps
    // Nec
    auto necBounds =
        std::make_unique<Acts::CylinderVolumeBounds>(imr, 300., 50.);

    auto necTransform = Acts::Transform3::Identity();
    necTransform.pretranslate(Acts::Vector3(0., 0., -250));
    auto necVolume = DetectorVolumeFactory::construct(
        portalGenerator, tContext, "Nec", necTransform, std::move(necBounds),
        detail::allPortals());

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
        std::move(pecInnerBounds), detail::allPortals());
    auto pecOuter = DetectorVolumeFactory::construct(
        portalGenerator, tContext, "PecOuter", pecTransform,
        std::move(pecOuterBounds), detail::allPortals());

    std::vector<std::shared_ptr<DetectorVolume>> pecVolumes = {pecInner,
                                                               pecOuter};
    auto pecContainer = connectDetectorVolumesInR(tContext, pecVolumes);

    auto overallContainer = connectContainersInZ(
        tContext, {necContainer, centralContainer, pecContainer});

    //  Add them togeter
    std::vector<std::shared_ptr<DetectorVolume>> dVolumes;
    dVolumes.push_back(innerPipe);
    dVolumes.push_back(necVolume);
    dVolumes.insert(dVolumes.end(), rVolumes.begin(), rVolumes.end());
    dVolumes.insert(dVolumes.end(), zVolumes.begin(), zVolumes.end());
    dVolumes.push_back(pecInner);
    dVolumes.push_back(pecOuter);

    auto detector = Detector::makeShared("DetectorFromProtoContainer", dVolumes,
                                         detail::tryAllVolumes());
  } // test with different innermost radii
}

BOOST_AUTO_TEST_SUITE_END()
