// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "ActsFatras/Digitization/PlanarSurfaceDrift.hpp"

#include <array>
#include <fstream>
#include <memory>
#include <string>

namespace ActsFatras {

BOOST_AUTO_TEST_SUITE(Digitization)

BOOST_AUTO_TEST_CASE(PlanarSurfaceDriftSimpleTests) {
  Acts::GeometryContext geoCtx;

  ActsFatras::PlanarSurfaceDrift psd;

  Acts::Vector3 cPosition = Acts::Vector3(10., 50., 12.);
  Acts::Vector3 cNormal = Acts::Vector3(1., 1., 1.).normalized();

  std::shared_ptr<Acts::PlaneSurface> planeSurface =
      Acts::CurvilinearSurface(cPosition, cNormal).planeSurface();

  double depletion = 0.250;

  // Nominal intersection
  Acts::Vector3 noDrift(0., 0., 0.);

  // Intersect surface at normal direction and no drift
  //
  // -> resulting segment must have entry & exit at (0,0) local coordinates
  auto [noDriftSegment, oSegment] =
      psd.toReadout(geoCtx, *planeSurface, depletion, cPosition, cNormal,
                    noDrift)
          .value();

  CHECK_CLOSE_ABS(noDriftSegment[0].x(), 0., Acts::s_epsilon);
  CHECK_CLOSE_ABS(noDriftSegment[0].y(), 0., Acts::s_epsilon);
  CHECK_CLOSE_ABS(noDriftSegment[1].x(), 0., Acts::s_epsilon);
  CHECK_CLOSE_ABS(noDriftSegment[1].y(), 0., Acts::s_epsilon);

  Acts::Vector3 particleDir = Acts::Vector3(2., 1., 1.).normalized();
  // Intersect surface at particleDirection != normal and no drift
  //
  // -> local segment must be symmetric around (0,0)
  auto [noDriftSegment1, dScale1] =
      psd.toReadout(geoCtx, *planeSurface, depletion, cPosition, particleDir,
                    noDrift)
          .value();

  CHECK_CLOSE_ABS(noDriftSegment1[0].x(), -noDriftSegment1[1].x(),
                  Acts::s_epsilon);
  CHECK_CLOSE_ABS(noDriftSegment1[0].y(), -noDriftSegment1[1].y(),
                  Acts::s_epsilon);
}

BOOST_AUTO_TEST_CASE(PlanarSurfaceDriftEnhancedTests) {
  // A really thick sensor :-)
  double thickness = 20.;

  Acts::Vector3 cPosition = Acts::Vector3(50., 50., 0.);

  Acts::Vector3 localZ = Acts::Vector3(1., 1., 0.).normalized();
  Acts::Vector3 localY = Acts::Vector3(0., 0., 1.);
  Acts::Vector3 localX = -localZ.cross(localY).normalized();

  Acts::RotationMatrix3 rotationMatrix;
  rotationMatrix.col(0) = localX;
  rotationMatrix.col(1) = localY;
  rotationMatrix.col(2) = localZ;

  auto entryTransform = Acts::Transform3(
      Acts::Translation3(cPosition - 0.5 * thickness * localZ) *
      rotationMatrix);
  auto centralTransform =
      Acts::Transform3(Acts::Translation3(cPosition) * rotationMatrix);
  auto exitTransform = Acts::Transform3(
      Acts::Translation3(cPosition + 0.5 * thickness * localZ) *
      rotationMatrix);

  // Create the entry and exit surface
  auto entrySurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      entryTransform, std::make_shared<Acts::RectangleBounds>(100., 100.));
  auto readoutSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      centralTransform, std::make_shared<Acts::RectangleBounds>(100., 100.));
  auto exitSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      exitTransform, std::make_shared<Acts::RectangleBounds>(100., 100.));

  std::vector<std::shared_ptr<Acts::Surface>> surfaces = {
      entrySurface, readoutSurface, exitSurface};

  std::vector<std::string> surfaceNames = {"EntrySurface", "ReadoutSurface",
                                           "ExitSurface"};

  Acts::Vector3 startingPoint = Acts::Vector3(0., 0., 0.);
  Acts::Vector3 particleDirection = Acts::Vector3(1.5, 0.3, 1.25).normalized();

  // Intersection positions
  std::vector<Acts::Vector3> intersectionPositions = {};

  // A test context
  Acts::GeometryContext tContext;

  // Original 3D segment output
  std::ofstream fiout;
  fiout.open("PlanarSurfaceDrift_segment_orig_3D.obj");

  for (auto [isf, sf] : Acts::enumerate(surfaces)) {
    // Intersect surface at particleDirection != normal and no drift
    auto sIntersection =
        sf->intersect(tContext, startingPoint, particleDirection).closest();
    if (!sIntersection.isValid()) {
      BOOST_FAIL("Intersection failed");
    }

    // Record & write
    intersectionPositions.push_back(sIntersection.position());
    fiout << "v " << sIntersection.position().x() << " "
          << sIntersection.position().y() << " " << sIntersection.position().z()
          << "\n";
    // Write out all the surfaces in separate files
    std::ofstream fsout;
    fsout.open("PlanarSurfaceDrift_" + surfaceNames[isf] + ".obj");

    auto polyHedron = sf->polyhedronRepresentation(tContext);
    for (const auto& vertex : polyHedron.vertices) {
      fsout << "v " << vertex.x() << " " << vertex.y() << " " << vertex.z()
            << "\n";
    }
    for (const auto& face : polyHedron.faces) {
      fsout << "f";
      for (const auto& vertex : face) {
        fsout << " " << vertex + 1;
      }
      fsout << "\n";
    }
    fsout.close();
  }

  fiout << "l 1 2\n";
  fiout << "l 2 3\n";
  fiout.close();

  // Readout grid output
  std::ofstream gout;
  gout.open("PlanarSurfaceDrift_readout_2D.obj");
  for (unsigned int ixy = 0; ixy < 40; ++ixy) {
    auto lx = cPosition - 95 * localY - 95 * localX + ixy * 5 * localX;
    auto rx = cPosition + 95 * localY - 95 * localX + ixy * 5 * localX;
    gout << "v " << lx.x() << " " << lx.y() << " " << lx.z() << "\n";
    gout << "v " << rx.x() << " " << rx.y() << " " << rx.z() << "\n";
    gout << "l " << ixy * 4 + 1 << " " << ixy * 4 + 2 << "\n";
    auto ly = cPosition - 95 * localX - 95 * localY + ixy * 5 * localY;
    auto ry = cPosition + 95 * localX - 95 * localY + ixy * 5 * localY;
    gout << "v " << ly.x() << " " << ly.y() << " " << ly.z() << "\n";
    gout << "v " << ry.x() << " " << ry.y() << " " << ry.z() << "\n";
    gout << "l " << ixy * 4 + 3 << " " << ixy * 4 + 4 << "\n";
  }
  gout.close();

  std::vector<Acts::Vector3> driftScenarios = {
      Acts::Vector3(0., 0., 1.),
      Acts::Vector3(0.8, 0.1, 1.).normalized(),
      Acts::Vector3(0.8, 0.1, -1.).normalized(),
  };

  // The drift module
  ActsFatras::PlanarSurfaceDrift psd;

  auto writeSegment = [&](const std::string& tag, const Acts::Vector3& p0,
                          const Acts::Vector3& p1) -> void {
    std::ofstream out;
    out.open(tag + ".obj");
    out << "v " << p0.x() << " " << p0.y() << " " << p0.z() << "\n";
    out << "v " << p1.x() << " " << p1.y() << " " << p1.z() << "\n";
    out << "l 1 2\n";
    out.close();
  };

  // Run through different drfit scenarios
  for (auto [id, drift] : Acts::enumerate(driftScenarios)) {
    std::string scenario =
        "PlanarSurfaceDrift_scenario_" + std::to_string(id) + "_";

    // Without lorentz angle drift
    auto tr =
        psd.toReadout(tContext, *readoutSurface, thickness,
                      intersectionPositions[1u], particleDirection, drift);
    BOOST_CHECK(tr.ok());

    auto [dSegment, oSegment] = tr.value();

    Acts::Vector3 oEntry = readoutSurface->transform(tContext) * oSegment[0];
    Acts::Vector3 oExit = readoutSurface->transform(tContext) * oSegment[1];
    writeSegment(scenario + "drifted_3D", oEntry, oExit);

    // Checking if entry and exit position are correct
    BOOST_CHECK(oEntry.isApprox(intersectionPositions[0u], Acts::s_epsilon));
    BOOST_CHECK(oExit.isApprox(intersectionPositions[2u], Acts::s_epsilon));

    Acts::Vector3 driftedEntry =
        readoutSurface->transform(tContext) *
        Acts::Vector3(dSegment[0].x(), dSegment[0].y(), 0.);
    Acts::Vector3 driftedExit =
        readoutSurface->transform(tContext) *
        Acts::Vector3(dSegment[1].x(), dSegment[1].y(), 0.);

    writeSegment(scenario + "_2D", driftedEntry, driftedExit);

    // Drift entry to readout
    double dts = 0.5 * thickness / drift.z();

    Acts::Vector3 drift3D =
        readoutSurface->transform(tContext).linear() * drift;
    writeSegment(scenario + "entry_to_readout_3D", intersectionPositions[0u],
                 intersectionPositions[0u] + dts * drift3D);

    BOOST_CHECK(Acts::Vector3(intersectionPositions[0u] + dts * drift3D)
                    .isApprox(driftedEntry, Acts::s_epsilon));

    writeSegment(scenario + "exit_to_readout_3D", intersectionPositions[2u],
                 intersectionPositions[2u] - dts * drift3D);

    BOOST_CHECK(Acts::Vector3(intersectionPositions[2u] - dts * drift3D)
                    .isApprox(driftedExit, Acts::s_epsilon));
  }

  // Special casing: particle normal to the surface
  auto trn = psd.toReadout(tContext, *readoutSurface, thickness,
                           intersectionPositions[1u], localZ, localZ);
  BOOST_CHECK(trn.ok());
  auto [dSegmentN, oSegmentN] = trn.value();
  BOOST_CHECK(dSegmentN[0].isApprox(dSegmentN[1], Acts::s_epsilon));
  CHECK_CLOSE_ABS(Acts::Vector3(oSegmentN[1] - oSegmentN[0]).norm(), thickness,
                  Acts::s_epsilon);

  // Error checking, particle along the xy plane
  auto trnErr = psd.toReadout(tContext, *readoutSurface, thickness,
                              intersectionPositions[1u], localX, localY);
  BOOST_CHECK(!trnErr.ok());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsFatras
