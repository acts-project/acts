// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/GridDetectorBuilder.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"

#include <memory>
#include <set>
#include <vector>

using namespace Acts;
using namespace Acts::Test;
using namespace Acts::Experimental;

GeometryContext tContext;
CylindricalTrackingGeometry cGeometry = CylindricalTrackingGeometry(tContext);

namespace {
/// Helper mehtod that allows to use the already existing testing
/// infrastructure with the new const-correct detector design
///
std::vector<std::shared_ptr<Acts::Surface>> unpackSurfaces(
    const std::vector<const Acts::Surface*>& surfaces) {
  std::vector<std::shared_ptr<Acts::Surface>> uSurfaces;
  uSurfaces.reserve(surfaces.size());
  for (const auto& s : surfaces) {
    Surface* ncs = const_cast<Surface*>(s);
    uSurfaces.push_back(ncs->getSharedPtr());
  }
  return uSurfaces;
}

/// @brief  Simple helper struct acting as a surface provider
struct SurfaceProvider {
  std::vector<std::shared_ptr<Acts::Surface>> surfaces;

  std::vector<std::shared_ptr<Acts::Surface>> operator()() { return surfaces; }
};

}  // namespace

BOOST_AUTO_TEST_SUITE(Detector)

// Test the creation of a ring like structure
BOOST_AUTO_TEST_CASE(GridDetectorBuilder_construction) {
  Acts::Experimental::GridDetectorBuilder::Config gdbConfig;

  // A non-valid buiilder: no binning provided
  BOOST_CHECK_THROW(Acts::Experimental::GridDetectorBuilder(
                        gdbConfig, Acts::getDefaultLogger("GridDetectorBuilder",
                                                          Logging::VERBOSE)),
                    std::invalid_argument);

  // A non-valid buiilder: only 2 binning description provided
  gdbConfig.binning = {
      {BinningData(BinningOption::open, BinningValue::binZ, 10, 0., 100.), 0},
      {BinningData(BinningOption::closed, BinningValue::binPhi, 10, -M_PI,
                   M_PI),
       0}};

  BOOST_CHECK_THROW(Acts::Experimental::GridDetectorBuilder(
                        gdbConfig, Acts::getDefaultLogger("GridDetectorBuilder",
                                                          Logging::VERBOSE)),
                    std::invalid_argument);

  // A non-valid buiilder: wrong binning description provided
  gdbConfig.binning = {
      {BinningData(BinningOption::open, BinningValue::binX, 10, 0., 100.), 0},
      {BinningData(BinningOption::open, BinningValue::binY, 10, 0., 100.), 0},
      {BinningData(BinningOption::open, BinningValue::binR, 10, 0., 100.), 0}};

  BOOST_CHECK_THROW(Acts::Experimental::GridDetectorBuilder(
                        gdbConfig, Acts::getDefaultLogger("GridDetectorBuilder",
                                                          Logging::VERBOSE)),
                    std::invalid_argument);

  // A valid builder - z,r,phi
  gdbConfig.binning = {
      {BinningData(BinningOption::open, BinningValue::binZ, 10, 0., 100.), 0},
      {BinningData(BinningOption::open, BinningValue::binR, 10, 0., 100.), 0},
      {BinningData(BinningOption::closed, BinningValue::binPhi, 10, -M_PI,
                   M_PI),
       0}};

  BOOST_CHECK_NO_THROW(Acts::Experimental::GridDetectorBuilder(
      gdbConfig,
      Acts::getDefaultLogger("GridDetectorBuilder", Logging::VERBOSE)));

  // Another valid builder - r,z,phi (unordere)
  // A valid builder - z,r,phi
  gdbConfig.binning = {
      {BinningData(BinningOption::open, BinningValue::binR, 10, 0., 100.), 0},
      {BinningData(BinningOption::open, BinningValue::binZ, 10, 0., 100.), 0},
      {BinningData(BinningOption::closed, BinningValue::binPhi, 10, -M_PI,
                   M_PI),
       0}};

  BOOST_CHECK_NO_THROW(Acts::Experimental::GridDetectorBuilder(
      gdbConfig,
      Acts::getDefaultLogger("GridDetectorBuilder", Logging::VERBOSE)));

  // And finally a valid builder - x,y,z
  gdbConfig.binning = {
      {BinningData(BinningOption::open, BinningValue::binX, 10, 0., 100.), 0},
      {BinningData(BinningOption::open, BinningValue::binY, 10, 0., 100.), 0},
      {BinningData(BinningOption::open, BinningValue::binZ, 10, 0., 100.), 0}};

  BOOST_CHECK_NO_THROW(Acts::Experimental::GridDetectorBuilder(
      gdbConfig,
      Acts::getDefaultLogger("GridDetectorBuilder", Logging::VERBOSE)));

  // Invalid again builder - x,y,z - but one is closed
  gdbConfig.binning = {
      {BinningData(BinningOption::open, BinningValue::binX, 10, 0., 100.), 0},
      {BinningData(BinningOption::closed, BinningValue::binY, 10, 0., 100.), 0},
      {BinningData(BinningOption::open, BinningValue::binZ, 10, 0., 100.), 0}};

  BOOST_CHECK_THROW(Acts::Experimental::GridDetectorBuilder(
                        gdbConfig, Acts::getDefaultLogger("GridDetectorBuilder",
                                                          Logging::VERBOSE)),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(GridDetectorBuilder_eqBeqBeqC_zrphi) {
  // Create a tracking geometry

  auto trackingGeometry = cGeometry();

  // Collect the surfaces
  struct SurfaceCollector {
    std::vector<const Acts::Surface*> surfaces;

    void operator()(const Acts::Surface* s) {
      if (s->geometryId().boundary() == 0u) {
        surfaces.push_back(s);
      }
    }
  };

  SurfaceCollector sc;
  trackingGeometry->visitSurfaces(sc);

  // Provide the surfaces
  SurfaceProvider sp;
  sp.surfaces = unpackSurfaces(sc.surfaces);

  int binsZ = 50;
  int binsR = 6;
  int binsPhi = 36;
  // The configuration
  Acts::Experimental::GridDetectorBuilder::Config gdbConfig;
  gdbConfig.name = "EqEqEq_zrphi";
  gdbConfig.surfaces = sp;
  gdbConfig.binning = {
      {BinningData(BinningOption::open, BinningValue::binZ, binsZ, -1100.,
                   1100.),
       0},
      {BinningData(BinningOption::open, BinningValue::binR, binsR, 0., 300.),
       0},
      {BinningData(BinningOption::closed, BinningValue::binPhi, binsPhi, -M_PI,
                   M_PI),
       0}};

  // The builder
  auto eeeBuilder = Acts::Experimental::GridDetectorBuilder(
      gdbConfig,
      Acts::getDefaultLogger("GridDetectorBuilder", Logging::VERBOSE));

  auto eeeDetector = eeeBuilder.construct(tContext);
  BOOST_CHECK(eeeDetector != nullptr);

  auto eeeDetectorVolumes = eeeDetector->volumes();
  BOOST_CHECK(eeeDetectorVolumes.size() ==
              std::size_t(binsZ * binsR * binsPhi));

  // Count the unqiue surfaces
  std::set<const Surface*> uniqueSurfaces;
  // Count the surface hits
  std::size_t surfaceHits = 0;

  // Check the number of portals
  //
  // Individual glueing setup for
  //
  // Starting point:
  // ->  binsZ * binsPhi inner volumes with 5 portal each, (binsR-1) * binsZ *
  // binsPhi outer volumes with 6 portal each Unconnected portals:
  // ->  binZ * binPhi * 5 + binZ * binPhi * (binR-1) * 6
  // After connection in phi: each volume loses one portal (circular):
  // -> binZ * binPhi * 4 + binZ * binPhi * (binR-1) * 5
  // After connection in r: all volumes in the outer rings lose one portal
  // -> binZ * binPhi * 4 + binZ * binPhi * (binR-1) * 4
  // After connection in z, (rBins * phiBins * zBins-1) lose one portal each
  // >>> binZ * binPhi * 4 + binZ * binPhi * (binR-1) * 4 - (rBins * phiBins *
  // zBins-1)
  std::size_t expectedPortals = binsZ * binsPhi * 4 +
                                binsZ * binsPhi * (binsR - 1) * 4 -
                                (binsR * binsPhi * (binsZ - 1));

  // We should expect 156 indivudal portal objects
  // Let's coutn them
  std::set<const Portal*> uniquePortals;
  for (const auto& edv : eeeDetectorVolumes) {
    for (const auto& p : edv->portals()) {
      uniquePortals.insert(p);
      for (const auto& s : edv->surfaces()) {
        uniqueSurfaces.insert(s);
      }
      surfaceHits += edv->surfaces().size();
    }
  }
  // Check if portal setup if fine
  BOOST_CHECK(uniquePortals.size() == expectedPortals);

  // Visualization
  ObjVisualization3D obj;
  for (auto& edv : eeeDetectorVolumes) {
    GeometryView3D::drawDetectorVolume(obj, *edv, tContext);
  }
  obj.write(gdbConfig.name);
}

BOOST_AUTO_TEST_SUITE_END()
