// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/GridDetectorVolumesBuilder.hpp"
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

/// @brief helper struct to collect the surface
struct SurfaceCollector {
  std::vector<const Acts::Surface*> surfaces = {};

  void operator()(const Acts::Surface* s) {
    if (s->geometryId().boundary() == 0u) {
      surfaces.push_back(s);
    }
  }
};

/// @brief  Simple helper struct acting as a surface provider
struct SurfaceProvider {
  std::vector<std::shared_ptr<Acts::Surface>> surfaces;

  std::vector<std::shared_ptr<Acts::Surface>> operator()() { return surfaces; }
};

/// @brief Helper method to run the test - base code
///
/// @param cfg is the configuration of the builder
/// @param extendsToZero is a flag to indicate if rMin == 0
/// @param runPortalChecks is the flag to run portal checks
void gridZRPhhiBase(
    const Acts::Experimental::GridDetectorVolumesBuilder::Config& cfg,
    bool extendsToZero = false, bool runPortalChecks = true) {
  std::size_t binsZ = cfg.binning[0].bins();
  std::size_t binsR = cfg.binning[1].bins();
  std::size_t binsPhi = cfg.binning[2].bins();

  // The builder
  auto eeeBuilder = Acts::Experimental::GridDetectorVolumesBuilder(
      cfg,
      Acts::getDefaultLogger("GridDetectorVolumesBuilder", Logging::VERBOSE));

  auto [volumes, portalContainer, roots] = eeeBuilder.construct(tContext);

  BOOST_CHECK(volumes.size() == std::size_t(binsZ * binsR * binsPhi));

  // Count the unqiue surfaces
  std::set<const Surface*> uniqueSurfaces;

  // Check the number of portals
  //
  // Individual glueing setup for CYLINDER Like volumes
  //
  // Attention, whether the innermost has 5 or 6 portals depends on rMin > 0
  //
  //
  // Starting point:
  // ->  binsZ * binsPhi inner volumes with 5/6 portal each, (binsR-1) * binsZ *
  // binsPhi outer volumes with 6 portal each Unconnected portals:
  // ->  binZ * binPhi * 5 + binZ * binPhi * (binR-1) * 6
  // After connection in phi: each volume loses one portal (circular):
  // -> binZ * binPhi * 4 + binZ * binPhi * (binR-1) * 5
  // After connection in r: all volumes in the outer rings lose one portal
  // -> binZ * binPhi * 4 + binZ * binPhi * (binR-1) * 4
  // After connection in z, (rBins * phiBins * zBins-1) lose one portal each
  // >>> binZ * binPhi * 4 + binZ * binPhi * (binR-1) * 4 - (rBins * phiBins *
  // zBins-1)
  //
  int innermostAddon = extendsToZero > 0. ? 0 : 1;
  //
  std::size_t expectedPortals = binsZ * binsPhi * (4 + innermostAddon) +
                                binsZ * binsPhi * (binsR - 1) * 4u -
                                (binsR * binsPhi * (binsZ - 1));

  // We should expect 156 indivudal portal objects
  // Let's coutn them
  std::set<const Portal*> uniquePortals;
  for (const auto& edv : volumes) {
    for (const auto& p : edv->portals()) {
      uniquePortals.insert(p);
      for (const auto& s : edv->surfaces()) {
        uniqueSurfaces.insert(s);
      }
    }
  }
  // Check if portal setup if fine
  if (runPortalChecks) {
    BOOST_CHECK(uniquePortals.size() == expectedPortals);
    // Check if the portal container is present
  }

  // Full portal definition
  BOOST_CHECK(portalContainer.size() == 6u);

  BOOST_CHECK(portalContainer[0u].size() == binsR * binsPhi);
  BOOST_CHECK(portalContainer[1u].size() == binsR * binsPhi);
  BOOST_CHECK(portalContainer[2u].size() == binsZ * binsPhi);
  // Check validity if inner one is present
  if (portalContainer[3].size() > 0u) {
    BOOST_CHECK(portalContainer[3u].size() == binsZ * binsPhi);
  }
  // Check validity if phi is not closed
  if (portalContainer[4u].size() > 0u) {
    BOOST_CHECK(portalContainer[4u].size() == binsZ * binsR);
  }
  if (portalContainer[5u].size() > 0u) {
    BOOST_CHECK(portalContainer[5u].size() == binsZ * binsR);
  }

  // Visualization
  ObjVisualization3D obj;
  for (auto& edv : volumes) {
    GeometryView3D::drawDetectorVolume(obj, *edv, tContext);
  }
  obj.write(cfg.name);
}

/// @brief Helper method to run the test - equidistant binning
///
/// @param detectorName is the name of the detector
/// @param binsZ is the number of bins in z
/// @param binsR is the number of bins in r
/// @param binsPhi is the number of bins in phi
/// @param rMin is the minimum radius
/// @param approximateCylinders is the flag to approximate cylinders
/// @param runPortalChecks is the flag to run portal checks
/// @param phiBoundary - special check for open phi setup
void gridZRPhhiTest(const std::string& detectorName, int binsZ = 50,
                    int binsR = 6, int binsPhi = 36, Acts::ActsScalar rMin = 0.,
                    bool approximateCylinders = false,
                    bool runPortalChecks = true,
                    const std::array<Acts::ActsScalar, 2>& phiRange = {-M_PI,
                                                                       M_PI}) {
  // Create a tracking geometry for reference conversion
  auto trackingGeometry = cGeometry();

  SurfaceCollector sc;
  trackingGeometry->visitSurfaces(sc);

  // Provide the surfaces
  SurfaceProvider sp;
  sp.surfaces = unpackSurfaces(sc.surfaces);

  // The configuration
  Acts::Experimental::GridDetectorVolumesBuilder::Config gdbConfig;
  gdbConfig.name = detectorName;
  gdbConfig.surfaces = sp;
  gdbConfig.binning = {
      ProtoBinning(BinningValue::binZ, Acts::detail::AxisBoundaryType::Bound,
                   -1100., 1100., binsZ, 0),
      ProtoBinning(BinningValue::binR, Acts::detail::AxisBoundaryType::Bound,
                   rMin, 300., binsR, 0),
      ProtoBinning(BinningValue::binPhi, Acts::detail::AxisBoundaryType::Closed,
                   phiRange[0u], phiRange[1u], binsPhi, 0),
  };

  gdbConfig.approximateCylinders = approximateCylinders;

  gridZRPhhiBase(gdbConfig, (rMin == 0.), runPortalChecks);
}

/// @brief Helper method to run the test - variable binning
///
/// @param detectorName is the name of the detector
/// @param edgesZ are the bin edges in z
/// @param edgesR are the bin edges in R
/// @param binsPhi is the number of bins in phi
/// @param approximateCylinders is the flag to approximate cylinders
/// @param innermostSignleBin is a flat to make the innermost bin a single bin
/// @param runPortalChecks is the flag to run portal checks
void gridZRPhhiTest(const std::string& detectorName,
                    const std::vector<Acts::ActsScalar>& edgesZ,
                    const std::vector<Acts::ActsScalar>& edgesR,
                    std::size_t binsPhi, bool approximateCylinders = false,
                    bool runPortalChecks = true) {
  // Create a tracking geometry for reference conversion
  auto trackingGeometry = cGeometry();

  SurfaceCollector sc;
  trackingGeometry->visitSurfaces(sc);

  // Provide the surfaces
  SurfaceProvider sp;
  sp.surfaces = unpackSurfaces(sc.surfaces);

  // The configuration
  Acts::Experimental::GridDetectorVolumesBuilder::Config gdbConfig;
  gdbConfig.name = detectorName;
  gdbConfig.surfaces = sp;
  gdbConfig.binning = {
      ProtoBinning(BinningValue::binZ, Acts::detail::AxisBoundaryType::Bound,
                   edgesZ, 0),
      ProtoBinning(BinningValue::binR, Acts::detail::AxisBoundaryType::Bound,
                   edgesR, 0),
      ProtoBinning(BinningValue::binPhi, Acts::detail::AxisBoundaryType::Closed,
                   -M_PI, M_PI, binsPhi, 0)};

  gdbConfig.approximateCylinders = approximateCylinders;

  gridZRPhhiBase(gdbConfig, (edgesR[0u] == 0.), runPortalChecks);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(Detector)

// Test the creation of a ring like structure
BOOST_AUTO_TEST_CASE(GridDetectorVolumesBuilder_construction) {
  Acts::Experimental::GridDetectorVolumesBuilder::Config gdbConfig;

  // A non-valid buiilder: no binning provided
  BOOST_CHECK_THROW(
      Acts::Experimental::GridDetectorVolumesBuilder(
          gdbConfig, Acts::getDefaultLogger("GridDetectorVolumesBuilder",
                                            Logging::VERBOSE)),
      std::invalid_argument);

  // A non-valid buiilder: only 2 binning description provided
  gdbConfig.binning = {
      ProtoBinning(BinningValue::binZ, Acts::detail::AxisBoundaryType::Bound,
                   0., 100., 10, 0),
      ProtoBinning(BinningValue::binPhi, Acts::detail::AxisBoundaryType::Closed,
                   -M_PI, M_PI, 6, 0)};

  BOOST_CHECK_THROW(
      Acts::Experimental::GridDetectorVolumesBuilder(
          gdbConfig, Acts::getDefaultLogger("GridDetectorVolumesBuilder",
                                            Logging::VERBOSE)),
      std::invalid_argument);

  // A non-valid buiilder: wrong binning description provided
  gdbConfig.binning = {
      ProtoBinning(BinningValue::binX, Acts::detail::AxisBoundaryType::Bound,
                   0., 100., 10, 0),
      ProtoBinning(BinningValue::binY, Acts::detail::AxisBoundaryType::Bound,
                   0., 100., 10, 0),
      ProtoBinning(BinningValue::binR, Acts::detail::AxisBoundaryType::Bound,
                   0., 100., 10, 0)};

  BOOST_CHECK_THROW(
      Acts::Experimental::GridDetectorVolumesBuilder(
          gdbConfig, Acts::getDefaultLogger("GridDetectorVolumesBuilder",
                                            Logging::VERBOSE)),
      std::invalid_argument);

  // A valid builder - z,r,phi
  gdbConfig.binning = {
      ProtoBinning(BinningValue::binZ, Acts::detail::AxisBoundaryType::Bound,
                   0., 100., 10, 0),
      ProtoBinning(BinningValue::binR, Acts::detail::AxisBoundaryType::Bound,
                   0., 100., 10, 0),
      ProtoBinning(BinningValue::binPhi, Acts::detail::AxisBoundaryType::Closed,
                   -M_PI, M_PI, 6, 0)};

  BOOST_CHECK_NO_THROW(Acts::Experimental::GridDetectorVolumesBuilder(
      gdbConfig,
      Acts::getDefaultLogger("GridDetectorVolumesBuilder", Logging::VERBOSE)));

  // Another valid builder - r,z,phi (unordere)
  // A valid builder - z,r,phi
  gdbConfig.binning = {
      ProtoBinning(BinningValue::binR, Acts::detail::AxisBoundaryType::Bound,
                   0., 100., 10, 0),
      ProtoBinning(BinningValue::binZ, Acts::detail::AxisBoundaryType::Bound,
                   0., 100., 10, 0),
      ProtoBinning(BinningValue::binPhi, Acts::detail::AxisBoundaryType::Closed,
                   -M_PI, M_PI, 6, 0)};

  BOOST_CHECK_NO_THROW(Acts::Experimental::GridDetectorVolumesBuilder(
      gdbConfig,
      Acts::getDefaultLogger("GridDetectorVolumesBuilder", Logging::VERBOSE)));

  // And finally a valid builder - x,y,z
  gdbConfig.binning = {
      ProtoBinning(BinningValue::binX, Acts::detail::AxisBoundaryType::Bound,
                   0., 100., 10, 0),
      ProtoBinning(BinningValue::binY, Acts::detail::AxisBoundaryType::Bound,
                   0., 100., 10, 0),
      ProtoBinning(BinningValue::binZ, Acts::detail::AxisBoundaryType::Bound,
                   0., 100., 10, 0)};

  BOOST_CHECK_NO_THROW(Acts::Experimental::GridDetectorVolumesBuilder(
      gdbConfig,
      Acts::getDefaultLogger("GridDetectorVolumesBuilder", Logging::VERBOSE)));

  // Invalid again builder - x,y,z - but one is closed
  gdbConfig.binning = {
      ProtoBinning(BinningValue::binX, Acts::detail::AxisBoundaryType::Bound,
                   0., 100., 10, 0),
      ProtoBinning(BinningValue::binY, Acts::detail::AxisBoundaryType::Closed,
                   0., 100., 10, 0),
      ProtoBinning(BinningValue::binZ, Acts::detail::AxisBoundaryType::Bound,
                   0., 100., 10, 0)};

  BOOST_CHECK_THROW(
      Acts::Experimental::GridDetectorVolumesBuilder(
          gdbConfig, Acts::getDefaultLogger("GridDetectorVolumesBuilder",
                                            Logging::VERBOSE)),
      std::invalid_argument);
}

// Test equidistant binning
BOOST_AUTO_TEST_CASE(GridDetectorVolumesBuilder_eqB_z) {
  BOOST_CHECK_NO_THROW(
      gridZRPhhiTest("EqEqEq_z_10", 10, 1, 1, 0., false, false));
}

BOOST_AUTO_TEST_CASE(GridDetectorVolumesBuilder_eqB_r) {
  BOOST_CHECK_NO_THROW(
      gridZRPhhiTest("EqEqEq_r_10", 1, 10, 1, 0., false, false));
}

BOOST_AUTO_TEST_CASE(GridDetectorVolumesBuilder_eqBeqB_zr) {
  BOOST_CHECK_NO_THROW(
      gridZRPhhiTest("EqEqEq_zr_10_20", 10, 20, 1, 0., false, false));
}

BOOST_AUTO_TEST_CASE(GridDetectorVolumesBuilder_eqBeqBeqC_zrphi_to0) {
  BOOST_CHECK_NO_THROW(gridZRPhhiTest("EqEqEq_zrphi_10_2_6", 10, 2, 6));
}

BOOST_AUTO_TEST_CASE(GridDetectorVolumesBuilder_eqBeqBeqC_zrphi) {
  BOOST_CHECK_NO_THROW(
      gridZRPhhiTest("EqEqEq_zrphi_19_8_36_rmin10", 10, 8, 36, 10.));
}

BOOST_AUTO_TEST_CASE(GridDetectorVolumesBuilder_eqBeqBeqB_zrphi_open) {
  BOOST_CHECK_NO_THROW(gridZRPhhiTest("EqEqEq_zrphi_19_8_36_rmin10_open", 10, 8,
                                      18, 10., false, false, {-1.8, 1.8}));
}

BOOST_AUTO_TEST_CASE(GridDetectorVolumesBuilder_eqBeqBeqC_zrphi_trapezoids) {
  BOOST_CHECK_NO_THROW(
      gridZRPhhiTest("EqEqEq_zrphi_trapezoids", 10, 8, 36, 10., true));

  // This should throw and exception - minimum 3 bins in phi are needed for
  // trapezoids
  BOOST_CHECK_THROW(
      gridZRPhhiTest("EqEqEq_zrphi_trapezoids_throw", 10, 8, 2, 10., true),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(
    GridDetectorVolumesBuilder_eqBeqBeqB_zrphi_open_trapezoids) {
  BOOST_CHECK_NO_THROW(
      gridZRPhhiTest("EqEqEq_zrphi_19_8_36_rmin10_open_trapezoids", 10, 8, 18,
                     10., true, false, {-1.8, 1.8}));
}

// Test variable binning
BOOST_AUTO_TEST_CASE(GridDetectorVolumesBuilder_varB_z) {
  BOOST_CHECK_NO_THROW(gridZRPhhiTest("Var_z_v", {-800, 200, 500, 1000},
                                      {10, 200}, 1, false, false));
}

// Test variable binning
BOOST_AUTO_TEST_CASE(GridDetectorVolumesBuilder_varB_r) {
  BOOST_CHECK_NO_THROW(
      gridZRPhhiTest("Var_r_v", {-800, 800}, {10, 180, 200}, 1, false, false));
}

BOOST_AUTO_TEST_SUITE_END()
