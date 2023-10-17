// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Navigation/NavigationDelegates.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <vector>

using namespace Acts;
using namespace Acts::Test;
using namespace Acts::Experimental;

GeometryContext tContext;
CylindricalTrackingGeometry cGeometry = CylindricalTrackingGeometry(tContext);

namespace {
/// Helper method that allows to use the already existing testing
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

}  // namespace

BOOST_AUTO_TEST_SUITE(Detector)

// Test the creation of a ring like structure
BOOST_AUTO_TEST_CASE(LayerStructureBuilder_creationRing) {
  // Detector store
  CylindricalTrackingGeometry::DetectorStore dStore;
  auto rSurfaces = cGeometry.surfacesRing(dStore, 6.4, 12.4, 36., 0.125, 0.,
                                          55., -800, 2., 22u);

  auto endcapSurfaces = std::make_shared<LayerStructureBuilder::SurfacesHolder>(
      unpackSurfaces(rSurfaces));
  // Configure the layer structure builder
  Acts::Experimental::LayerStructureBuilder::Config lsConfig;
  lsConfig.auxiliary = "*** Endcap with 22 surfaces ***";
  lsConfig.surfacesProvider = endcapSurfaces;
  lsConfig.binnings = {ProtoBinning(Acts::binPhi,
                                    Acts::detail::AxisBoundaryType::Closed,
                                    -M_PI, M_PI, 22u, 1u)};

  auto endcapBuilder = Acts::Experimental::LayerStructureBuilder(
      lsConfig, Acts::getDefaultLogger("EndcapBuilder", Logging::VERBOSE));

  auto [surfaces0, volumes0, surfacesUpdator0, volumeUpdator0] =
      endcapBuilder.construct(tContext);

  BOOST_CHECK(surfaces0.size() == 22u);
  BOOST_CHECK(surfacesUpdator0.connected());
  BOOST_CHECK(volumes0.empty());
  BOOST_CHECK(volumeUpdator0.connected());

  using LayerSupport = Acts::Experimental::ProtoSupport;

  lsConfig.auxiliary = "*** Endcap with 22 surfaces + 1 support disc ***";
  lsConfig.supports = {LayerSupport{
      {15., 10., 10., 0., 0.}, Acts::Surface::SurfaceType::Disc, {binZ, binR}}};
  endcapBuilder = Acts::Experimental::LayerStructureBuilder(
      lsConfig, Acts::getDefaultLogger("EndcapBuilder", Logging::VERBOSE));

  auto [surfaces1, volumes1, surfacesUpdator1, volumeUpdator1] =
      endcapBuilder.construct(tContext);

  BOOST_CHECK(surfaces1.size() == 22u + 1u);
  BOOST_CHECK(surfaces1.back()->type() == Acts::Surface::SurfaceType::Disc);
  BOOST_CHECK(surfacesUpdator1.connected());
  BOOST_CHECK(volumes1.empty());
  BOOST_CHECK(volumeUpdator1.connected());

  lsConfig.auxiliary =
      "*** Endcap with 22 surfaces + 1 support -> split into 11 planes ***";
  lsConfig.supports = {LayerSupport{{15., 10., 10., 0., 0.},
                                    Acts::Surface::SurfaceType::Disc,
                                    {binZ, binR},
                                    11u}};
  endcapBuilder = Acts::Experimental::LayerStructureBuilder(
      lsConfig, Acts::getDefaultLogger("EndcapBuilder", Logging::VERBOSE));

  auto [surfaces2, volumes2, surfacesUpdator2, volumeUpdator2] =
      endcapBuilder.construct(tContext);

  BOOST_CHECK(surfaces2.size() == 22u + 11u);
  BOOST_CHECK(surfaces2.back()->type() == Acts::Surface::SurfaceType::Plane);
  BOOST_CHECK(surfacesUpdator2.connected());
  BOOST_CHECK(volumes2.empty());
  BOOST_CHECK(volumeUpdator2.connected());
}

// Test the creation of a cylindrical structure
BOOST_AUTO_TEST_CASE(LayerStructureBuilder_creationCylinder) {
  CylindricalTrackingGeometry::DetectorStore dStore;
  auto cSurfaces = cGeometry.surfacesCylinder(dStore, 8.4, 36., 0.15, 0.145, 72,
                                              3., 2., {32u, 14u});

  auto barrelSurfaces = std::make_shared<LayerStructureBuilder::SurfacesHolder>(
      unpackSurfaces(cSurfaces));

  // Configure the layer structure builder
  Acts::Experimental::LayerStructureBuilder::Config lsConfig;
  lsConfig.auxiliary = "*** Barrel with 448 surfaces ***";
  lsConfig.surfacesProvider = barrelSurfaces;
  lsConfig.binnings = {Acts::Experimental::ProtoBinning{
                           Acts::binZ, Acts::detail::AxisBoundaryType::Bound,
                           -480., 480., 14u, 1u},
                       Acts::Experimental::ProtoBinning(
                           Acts::binPhi, Acts::detail::AxisBoundaryType::Closed,
                           -M_PI, M_PI, 32u, 1u)};

  auto barrelBuilder = Acts::Experimental::LayerStructureBuilder(
      lsConfig, Acts::getDefaultLogger("BarrelBuilder", Logging::VERBOSE));

  auto [surfaces0, volumes0, surfacesUpdator0, volumeUpdator0] =
      barrelBuilder.construct(tContext);

  BOOST_CHECK(surfaces0.size() == 448u);
  BOOST_CHECK(surfacesUpdator0.connected());
  BOOST_CHECK(volumes0.empty());
  BOOST_CHECK(volumeUpdator0.connected());

  using LayerSupport = Acts::Experimental::ProtoSupport;

  lsConfig.auxiliary = "*** Barrel with 448 surfaces + 1 support cylinder ***";
  lsConfig.supports = {LayerSupport{{15., 10., 10., 0., 0.},
                                    Acts::Surface::SurfaceType::Cylinder,
                                    {binZ, binR}}};

  barrelBuilder = Acts::Experimental::LayerStructureBuilder(
      lsConfig, Acts::getDefaultLogger("BarrelBuilder", Logging::VERBOSE));

  auto [surfaces1, volumes1, surfacesUpdator1, volumeUpdator1] =
      barrelBuilder.construct(tContext);

  BOOST_CHECK(surfaces1.size() == 448u + 1u);
  BOOST_CHECK(surfacesUpdator1.connected());
  BOOST_CHECK(volumes1.empty());
  BOOST_CHECK(volumeUpdator1.connected());

  lsConfig.auxiliary =
      "*** Barrel with 448 surfaces + 1 support -> split into 32 planes ***";
  lsConfig.supports = {LayerSupport{{15., 10., 10., 0., 0.},
                                    Acts::Surface::SurfaceType::Cylinder,
                                    {binZ, binR},
                                    32u}};

  barrelBuilder = Acts::Experimental::LayerStructureBuilder(
      lsConfig, Acts::getDefaultLogger("BarrelBuilder", Logging::VERBOSE));

  auto [surfaces2, volumes2, surfacesUpdator2, volumeUpdator2] =
      barrelBuilder.construct(tContext);

  BOOST_CHECK(surfaces2.size() == 448u + 32u);
  BOOST_CHECK(surfacesUpdator2.connected());
  BOOST_CHECK(volumes2.empty());
  BOOST_CHECK(volumeUpdator2.connected());
}

BOOST_AUTO_TEST_SUITE_END()
