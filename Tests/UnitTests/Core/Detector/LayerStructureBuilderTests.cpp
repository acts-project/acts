// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <memory>
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

// Test the creation of the Ring
BOOST_AUTO_TEST_CASE(LayerStructureBuilder_creationRing) {
  // Detector store
  CylindricalTrackingGeometry::DetectorStore dStore;
  auto rSurfaces = cGeometry.surfacesRing(dStore, 6.4, 12.4, 36., 0.125, 0.,
                                          55., -800, 2., 22u);

  SurfaceProvider endcap{unpackSurfaces(rSurfaces)};

  using LayerBinning = Acts::Experimental::LayerStructureBuilder::Binning;

  // Configure the layer structure builder
  Acts::Experimental::LayerStructureBuilder::Config lsConfig;
  lsConfig.auxilliary = "*** Endcap with 22 surfaces ***";
  lsConfig.surfaces = endcap;
  lsConfig.binnings = {LayerBinning{
      Acts::BinningData(Acts::closed, Acts::binPhi, 22u, -M_PI, M_PI), 1u}};

  auto endcapBuilder = Acts::Experimental::LayerStructureBuilder(
      lsConfig, Acts::getDefaultLogger("EndcapBuilder", Logging::VERBOSE));

  auto [surfaces0, volumes0, surfacesUpdator0, volumeUpdator0] =
      endcapBuilder.construct(tContext);

  BOOST_CHECK(surfaces0.size() == 22u);
  BOOST_CHECK(surfacesUpdator0.connected());
  BOOST_CHECK(volumes0.empty());
  BOOST_CHECK(volumeUpdator0.connected());

  using LayerSupport = Acts::Experimental::LayerStructureBuilder::Support;

  lsConfig.auxilliary = "*** Endcap with 22 surfaces + 1 support disc ***";
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

  lsConfig.auxilliary =
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

// Test the creation of the Cylinder
BOOST_AUTO_TEST_CASE(LayerStructureKDT_creationCylinder) {
  CylindricalTrackingGeometry::DetectorStore dStore;
  auto cSurfaces = cGeometry.surfacesCylinder(dStore, 8.4, 36., 0.15, 0.145, 72,
                                              3., 2., {32u, 14u});

  SurfaceProvider barrel{unpackSurfaces(cSurfaces)};

  using LayerBinning = Acts::Experimental::LayerStructureBuilder::Binning;

  // Configure the layer structure builder
  Acts::Experimental::LayerStructureBuilder::Config lsConfig;
  lsConfig.auxilliary = "*** Barrel with 448 surfaces ***";
  lsConfig.surfaces = barrel;
  lsConfig.binnings = {
      LayerBinning{Acts::BinningData(Acts::open, Acts::binZ, 14u, -480., 480.),
                   1u},
      LayerBinning{
          Acts::BinningData(Acts::closed, Acts::binPhi, 32u, -M_PI, M_PI), 1u}};

  auto barrelBuilder = Acts::Experimental::LayerStructureBuilder(
      lsConfig, Acts::getDefaultLogger("BarrelBuilder", Logging::VERBOSE));

  auto [surfaces0, volumes0, surfacesUpdator0, volumeUpdator0] =
      barrelBuilder.construct(tContext);

  BOOST_CHECK(surfaces0.size() == 448u);
  BOOST_CHECK(surfacesUpdator0.connected());
  BOOST_CHECK(volumes0.empty());
  BOOST_CHECK(volumeUpdator0.connected());

  using LayerSupport = Acts::Experimental::LayerStructureBuilder::Support;

  lsConfig.auxilliary = "*** Barrel with 448 surfaces + 1 support cylinder ***";
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

  lsConfig.auxilliary =
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
