// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Experimental/CylindricalVolumeBuilder.hpp"
#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/GeometricExtent.hpp"
#include "Acts/Experimental/InternalBlueprint.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"

#include <exception>
#include <memory>
#include <vector>

#include "GeometryHelper.hpp"

namespace Acts {

namespace Test {

namespace {

/// Helper method to create a display file
///
/// @param dv the detector volume
/// @param vcfg the volume view configuration
/// @param scfg the surface view configuration
///
void display(const DetectorVolume& dv,
             const ViewConfig vcfg = ViewConfig({0, 0, 100}),
             const ViewConfig scfg = ViewConfig({100, 0, 0})) {
  GeometryContext gctx;

  Acts::ObjVisualization3D objVis;
  for (const auto& psf : dv.portals()) {
    Acts::GeometryView3D::drawSurface(objVis, psf->surfaceRepresentation(),
                                      gctx, Transform3::Identity(), vcfg);
  }
  objVis.write(dv.name() + std::string(".obj"));

  if (not dv.surfaces().empty()) {
    Acts::ObjVisualization3D objVisSen;
    for (auto sf : dv.surfaces()) {
      Acts::GeometryView3D::drawSurface(objVisSen, *sf, gctx,
                                        Transform3::Identity(), scfg);
    }
    objVisSen.write(dv.name() + std::string("_sensitives.obj"));
  }
}

}  // namespace

BOOST_AUTO_TEST_SUITE(Experimental)

BOOST_AUTO_TEST_CASE(CylindricalVolumeBuilder_) {
  GeometryContext gctx;

  CylindricalVolumeBuilder<> cvb;

  // Check invalid input exceptions
  GeometricExtent unbound{};
  BOOST_CHECK_THROW(cvb(gctx), std::invalid_argument);

  // Check exception if more than one internal blue print is given
  cvb.internalBlueprints = {InternalBlueprint(gctx), InternalBlueprint(gctx)};
  BOOST_CHECK_THROW(cvb(gctx), std::invalid_argument);

  // Volumes without surfaces ---------------

  // A simple solid
  GeometricExtent solidExtent;
  solidExtent.set(binR, 0., 10.);
  solidExtent.set(binZ, -10., 10.);

  cvb.volumeExtent = solidExtent;
  cvb.name = "solid";
  auto solids = cvb(gctx);
  BOOST_CHECK(solids.size() == 1);
  auto solid = solids[0];
  BOOST_CHECK(solid->name() == "solid");
  display(*solid, ViewConfig({0, 0, 100}));

  // A shifted solid
  GeometricExtent shiftedSolidExtent;
  shiftedSolidExtent.set(binR, 0., 10.);
  shiftedSolidExtent.set(binZ, -40., -20.);
  cvb.volumeExtent = shiftedSolidExtent;
  cvb.name = "shifted_solid";
  auto shiftedSolids = cvb(gctx);
  BOOST_CHECK(shiftedSolids.size() == 1);
  auto shiftedSolid = shiftedSolids[0];
  BOOST_CHECK(shiftedSolid->name() == "shifted_solid");
  display(*shiftedSolid, ViewConfig({0, 0, 40}));

  // A shifted tube
  GeometricExtent shiftedTubeExtent;
  shiftedTubeExtent.set(binR, 2., 10.);
  shiftedTubeExtent.set(binZ, 20., 40.);
  cvb.volumeExtent = shiftedTubeExtent;
  cvb.name = "shifted_tube";
  auto shiftedTubes = cvb(gctx);
  BOOST_CHECK(shiftedTubes.size() == 1);
  auto shiftedTube = shiftedTubes[0];
  BOOST_CHECK(shiftedTube->name() == "shifted_tube");
  display(*shiftedTube, ViewConfig({0, 0, 160}));

  // A shifted sectoral tube
  GeometricExtent shiftedTubeSectorExtent;
  shiftedTubeSectorExtent.set(binR, 2., 10.);
  shiftedTubeSectorExtent.set(binZ, 50., 70.);
  shiftedTubeSectorExtent.set(binPhi, -0.85, 0.85);
  cvb.volumeExtent = shiftedTubeSectorExtent;
  cvb.name = "shifted_tube_sector";
  auto shiftedTubeSectors = cvb(gctx);
  BOOST_CHECK(shiftedTubeSectors.size() == 1);
  auto shiftedTubeSector = shiftedTubeSectors[0];
  BOOST_CHECK(shiftedTubeSector->name() == "shifted_tube_sector");
  display(*shiftedTubeSector, ViewConfig({0, 0, 220}));

  // Volumes with surfaces ---------------

  auto layer = surfacesCylinder(8.4, 36., 0.145, 32., 0.5, 5., {16, 14});
  InternalBlueprint ilb(gctx, layer, AllInternalSurfaces{4},
                        {{binZ, {0., 0.}}, {binR, {0., 0.}}}, "barrel_layer");
  cvb.internalBlueprints = {ilb};
  cvb.name = "barrel_layer_tight";
  auto barrelLayerTightVec = cvb(gctx);

  BOOST_CHECK(barrelLayerTightVec.size() == 1);
  auto barrelLayerTight = barrelLayerTightVec[0];

  display(*barrelLayerTight, ViewConfig({0, 0, 220}), ViewConfig({100, 0, 0}));

  // Build with envelope
  ExtentEnvelope envelopes = zeroEnvelopes;
  envelopes[Acts::binR] = {5., 5.};
  envelopes[Acts::binZ] = {5., 5.};

  cvb.volumeExtent = GeometricExtent{envelopes};
  cvb.name = "barrel_layer_envelope";
  auto barrelLayerEnvVec = cvb(gctx);

  BOOST_CHECK(barrelLayerEnvVec.size() == 1);
  auto barrelLayerEnv = barrelLayerEnvVec[0];

  display(*barrelLayerEnv, ViewConfig({0, 0, 220}), ViewConfig({100, 0, 0}));

  // Build with an external size
  GeometricExtent layerTubeExtent;
  layerTubeExtent.set(binR, 5., 100.);
  layerTubeExtent.set(binZ, -1000., 1000.);

  cvb.volumeExtent = layerTubeExtent;
  cvb.name = "barrel_layer_external";
  auto barrelLayerTubeVec = cvb(gctx);

  BOOST_CHECK(barrelLayerTubeVec.size() == 1);
  auto barrelLayerTube = barrelLayerTubeVec[0];

  display(*barrelLayerTube, ViewConfig({0, 0, 220}), ViewConfig({100, 0, 0}));

  // Build with an invalid external size - catch exception
  GeometricExtent layerTubeInvalid;
  layerTubeInvalid.set(binR, 5., 100.);
  layerTubeInvalid.set(binZ, -100., 100.);
  cvb.volumeExtent = layerTubeInvalid;

  BOOST_CHECK_THROW(cvb(gctx),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(CylindricalVolumeBuilderR_) {
  GeometryContext gctx;

  CylindricalVolumeBuilder<binR> cvbr;

  auto layer0 = surfacesCylinder(8.4, 36., 0.145, 32., 0.5, 5., {16, 14});
  auto layer1 = surfacesCylinder(8.4, 36., 0.145, 72., 2., 5., {32, 14});

  InternalBlueprint ilb0(gctx, layer0, AllInternalSurfaces{4},
                         {{binZ, {5., 5.}}, {binR, {2., 2.}}}, "barrel_layer0");

  InternalBlueprint ilb1(gctx, layer1, AllInternalSurfaces{4},
                         {{binZ, {10., 10.}}, {binR, {2., 2.}}},
                         "barrel_layer1");

  // Test excpetion throwing
  // - no internal sturcture given
  BOOST_CHECK_THROW(cvbr(gctx), std::invalid_argument);

  // - given extent is too small
  cvbr.volumeExtent.set(binR, 70., 71.);
  BOOST_CHECK_THROW(cvbr(gctx), std::invalid_argument);

  // No extent, only the volumes are built with one gap volume
  cvbr.internalBlueprints = {ilb0, ilb1};
  cvbr.name = "barrel_unbound";
  auto volumesInR = cvbr(gctx);
  BOOST_CHECK(volumesInR.size() == 3u);
  for (const auto& v : volumesInR) {
    display(*v, ViewConfig({0, 70, 0}));
  }

  // This should constrain only in R, but leave Z
  cvbr.volumeExtent = GeometricExtent{};
  cvbr.volumeExtent.set(binR, 5., 100.);
  cvbr.name = "barrel_r_extent";

  volumesInR = cvbr(gctx);
  BOOST_CHECK(volumesInR.size() == 5u);

  for (const auto& v : volumesInR) {
    display(*v, ViewConfig({0, 100, 0}));
  }
}

BOOST_AUTO_TEST_CASE(CylindricalVolumeBuilderZ_) {
  GeometryContext gctx;

  CylindricalVolumeBuilder<binZ> cvbz;
  auto disc0 = surfacesRing(4.4, 12.4, 32, 0., 44., -600., 1., 20);
  auto disc1 = surfacesRing(4.4, 12.4, 32, 0., 44., -550., 1., 20);
  auto disc2 = surfacesRing(4.4, 12.4, 32, 0., 44., -500., 1., 20);

  // This should constrain only in R, but leave Z
  GeometricExtent extentInRZ;
  extentInRZ.set(binR, 5., 100.);
  extentInRZ.set(binZ, -650, -490.);
  cvbz.volumeExtent = extentInRZ;

  InternalBlueprint ild0(gctx, disc0, AllInternalSurfaces{4},
                         {{binZ, {1., 1.}}, {binR, {2., 2.}}}, "barrel_disc0");

  InternalBlueprint ild1(gctx, disc1, AllInternalSurfaces{4},
                         {{binZ, {1., 1.}}, {binR, {2., 2.}}}, "barrel_disc1");

  InternalBlueprint ild2(gctx, disc2, AllInternalSurfaces{4},
                         {{binZ, {1., 1.}}, {binR, {2., 2.}}}, "barrel_disc2");

  cvbz.internalBlueprints = {ild0, ild1, ild2};
  cvbz.name = "endcap_rz_extent";

  auto volumesInZ = cvbz(gctx);
  BOOST_CHECK(volumesInZ.size() == 7u);

  for (const auto& v : volumesInZ) {
    display(*v, ViewConfig({0, 40, 0}));
  }
}

BOOST_AUTO_TEST_CASE(CylindricalVolumeBuilderPhi_) {
  GeometryContext gctx;

  CylindricalVolumeBuilder<binPhi> cvbp;
  auto phi0 = surfacesPhiSector(40., 50., 1.4, 10., 0.4, 0.1, 0.3, 3);
  auto phi1 = surfacesPhiSector(40., 50., 1.4, 10., 0.4, 0.6, 0.8, 3);
  auto phi2 = surfacesPhiSector(40., 50., 1.4, 10., 0.4, 1.1, 1.3, 3);

  InternalBlueprint ilp0(
      gctx, phi0, AllInternalSurfaces{6},
      {{binPhi, {0.02, 0.02}}, {binR, {2., 2.}}, {binZ, {5., 5.}}},
      "barrel_phi0");

  InternalBlueprint ilp1(
      gctx, phi1, AllInternalSurfaces{6},
      {{binPhi, {0.02, 0.02}}, {binR, {2., 2.}}, {binZ, {5., 5.}}},
      "barrel_phi1");

  InternalBlueprint ilp2(
      gctx, phi2, AllInternalSurfaces{6},
      {{binPhi, {0.02, 0.02}}, {binR, {2., 2.}}, {binZ, {5., 5.}}},
      "barrel_phi2");

  cvbp.internalBlueprints = {ilp0, ilp1, ilp2};
  cvbp.name = "sectoral_layer";

  auto volumesInPhi = cvbp(gctx);
  BOOST_CHECK(volumesInPhi.size() == 5u);

  for (const auto& v : volumesInPhi) {
    display(*v, ViewConfig({120, 120, 0}));
  }

  // Make a full sector of blueprints ... with boundary at -M_PI
  unsigned int nPhi = 10;
  std::vector<InternalBlueprint> ilps;
  ActsScalar phiStep = M_PI / nPhi;
  for (unsigned int iphi = 0; iphi < nPhi; ++iphi) {
    unsigned int cphi = 2 * iphi;
    auto phis =
        surfacesPhiSector(40., 50., 1.4, 10., 0.4, -M_PI + (cphi)*phiStep,
                          -M_PI + (cphi + 1) * phiStep, 3);
    ilps.push_back(InternalBlueprint(
        gctx, phis, AllInternalSurfaces{6},
        {{binPhi, {0.02, 0.02}}, {binR, {2., 2.}}, {binZ, {5., 5.}}},
        "barrel_sector"));
  }

  GeometricExtent fullPhiExtent;
  fullPhiExtent.set(binPhi, -M_PI, M_PI);

  cvbp.volumeExtent = fullPhiExtent;
  cvbp.internalBlueprints = ilps;
  cvbp.name = "sectoral_full_layer";

  volumesInPhi = cvbp(gctx);

  for (const auto& v : volumesInPhi) {
    display(*v, ViewConfig({180, 180, 0}));
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts
