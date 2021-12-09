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
#include "Acts/Experimental/CylindricalVolumeBuilders.hpp"
#include "Acts/Experimental/DetectorVolume.hpp"
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

  CylindricalVolumeBuilder cvb;

  // Check invalid input exceptions
  Extent unbound(false);
  BOOST_CHECK_THROW(cvb(gctx, {}, unbound, "unbound"), std::invalid_argument);

  // Check exception if more than one internal blue print is given
  std::vector<InternalBlueprint> ibps = {InternalBlueprint(gctx),
                                         InternalBlueprint(gctx)};
  BOOST_CHECK_THROW(cvb(gctx, ibps, Extent(), "blueprints"),
                    std::invalid_argument);

  // Volumes without surfaces

  // A simple solid
  Extent solidExtent;
  solidExtent.ranges[binR].second = 10.;
  solidExtent.ranges[binZ].first = -10.;
  solidExtent.ranges[binZ].second = 10.;

  auto solids = cvb(gctx, {}, solidExtent, "solid");
  BOOST_CHECK(solids.size() == 1);
  auto solid = solids[0];
  BOOST_CHECK(solid->name() == "solid");
  display(*solid, ViewConfig({0, 0, 100}));

  // A shifted solid
  Extent shiftedSolidExtent;
  shiftedSolidExtent.ranges[binR].second = 10.;
  shiftedSolidExtent.ranges[binZ].first = -40.;
  shiftedSolidExtent.ranges[binZ].second = -20.;

  auto shiftedSolids = cvb(gctx, {}, shiftedSolidExtent, "shifed_solid");
  BOOST_CHECK(shiftedSolids.size() == 1);
  auto shiftedSolid = shiftedSolids[0];
  BOOST_CHECK(shiftedSolid->name() == "shifed_solid");
  display(*shiftedSolid, ViewConfig({0, 0, 40}));

  // A shifted tube
  Extent shiftedTubeExtent;
  shiftedTubeExtent.ranges[binR].first = 2.;
  shiftedTubeExtent.ranges[binR].second = 10.;
  shiftedTubeExtent.ranges[binZ].first = 20.;
  shiftedTubeExtent.ranges[binZ].second = 40.;

  auto shiftedTubes = cvb(gctx, {}, shiftedTubeExtent, "shifed_tube");
  BOOST_CHECK(shiftedTubes.size() == 1);
  auto shiftedTube = shiftedTubes[0];
  BOOST_CHECK(shiftedTube->name() == "shifed_tube");
  display(*shiftedTube, ViewConfig({0, 0, 160}));

  // A shifted sectoral tube
  Extent shiftedTubeSectorExtent;

  shiftedTubeSectorExtent.ranges[binR].first = 2.;
  shiftedTubeSectorExtent.ranges[binR].second = 10.;
  shiftedTubeSectorExtent.ranges[binZ].first = 50.;
  shiftedTubeSectorExtent.ranges[binZ].second = 70.;
  shiftedTubeSectorExtent.ranges[binPhi].first = -0.85;
  shiftedTubeSectorExtent.ranges[binPhi].second = 0.85;

  auto shiftedTubeSectors =
      cvb(gctx, {}, shiftedTubeSectorExtent, "shifed_tube_sector");
  BOOST_CHECK(shiftedTubeSectors.size() == 1);
  auto shiftedTubeSector = shiftedTubeSectors[0];
  BOOST_CHECK(shiftedTubeSector->name() == "shifed_tube_sector");
  display(*shiftedTubeSector, ViewConfig({0, 0, 220}));

  
  // Volume without surfaces
  auto layer = surfacesCylinder(8.4, 36., 0.145, 32., 0.5, 5., {16, 14});
  InternalBlueprint ilb(gctx, layer);
  auto barrelLayerTightVec =
      cvb(gctx, {ilb}, Extent{false}, "barrel_layer_tight");

  BOOST_CHECK(barrelLayerTightVec.size() == 1);
  auto barrelLayerTight = barrelLayerTightVec[0];

  display(*barrelLayerTight, ViewConfig({0, 0, 220}), ViewConfig({100, 0, 0}));

  /// Build with envelope
  auto envelope = std::vector<std::pair<ActsScalar, ActsScalar>>(
      size_t(binValues), {0., 0.});
  envelope[Acts::binR] = { 5., 5.};

  InternalBlueprint ilbEnv(gctx, layer, AllSurfaces{}, Extent{}, envelope);

  auto barrelLayerEnvVec =
      cvb(gctx, {ilbEnv}, Extent{false}, "barrel_layer_envelope");

  BOOST_CHECK(barrelLayerEnvVec.size() == 1);
  auto barrelLayerEnv = barrelLayerEnvVec[0];

  display(*barrelLayerEnv, ViewConfig({0, 0, 220}), ViewConfig({100, 0, 0}));


}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts
