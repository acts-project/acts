// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Detector/Blueprint.hpp"
#include "Acts/Detector/CylindricalContainerBuilder.hpp"
#include "Acts/Detector/DetectorBuilder.hpp"
#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Detector/IndexedRootVolumeFinderBuilder.hpp"
#include "Acts/Detector/detail/BlueprintDrawer.hpp"
#include "Acts/Detector/detail/BlueprintHelper.hpp"
#include "Acts/Detector/interface/IInternalStructureBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Utilities/BinningData.hpp"

#include <fstream>

using namespace Acts;

namespace ActsTests {

template <typename surface_type>
class SurfaceBuilder : public Experimental::IInternalStructureBuilder {
 public:
  SurfaceBuilder(const Transform3& trf, double p0, double p1)
      : m_surface(Surface::makeShared<surface_type>(trf, p0, p1)) {}
  /// Conrstruct and return the internal structure creation
  ///
  /// @param gctx the geometry context at the creation of the internal structure
  ///
  /// @return a consistent set of detector volume internals
  Experimental::InternalStructure construct(
      [[maybe_unused]] const GeometryContext& gctx) const final {
    // Trivialities first: internal volumes
    std::vector<std::shared_ptr<Experimental::DetectorVolume>> internalVolumes =
        {};
    Experimental::ExternalNavigationDelegate internalVolumeUpdater =
        Experimental::tryNoVolumes();

    // Retrieve the layer surfaces
    Experimental::InternalNavigationDelegate internalCandidatesUpdater =
        Experimental::tryAllPortalsAndSurfaces();

    // Return the internal structure
    return Experimental::InternalStructure{{m_surface},
                                           internalVolumes,
                                           std::move(internalCandidatesUpdater),
                                           std::move(internalVolumeUpdater)};
  }

 private:
  std::shared_ptr<Surface> m_surface;
};

BOOST_AUTO_TEST_SUITE(DetectorSuite)

BOOST_AUTO_TEST_CASE(CylindricalDetectorFromBlueprintTest) {
  GeometryContext tContext;

  // This tests shows how to careate cylindrical detector from a detector
  // blueprint.
  //
  // In general, the blueprint (lines below) is generated through reading in
  // or by parsing the geometry model (DD4heo, TGeo, Geant4, etc.). For
  // testing purpose, let us create the blueprint manually.
  //

  // Blueprint starts here ----------------

  // Detector dimensions
  double detectorIr = 0.;
  double detectorOr = 120.;
  double detectorHz = 400.;

  // Beam pipe
  double beamPipeOr = 20.;

  // Pixel system
  double pixelIr = 25;
  double pixelOr = 115;
  double pixelEcHz = 50;
  double pixelEcLayerHz = 10;

  // Create  root node
  std::vector<AxisDirection> detectorBinning = {AxisDirection::AxisR};
  std::vector<double> detectorBoundaries = {detectorIr, detectorOr, detectorHz};

  // The root node - detector
  auto detectorBpr = std::make_unique<Experimental::Gen2Blueprint::Node>(
      "detector", Transform3::Identity(), VolumeBounds::eCylinder,
      detectorBoundaries, detectorBinning);

  // The beam pipe
  std::vector<double> beamPipeBoundaries = {detectorIr, beamPipeOr, detectorHz};

  auto beamPipeStructure = std::make_shared<SurfaceBuilder<CylinderSurface>>(
      Transform3::Identity(), 18, 0.99 * detectorHz);
  auto beamPipe = std::make_unique<Experimental::Gen2Blueprint::Node>(
      "beam_pipe", Transform3::Identity(), VolumeBounds::eCylinder,
      beamPipeBoundaries, beamPipeStructure);
  detectorBpr->add(std::move(beamPipe));

  // A pixel system
  std::vector<double> pixelBoundaries = {pixelIr, pixelOr, detectorHz};
  std::vector<AxisDirection> pixelBinning = {AxisDirection::AxisZ};
  auto pixel = std::make_unique<Experimental::Gen2Blueprint::Node>(
      "pixel", Transform3::Identity(), VolumeBounds::eCylinder, pixelBoundaries,
      pixelBinning);

  // Nec: Small differences to check if the adjustments are made
  std::vector<double> pixelEcBoundaries = {pixelIr, pixelOr - 5., pixelEcHz};
  std::vector<AxisDirection> pixelEcBinning = {AxisDirection::AxisZ};

  Transform3 pixelNecTransform =
      Transform3::Identity() * Translation3(0., 0., -detectorHz + pixelEcHz);

  auto pixelNec = std::make_unique<Experimental::Gen2Blueprint::Node>(
      "pixel_nec", pixelNecTransform, VolumeBounds::eCylinder,
      pixelEcBoundaries, pixelEcBinning);

  // Add a single encap layer
  std::vector<double> pixelNecBoundaries = {pixelIr + 2, pixelOr - 7.,
                                            pixelEcLayerHz};

  auto pixelNecLayerStructure = std::make_shared<SurfaceBuilder<DiscSurface>>(
      pixelNecTransform, pixelIr + 10., pixelOr - 10.);

  auto pixelNecLayer = std::make_unique<Experimental::Gen2Blueprint::Node>(
      "pixel_nec_layer", pixelNecTransform, VolumeBounds::eCylinder,
      pixelNecBoundaries, pixelNecLayerStructure);

  pixelNec->add(std::move(pixelNecLayer));

  // Barrel
  std::vector<double> pixelBarrelBoundaries = {pixelIr + 1, pixelOr - 1.,
                                               detectorHz - 2 * pixelEcHz};
  std::vector<AxisDirection> pixelBarrelBinning = {AxisDirection::AxisR};

  auto pixelBarrel = std::make_unique<Experimental::Gen2Blueprint::Node>(
      "pixel_barrel", Transform3::Identity(), VolumeBounds::eCylinder,
      pixelBarrelBoundaries, pixelBarrelBinning);

  auto pixelBarrelL0Structure =
      std::make_shared<SurfaceBuilder<CylinderSurface>>(
          Transform3::Identity(), 62.5, detectorHz - 2 * pixelEcHz - 10.);
  std::vector<double> pixelBarrelL0Boundaries = {60, 65.,
                                                 detectorHz - 2 * pixelEcHz};
  auto pixelBarrelL0 = std::make_unique<Experimental::Gen2Blueprint::Node>(
      "pixel_barrel_l0", Transform3::Identity(), VolumeBounds::eCylinder,
      pixelBarrelL0Boundaries, pixelBarrelL0Structure);

  auto pixelBarrelL1Structure =
      std::make_shared<SurfaceBuilder<CylinderSurface>>(
          Transform3::Identity(), 102.5, detectorHz - 2 * pixelEcHz - 10.);

  std::vector<double> pixelBarrelL1Boundaries = {100, 105.,
                                                 detectorHz - 2 * pixelEcHz};
  auto pixelBarrelL1 = std::make_unique<Experimental::Gen2Blueprint::Node>(
      "pixel_barrel_l1", Transform3::Identity(), VolumeBounds::eCylinder,
      pixelBarrelL1Boundaries, pixelBarrelL1Structure);
  pixelBarrel->add(std::move(pixelBarrelL0));
  pixelBarrel->add(std::move(pixelBarrelL1));

  Transform3 pixelPecTransform =
      Transform3::Identity() * Translation3(0., 0., detectorHz - pixelEcHz);

  auto pixelPec = std::make_unique<Experimental::Gen2Blueprint::Node>(
      "pixel_pec", pixelPecTransform, VolumeBounds::eCylinder,
      pixelEcBoundaries, pixelEcBinning);

  std::vector<double> pixelPecBoundaries = {pixelIr + 2, pixelOr - 7., 10.};

  auto pixelPecLayerStructure = std::make_shared<SurfaceBuilder<DiscSurface>>(
      pixelPecTransform, pixelIr + 10., pixelOr - 10.);

  auto pixelPecLayer = std::make_unique<Experimental::Gen2Blueprint::Node>(
      "pixel_pec_layer", pixelPecTransform, VolumeBounds::eCylinder,
      pixelPecBoundaries, pixelPecLayerStructure);

  pixelPec->add(std::move(pixelPecLayer));

  // Adding pixel
  pixel->add(std::move(pixelNec));
  pixel->add(std::move(pixelPec));
  pixel->add(std::move(pixelBarrel));

  detectorBpr->add(std::move(pixel));

  // An Indexed volume finder will be attached
  std::vector<AxisDirection> rootVolumeBinning = {AxisDirection::AxisZ,
                                                  AxisDirection::AxisR};
  detectorBpr->rootVolumeFinderBuilder =
      std::make_shared<Experimental::IndexedRootVolumeFinderBuilder>(
          rootVolumeBinning);

  // A geo ID generator
  detectorBpr->geoIdGenerator =
      std::make_shared<Experimental::GeometryIdGenerator>(
          Experimental::GeometryIdGenerator::Config{},
          getDefaultLogger("RecursiveIdGenerator", Logging::VERBOSE));

  // Complete and fill gaps
  Experimental::detail::BlueprintHelper::fillGaps(*detectorBpr);

  std::fstream fs("cylindrical_detector_blueprint.dot", std::ios::out);
  Experimental::detail::BlueprintDrawer::dotStream(fs, *detectorBpr);
  fs.close();

  // ----------------------------- end of blueprint

  // Create a Cylindrical detector builder from this blueprint
  auto detectorBuilder =
      std::make_shared<Experimental::CylindricalContainerBuilder>(
          *detectorBpr, Logging::VERBOSE);

  // Detector builder
  Experimental::DetectorBuilder::Config dCfg;
  dCfg.auxiliary =
      "*** Test : auto generated cylindrical detector builder  ***";
  dCfg.name = "Cylindrical detector from blueprint";
  dCfg.builder = detectorBuilder;
  dCfg.geoIdGenerator = detectorBpr->geoIdGenerator;

  auto detector = Experimental::DetectorBuilder(dCfg).construct(tContext);

  BOOST_REQUIRE_NE(detector, nullptr);

  // There should be 14 volumes, and they should be built in order
  // beam_pipe
  // detector_gap_0
  // pixel_nec_gap_0
  // pixel_nec_layer
  // pixel_nec_gap_1
  // pixel_barrel_gap_0
  // pixel_barrel_l0
  // pixel_barrel_gap_1
  // pixel_barrel_l1
  // pixel_barrel_gap_2
  // pixel_pec_gap_0
  // pixel_pec_layer
  // pixel_pec_gap_1
  // detector_gap_1
  BOOST_CHECK_EQUAL(detector->volumes().size(), 14u);
  BOOST_CHECK_EQUAL(detector->volumes()[0]->name(), "beam_pipe");
  BOOST_CHECK_EQUAL(detector->volumes()[1]->name(), "detector_gap_0");
  BOOST_CHECK_EQUAL(detector->volumes()[2]->name(), "pixel_nec_gap_0");
  BOOST_CHECK_EQUAL(detector->volumes()[3]->name(), "pixel_nec_layer");
  BOOST_CHECK_EQUAL(detector->volumes()[4]->name(), "pixel_nec_gap_1");
  BOOST_CHECK_EQUAL(detector->volumes()[5]->name(), "pixel_barrel_gap_0");
  BOOST_CHECK_EQUAL(detector->volumes()[6]->name(), "pixel_barrel_l0");
  BOOST_CHECK_EQUAL(detector->volumes()[7]->name(), "pixel_barrel_gap_1");
  BOOST_CHECK_EQUAL(detector->volumes()[8]->name(), "pixel_barrel_l1");
  BOOST_CHECK_EQUAL(detector->volumes()[9]->name(), "pixel_barrel_gap_2");
  BOOST_CHECK_EQUAL(detector->volumes()[10]->name(), "pixel_pec_gap_0");
  BOOST_CHECK_EQUAL(detector->volumes()[11]->name(), "pixel_pec_layer");
  BOOST_CHECK_EQUAL(detector->volumes()[12]->name(), "pixel_pec_gap_1");
  BOOST_CHECK_EQUAL(detector->volumes()[13]->name(), "detector_gap_1");
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
