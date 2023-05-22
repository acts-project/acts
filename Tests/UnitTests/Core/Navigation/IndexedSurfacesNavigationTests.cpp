// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Detector/detail/IndexedGridFiller.hpp"
#include "Acts/Detector/detail/PortalHelper.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/MultiWireLayerUpdators.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Navigation/NavigationStateFillers.hpp"
#include "Acts/Navigation/NavigationStateUpdators.hpp"
#include "Acts/Navigation/NextNavigator.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/SurfaceCollector.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/IVisualization3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "ActsExamples/Geant4Detector/Geant4Detector.hpp"
#include "ActsExamples/MuonSpectrometerMockupDetector/MockupSectorBuilder.hpp"

#include <array>
#include <memory>
#include <vector>

using namespace Acts;
using namespace ActsExamples;
using namespace Acts::Experimental;

// A test context
Acts::GeometryContext tgContext;
Acts::MagneticFieldContext mfContext;

BOOST_AUTO_TEST_SUITE(Experimental)

// This test checks whether the surface candidate updator works as expected
// giving the expected number of surfaces as candidates

BOOST_AUTO_TEST_CASE(MultiWireLayerUpdator) {
  // create the mockup geometry for one sector

  auto mockup_config = MockupSectorBuilder::Config();

  auto mockup_chamberConfig_inner = MockupSectorBuilder::ChamberConfig();
  auto mockup_chamberConfig_middle = MockupSectorBuilder::ChamberConfig();
  auto mockup_chamberConfig_outer = MockupSectorBuilder::ChamberConfig();

  mockup_config.gdmlPath =
      " ../../../../acts/Examples/Detectors/MuonSpectrometerMockupDetector/"
      "MuonChamber.gdml";
  mockup_config.NumberOfSectors = 1;

  mockup_chamberConfig_inner.name = "Inner_Detector_Chamber";
  mockup_chamberConfig_inner.SensitiveNames = {"Inner_Skin"};
  mockup_chamberConfig_inner.PassiveNames = {"xx"};

  mockup_chamberConfig_middle.name = "Middle_Detector_Chamber";
  mockup_chamberConfig_middle.SensitiveNames = {"Middle_Skin"};
  mockup_chamberConfig_middle.PassiveNames = {"xx"};

  mockup_chamberConfig_outer.name = "Outer_Detector_Chamber";
  mockup_chamberConfig_outer.SensitiveNames = {"Outer_Skin"};
  mockup_chamberConfig_outer.PassiveNames = {"xx"};

  MockupSectorBuilder mockup_builder(mockup_config);

  auto detectorVolume_inner_chamber =
      mockup_builder.buildChamber(mockup_chamberConfig_inner);

  auto detectorVolume_middle_chamber =
      mockup_builder.buildChamber(mockup_chamberConfig_middle);

  auto detectorVolume_outer_chamber =
      mockup_builder.buildChamber(mockup_chamberConfig_outer);

  std::cout << "here" << std::endl;

  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
      detector_volumes = {};

  detector_volumes.push_back(detectorVolume_inner_chamber);
  detector_volumes.push_back(detectorVolume_middle_chamber);
  detector_volumes.push_back(detectorVolume_outer_chamber);

  auto detectorVolume_sector = mockup_builder.buildSector(detector_volumes);

  auto detector_sector = Acts::Experimental::Detector::makeShared(
      "Detector", {detectorVolume_sector},
      Acts::Experimental::tryRootVolumes());

  // Test if the surface candidates are as expected
  Acts::Experimental::NavigationState nState;
  nState.currentDetector = detector_sector.get();
  nState.position = Vector3(0., 0., 0.);
  nState.direction = Vector3(0., 1., 0.);

  // assign the detector volume of the multilayer of the first chamber to the
  // navigation state

  for (const auto& inner_vol : detector_volumes[0]->volumes()) {
    nState.position = {
        0., inner_vol->center().y() - inner_vol->volumeBounds().values()[1],
        0.};
    nState.currentVolume = inner_vol;
    nState.currentVolume->updateNavigationState(tgContext, nState);

    // the 12 straw surfaces (three for each layer) and the 6 protals of the
    // multiwire layer detector volume
    BOOST_CHECK(nState.surfaceCandidates.size() ==
                nState.currentVolume->portals().size() + 12u);
  }

  nState.surfaceCandidates.clear();

  for (const auto& inner_vol : detector_volumes[1]->volumes()) {
    nState.position = {
        0., inner_vol->center().y() - inner_vol->volumeBounds().values()[1],
        0.};
    nState.currentVolume = inner_vol;
    nState.currentVolume->updateNavigationState(tgContext, nState);

    // the 9 straw surfaces (three for each layer) and the 6 protals of the
    // multiwire layer detector volume
    BOOST_CHECK(nState.surfaceCandidates.size() ==
                nState.currentVolume->portals().size() + 9u);
  }

  nState.surfaceCandidates.clear();

  for (const auto& inner_vol : detector_volumes[2]->volumes()) {
    nState.position = {
        0., inner_vol->center().y() - inner_vol->volumeBounds().values()[1],
        0.};
    nState.currentVolume = inner_vol;
    nState.currentVolume->updateNavigationState(tgContext, nState);

    // the 9 straw surfaces (three for each layer) and the 6 protals of the
    // multiwire layer detector volume
    BOOST_CHECK(nState.surfaceCandidates.size() ==
                nState.currentVolume->portals().size() + 9u);
  }
}

BOOST_AUTO_TEST_SUITE_END()
