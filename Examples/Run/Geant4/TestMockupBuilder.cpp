// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/Geant4Detector/Geant4Detector.hpp"
#include "ActsExamples/MuonSpectrometerMockupDetector/MockupSectorBuilder.hpp"

#include <memory>
#include <string>
#include <vector>

namespace Acts {
namespace Experimental {
class DetectorVolume;
}  // namespace Experimental
}  // namespace Acts

using namespace Acts;
using namespace ActsExamples;

struct StrawSelector {
  /// Call operator
  /// @param sf The input surface to be checked
  bool operator()(const Surface& sf) const {
    return (sf.type() == Surface::Straw);
  }
};

int main() {
  auto mockup_config = MockupSectorBuilder::Config();

  auto mockup_chamberConfig_inner = MockupSectorBuilder::ChamberConfig();
  auto mockup_chamberConfig_middle = MockupSectorBuilder::ChamberConfig();
  auto mockup_chamberConfig_outer = MockupSectorBuilder::ChamberConfig();

  mockup_config.gdmlPath =
      " ../../../../acts/Examples/Detectors/MuonSpectrometerMockupDetector/"
      "MuonChamber.gdml";
  mockup_config.NumberOfSectors = 8;

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

  GeometryContext gctx;
  auto detectorVolume_inner_chamber =
      mockup_builder.buildChamber(mockup_chamberConfig_inner);

  auto detectorVolume_middle_chamber =
      mockup_builder.buildChamber(mockup_chamberConfig_middle);

  auto detectorVolume_outer_chamber =
      mockup_builder.buildChamber(mockup_chamberConfig_outer);

  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
      detector_volumes = {detectorVolume_inner_chamber,
                          detectorVolume_middle_chamber,
                          detectorVolume_outer_chamber};

  auto detectorVolume_sector = mockup_builder.buildSector(detector_volumes);

  mockup_builder.drawSector(detectorVolume_sector, "sector_test.obj");

  auto detector_sector = Acts::Experimental::Detector::makeShared(
      "Detector", {detectorVolume_sector},
      Acts::Experimental::tryRootVolumes());
}
