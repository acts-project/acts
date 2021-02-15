// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/TGeoDetector/TGeoDetector.hpp"

#include "SeedingExample.hpp"

int main(int argc, char* argv[]) {
  /// selection must not have duplicates.
  std::vector<Acts::GeometryIdentifier> layersForSeeding = {
      // ITk barrel layers
      // the selection intentionally contains duplicates to demonstrate the
      // automatic selection normalization. setting only the volume already
      // selects all layers within it. the explicit layers in the selection
      // should have no effect.
 // ITk barrel
 Acts::GeometryIdentifier().setVolume(9),  // L0, L1
 Acts::GeometryIdentifier().setVolume(16), // L2,L3,L4
 // ITk negative endcap
 Acts::GeometryIdentifier().setVolume(8),  // R0, R1
 Acts::GeometryIdentifier().setVolume(13), // R2,
 Acts::GeometryIdentifier().setVolume(14), // R3,
 Acts::GeometryIdentifier().setVolume(15), // R4,
 // ITk positive endcap
 Acts::GeometryIdentifier().setVolume(10),  // R0, R1
 Acts::GeometryIdentifier().setVolume(18), // R2,
 Acts::GeometryIdentifier().setVolume(19), // R3,
 Acts::GeometryIdentifier().setVolume(20), // R4,

  };

  return runSeedingExample(argc, argv, std::make_shared<TGeoDetector>(),
                           layersForSeeding);
}
