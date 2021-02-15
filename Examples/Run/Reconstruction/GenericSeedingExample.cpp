// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/GenericDetector/GenericDetector.hpp"

#include "SeedingExample.hpp"

int main(int argc, char* argv[]) {
  /// selection must not have duplicates.
  std::vector<Acts::GeometryIdentifier> layersForSeeding = {
      // generic detector barrel layers
      // the selection intentionally contains duplicates to demonstrate the
      // automatic selection normalization. setting only the volume already
      // selects all layers within it. the explicit layers in the selection
      // should have no effect.
      Acts::GeometryIdentifier().setVolume(8),
      Acts::GeometryIdentifier().setVolume(8).setLayer(2),
      Acts::GeometryIdentifier().setVolume(8).setLayer(4),
      Acts::GeometryIdentifier().setVolume(8).setLayer(6),
      Acts::GeometryIdentifier().setVolume(8).setLayer(8),
      // generic detector positive endcap layers
      Acts::GeometryIdentifier().setVolume(9).setLayer(2),
      Acts::GeometryIdentifier().setVolume(9).setLayer(4),
      Acts::GeometryIdentifier().setVolume(9).setLayer(6),
      Acts::GeometryIdentifier().setVolume(9).setLayer(8),
      // generic detector negative endcap layers
      Acts::GeometryIdentifier().setVolume(7).setLayer(14),
      Acts::GeometryIdentifier().setVolume(7).setLayer(12),
      Acts::GeometryIdentifier().setVolume(7).setLayer(10),
      Acts::GeometryIdentifier().setVolume(7).setLayer(8),
  };

  return runSeedingExample(argc, argv, std::make_shared<GenericDetector>(),
                           layersForSeeding);
}
