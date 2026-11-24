// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "ActsPlugins/Root/RootMaterialTrackIo.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <stdexcept>

#include "TTree.h"

using namespace Acts;
using namespace Acts::UnitLiterals;
using namespace ActsPlugins;

namespace {

/// A track with two steps through silicon along the x axis
RecordedMaterialTrack makeTrack() {
  RecordedMaterialTrack track;
  track.first = {Vector3(0., 0., 0.), Vector3(1., 0., 0.)};
  Material silicon = Material::fromMassDensity(93.7_mm, 465.2_mm, 28.0855, 14.,
                                               2.329_g / 1_cm3);
  for (double x : {10_mm, 20_mm}) {
    MaterialInteraction interaction;
    interaction.position = Vector3(x, 0., 0.);
    interaction.direction = Vector3(1., 0., 0.);
    interaction.materialSlab = MaterialSlab(silicon, 1_mm);
    interaction.pathCorrection = 1_mm;
    track.second.materialInteractions.push_back(interaction);
  }
  return track;
}

}  // namespace

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(RootSuite)

BOOST_AUTO_TEST_CASE(RootMaterialTrackIoWriteWithoutContext) {
  RootMaterialTrackIo::Config cfg;
  cfg.prePostStepInfo = true;
  cfg.recalculateTotals = true;

  TTree tree("material_tracks", "material_tracks");
  tree.SetDirectory(nullptr);
  RootMaterialTrackIo writer(cfg);
  writer.connectForWrite(tree);

  const RecordedMaterialTrack track = makeTrack();
  writer.write(42, track);
  tree.Fill();
  BOOST_CHECK_EQUAL(tree.GetEntries(), 1);

  RootMaterialTrackIo reader(cfg);
  reader.connectForRead(tree);
  tree.GetEntry(0);
  BOOST_CHECK_EQUAL(reader.eventId(), 42u);

  RecordedMaterialTrack readBack = reader.read();
  CHECK_CLOSE_ABS(readBack.first.second, track.first.second, 1e-6);
  BOOST_REQUIRE_EQUAL(readBack.second.materialInteractions.size(), 2u);
  for (std::size_t i = 0; i < 2; ++i) {
    const auto& written = track.second.materialInteractions[i];
    const auto& read = readBack.second.materialInteractions[i];
    CHECK_CLOSE_ABS(read.position, written.position, 1e-4);
    CHECK_CLOSE_REL(read.materialSlab.thickness(),
                    written.materialSlab.thickness(), 1e-6);
    CHECK_CLOSE_REL(read.materialSlab.material().X0(),
                    written.materialSlab.material().X0(), 1e-6);
  }
  // The totals are recalculated from the steps
  CHECK_CLOSE_REL(readBack.second.materialInX0, 2_mm / 93.7_mm, 1e-5);
}

BOOST_AUTO_TEST_CASE(RootMaterialTrackIoWriteWithoutContextSurfaceInfo) {
  RootMaterialTrackIo::Config cfg;
  cfg.surfaceInfo = true;

  TTree tree("material_tracks", "material_tracks");
  tree.SetDirectory(nullptr);
  RootMaterialTrackIo writer(cfg);
  writer.connectForWrite(tree);

  // Surface information needs a geometry context
  BOOST_CHECK_THROW(writer.write(0, makeTrack()), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
