// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/Detector.hpp"
#include "Acts/Geometry/DetectorBuilder.hpp"
#include "Acts/Geometry/DetectorVolume.hpp"
#include "Acts/Geometry/ProtoDetector.hpp"
#include "Acts/Geometry/ProtoVolumeConverter.hpp"
#include "Acts/Utilities/Logger.hpp"

Acts::GeometryContext tContext;

BOOST_AUTO_TEST_SUITE(Experimental)

BOOST_AUTO_TEST_CASE(DetectorBuilderSimpleVolume) {
  using namespace Acts::Experimental;

  Acts::ProtoVolume cylinder;
  cylinder.name = "cylinder";
  cylinder.extent.set(Acts::binR, 0., 200.);
  cylinder.extent.set(Acts::binZ, -200., 200.);
  cylinder.blockBuilder = SingleBlockBuilder<>{cylinder};

  Acts::ProtoDetector protoDetector;
  protoDetector.name = "cylinder-detector";
  protoDetector.worldVolume = cylinder;

  DetectorBuilder::Config dCfg{Acts::Logging::VERBOSE};
  DetectorBuilder dBuilder(dCfg);

  auto detector = dBuilder.construct(tContext, protoDetector);

  // Check that we have one volume
  BOOST_CHECK(detector->name() == "cylinder-detector");
  BOOST_CHECK(detector->volumes().size() == 1u);
}

BOOST_AUTO_TEST_CASE(DetectorBuilderContainerVolume) {
  using namespace Acts::Experimental;

  Acts::ProtoVolume nec;
  nec.name = "nec";
  nec.extent.set(Acts::binR, 0., 200.);
  nec.extent.set(Acts::binZ, -600, -200.);
  nec.blockBuilder = SingleBlockBuilder<>{nec};

  Acts::ProtoVolume barrel;
  barrel.name = "barrel";
  barrel.extent.set(Acts::binR, 0., 200.);
  barrel.extent.set(Acts::binZ, -200., 200.);
  barrel.blockBuilder = SingleBlockBuilder<>{barrel};

  Acts::ProtoVolume pec;
  pec.name = "pec";
  pec.extent.set(Acts::binR, 0., 200.);
  pec.extent.set(Acts::binZ, 200, 600.);
  pec.blockBuilder = SingleBlockBuilder<>{pec};

  Acts::ProtoVolume system;
  system.name = "nec-barrel-pec";

  system.container = Acts::ProtoVolume::ContainerStructure{
      {nec, barrel, pec},
      {Acts::BinningData(Acts::open, Acts::binZ, {0., 1.})}};

  system.blockBuilder = ContainerBlockBuilder{system};

  Acts::ProtoDetector protoDetector;
  protoDetector.name = "cylinder-detector";
  protoDetector.worldVolume = system;

  Acts::Experimental::DetectorBuilder::Config dCfg{Acts::Logging::VERBOSE};
  Acts::Experimental::DetectorBuilder dBuilder(dCfg);

  auto detector = dBuilder.construct(tContext, protoDetector);

  // Check that we have one volume
  BOOST_CHECK(detector->name() == "cylinder-detector");
  BOOST_CHECK(detector->volumes().size() == 3u);

  Acts::ProtoVolume beamPipe;
  beamPipe.name = "beampipe";
  beamPipe.extent.set(Acts::binR, 0., 20.);
  beamPipe.extent.set(Acts::binZ, -600, 600.);
  beamPipe.blockBuilder = SingleBlockBuilder<>{beamPipe};
  // Slightly changed with beam pipe
  nec.extent.set(Acts::binR, 20., 200.);
  nec.blockBuilder = SingleBlockBuilder<>{nec};
  barrel.extent.set(Acts::binR, 20., 200.);
  barrel.blockBuilder = SingleBlockBuilder<>{barrel};
  pec.extent.set(Acts::binR, 20., 200.);
  pec.blockBuilder = SingleBlockBuilder<>{pec};
  // Re-assign the system
  system.container = Acts::ProtoVolume::ContainerStructure{
      {nec, barrel, pec},
      {Acts::BinningData(Acts::open, Acts::binZ, {0., 1.})}};
  system.blockBuilder = ContainerBlockBuilder{system};
  // System with beam pipe
  Acts::ProtoVolume bpsystem;
  bpsystem.name = "beampipe-nec-barrel-pec";

  bpsystem.container = Acts::ProtoVolume::ContainerStructure{
      {beamPipe, system},
      {Acts::BinningData(Acts::open, Acts::binR, {0., 1.})}};
  bpsystem.blockBuilder = ContainerBlockBuilder{bpsystem};

  protoDetector.name = "cylinder-detector-beam-pipe";
  protoDetector.worldVolume = bpsystem;

  detector = dBuilder.construct(tContext, protoDetector);
  BOOST_CHECK(detector->name() == "cylinder-detector-beam-pipe");
  BOOST_CHECK(detector->volumes().size() == 4u);
}

BOOST_AUTO_TEST_SUITE_END()