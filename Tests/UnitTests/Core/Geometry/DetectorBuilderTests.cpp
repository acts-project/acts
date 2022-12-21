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

  Acts::ProtoDetector cylinderDetector;
  cylinderDetector.name = "cylinder-detector";
  cylinderDetector.worldVolume = cylinder;

  DetectorBuilder::Config dCfg{cylinderDetector, Acts::Logging::VERBOSE};
  DetectorBuilder dBuilder(dCfg);

  auto detector = dBuilder.construct(tContext);

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
  system.constituentVolumes = {nec, barrel, pec};
  system.constituentBinning = {
      Acts::BinningData(Acts::open, Acts::binZ, {0., 1.})};
  system.blockBuilder = ContainerBlockBuilder{system};

  Acts::ProtoDetector cylinderDetector;
  cylinderDetector.name = "cylinder-detector";
  cylinderDetector.worldVolume = system;

  Acts::Experimental::DetectorBuilder::Config dCfg{cylinderDetector,
                                                   Acts::Logging::VERBOSE};
  Acts::Experimental::DetectorBuilder dBuilder(dCfg);

  auto detector = dBuilder.construct(tContext);

  // Check that we have one volume
  BOOST_CHECK(detector->name() == "cylinder-detector");
  BOOST_CHECK(detector->volumes().size() == 3u);
}

BOOST_AUTO_TEST_SUITE_END()
