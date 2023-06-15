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
#include "Acts/Detector/DetectorBuilder.hpp"
#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Detector/interface/IDetectorComponentBuilder.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <stdexcept>

class CompBuilder final : public Acts::Experimental::IDetectorComponentBuilder {
 public:
  Acts::Experimental::DetectorComponent construct(
      const Acts::GeometryContext& gctx) const final {
    auto bounds = std::make_unique<Acts::CuboidVolumeBounds>(10., 10., 10.);
    // Construct the DetectorVolume
    auto dVolume = Acts::Experimental::DetectorVolumeFactory::construct(
        Acts::Experimental::defaultPortalGenerator(), gctx, "TestVolume",
        Acts::Transform3::Identity(), std::move(bounds),
        Acts::Experimental::tryAllPortals());

    // Fill the portal container
    Acts::Experimental::DetectorComponent::PortalContainer portalContainer;
    for (auto [ip, p] : Acts::enumerate(dVolume->portalPtrs())) {
      portalContainer[ip] = p;
    }

    return Acts::Experimental::DetectorComponent{
        {dVolume},
        portalContainer,
        {{dVolume}, Acts::Experimental::tryRootVolumes()}};
  }
};

BOOST_AUTO_TEST_SUITE(Detector)

BOOST_AUTO_TEST_CASE(DetectorBuilder_Misconfigured) {
  // Detector builder
  Acts::Experimental::DetectorBuilder::Config dCfg;
  dCfg.auxilliary = "*** Test X * Misconfigued ***";
  dCfg.name = "EmptyCylinder";
  dCfg.builder = nullptr;

  BOOST_CHECK_THROW(auto a = Acts::Experimental::DetectorBuilder(dCfg),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(DetectorBuilder_test) {
  Acts::GeometryContext tContext;

  // Detector builder
  Acts::Experimental::DetectorBuilder::Config dCfg;
  dCfg.auxilliary = "*** Test : Detector ***";
  dCfg.name = "TestDetector";
  dCfg.builder = std::make_shared<CompBuilder>();

  BOOST_CHECK_NO_THROW(auto a = Acts::Experimental::DetectorBuilder(dCfg));

  auto detector = Acts::Experimental::DetectorBuilder(dCfg).construct(tContext);

  BOOST_CHECK_EQUAL(detector->name(), "TestDetector");
  BOOST_CHECK(detector->rootVolumes().size() == 1);
}

BOOST_AUTO_TEST_SUITE_END()
