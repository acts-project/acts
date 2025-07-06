// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/MultiWireStructureBuilder.hpp"
#include "Acts/Geometry/MultiWireVolumeBuilder.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Navigation/NavigationStateFillers.hpp"
#include "Acts/Navigation/NavigationStateUpdaters.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <fstream>
#include <memory>
#include <numbers>
#include <string>
#include <vector>

using namespace Acts;
using namespace Acts::Experimental;
using namespace Acts::detail;
using namespace Acts::Test;

GeometryContext tContext;
constexpr std::size_t nSurfacesX = 15;
constexpr std::size_t nSurfacesY = 4;
constexpr double radius = 15.0;
constexpr double halfZ = 250.0;

BOOST_AUTO_TEST_SUITE(Experimental)
auto logger = getDefaultLogger("MultiWireNavigationTests", Logging::VERBOSE);

// a function that constructs and returns detector elements for straw surfaces
void generateStrawSurfaces(
    std::size_t nStraws, std::size_t nLayers, double radius, double halfZ,
    std::vector<std::shared_ptr<Surface>>& strawSurfaces,
    std::vector<std::unique_ptr<DetectorElementStub>>& elements) {
  // The transform of the 1st surface
  Vector3 ipos = {-0.5 * nStraws * 2 * radius + radius,
                  -0.5 * nLayers * 2 * radius + radius, 0.};

  Vector3 pos = ipos;
  auto lineBounds = std::make_shared<LineBounds>(radius, halfZ);
  int id = 1;

  // Generate the surfaces
  for (std::size_t i = 0; i < nLayers; i++) {
    for (std::size_t j = 0; j < nStraws; j++) {
      pos.x() = ipos.x() + 2 * j * radius;

      auto& element =
          elements.emplace_back(std::make_unique<DetectorElementStub>(
              Transform3(Translation3(pos)), lineBounds, 0));

      element->surface().assignGeometryId(GeometryIdentifier(id++));

      element->surface().assignDetectorElement(*element);

      strawSurfaces.push_back(element->surface().getSharedPtr());
    }

    pos.y() = ipos.y() + 2 * (i + 1) * radius;
  }
}

BOOST_AUTO_TEST_CASE(Navigation_in_Indexed_Surfaces) {
  std::vector<std::shared_ptr<Surface>> strawSurfaces = {};
  std::vector<std::unique_ptr<DetectorElementStub>> detElements = {};

  generateStrawSurfaces(nSurfacesX, nSurfacesY, radius, halfZ, strawSurfaces,
                        detElements);

  std::vector<double> vBounds = {0.5 * nSurfacesX * 2 * radius,
                                 0.5 * nSurfacesX * 2 * radius,
                                 0.5 * nSurfacesY * 2 * radius, halfZ};

  MultiWireStructureBuilder::Config mlCfg;
  mlCfg.name = "Multi_Layer_With_Wires";
  mlCfg.mlSurfaces = strawSurfaces;

  mlCfg.mlBinning = {
      {DirectedProtoAxis(AxisDirection::AxisX, AxisBoundaryType::Bound,
                         -vBounds[0], vBounds[0], nSurfacesX),
       1u},
      {DirectedProtoAxis(AxisDirection::AxisY, AxisBoundaryType::Bound,
                         -vBounds[2], vBounds[2], nSurfacesY),
       0u}};
  auto boundsPtr = std::make_unique<TrapezoidVolumeBounds>(
      vBounds[0], vBounds[1], vBounds[2], vBounds[3]);
  mlCfg.mlBounds = boundsPtr->values();

  MultiWireStructureBuilder mlBuilder(mlCfg);
  auto [volumes, portals, roots] = mlBuilder.construct(tContext);

  Acts::Experimental::NavigationState nState;
  nState.position = Acts::Vector3(0., -59., 0.);
  nState.direction = Acts::Vector3(0., 1., 0.);

  nState.currentVolume = volumes.front().get();
  nState.currentVolume->updateNavigationState(tContext, nState);

  // check the surface candidates after update (12 surfaces + 6 portals but only
  // 4 surfaces are reachable (one of each layer and one portal)
  BOOST_CHECK_EQUAL(nState.surfaceCandidates.size(), 5u);
}

// This tests the multilayer navigation policy for gen3 geometry interface
BOOST_AUTO_TEST_CASE(MultiLayer_NavigationPolicy) {
  // Create the surfaces
  std::vector<std::shared_ptr<Surface>> strawSurfaces = {};
  std::vector<std::unique_ptr<DetectorElementStub>> detElements = {};

  generateStrawSurfaces(nSurfacesX, nSurfacesY, radius, halfZ, strawSurfaces,
                        detElements);

  std::vector<double> vBounds = {0.5 * nSurfacesX * 2 * radius,
                                 0.5 * nSurfacesX * 2 * radius,
                                 0.5 * nSurfacesY * 2 * radius, halfZ};
  // Create the multi wire volume (tracking volume this time!)
  MultiWireVolumeBuilder::Config mwCfg;
  mwCfg.name = "MultiWireVolume";
  mwCfg.mlSurfaces = strawSurfaces;
  mwCfg.binning = {
      {DirectedProtoAxis(AxisDirection::AxisX, AxisBoundaryType::Bound,
                         -vBounds[0], vBounds[0], nSurfacesX),
       1u},
      {DirectedProtoAxis(AxisDirection::AxisY, AxisBoundaryType::Bound,
                         -vBounds[2], vBounds[2], nSurfacesY),
       0u}};
  auto boundsPtr = std::make_shared<Acts::TrapezoidVolumeBounds>(
      vBounds[0], vBounds[1], vBounds[2], vBounds[3]);
  mwCfg.bounds = boundsPtr;
  mwCfg.transform = Transform3(Translation3(Vector3(0., 0., 0.)));

  // Build the volume
  MultiWireVolumeBuilder mwBuilder(mwCfg);
  std::unique_ptr<Acts::TrackingVolume> volume =
      mwBuilder.buildVolume(tContext);

  // Check the volume
  // we do not except any children volumes
  BOOST_CHECK(volume->volumes().empty());

  // we expect 15*4 = 60 surfaces
  BOOST_CHECK(volume->surfaces().size() == 60u);

  BOOST_CHECK(volume->portals().size() == 6u);

  // check the navigation policy
  NavigationStream main;
  AppendOnlyNavigationStream stream{main};
  Vector3 startPos = {0., -59., 0.};
  Vector3 startDir = {0., 1., 0.};
  NavigationArguments args{startPos, startDir};

  volume->initializeNavigationCandidates(args, stream, *logger);

  // we expect 18 candidates (12 surfaces + 6 portals)
  BOOST_CHECK_EQUAL(main.candidates().size(), 18u);

  // check if we have duplicated surface candidates
  auto it = std::unique(main.candidates().begin(), main.candidates().end(),
                        [](const auto& lhs, const auto& rhs) {
                          return lhs.surface() == rhs.surface();
                        });
  BOOST_CHECK(it == main.candidates().end());

  // try with a different direction
  double angle = std::numbers::pi / 4.;
  startDir = {std::cos(angle), std::sin(angle), 0.};
  args.direction = startDir;
  // clear the candidates and re initialize with new arguments
  main.candidates().clear();
  volume->initializeNavigationCandidates(args, stream, *logger);
  // we expect 18 candidates (12 surfaces + 6 portals)
  BOOST_CHECK_EQUAL(main.candidates().size(), 18u);
}

BOOST_AUTO_TEST_SUITE_END()
