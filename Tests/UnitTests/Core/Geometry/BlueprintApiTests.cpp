// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/LayerBlueprintNode.hpp"
#include "Acts/Geometry/MaterialDesignatorBlueprintNode.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/PortalDesignatorBlueprintNode.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Geometry/VolumeResizeStrategy.hpp"
#include "Acts/Material/MergedMaterialMarker.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "ActsTests/CommonHelpers/DetectorElementStub.hpp"

#include <fstream>
#include <random>
#include <vector>

using namespace Acts;
using namespace UnitLiterals;

using Experimental::Blueprint;
using Experimental::LayerBlueprintNode;
using Experimental::MaterialDesignatorBlueprintNode;

namespace ActsTests {

auto logger = getDefaultLogger("UnitTests", Logging::DEBUG);

auto gctx = GeometryContext::dangerouslyDefaultConstruct();

inline std::vector<std::shared_ptr<Surface>> makeFanLayer(
    const Transform3& base,
    std::vector<std::unique_ptr<SurfacePlacementBase>>& elements,
    double r = 300_mm, std::size_t nSensors = 8, double thickness = 0) {
  auto recBounds = std::make_shared<RectangleBounds>(40_mm, 60_mm);

  double deltaPhi = 2 * std::numbers::pi / nSensors;
  std::vector<std::shared_ptr<Surface>> surfaces;
  for (std::size_t i = 0; i < nSensors; i++) {
    // Create a fan of sensors

    Transform3 trf = base * AngleAxis3{deltaPhi * i, Vector3::UnitZ()} *
                     Translation3(Vector3::UnitX() * r);

    if (i % 2 == 0) {
      trf = trf * Translation3{Vector3::UnitZ() * 5_mm};
    }

    auto& element = elements.emplace_back(
        std::make_unique<DetectorElementStub>(trf, recBounds, thickness));

    element->surface().assignSurfacePlacement(*element);

    surfaces.push_back(element->surface().getSharedPtr());
  }
  return surfaces;
}

inline std::vector<std::shared_ptr<Surface>> makeBarrelLayer(
    const Transform3& base,
    std::vector<std::unique_ptr<SurfacePlacementBase>>& elements,
    double r = 300_mm, std::size_t nStaves = 10, int nSensorsPerStave = 8,
    double thickness = 0, double hlPhi = 40_mm, double hlZ = 60_mm) {
  auto recBounds = std::make_shared<RectangleBounds>(hlPhi, hlZ);

  double deltaPhi = 2 * std::numbers::pi / nStaves;
  std::vector<std::shared_ptr<Surface>> surfaces;

  for (std::size_t istave = 0; istave < nStaves; istave++) {
    for (int isensor = -nSensorsPerStave; isensor <= nSensorsPerStave;
         isensor++) {
      double z = isensor * (2 * hlZ + 5_mm);

      Transform3 trf = base * Translation3(Vector3::UnitZ() * z) *
                       AngleAxis3{deltaPhi * istave, Vector3::UnitZ()} *
                       Translation3(Vector3::UnitX() * r) *
                       AngleAxis3{10_degree, Vector3::UnitZ()} *
                       AngleAxis3{90_degree, Vector3::UnitY()} *
                       AngleAxis3{90_degree, Vector3::UnitZ()};
      auto& element = elements.emplace_back(
          std::make_unique<DetectorElementStub>(trf, recBounds, thickness));
      element->surface().assignSurfacePlacement(*element);
      surfaces.push_back(element->surface().getSharedPtr());
    }
  }

  return surfaces;
}

}  // namespace ActsTests

using namespace ActsTests;

BOOST_AUTO_TEST_SUITE(GeometrySuite);

void pseudoNavigation(const TrackingGeometry& trackingGeometry,
                      Vector3 position, const Vector3& direction,
                      std::ostream& csv, std::size_t run,
                      std::size_t substepsPerCm, const Logger& logger) {
  ACTS_VERBOSE("start navigation " << run);
  ACTS_VERBOSE("dir: " << direction.transpose());
  ACTS_VERBOSE(direction.norm());

  std::mt19937 rng{static_cast<unsigned int>(run)};
  std::uniform_real_distribution<> dist{0.01, 0.99};

  const auto* volume = trackingGeometry.lowestTrackingVolume(gctx, position);
  BOOST_REQUIRE_NE(volume, nullptr);
  ACTS_VERBOSE(volume->volumeName());

  NavigationStream main;
  const TrackingVolume* currentVolume = volume;

  csv << run << "," << position[0] << "," << position[1] << "," << position[2];
  csv << "," << volume->geometryId().volume();
  csv << "," << volume->geometryId().boundary();
  csv << "," << volume->geometryId().sensitive();
  csv << std::endl;

  ACTS_VERBOSE("start pseudo navigation");

  for (std::size_t i = 0; i < 100; i++) {
    main = NavigationStream{};
    AppendOnlyNavigationStream stream{main};

    NavigationArguments navArgs{.position = position, .direction = direction};
    NavigationPolicyStateManager stateManager;
    currentVolume->navigationPolicy()->createState(gctx, navArgs, stateManager,
                                                   logger);
    auto policyState = stateManager.currentState();
    currentVolume->initializeNavigationCandidates(gctx, navArgs, policyState,
                                                  stream, logger);

    ACTS_VERBOSE(main.candidates().size() << " candidates");

    for (const auto& candidate : main.candidates()) {
      ACTS_VERBOSE(" -> " << candidate.surface().geometryId());
      ACTS_VERBOSE("    " << candidate.surface().toStream(gctx));
    }

    ACTS_VERBOSE("initializing candidates");
    main.initialize(gctx, {position, direction}, BoundaryTolerance::None());

    ACTS_VERBOSE(main.candidates().size() << " candidates remaining");

    for (const auto& candidate : main.candidates()) {
      ACTS_VERBOSE(" -> " << candidate.surface().geometryId());
      ACTS_VERBOSE("    " << candidate.surface().toStream(gctx));
    }

    if (main.currentCandidate().surface().isOnSurface(gctx, position,
                                                      direction)) {
      ACTS_VERBOSE("Already on surface at initialization, skipping candidate");

      auto id = main.currentCandidate().surface().geometryId();
      csv << run << "," << position[0] << "," << position[1] << ","
          << position[2];
      csv << "," << id.volume();
      csv << "," << id.boundary();
      csv << "," << id.sensitive();
      csv << std::endl;
      if (!main.switchToNextCandidate()) {
        ACTS_WARNING("candidates exhausted unexpectedly");
        break;
      }
    }

    auto writeIntersection = [&](const Vector3& pos, const Surface& surface) {
      csv << run << "," << pos[0] << "," << pos[1] << "," << pos[2];
      csv << "," << surface.geometryId().volume();
      csv << "," << surface.geometryId().boundary();
      csv << "," << surface.geometryId().sensitive();
      csv << std::endl;
    };

    bool terminated = false;
    while (main.remainingCandidates() > 0) {
      const auto& candidate = main.currentCandidate();

      ACTS_VERBOSE(candidate.position().transpose());

      ACTS_VERBOSE("moving to position: " << position.transpose() << " (r="
                                          << VectorHelpers::perp(position)
                                          << ")");

      Vector3 delta = candidate.position() - position;

      std::size_t substeps =
          std::max(1l, std::lround(delta.norm() / 10_cm * substepsPerCm));

      for (std::size_t j = 0; j < substeps; j++) {
        // position += delta / (substeps + 1);
        Vector3 subpos = position + dist(rng) * delta;
        csv << run << "," << subpos[0] << "," << subpos[1] << "," << subpos[2];
        csv << "," << currentVolume->geometryId().volume();
        csv << ",0,0";  // zero boundary and sensitive ids
        csv << std::endl;
      }

      position = candidate.position();
      ACTS_VERBOSE("                 -> "
                   << position.transpose()
                   << " (r=" << VectorHelpers::perp(position) << ")");

      writeIntersection(position, candidate.surface());

      if (candidate.isPortalTarget()) {
        ACTS_VERBOSE("On portal: " << candidate.surface().toStream(gctx));
        currentVolume =
            candidate.portal().resolveVolume(gctx, position, direction).value();

        if (currentVolume == nullptr) {
          ACTS_VERBOSE("switched to nullptr -> we're done");
          terminated = true;
        }
        break;

      } else {
        ACTS_VERBOSE("Not on portal");
      }

      main.switchToNextCandidate();
    }

    if (terminated) {
      ACTS_VERBOSE("Terminate pseudo navigation");
      break;
    }

    ACTS_VERBOSE("switched to " << currentVolume->volumeName());

    ACTS_VERBOSE("-----");
  }
}

BOOST_AUTO_TEST_CASE(NodeApiTestContainers) {
  // Transform3 base{AngleAxis3{30_degree, Vector3{1, 0, 0}}};
  Transform3 base{Transform3::Identity()};

  std::vector<std::unique_ptr<SurfacePlacementBase>> detectorElements;
  auto makeFan = [&](const Transform3& layerBase, auto&&..., double r,
                     std::size_t nSensors, double thickness) {
    return makeFanLayer(layerBase, detectorElements, r, nSensors, thickness);
  };

  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {0_mm, 20_mm};
  auto root = std::make_unique<Blueprint>(cfg);

  root->addMaterial("GlobalMaterial", [&](MaterialDesignatorBlueprintNode&
                                              mat) {
    using enum AxisDirection;
    using enum AxisBoundaryType;
    using enum CylinderVolumeBounds::Face;

    // Configure cylinder faces with proper binning
    mat.configureFace(OuterCylinder, {AxisRPhi, Bound, 20}, {AxisZ, Bound, 20});
    mat.configureFace(NegativeDisc, {AxisR, Bound, 15}, {AxisPhi, Bound, 25});
    mat.configureFace(PositiveDisc, {AxisR, Bound, 15}, {AxisPhi, Bound, 25});

    mat.addCylinderContainer("Detector", AxisDirection::AxisR, [&](auto& det) {
      det.addCylinderContainer("Pixel", AxisDirection::AxisZ, [&](auto& cyl) {
        cyl.setAttachmentStrategy(VolumeAttachmentStrategy::Gap)
            .setResizeStrategy(VolumeResizeStrategy::Gap);

        cyl.addCylinderContainer(
            "PixelNegativeEndcap", AxisDirection::AxisZ, [&](auto& ec) {
              ec.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);

              auto makeLayer = [&](const Transform3& trf, auto& layer) {
                std::vector<std::shared_ptr<Surface>> surfaces;
                auto layerSurfaces = makeFan(trf, 300_mm, 10, 2_mm);
                std::copy(layerSurfaces.begin(), layerSurfaces.end(),
                          std::back_inserter(surfaces));
                layerSurfaces = makeFan(trf, 500_mm, 16, 2_mm);
                std::copy(layerSurfaces.begin(), layerSurfaces.end(),
                          std::back_inserter(surfaces));

                layer.setSurfaces(surfaces)
                    .setLayerType(LayerBlueprintNode::LayerType::Disc)
                    .setEnvelope(ExtentEnvelope{{
                        .z = {5_mm, 5_mm},
                        .r = {10_mm, 20_mm},
                    }})
                    .setTransform(base);
              };

              ec.addLayer("PixelNeg1", [&](auto& layer) {
                makeLayer(base * Translation3{Vector3{0, 0, -700_mm}}, layer);
              });

              ec.addLayer("PixelNeg2", [&](auto& layer) {
                makeLayer(base * Translation3{Vector3{0, 0, -500_mm}}, layer);
              });
            });

        cyl.addCylinderContainer(
            "PixelBarrel", AxisDirection::AxisR, [&](auto& brl) {
              brl.setAttachmentStrategy(VolumeAttachmentStrategy::Gap)
                  .setResizeStrategy(VolumeResizeStrategy::Gap);

              auto makeLayer = [&](const std::string& name, double r,
                                   std::size_t nStaves, int nSensorsPerStave) {
                brl.addLayer(name, [&](auto& layer) {
                  std::vector<std::shared_ptr<Surface>> surfaces =
                      makeBarrelLayer(base, detectorElements, r, nStaves,
                                      nSensorsPerStave, 2.5_mm, 10_mm, 20_mm);

                  layer.setSurfaces(surfaces)
                      .setLayerType(LayerBlueprintNode::LayerType::Cylinder)
                      .setEnvelope(ExtentEnvelope{{
                          .z = {5_mm, 5_mm},
                          .r = {1_mm, 1_mm},
                      }})
                      .setTransform(base);
                });
              };

              makeLayer("PixelLayer0", 30_mm, 18, 5);
              makeLayer("PixelLayer1", 90_mm, 30, 6);

              brl.addStaticVolume(base,
                                  std::make_shared<CylinderVolumeBounds>(
                                      100_mm, 110_mm, 250_mm),
                                  "PixelSupport");

              makeLayer("PixelLayer2", 150_mm, 40, 7);
              makeLayer("PixelLayer3", 250_mm, 70, 8);
            });

        auto& ec =
            cyl.addCylinderContainer("PixelPosWrapper", AxisDirection::AxisR);
        ec.setResizeStrategy(VolumeResizeStrategy::Gap);
        ec.addStaticVolume(std::make_unique<TrackingVolume>(
            base * Translation3{Vector3{0, 0, 600_mm}},
            std::make_shared<CylinderVolumeBounds>(150_mm, 390_mm, 200_mm),
            "PixelPositiveEndcap"));
      });

      det.addStaticVolume(
          base, std::make_shared<CylinderVolumeBounds>(0_mm, 23_mm, 1000_mm),
          "BeamPipe");
    });
  });

  std::ofstream dot{"api_test_container.dot"};
  root->graphviz(dot);

  auto trackingGeometry = root->construct({}, gctx, *logger);

  BOOST_REQUIRE(trackingGeometry);
  BOOST_CHECK(trackingGeometry->geometryVersion() ==
              TrackingGeometry::GeometryVersion::Gen3);

  trackingGeometry->visitVolumes([&](const TrackingVolume* volume) {
    std::cout << volume->volumeName() << std::endl;
    std::cout << " -> id: " << volume->geometryId() << std::endl;
    std::cout << " -> " << volume->portals().size() << " portals" << std::endl;
  });

  ObjVisualization3D vis;

  trackingGeometry->visualize(vis, gctx, {}, {});

  vis.write("api_test_container.obj");

  Vector3 position = Vector3::Zero();
  std::ofstream csv{"api_test_container.csv"};
  csv << "x,y,z,volume,boundary,sensitive" << std::endl;

  std::mt19937 rnd{42};

  std::uniform_real_distribution<> dist{-1, 1};

  double etaWidth = 3;
  double thetaMin = 2 * std::atan(std::exp(-etaWidth));
  double thetaMax = 2 * std::atan(std::exp(etaWidth));
  std::uniform_real_distribution<> thetaDist{thetaMin, thetaMax};

  using namespace UnitLiterals;

  for (std::size_t i = 0; i < 5000; i++) {
    double theta = thetaDist(rnd);
    double phi = 2 * std::numbers::pi * dist(rnd);

    Vector3 direction;
    direction[0] = std::sin(theta) * std::cos(phi);
    direction[1] = std::sin(theta) * std::sin(phi);
    direction[2] = std::cos(theta);

    pseudoNavigation(*trackingGeometry, position, direction, csv, i, 2,
                     *logger->clone(std::nullopt, Logging::DEBUG));
  }
}

BOOST_AUTO_TEST_CASE(NodeApiTestCuboid) {
  Transform3 base{Transform3::Identity()};

  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {0_mm, 20_mm};
  auto root = std::make_unique<Blueprint>(cfg);

  root->addMaterial("GlobalMaterial", [&](MaterialDesignatorBlueprintNode&
                                              mat) {
    using enum AxisDirection;
    using enum AxisBoundaryType;
    using enum CuboidVolumeBounds::Face;

    // Configure valid axis combinations for each face type
    mat.configureFace(NegativeXFace, {AxisX, Bound, 20}, {AxisY, Bound, 20});
    mat.configureFace(PositiveXFace, {AxisX, Bound, 20}, {AxisY, Bound, 20});
    mat.configureFace(NegativeYFace, {AxisX, Bound, 15}, {AxisY, Bound, 25});
    mat.configureFace(PositiveYFace, {AxisX, Bound, 15}, {AxisY, Bound, 25});
    mat.configureFace(NegativeZFace, {AxisX, Bound, 15}, {AxisY, Bound, 25});
    mat.configureFace(PositiveZFace, {AxisX, Bound, 15}, {AxisY, Bound, 25});

    mat.addStaticVolume(
        base, std::make_shared<CuboidVolumeBounds>(100_mm, 100_mm, 100_mm),
        "TestVolume");
  });

  auto trackingGeometry = root->construct({}, gctx, *logger);
  BOOST_REQUIRE(trackingGeometry);
  BOOST_CHECK(trackingGeometry->geometryVersion() ==
              TrackingGeometry::GeometryVersion::Gen3);
}

// Reproduces the "material on a merged portal" failure: material is designated
// on a portal face that is subsequently merged during container stacking. The
// material designator wraps a child of a z-stacking container; since stacking
// in z merges the OuterCylinder portals of all children, the designated face
// cannot survive the merge.
//
// This is detected early by the container node, before the stack shell is
// built, producing a node-scoped error. The deeper shell-level reporting is
// exercised directly in the CylinderPortalShell tests.
BOOST_AUTO_TEST_CASE(MaterialOnMergedPortalThrows) {
  Transform3 base{Transform3::Identity()};

  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {0_mm, 20_mm};
  auto root = std::make_unique<Blueprint>(cfg);

  root->addCylinderContainer("Stack", AxisDirection::AxisZ, [&](auto& stack) {
    using enum AxisDirection;
    using enum AxisBoundaryType;
    using enum CylinderVolumeBounds::Face;

    // First child: a static volume whose OuterCylinder face is given material.
    // This is the face that the parent z-stack will try to merge.
    stack.addMaterial("Material", [&](auto& mat) {
      mat.configureFace(OuterCylinder, {AxisRPhi, Bound, 20},
                        {AxisZ, Bound, 20});
      mat.addStaticVolume(
          base * Translation3{Vector3{0, 0, -200_mm}},
          std::make_shared<CylinderVolumeBounds>(0_mm, 100_mm, 100_mm),
          "VolumeA");
    });

    // Second child: plain static volume on the other side in z.
    stack.addStaticVolume(
        base * Translation3{Vector3{0, 0, 200_mm}},
        std::make_shared<CylinderVolumeBounds>(0_mm, 100_mm, 100_mm),
        "VolumeB");
  });

  // The early-detection error should name the offending face, the shell
  // involved, and explain that material was placed on a merged face.
  bool thrown = false;
  {
    Logging::ScopedFailureThreshold threshold{Logging::Level::FATAL};
    try {
      root->construct({}, gctx, *logger);
    } catch (const PortalMergingException& e) {
      thrown = true;
      std::string msg = e.what();
      BOOST_CHECK(msg.find("OuterCylinder") != std::string::npos);
      BOOST_CHECK(msg.find("VolumeA") != std::string::npos);
      BOOST_CHECK(msg.find("material") != std::string::npos);
    }
  }
  BOOST_CHECK(thrown);
}

BOOST_AUTO_TEST_CASE(MaterialOnMergedPortalKeepGoing) {
  // Same blueprint as MaterialOnMergedPortalThrows, but constructed with the
  // keep-going option. Construction must succeed, and the merged outer cylinder
  // surface must carry a MergedMaterialMarker.
  Transform3 base{Transform3::Identity()};

  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {0_mm, 20_mm};
  auto root = std::make_unique<Blueprint>(cfg);

  root->addCylinderContainer("Stack", AxisDirection::AxisZ, [&](auto& stack) {
    using enum AxisDirection;
    using enum AxisBoundaryType;
    using enum CylinderVolumeBounds::Face;

    stack.addMaterial("Material", [&](auto& mat) {
      mat.configureFace(OuterCylinder, {AxisRPhi, Bound, 20},
                        {AxisZ, Bound, 20});
      mat.addStaticVolume(
          base * Translation3{Vector3{0, 0, -200_mm}},
          std::make_shared<CylinderVolumeBounds>(0_mm, 100_mm, 100_mm),
          "VolumeA");
    });

    stack.addStaticVolume(
        base * Translation3{Vector3{0, 0, 200_mm}},
        std::make_shared<CylinderVolumeBounds>(0_mm, 100_mm, 100_mm),
        "VolumeB");
  });

  Experimental::BlueprintOptions options;
  options.keepGoingOnMaterialMergeFailure = true;

  std::unique_ptr<const TrackingGeometry> trackingGeometry;
  {
    Logging::ScopedFailureThreshold threshold{Logging::Level::FATAL};
    BOOST_REQUIRE_NO_THROW(trackingGeometry =
                               root->construct(options, gctx, *logger));
  }
  BOOST_REQUIRE(trackingGeometry != nullptr);

  std::size_t markerCount = 0;
  trackingGeometry->visitSurfaces(
      [&](const Surface* surface) {
        if (surface != nullptr && surface->surfaceMaterial() != nullptr &&
            dynamic_cast<const MergedMaterialMarker*>(
                surface->surfaceMaterial()) != nullptr) {
          ++markerCount;
        }
      },
      false);
  BOOST_CHECK_GE(markerCount, 1u);
}

BOOST_AUTO_TEST_CASE(MaterialOnMergedPortalKeepGoingSingleChildFalseWarning) {
  // Reproducer for a false-positive in the early material-clash check:
  // a single-child AxisZ stack must NOT trigger a merge-material error because
  // no portal merge actually happens.  With keepGoingOnMaterialMergeFailure at
  // its default (false / strict), construction must succeed without throwing.
  Transform3 base{Transform3::Identity()};

  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {0_mm, 20_mm};
  auto root = std::make_unique<Blueprint>(cfg);

  root->addCylinderContainer("Stack", AxisDirection::AxisZ, [&](auto& stack) {
    using enum AxisDirection;
    using enum AxisBoundaryType;
    using enum CylinderVolumeBounds::Face;

    stack.addMaterial("Material", [&](auto& mat) {
      mat.configureFace(OuterCylinder, {AxisRPhi, Bound, 20},
                        {AxisZ, Bound, 20});
      mat.addStaticVolume(
          base, std::make_shared<CylinderVolumeBounds>(0_mm, 100_mm, 100_mm),
          "VolumeA");
    });
  });

  // Use strict mode (keepGoingOnMaterialMergeFailure = false by default).
  // The single-child stack must not be mistaken for a real merge, so no
  // exception should be thrown.
  Experimental::BlueprintOptions options;

  std::unique_ptr<const TrackingGeometry> trackingGeometry;
  BOOST_REQUIRE_NO_THROW(trackingGeometry =
                             root->construct(options, gctx, *logger));
  BOOST_REQUIRE(trackingGeometry != nullptr);

  std::size_t markerCount = 0;
  trackingGeometry->visitSurfaces(
      [&](const Surface* surface) {
        if (surface != nullptr && surface->surfaceMaterial() != nullptr &&
            dynamic_cast<const MergedMaterialMarker*>(
                surface->surfaceMaterial()) != nullptr) {
          ++markerCount;
        }
      },
      false);
  BOOST_CHECK_EQUAL(markerCount, 0u);
}

// Tag the fused face between two z-stacked volumes and look the portal back up
// from the final geometry. The tagged portal is shared by both volumes.
BOOST_AUTO_TEST_CASE(PortalTagLookup) {
  using Experimental::PortalDesignatorBlueprintNode;
  Transform3 base{Transform3::Identity()};

  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {0_mm, 20_mm};
  auto root = std::make_unique<Blueprint>(cfg);

  root->addCylinderContainer("Stack", AxisDirection::AxisZ, [&](auto& stack) {
    using enum CylinderVolumeBounds::Face;

    // Tag VolumeA's PositiveDisc, the face fused with VolumeB's NegativeDisc.
    stack.addPortalDesignator("Tags", [&](auto& tags) {
      tags.tagFace(PositiveDisc, "tracker_calo_boundary");
      tags.addStaticVolume(
          base * Translation3{Vector3{0, 0, -200_mm}},
          std::make_shared<CylinderVolumeBounds>(0_mm, 100_mm, 100_mm),
          "VolumeA");
    });

    stack.addStaticVolume(
        base * Translation3{Vector3{0, 0, 200_mm}},
        std::make_shared<CylinderVolumeBounds>(0_mm, 100_mm, 100_mm),
        "VolumeB");
  });

  std::unique_ptr<const TrackingGeometry> trackingGeometry =
      root->construct({}, gctx, *logger);
  BOOST_REQUIRE(trackingGeometry != nullptr);

  const Portal* portal = trackingGeometry->findPortal("tracker_calo_boundary");
  BOOST_REQUIRE(portal != nullptr);

  BOOST_CHECK(trackingGeometry->findPortal("does_not_exist") == nullptr);

  // The tagged portal must be the shared portal reachable from both volumes.
  auto containsPortal = [&](const TrackingVolume& volume) {
    for (const auto& p : volume.portals()) {
      if (&p == portal) {
        return true;
      }
    }
    return false;
  };

  const TrackingVolume* volumeA = nullptr;
  const TrackingVolume* volumeB = nullptr;
  trackingGeometry->apply([&](const TrackingVolume& volume) {
    if (volume.volumeName() == "VolumeA") {
      volumeA = &volume;
    } else if (volume.volumeName() == "VolumeB") {
      volumeB = &volume;
    }
  });

  BOOST_REQUIRE(volumeA != nullptr);
  BOOST_REQUIRE(volumeB != nullptr);
  BOOST_CHECK(containsPortal(*volumeA));
  BOOST_CHECK(containsPortal(*volumeB));
}

// Tagging two distinct (non-fused) portals with the same label must be detected
// as a collision when the geometry is closed.
BOOST_AUTO_TEST_CASE(PortalTagDuplicateThrows) {
  using Experimental::PortalDesignatorBlueprintNode;
  Transform3 base{Transform3::Identity()};

  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {0_mm, 20_mm};
  auto root = std::make_unique<Blueprint>(cfg);

  root->addCylinderContainer("Stack", AxisDirection::AxisZ, [&](auto& stack) {
    using enum CylinderVolumeBounds::Face;

    // Tag VolumeA's NegativeDisc (an outer boundary, not fused).
    stack.addPortalDesignator("TagsA", [&](auto& tags) {
      tags.tagFace(NegativeDisc, "dup");
      tags.addStaticVolume(
          base * Translation3{Vector3{0, 0, -200_mm}},
          std::make_shared<CylinderVolumeBounds>(0_mm, 100_mm, 100_mm),
          "VolumeA");
    });

    // Tag VolumeB's PositiveDisc (a different outer boundary) with the same
    // tag.
    stack.addPortalDesignator("TagsB", [&](auto& tags) {
      tags.tagFace(PositiveDisc, "dup");
      tags.addStaticVolume(
          base * Translation3{Vector3{0, 0, 200_mm}},
          std::make_shared<CylinderVolumeBounds>(0_mm, 100_mm, 100_mm),
          "VolumeB");
    });
  });

  Logging::ScopedFailureThreshold threshold{Logging::Level::FATAL};
  BOOST_CHECK_THROW(root->construct({}, gctx, *logger), std::invalid_argument);
}

// Same as PortalTagLookup, but for a cuboid x-stack: VolumeA's PositiveXFace is
// fused with VolumeB's NegativeXFace.
BOOST_AUTO_TEST_CASE(PortalTagLookupCuboid) {
  using Experimental::PortalDesignatorBlueprintNode;
  Transform3 base{Transform3::Identity()};

  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisX] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisY] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  auto root = std::make_unique<Blueprint>(cfg);

  root->addCuboidContainer("Stack", AxisDirection::AxisX, [&](auto& stack) {
    using enum CuboidVolumeBounds::Face;

    stack.addPortalDesignator("Tags", [&](auto& tags) {
      tags.tagFace(PositiveXFace, "cuboid_boundary");
      tags.addStaticVolume(
          base * Translation3{Vector3{-200_mm, 0, 0}},
          std::make_shared<CuboidVolumeBounds>(100_mm, 100_mm, 100_mm),
          "VolumeA");
    });

    stack.addStaticVolume(
        base * Translation3{Vector3{200_mm, 0, 0}},
        std::make_shared<CuboidVolumeBounds>(100_mm, 100_mm, 100_mm),
        "VolumeB");
  });

  std::unique_ptr<const TrackingGeometry> trackingGeometry =
      root->construct({}, gctx, *logger);
  BOOST_REQUIRE(trackingGeometry != nullptr);

  const Portal* portal = trackingGeometry->findPortal("cuboid_boundary");
  BOOST_REQUIRE(portal != nullptr);

  auto containsPortal = [&](const TrackingVolume& volume) {
    for (const auto& p : volume.portals()) {
      if (&p == portal) {
        return true;
      }
    }
    return false;
  };

  const TrackingVolume* volumeA = nullptr;
  const TrackingVolume* volumeB = nullptr;
  trackingGeometry->apply([&](const TrackingVolume& volume) {
    if (volume.volumeName() == "VolumeA") {
      volumeA = &volume;
    } else if (volume.volumeName() == "VolumeB") {
      volumeB = &volume;
    }
  });

  BOOST_REQUIRE(volumeA != nullptr);
  BOOST_REQUIRE(volumeB != nullptr);
  BOOST_CHECK(containsPortal(*volumeA));
  BOOST_CHECK(containsPortal(*volumeB));
}

BOOST_AUTO_TEST_SUITE_END();
