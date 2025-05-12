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

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeStack.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/LayerBlueprintNode.hpp"
#include "Acts/Geometry/MaterialDesignatorBlueprintNode.hpp"
#include "Acts/Geometry/ProcessorBlueprintNode.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"

#include <fstream>
#include <memory>
#include <stdexcept>
#include <vector>

using namespace Acts::UnitLiterals;

using Acts::Experimental::Blueprint;
using Acts::Experimental::BlueprintNode;
using Acts::Experimental::BlueprintOptions;
using Acts::Experimental::LayerBlueprintNode;
using Acts::Experimental::MaterialDesignatorBlueprintNode;
using Acts::Experimental::StaticBlueprintNode;

namespace Acts::Test {

auto logger = Acts::getDefaultLogger("UnitTests", Acts::Logging::INFO);

GeometryContext gctx;

namespace {

auto nameLookup(const TrackingGeometry& geo) {
  return [&](const std::string& name) -> const TrackingVolume& {
    const TrackingVolume* volume = nullptr;

    geo.visitVolumes([&](const TrackingVolume* v) {
      if (v->volumeName() == name) {
        volume = v;
      }
    });

    if (volume == nullptr) {
      throw std::runtime_error("Volume not found: " + name);
    }
    return *volume;
  };
}

std::size_t countVolumes(const TrackingGeometry& geo) {
  std::size_t nVolumes = 0;
  geo.visitVolumes([&](const TrackingVolume* /*volume*/) { nVolumes++; });
  return nVolumes;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(Geometry);

BOOST_AUTO_TEST_SUITE(BlueprintNodeTest);

BOOST_AUTO_TEST_CASE(InvalidRoot) {
  Logging::ScopedFailureThreshold threshold{Logging::Level::FATAL};

  Blueprint::Config cfg;
  Blueprint root{cfg};
  BOOST_CHECK_THROW(root.construct({}, gctx, *logger), std::logic_error);

  // More than one child is also invalid
  auto cylBounds = std::make_shared<CylinderVolumeBounds>(10_mm, 20_mm, 100_mm);
  root.addChild(
      std::make_unique<StaticBlueprintNode>(std::make_unique<TrackingVolume>(
          Transform3::Identity(), cylBounds, "child1")));
  root.addChild(
      std::make_unique<StaticBlueprintNode>(std::make_unique<TrackingVolume>(
          Transform3::Identity(), cylBounds, "child2")));

  BOOST_CHECK_THROW(root.construct({}, gctx, *logger), std::logic_error);
}

class DummyNode : public BlueprintNode {
 public:
  explicit DummyNode(const std::string& name) : m_name(name) {}

  const std::string& name() const override { return m_name; }

  Volume& build(const BlueprintOptions& /*options*/,
                const GeometryContext& /*gctx*/,
                const Acts::Logger& /*logger*/) override {
    throw std::logic_error("Not implemented");
  }

  PortalShellBase& connect(const BlueprintOptions& /*options*/,
                           const GeometryContext& /*gctx*/,
                           const Logger& /*logger */) override {
    throw std::logic_error("Not implemented");
  }

  void finalize(const BlueprintOptions& /*options*/,
                const GeometryContext& /*gctx*/, TrackingVolume& /*parent*/,
                const Logger& /*logger*/) override {
    throw std::logic_error("Not implemented");
  }

 private:
  std::string m_name;
};

BOOST_AUTO_TEST_CASE(RootCannotBeChild) {
  auto node = std::make_shared<DummyNode>("node");
  Blueprint::Config cfg;
  auto root = std::make_shared<Blueprint>(cfg);

  BOOST_CHECK_THROW(node->addChild(root), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(AddChildInvalid) {
  auto node = std::make_shared<DummyNode>("node");

  // Add self
  BOOST_CHECK_THROW(node->addChild(node), std::invalid_argument);

  // Add nullptr
  BOOST_CHECK_THROW(node->addChild(nullptr), std::invalid_argument);

  auto nodeB = std::make_shared<DummyNode>("nodeB");
  auto nodeC = std::make_shared<DummyNode>("nodeC");

  node->addChild(nodeB);
  nodeB->addChild(nodeC);
  BOOST_CHECK_THROW(nodeC->addChild(node), std::invalid_argument);

  // already has parent, can't be added as a child anywhere else
  BOOST_CHECK_THROW(node->addChild(nodeC), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(Depth) {
  auto node1 = std::make_shared<DummyNode>("node1");
  auto node2 = std::make_shared<DummyNode>("node2");
  auto node3 = std::make_shared<DummyNode>("node3");

  BOOST_CHECK_EQUAL(node1->depth(), 0);
  BOOST_CHECK_EQUAL(node2->depth(), 0);
  BOOST_CHECK_EQUAL(node3->depth(), 0);

  node2->addChild(node3);
  BOOST_CHECK_EQUAL(node2->depth(), 0);
  BOOST_CHECK_EQUAL(node3->depth(), 1);

  node1->addChild(node2);
  BOOST_CHECK_EQUAL(node1->depth(), 0);
  BOOST_CHECK_EQUAL(node2->depth(), 1);
  BOOST_CHECK_EQUAL(node3->depth(), 2);
}

BOOST_AUTO_TEST_CASE(Static) {
  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {1_mm, 2_mm};
  Blueprint root{cfg};

  double hlZ = 30_mm;
  auto cylBounds = std::make_shared<CylinderVolumeBounds>(10_mm, 20_mm, hlZ);
  auto cyl = std::make_unique<TrackingVolume>(Transform3::Identity(), cylBounds,
                                              "child");

  root.addStaticVolume(std::move(cyl));

  BOOST_CHECK_EQUAL(root.children().size(), 1);

  auto tGeometry = root.construct({}, gctx, *logger);

  BOOST_REQUIRE(tGeometry);

  BOOST_CHECK(tGeometry->geometryVersion() ==
              TrackingGeometry::GeometryVersion::Gen3);

  BOOST_CHECK_EQUAL(tGeometry->highestTrackingVolume()->volumes().size(), 1);

  BOOST_CHECK_EQUAL(countVolumes(*tGeometry), 2);

  auto lookup = nameLookup(*tGeometry);
  auto actCyl =
      dynamic_cast<const CylinderVolumeBounds&>(lookup("child").volumeBounds());
  // Size as given
  BOOST_CHECK_EQUAL(actCyl.get(CylinderVolumeBounds::eMinR), 10_mm);
  BOOST_CHECK_EQUAL(actCyl.get(CylinderVolumeBounds::eMaxR), 20_mm);
  BOOST_CHECK_EQUAL(actCyl.get(CylinderVolumeBounds::eHalfLengthZ), hlZ);

  auto worldCyl =
      dynamic_cast<const CylinderVolumeBounds&>(lookup("World").volumeBounds());
  BOOST_CHECK_EQUAL(worldCyl.get(CylinderVolumeBounds::eMinR), 9_mm);
  BOOST_CHECK_EQUAL(worldCyl.get(CylinderVolumeBounds::eMaxR), 22_mm);
  BOOST_CHECK_EQUAL(worldCyl.get(CylinderVolumeBounds::eHalfLengthZ),
                    hlZ + 20_mm);

  BOOST_CHECK_EQUAL(lookup("World").portals().size(), 8);
}

BOOST_AUTO_TEST_CASE(CylinderContainer) {
  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {2_mm, 20_mm};
  auto root = std::make_unique<Blueprint>(cfg);

  auto& cyl = root->addCylinderContainer("Container", AxisDirection::AxisZ);
  cyl.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);

  double z0 = -200_mm;
  double hlZ = 30_mm;
  auto cylBounds = std::make_shared<CylinderVolumeBounds>(10_mm, 20_mm, hlZ);
  for (std::size_t i = 0; i < 3; i++) {
    auto childCyl = std::make_unique<TrackingVolume>(
        Transform3::Identity() *
            Translation3{Vector3{0, 0, z0 + i * 2 * hlZ * 1.2}},
        cylBounds, "child" + std::to_string(i));
    cyl.addStaticVolume(std::move(childCyl));
  }

  auto tGeometry = root->construct({}, gctx, *logger);
  BOOST_CHECK(tGeometry->geometryVersion() ==
              TrackingGeometry::GeometryVersion::Gen3);

  // 4 real volumes + 2 gaps
  BOOST_CHECK_EQUAL(countVolumes(*tGeometry), 6);

  auto lookup = nameLookup(*tGeometry);
  auto worldCyl =
      dynamic_cast<const CylinderVolumeBounds&>(lookup("World").volumeBounds());
  BOOST_CHECK_EQUAL(worldCyl.get(CylinderVolumeBounds::eMinR), 8_mm);
  BOOST_CHECK_EQUAL(worldCyl.get(CylinderVolumeBounds::eMaxR), 40_mm);
  BOOST_CHECK_EQUAL(worldCyl.get(CylinderVolumeBounds::eHalfLengthZ), 122_mm);

  BOOST_CHECK_EQUAL(lookup("World").portals().size(), 8);

  for (std::size_t i = 0; i < 3; i++) {
    auto actCyl = dynamic_cast<const CylinderVolumeBounds&>(
        lookup("child" + std::to_string(i)).volumeBounds());
    BOOST_CHECK_EQUAL(actCyl.get(CylinderVolumeBounds::eMinR), 10_mm);
    BOOST_CHECK_EQUAL(actCyl.get(CylinderVolumeBounds::eMaxR), 20_mm);
    BOOST_CHECK_EQUAL(actCyl.get(CylinderVolumeBounds::eHalfLengthZ), hlZ);
  }

  for (std::size_t i = 0; i < 2; i++) {
    auto gapCyl = dynamic_cast<const CylinderVolumeBounds&>(
        lookup("Container::Gap" + std::to_string(i + 1)).volumeBounds());
    BOOST_CHECK_EQUAL(gapCyl.get(CylinderVolumeBounds::eMinR), 10_mm);
    BOOST_CHECK_EQUAL(gapCyl.get(CylinderVolumeBounds::eMaxR), 20_mm);
    BOOST_CHECK_EQUAL(gapCyl.get(CylinderVolumeBounds::eHalfLengthZ), 6_mm);
  }
}

BOOST_AUTO_TEST_CASE(Confined) {
  Transform3 base{Transform3::Identity()};

  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {2_mm, 20_mm};
  auto root = std::make_unique<Blueprint>(cfg);

  root->addStaticVolume(
      base, std::make_shared<CylinderVolumeBounds>(50_mm, 400_mm, 1000_mm),
      "PixelWrapper", [&](auto& wrap) {
        double rMin = 100_mm;
        double rMax = 350_mm;
        double hlZ = 100_mm;

        wrap.addStaticVolume(
            base * Translation3{Vector3{0, 0, -600_mm}},
            std::make_shared<CylinderVolumeBounds>(rMin, rMax, hlZ),
            "PixelNeg1");

        wrap.addStaticVolume(
            base * Translation3{Vector3{0, 0, -200_mm}},
            std::make_shared<CylinderVolumeBounds>(rMin, rMax, hlZ),
            "PixelNeg2");

        wrap.addStaticVolume(
            base * Translation3{Vector3{0, 0, 200_mm}},
            std::make_shared<CylinderVolumeBounds>(rMin, rMax, hlZ),
            "PixelPos1");

        wrap.addStaticVolume(
            base * Translation3{Vector3{0, 0, 600_mm}},
            std::make_shared<CylinderVolumeBounds>(rMin, rMax, hlZ),
            "PixelPos2");
      });

  auto trackingGeometry = root->construct({}, gctx, *logger);

  // overall dimensions are the wrapper volume + envelope
  auto lookup = nameLookup(*trackingGeometry);
  auto worldCyl =
      dynamic_cast<const CylinderVolumeBounds&>(lookup("World").volumeBounds());
  BOOST_CHECK_EQUAL(worldCyl.get(CylinderVolumeBounds::eMinR), 48_mm);
  BOOST_CHECK_EQUAL(worldCyl.get(CylinderVolumeBounds::eMaxR), 420_mm);
  BOOST_CHECK_EQUAL(worldCyl.get(CylinderVolumeBounds::eHalfLengthZ), 1020_mm);

  // 4 outer portals and 4 inner
  BOOST_CHECK_EQUAL(lookup("World").portals().size(), 8);
  BOOST_CHECK_EQUAL(lookup("World").volumes().size(), 1);

  auto wrapperCyl = dynamic_cast<const CylinderVolumeBounds&>(
      lookup("PixelWrapper").volumeBounds());
  BOOST_CHECK_EQUAL(wrapperCyl.get(CylinderVolumeBounds::eMinR), 50_mm);
  BOOST_CHECK_EQUAL(wrapperCyl.get(CylinderVolumeBounds::eMaxR), 400_mm);
  BOOST_CHECK_EQUAL(wrapperCyl.get(CylinderVolumeBounds::eHalfLengthZ),
                    1000_mm);
  BOOST_CHECK_EQUAL(lookup("PixelWrapper").portals().size(), 4 + 4 * 4);
  BOOST_CHECK_EQUAL(lookup("PixelWrapper").volumes().size(), 4);

  for (const auto& name :
       {"PixelNeg1", "PixelNeg2", "PixelPos1", "PixelPos2"}) {
    auto actCyl =
        dynamic_cast<const CylinderVolumeBounds&>(lookup(name).volumeBounds());
    BOOST_CHECK_EQUAL(actCyl.get(CylinderVolumeBounds::eMinR), 100_mm);
    BOOST_CHECK_EQUAL(actCyl.get(CylinderVolumeBounds::eMaxR), 350_mm);
    BOOST_CHECK_EQUAL(actCyl.get(CylinderVolumeBounds::eHalfLengthZ), 100_mm);
    BOOST_CHECK_EQUAL(lookup(name).portals().size(), 4);
  }
}

BOOST_AUTO_TEST_CASE(DiscLayer) {
  double yrot = 45_degree;
  Transform3 base = Transform3::Identity() * AngleAxis3{yrot, Vector3::UnitY()};

  std::vector<std::shared_ptr<Surface>> surfaces;
  std::vector<std::unique_ptr<DetectorElementBase>> elements;
  double r = 300_mm;
  std::size_t nSensors = 8;
  double thickness = 2.5_mm;
  auto recBounds = std::make_shared<RectangleBounds>(40_mm, 60_mm);

  double deltaPhi = 2 * std::numbers::pi / nSensors;
  for (std::size_t i = 0; i < nSensors; i++) {
    // Create a fan of sensors

    Transform3 trf = base * AngleAxis3{deltaPhi * i, Vector3::UnitZ()} *
                     Translation3(Vector3::UnitX() * r);

    if (i % 2 == 0) {
      trf = trf * Translation3{Vector3::UnitZ() * 5_mm};
    }

    auto& element = elements.emplace_back(
        std::make_unique<DetectorElementStub>(trf, recBounds, thickness));

    element->surface().assignDetectorElement(*element);

    surfaces.push_back(element->surface().getSharedPtr());
  }

  std::function<void(LayerBlueprintNode&)> withSurfaces =
      [&surfaces, &base](LayerBlueprintNode& layer) {
        layer.setSurfaces(surfaces)
            .setLayerType(LayerBlueprintNode::LayerType::Disc)
            .setEnvelope(ExtentEnvelope{{
                .z = {0.1_mm, 0.1_mm},
                .r = {1_mm, 1_mm},
            }})
            .setTransform(base);
      };

  std::function<void(LayerBlueprintNode&)> withProtoLayer =
      [&surfaces, &base](LayerBlueprintNode& layer) {
        MutableProtoLayer protoLayer{gctx, surfaces, base.inverse()};
        layer.setProtoLayer(protoLayer)
            .setLayerType(LayerBlueprintNode::LayerType::Disc)
            .setEnvelope(ExtentEnvelope{{
                .z = {0.1_mm, 0.1_mm},
                .r = {1_mm, 1_mm},
            }});
      };

  for (const auto& func : {withSurfaces, withProtoLayer}) {
    Blueprint root{{.envelope = ExtentEnvelope{{
                        .z = {2_mm, 2_mm},
                        .r = {3_mm, 5_mm},
                    }}}};

    root.addLayer("Layer0", [&](auto& layer) { func(layer); });

    auto trackingGeometry = root.construct({}, gctx, *logger);

    std::size_t nSurfaces = 0;

    trackingGeometry->visitSurfaces([&](const Surface* surface) {
      if (surface->associatedDetectorElement() != nullptr) {
        nSurfaces++;
      }
    });

    BOOST_CHECK_EQUAL(nSurfaces, surfaces.size());
    BOOST_CHECK_EQUAL(countVolumes(*trackingGeometry), 2);
    auto lookup = nameLookup(*trackingGeometry);
    auto layerCyl = dynamic_cast<const CylinderVolumeBounds&>(
        lookup("Layer0").volumeBounds());
    BOOST_CHECK_CLOSE(layerCyl.get(CylinderVolumeBounds::eMinR), 258.9999999_mm,
                      1e-6);
    BOOST_CHECK_CLOSE(layerCyl.get(CylinderVolumeBounds::eMaxR),
                      346.25353003_mm, 1e-6);
    BOOST_CHECK_CLOSE(layerCyl.get(CylinderVolumeBounds::eHalfLengthZ), 3.85_mm,
                      1e-6);
  }
}

BOOST_AUTO_TEST_CASE(CylinderLayer) {
  double yrot = 0_degree;
  Transform3 base = Transform3::Identity() * AngleAxis3{yrot, Vector3::UnitY()};

  std::vector<std::shared_ptr<Surface>> surfaces;
  std::vector<std::unique_ptr<DetectorElementBase>> elements;

  double r = 300_mm;
  std::size_t nStaves = 10;
  int nSensorsPerStave = 8;
  double thickness = 0;
  double hlPhi = 40_mm;
  double hlZ = 60_mm;
  auto recBounds = std::make_shared<RectangleBounds>(hlPhi, hlZ);

  double deltaPhi = 2 * std::numbers::pi / nStaves;

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
      element->surface().assignDetectorElement(*element);
      surfaces.push_back(element->surface().getSharedPtr());
    }
  }

  std::function<void(LayerBlueprintNode&)> withSurfaces =
      [&surfaces, &base](LayerBlueprintNode& layer) {
        layer.setSurfaces(surfaces)
            .setLayerType(LayerBlueprintNode::LayerType::Cylinder)
            .setEnvelope(ExtentEnvelope{{
                .z = {10_mm, 10_mm},
                .r = {20_mm, 10_mm},
            }})
            .setTransform(base);
      };

  std::function<void(LayerBlueprintNode&)> withProtoLayer =
      [&surfaces, &base](LayerBlueprintNode& layer) {
        MutableProtoLayer protoLayer{gctx, surfaces, base.inverse()};
        layer.setProtoLayer(protoLayer)
            .setLayerType(LayerBlueprintNode::LayerType::Cylinder)
            .setEnvelope(ExtentEnvelope{{
                .z = {10_mm, 10_mm},
                .r = {20_mm, 10_mm},
            }});
      };

  for (const auto& func : {withSurfaces, withProtoLayer}) {
    Blueprint root{{.envelope = ExtentEnvelope{{
                        .z = {2_mm, 2_mm},
                        .r = {3_mm, 5_mm},
                    }}}};

    root.addLayer("Layer0", [&](auto& layer) { func(layer); });

    auto trackingGeometry = root.construct({}, gctx, *logger);

    std::size_t nSurfaces = 0;

    trackingGeometry->visitSurfaces([&](const Surface* surface) {
      if (surface->associatedDetectorElement() != nullptr) {
        nSurfaces++;
      }
    });

    BOOST_CHECK_EQUAL(nSurfaces, surfaces.size());
    BOOST_CHECK_EQUAL(countVolumes(*trackingGeometry), 2);
    auto lookup = nameLookup(*trackingGeometry);
    auto layerCyl = dynamic_cast<const CylinderVolumeBounds&>(
        lookup("Layer0").volumeBounds());
    BOOST_CHECK_EQUAL(lookup("Layer0").portals().size(), 4);
    BOOST_CHECK_CLOSE(layerCyl.get(CylinderVolumeBounds::eMinR), 275.6897761_mm,
                      1e-6);
    BOOST_CHECK_CLOSE(layerCyl.get(CylinderVolumeBounds::eMaxR), 319.4633358_mm,
                      1e-6);
    BOOST_CHECK_CLOSE(layerCyl.get(CylinderVolumeBounds::eHalfLengthZ), 1070_mm,
                      1e-6);
  }
}

BOOST_AUTO_TEST_CASE(Material) {
  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {1_mm, 2_mm};
  Blueprint root{cfg};

  double hlZ = 30_mm;
  auto cylBounds = std::make_shared<CylinderVolumeBounds>(10_mm, 20_mm, hlZ);
  auto cyl = std::make_unique<TrackingVolume>(Transform3::Identity(), cylBounds,
                                              "child");

  using enum AxisDirection;
  using enum CylinderVolumeBounds::Face;
  using enum AxisBoundaryType;

  root.addMaterial("Material", [&](auto& mat) {
    mat.configureFace(NegativeDisc, {AxisR, Bound, 5}, {AxisPhi, Bound, 10});
    mat.configureFace(PositiveDisc, {AxisR, Bound, 15}, {AxisPhi, Bound, 20});
    mat.configureFace(OuterCylinder, {AxisRPhi, Bound, 25}, {AxisZ, Bound, 30});

    mat.addStaticVolume(std::move(cyl));
  });

  auto trackingGeometry = root.construct({}, gctx, *logger);

  BOOST_CHECK_EQUAL(countVolumes(*trackingGeometry), 2);
  auto lookup = nameLookup(*trackingGeometry);
  auto& child = lookup("child");

  // Check negative disc material
  const auto* negDisc = child.portals()
                            .at(static_cast<std::size_t>(NegativeDisc))
                            .surface()
                            .surfaceMaterial();
  BOOST_REQUIRE_NE(negDisc, nullptr);
  const auto& negDiscMat =
      dynamic_cast<const ProtoGridSurfaceMaterial&>(*negDisc);
  // Check positive disc material
  const auto* posDisc = child.portals()
                            .at(static_cast<std::size_t>(PositiveDisc))
                            .surface()
                            .surfaceMaterial();
  BOOST_REQUIRE_NE(posDisc, nullptr);
  const auto& posDiscMat =
      dynamic_cast<const ProtoGridSurfaceMaterial&>(*posDisc);

  BOOST_CHECK_EQUAL(negDiscMat.binning().at(0).getAxis().getNBins(), 5);
  BOOST_CHECK_EQUAL(negDiscMat.binning().at(1).getAxis().getNBins(), 10);
  BOOST_CHECK_EQUAL(posDiscMat.binning().at(0).getAxis().getNBins(), 15);
  BOOST_CHECK_EQUAL(posDiscMat.binning().at(1).getAxis().getNBins(), 20);

  // Check outer cylinder material
  const auto* outerCyl = child.portals()
                             .at(static_cast<std::size_t>(OuterCylinder))
                             .surface()
                             .surfaceMaterial();
  BOOST_REQUIRE_NE(outerCyl, nullptr);
  const auto& outerCylMat =
      dynamic_cast<const ProtoGridSurfaceMaterial&>(*outerCyl);
  BOOST_CHECK_EQUAL(outerCylMat.binning().at(0).getAxis().getNBins(), 25);
  BOOST_CHECK_EQUAL(outerCylMat.binning().at(1).getAxis().getNBins(), 30);

  // Check that other faces have no material
  for (std::size_t i = 0; i < child.portals().size(); i++) {
    if (i != static_cast<std::size_t>(NegativeDisc) &&
        i != static_cast<std::size_t>(PositiveDisc) &&
        i != static_cast<std::size_t>(OuterCylinder)) {
      BOOST_CHECK_EQUAL(child.portals().at(i).surface().surfaceMaterial(),
                        nullptr);
    }
  }
}

BOOST_AUTO_TEST_CASE(MaterialInvalidAxisDirections) {
  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {1_mm, 2_mm};
  Blueprint root{cfg};

  using enum AxisDirection;
  using enum AxisBoundaryType;

  // Test invalid axis direction combinations for cylinder faces
  BOOST_CHECK_THROW(
      root.addMaterial("Material",
                       [&](auto& mat) {
                         mat.configureFace(
                             CylinderVolumeBounds::Face::NegativeDisc,
                             {AxisZ, Bound, 5}, {AxisPhi, Bound, 10});
                       }),
      std::invalid_argument);

  BOOST_CHECK_THROW(
      root.addMaterial("Material",
                       [&](auto& mat) {
                         mat.configureFace(
                             CylinderVolumeBounds::Face::OuterCylinder,
                             {AxisR, Bound, 5}, {AxisR, Bound, 10});
                       }),
      std::invalid_argument);

  // Test invalid axis direction combinations for cuboid faces
  BOOST_CHECK_THROW(
      root.addMaterial("Material",
                       [&](auto& mat) {
                         mat.configureFace(
                             CuboidVolumeBounds::Face::NegativeXFace,
                             {AxisX, Bound, 5}, {AxisZ, Bound, 10});
                       }),
      std::invalid_argument);

  BOOST_CHECK_THROW(
      root.addMaterial("Material",
                       [&](auto& mat) {
                         mat.configureFace(
                             CuboidVolumeBounds::Face::PositiveYFace,
                             {AxisY, Bound, 5}, {AxisX, Bound, 10});
                       }),
      std::invalid_argument);

  BOOST_CHECK_THROW(
      root.addMaterial("Material",
                       [&](auto& mat) {
                         mat.configureFace(
                             CuboidVolumeBounds::Face::NegativeZFace,
                             {AxisZ, Bound, 5}, {AxisY, Bound, 10});
                       }),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(MaterialMixedVolumeTypes) {
  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {1_mm, 2_mm};
  Blueprint root{cfg};

  using enum AxisDirection;
  using enum AxisBoundaryType;

  // Configure for cylinder first, then try to add cuboid - should throw
  BOOST_CHECK_THROW(
      root.addMaterial(
          "Material",
          [&](auto& mat) {
            mat.configureFace(CylinderVolumeBounds::Face::NegativeDisc,
                              {AxisR, Bound, 5}, {AxisPhi, Bound, 10});
            mat.configureFace(CuboidVolumeBounds::Face::NegativeXFace,
                              {AxisX, Bound, 5}, {AxisY, Bound, 10});
          }),
      std::runtime_error);

  // Configure for cuboid first, then try to add cylinder - should throw
  BOOST_CHECK_THROW(
      root.addMaterial(
          "Material",
          [&](auto& mat) {
            mat.configureFace(CuboidVolumeBounds::Face::NegativeXFace,
                              {AxisX, Bound, 5}, {AxisY, Bound, 10});
            mat.configureFace(CylinderVolumeBounds::Face::NegativeDisc,
                              {AxisR, Bound, 5}, {AxisPhi, Bound, 10});
          }),
      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(MaterialCuboid) {
  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {1_mm, 2_mm};
  Blueprint root{cfg};

  using enum AxisDirection;
  using enum AxisBoundaryType;
  using enum CuboidVolumeBounds::Face;

  double hlX = 30_mm;
  double hlY = 40_mm;
  double hlZ = 50_mm;
  auto cuboidBounds = std::make_shared<CuboidVolumeBounds>(hlX, hlY, hlZ);
  auto cuboid = std::make_unique<TrackingVolume>(Transform3::Identity(),
                                                 cuboidBounds, "child");

  auto mat = std::make_shared<MaterialDesignatorBlueprintNode>("Material");

  // Configure material for different faces with different binning
  mat->configureFace(NegativeXFace, {AxisX, Bound, 5}, {AxisY, Bound, 10});
  mat->configureFace(PositiveXFace, {AxisX, Bound, 15}, {AxisY, Bound, 20});
  mat->configureFace(NegativeYFace, {AxisX, Bound, 25}, {AxisY, Bound, 30});
  mat->configureFace(PositiveYFace, {AxisX, Bound, 35}, {AxisY, Bound, 40});
  mat->configureFace(NegativeZFace, {AxisX, Bound, 45}, {AxisY, Bound, 50});
  mat->configureFace(PositiveZFace, {AxisX, Bound, 55}, {AxisY, Bound, 60});

  mat->addChild(std::make_shared<StaticBlueprintNode>(std::move(cuboid)));

  root.addChild(mat);

  auto trackingGeometry =
      root.construct({}, gctx, *logger->clone(std::nullopt, Logging::VERBOSE));

  BOOST_REQUIRE(trackingGeometry);

  auto lookup = nameLookup(*trackingGeometry);
  auto& child = lookup("child");

  // Check that material is attached to all faces
  for (std::size_t i = 0; i < child.portals().size(); i++) {
    const auto* material = child.portals().at(i).surface().surfaceMaterial();
    BOOST_REQUIRE_NE(material, nullptr);

    const auto& gridMaterial =
        dynamic_cast<const ProtoGridSurfaceMaterial&>(*material);

    // Check binning based on face
    CuboidVolumeBounds::Face face = static_cast<CuboidVolumeBounds::Face>(i);
    switch (face) {
      case NegativeXFace:
        BOOST_CHECK_EQUAL(gridMaterial.binning().at(0).getAxis().getNBins(), 5);
        BOOST_CHECK_EQUAL(gridMaterial.binning().at(1).getAxis().getNBins(),
                          10);
        break;
      case PositiveXFace:
        BOOST_CHECK_EQUAL(gridMaterial.binning().at(0).getAxis().getNBins(),
                          15);
        BOOST_CHECK_EQUAL(gridMaterial.binning().at(1).getAxis().getNBins(),
                          20);
        break;
      case NegativeYFace:
        BOOST_CHECK_EQUAL(gridMaterial.binning().at(0).getAxis().getNBins(),
                          25);
        BOOST_CHECK_EQUAL(gridMaterial.binning().at(1).getAxis().getNBins(),
                          30);
        break;
      case PositiveYFace:
        BOOST_CHECK_EQUAL(gridMaterial.binning().at(0).getAxis().getNBins(),
                          35);
        BOOST_CHECK_EQUAL(gridMaterial.binning().at(1).getAxis().getNBins(),
                          40);
        break;
      case NegativeZFace:
        BOOST_CHECK_EQUAL(gridMaterial.binning().at(0).getAxis().getNBins(),
                          45);
        BOOST_CHECK_EQUAL(gridMaterial.binning().at(1).getAxis().getNBins(),
                          50);
        break;
      case PositiveZFace:
        BOOST_CHECK_EQUAL(gridMaterial.binning().at(0).getAxis().getNBins(),
                          55);
        BOOST_CHECK_EQUAL(gridMaterial.binning().at(1).getAxis().getNBins(),
                          60);
        break;
    }
  }
}

BOOST_AUTO_TEST_CASE(HomogeneousMaterialCylinder) {
  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {1_mm, 2_mm};
  Blueprint root{cfg};

  double hlZ = 30_mm;
  auto cylBounds = std::make_shared<CylinderVolumeBounds>(10_mm, 20_mm, hlZ);
  auto cyl = std::make_unique<TrackingVolume>(Transform3::Identity(), cylBounds,
                                              "child");

  using enum CylinderVolumeBounds::Face;

  // Create some homogeneous materials with different properties
  auto testMaterial = Acts::Material::fromMolarDensity(
      9.370_cm, 46.52_cm, 28.0855, 14, (2.329 / 28.0855) * 1_mol / 1_cm3);

  auto negDiscMat = std::make_shared<HomogeneousSurfaceMaterial>(
      MaterialSlab(testMaterial, 0.1_mm));
  auto posDiscMat = std::make_shared<HomogeneousSurfaceMaterial>(
      MaterialSlab(testMaterial, 0.2_mm));
  auto outerCylMat = std::make_shared<HomogeneousSurfaceMaterial>(
      MaterialSlab(testMaterial, 0.3_mm));

  root.addMaterial("Material", [&](auto& mat) {
    mat.configureFace(NegativeDisc, negDiscMat);
    mat.configureFace(PositiveDisc, posDiscMat);
    mat.configureFace(OuterCylinder, outerCylMat);

    mat.addStaticVolume(std::move(cyl));
  });

  auto trackingGeometry = root.construct({}, gctx, *logger);

  BOOST_CHECK_EQUAL(countVolumes(*trackingGeometry), 2);
  auto lookup = nameLookup(*trackingGeometry);
  auto& child = lookup("child");

  // Check negative disc material
  const auto* negDisc = child.portals()
                            .at(static_cast<std::size_t>(NegativeDisc))
                            .surface()
                            .surfaceMaterial();
  BOOST_REQUIRE_NE(negDisc, nullptr);
  BOOST_CHECK_EQUAL(negDisc, negDiscMat.get());

  // Check positive disc material
  const auto* posDisc = child.portals()
                            .at(static_cast<std::size_t>(PositiveDisc))
                            .surface()
                            .surfaceMaterial();
  BOOST_REQUIRE_NE(posDisc, nullptr);
  BOOST_CHECK_EQUAL(posDisc, posDiscMat.get());

  // Check outer cylinder material
  const auto* outerCyl = child.portals()
                             .at(static_cast<std::size_t>(OuterCylinder))
                             .surface()
                             .surfaceMaterial();
  BOOST_REQUIRE_NE(outerCyl, nullptr);
  BOOST_CHECK_EQUAL(outerCyl, outerCylMat.get());

  // Check that other faces have no material
  for (std::size_t i = 0; i < child.portals().size(); i++) {
    if (i != static_cast<std::size_t>(NegativeDisc) &&
        i != static_cast<std::size_t>(PositiveDisc) &&
        i != static_cast<std::size_t>(OuterCylinder)) {
      BOOST_CHECK_EQUAL(child.portals().at(i).surface().surfaceMaterial(),
                        nullptr);
    }
  }
}

BOOST_AUTO_TEST_CASE(HomogeneousMaterialCuboid) {
  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {1_mm, 2_mm};
  Blueprint root{cfg};

  using enum CuboidVolumeBounds::Face;

  double hlX = 30_mm;
  double hlY = 40_mm;
  double hlZ = 50_mm;
  auto cuboidBounds = std::make_shared<CuboidVolumeBounds>(hlX, hlY, hlZ);
  auto cuboid = std::make_unique<TrackingVolume>(Transform3::Identity(),
                                                 cuboidBounds, "child");

  // Create different homogeneous materials for each face
  auto testMaterial = Acts::Material::fromMolarDensity(
      9.370_cm, 46.52_cm, 28.0855, 14, (2.329 / 28.0855) * 1_mol / 1_cm3);

  auto negXMat = std::make_shared<HomogeneousSurfaceMaterial>(
      MaterialSlab(testMaterial, 0.1_mm));
  auto posXMat = std::make_shared<HomogeneousSurfaceMaterial>(
      MaterialSlab(testMaterial, 0.2_mm));
  auto negYMat = std::make_shared<HomogeneousSurfaceMaterial>(
      MaterialSlab(testMaterial, 0.3_mm));
  auto posYMat = std::make_shared<HomogeneousSurfaceMaterial>(
      MaterialSlab(testMaterial, 0.4_mm));
  auto negZMat = std::make_shared<HomogeneousSurfaceMaterial>(
      MaterialSlab(testMaterial, 0.5_mm));
  auto posZMat = std::make_shared<HomogeneousSurfaceMaterial>(
      MaterialSlab(testMaterial, 0.6_mm));

  root.addMaterial("Material", [&](auto& mat) {
    mat.configureFace(NegativeXFace, negXMat);
    mat.configureFace(PositiveXFace, posXMat);
    mat.configureFace(NegativeYFace, negYMat);
    mat.configureFace(PositiveYFace, posYMat);
    mat.configureFace(NegativeZFace, negZMat);
    mat.configureFace(PositiveZFace, posZMat);

    mat.addStaticVolume(std::move(cuboid));
  });

  auto trackingGeometry = root.construct({}, gctx, *logger);

  BOOST_REQUIRE(trackingGeometry);

  auto lookup = nameLookup(*trackingGeometry);
  auto& child = lookup("child");

  // Check that material is attached to all faces with correct properties
  for (std::size_t i = 0; i < child.portals().size(); i++) {
    const auto* material = child.portals().at(i).surface().surfaceMaterial();
    BOOST_REQUIRE_NE(material, nullptr);

    const auto* homMaterial =
        dynamic_cast<const HomogeneousSurfaceMaterial*>(material);
    BOOST_REQUIRE_NE(homMaterial, nullptr);

    // Check thickness based on face
    CuboidVolumeBounds::Face face = static_cast<CuboidVolumeBounds::Face>(i);
    switch (face) {
      case NegativeXFace:
        BOOST_CHECK_EQUAL(homMaterial, negXMat.get());
        break;
      case PositiveXFace:
        BOOST_CHECK_EQUAL(homMaterial, posXMat.get());
        break;
      case NegativeYFace:
        BOOST_CHECK_EQUAL(homMaterial, negYMat.get());
        break;
      case PositiveYFace:
        BOOST_CHECK_EQUAL(homMaterial, posYMat.get());
        break;
      case NegativeZFace:
        BOOST_CHECK_EQUAL(homMaterial, negZMat.get());
        break;
      case PositiveZFace:
        BOOST_CHECK_EQUAL(homMaterial, posZMat.get());
        break;
    }
  }
}

BOOST_AUTO_TEST_CASE(HomogeneousMaterialMixedVolumeTypes) {
  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {1_mm, 2_mm};
  Blueprint root{cfg};

  auto testMaterial = Acts::Material::fromMolarDensity(
      9.370_cm, 46.52_cm, 28.0855, 14, (2.329 / 28.0855) * 1_mol / 1_cm3);

  auto material = std::make_shared<HomogeneousSurfaceMaterial>(
      MaterialSlab(testMaterial, 0.1_mm));

  // Configure for cylinder first, then try to add cuboid - should throw
  BOOST_CHECK_THROW(
      root.addMaterial("Material",
                       [&](auto& mat) {
                         mat.configureFace(
                             CylinderVolumeBounds::Face::NegativeDisc,
                             material);
                         mat.configureFace(
                             CuboidVolumeBounds::Face::NegativeXFace, material);
                       }),
      std::runtime_error);

  // Configure for cuboid first, then try to add cylinder - should throw
  BOOST_CHECK_THROW(
      root.addMaterial("Material",
                       [&](auto& mat) {
                         mat.configureFace(
                             CuboidVolumeBounds::Face::NegativeXFace, material);
                         mat.configureFace(
                             CylinderVolumeBounds::Face::NegativeDisc,
                             material);
                       }),
      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(LayerCenterOfGravity) {
  // Test disc layer with center of gravity disabled
  {
    double yrot = 45_degree;
    Transform3 base =
        Transform3::Identity() * AngleAxis3{yrot, Vector3::UnitY()};

    std::vector<std::shared_ptr<Surface>> surfaces;
    std::vector<std::unique_ptr<DetectorElementBase>> elements;
    double r = 300_mm;
    std::size_t nSensors = 8;
    double thickness = 2.5_mm;
    auto recBounds = std::make_shared<RectangleBounds>(40_mm, 60_mm);

    double deltaPhi = 2 * std::numbers::pi / nSensors;
    for (std::size_t i = 0; i < nSensors; i++) {
      Transform3 trf = base * AngleAxis3{deltaPhi * i, Vector3::UnitZ()} *
                       Translation3(Vector3::UnitX() * r);

      if (i % 2 == 0) {
        trf = trf * Translation3{Vector3::UnitZ() * 5_mm};
      }

      auto& element = elements.emplace_back(
          std::make_unique<DetectorElementStub>(trf, recBounds, thickness));

      element->surface().assignDetectorElement(*element);
      surfaces.push_back(element->surface().getSharedPtr());
    }

    Blueprint root{{.envelope = ExtentEnvelope{{
                        .z = {2_mm, 2_mm},
                        .r = {3_mm, 5_mm},
                    }}}};

    root.addLayer("Layer0", [&](auto& layer) {
      layer.setSurfaces(surfaces)
          .setLayerType(LayerBlueprintNode::LayerType::Disc)
          .setEnvelope(ExtentEnvelope{{
              .z = {0.1_mm, 0.1_mm},
              .r = {1_mm, 1_mm},
          }})
          .setTransform(base)
          .setUseCenterOfGravity(false, false, false);  // Disable all axes
    });

    auto trackingGeometry = root.construct({}, gctx, *logger);
    auto lookup = nameLookup(*trackingGeometry);
    auto layerCyl = dynamic_cast<const CylinderVolumeBounds&>(
        lookup("Layer0").volumeBounds());

    // With center of gravity disabled, the layer should be at the origin
    BOOST_CHECK_CLOSE(layerCyl.get(CylinderVolumeBounds::eMinR), 258.9999999_mm,
                      1e-6);
    BOOST_CHECK_CLOSE(layerCyl.get(CylinderVolumeBounds::eMaxR),
                      346.25353003_mm, 1e-6);
    BOOST_CHECK_CLOSE(layerCyl.get(CylinderVolumeBounds::eHalfLengthZ), 3.85_mm,
                      1e-6);
  }

  // Test cylinder layer with center of gravity disabled
  {
    double yrot = 0_degree;
    Transform3 base =
        Transform3::Identity() * AngleAxis3{yrot, Vector3::UnitY()};

    std::vector<std::shared_ptr<Surface>> surfaces;
    std::vector<std::unique_ptr<DetectorElementBase>> elements;

    double r = 300_mm;
    std::size_t nStaves = 10;
    int nSensorsPerStave = 8;
    double thickness = 0;
    double hlPhi = 40_mm;
    double hlZ = 60_mm;
    auto recBounds = std::make_shared<RectangleBounds>(hlPhi, hlZ);

    double deltaPhi = 2 * std::numbers::pi / nStaves;

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
        element->surface().assignDetectorElement(*element);
        surfaces.push_back(element->surface().getSharedPtr());
      }
    }

    Blueprint root{{.envelope = ExtentEnvelope{{
                        .z = {2_mm, 2_mm},
                        .r = {3_mm, 5_mm},
                    }}}};

    root.addLayer("Layer0", [&](auto& layer) {
      layer.setSurfaces(surfaces)
          .setLayerType(LayerBlueprintNode::LayerType::Cylinder)
          .setEnvelope(ExtentEnvelope{{
              .z = {10_mm, 10_mm},
              .r = {20_mm, 10_mm},
          }})
          .setTransform(base)
          .setUseCenterOfGravity(false, false, false);  // Disable all axes
    });

    auto trackingGeometry = root.construct({}, gctx, *logger);
    auto lookup = nameLookup(*trackingGeometry);
    auto layerCyl = dynamic_cast<const CylinderVolumeBounds&>(
        lookup("Layer0").volumeBounds());

    // With center of gravity disabled, the layer should be at the origin
    BOOST_CHECK_EQUAL(lookup("Layer0").portals().size(), 4);
    BOOST_CHECK_CLOSE(layerCyl.get(CylinderVolumeBounds::eMinR), 275.6897761_mm,
                      1e-6);
    BOOST_CHECK_CLOSE(layerCyl.get(CylinderVolumeBounds::eMaxR), 319.4633358_mm,
                      1e-6);
    BOOST_CHECK_CLOSE(layerCyl.get(CylinderVolumeBounds::eHalfLengthZ), 1070_mm,
                      1e-6);
  }
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE(ProcessorBlueprintNodeTest);

BOOST_AUTO_TEST_CASE(BasicFunctionality) {
  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {1_mm, 2_mm};
  Blueprint root{cfg};

  double hlZ = 30_mm;
  auto cylBounds = std::make_shared<CylinderVolumeBounds>(10_mm, 20_mm, hlZ);
  auto cyl = std::make_unique<TrackingVolume>(Transform3::Identity(), cylBounds,
                                              "child");

  bool buildCalled = false;
  bool connectCalled = false;
  bool finalizeCalled = false;

  root.withProcessor([&](auto& processor) {
    processor.setName("Processor")
        .onBuild([&](Volume& vol) -> Volume& {
          buildCalled = true;
          return vol;
        })
        .onConnect([&](const BlueprintOptions& /*opts*/,
                       const GeometryContext& /*gctx*/,
                       const Logger& /*logger*/,
                       PortalShellBase& shell) -> PortalShellBase& {
          connectCalled = true;
          return shell;
        })
        .onFinalize([&](const BlueprintOptions& /*opts*/,
                        const GeometryContext& /*gctx*/,
                        TrackingVolume& /*parent*/,
                        const Logger& /*logger*/) { finalizeCalled = true; });
    // Add a static node as child
    processor.addChild(std::make_shared<StaticBlueprintNode>(std::move(cyl)));
  });

  auto trackingGeometry = root.construct({}, gctx, *logger);

  BOOST_CHECK(buildCalled);
  BOOST_CHECK(connectCalled);
  BOOST_CHECK(finalizeCalled);

  auto lookup = nameLookup(*trackingGeometry);
  BOOST_CHECK_NO_THROW(lookup("child"));
}

BOOST_AUTO_TEST_CASE(ErrorHandling) {
  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {1_mm, 2_mm};
  Blueprint root{cfg};

  auto& processor = root.withProcessor();
  processor.setName("Processor");

  // Test error when no child is added
  BOOST_CHECK_THROW(root.construct({}, gctx, *logger), std::runtime_error);

  // Test error when multiple children are added
  double hlZ = 30_mm;
  auto cylBounds = std::make_shared<CylinderVolumeBounds>(10_mm, 20_mm, hlZ);
  auto cyl1 = std::make_unique<TrackingVolume>(Transform3::Identity(),
                                               cylBounds, "child1");
  auto cyl2 = std::make_unique<TrackingVolume>(Transform3::Identity(),
                                               cylBounds, "child2");

  processor.addChild(std::make_shared<StaticBlueprintNode>(std::move(cyl1)));
  processor.addChild(std::make_shared<StaticBlueprintNode>(std::move(cyl2)));

  BOOST_CHECK_THROW(root.construct({}, gctx, *logger), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(VolumeModification) {
  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {1_mm, 2_mm};
  Blueprint root{cfg};

  double hlZ = 30_mm;
  auto cylBounds = std::make_shared<CylinderVolumeBounds>(10_mm, 20_mm, hlZ);
  auto cyl = std::make_unique<TrackingVolume>(Transform3::Identity(), cylBounds,
                                              "child");

  auto processor = std::make_shared<Experimental::ProcessorBlueprintNode>();
  processor->setName("Processor").onBuild([](Volume& vol) -> Volume& {
    auto& tVol = dynamic_cast<TrackingVolume&>(vol);
    tVol.setVolumeName("modified_child");
    return vol;
  });

  // Add a static node as child
  processor->addChild(std::make_shared<StaticBlueprintNode>(std::move(cyl)));

  root.addChild(processor);

  auto trackingGeometry = root.construct({}, gctx, *logger);

  auto lookup = nameLookup(*trackingGeometry);
  BOOST_CHECK_THROW(lookup("child"), std::runtime_error);
  BOOST_CHECK_NO_THROW(lookup("modified_child"));
}

BOOST_AUTO_TEST_CASE(StreamOperator) {
  auto processor = std::make_shared<Experimental::ProcessorBlueprintNode>();
  processor->setName("TestProcessor");

  std::stringstream ss;
  ss << *processor;
  BOOST_CHECK_EQUAL(ss.str(), "ProcessorBlueprintNode(TestProcessor)");
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE_END();

}  // namespace Acts::Test
