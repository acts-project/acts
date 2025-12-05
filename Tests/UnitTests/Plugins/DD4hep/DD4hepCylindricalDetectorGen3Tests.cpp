// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifierBlueprintNode.hpp"
#include "Acts/Geometry/LayerBlueprintNode.hpp"
#include "Acts/Geometry/MaterialDesignatorBlueprintNode.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/TransformRange.hpp"
#include "ActsPlugins/DD4hep/DD4hepConversionHelpers.hpp"
#include "ActsPlugins/DD4hep/DD4hepDetectorElement.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "ActsTests/CommonHelpers/TemporaryDirectory.hpp"

#include <format>
#include <fstream>
#include <string>

#include <DD4hep/DetElement.h>
#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/Detector.h>
#include <XML/Utilities.h>
#include <XMLFragments.hpp>

#include "DD4hepTestsHelper.hpp"

using namespace Acts;
using namespace ActsPlugins;

namespace ActsTests {

GeometryContext tContext;
CylindricalTrackingGeometry cGeometry = CylindricalTrackingGeometry(tContext);

const char* beampipe_head_xml =
    R""""(
    <detectors>
        <detector id="0" name="BeamPipe" type="BarrelDetector">
            <type_flags type="DetType_TRACKER + DetType_BEAMPIPE"/>
            <layers>
                <layer id="0"  name="BeamPipeLayer">
)"""";

const char* nec_head_xml =
    R""""(
        <detector id="1" name="PixelNegativeEndcap" type="EndcapDetector" readout="PixelReadout">
            <type_flags type="DetType_TRACKER + DetType_BARREL"/>
)"""";

const char* barrel_head_xml =
    R""""(
        <detector id="2" name="PixelBarrel" type="BarrelDetector" readout="PixelReadout">
            <type_flags type="DetType_TRACKER + DetType_BARREL"/>
)"""";

const char* pec_head_xml =
    R""""(
        <detector id="3" name="PixelPositiveEndcap" type="EndcapDetector" readout="PixelReadout">
            <type_flags type="DetType_TRACKER + DetType_BARREL"/>
)"""";

const char* plugin_xml =
    R"""(<plugins>
    <plugin name="DD4hep_ParametersPlugin">
      <argument value="world"/>
      <argument value="acts_world: bool = true"/>
      <argument value="acts_world_type: int = 3"/>
      <argument value="acts_world_bvalues_n: int = 3"/>
      <argument value="acts_world_bvalues_0: double = 0"/>
      <argument value="acts_world_bvalues_1: double = 150*mm"/>
      <argument value="acts_world_bvalues_2: double = 1000*mm"/>
      <argument value="acts_world_binning: str = r"/>
      <argument value="acts_world_geo_id: str = incremental"/>
      <argument value="acts_world_root_volume_finder: str = indexed"/>
    </plugin>
    <plugin name="DD4hep_ParametersPlugin">
      <argument value="BeamPipe"/>
      <argument value="acts_volume: bool = true"/>
      <argument value="acts_volume_type: int = 3"/>
      <argument value="acts_volume_bvalues_n: int = 3"/>
      <argument value="acts_volume_bvalues_0: double = 0"/>
      <argument value="acts_volume_bvalues_1: double = 20*mm"/>
      <argument value="acts_volume_bvalues_2: double = 1000*mm"/>
      <argument value="acts_volume_internals: bool = true"/>
      <argument value="acts_volume_internals_type: str = layer"/>
    </plugin>
    <plugin name="DD4hep_ParametersPlugin">
      <argument value="Pixel"/>
      <argument value="acts_container: bool = true"/>
      <argument value="acts_container_type: int = 3"/>
      <argument value="acts_container_bvalues_n: int = 3"/>
      <argument value="acts_container_bvalues_0: double = 20*mm"/>
      <argument value="acts_container_bvalues_1: double = 150*mm"/>
      <argument value="acts_container_bvalues_2: double = 1000*mm"/>
      <argument value="acts_container_binning: str = z"/>
    </plugin>
    <plugin name="DD4hep_ParametersPlugin">
      <argument value="PixelNegativeEndcap"/>
      <argument value="acts_container: bool = true"/>
      <argument value="acts_container_type: int = 3"/>
      <argument value="acts_container_bvalues_n: int = 3"/>
      <argument value="acts_container_bvalues_0: double = 20*mm"/>
      <argument value="acts_container_bvalues_1: double = 150*mm"/>
      <argument value="acts_container_bvalues_2: double = 200*mm"/>
      <argument value="acts_container_binning: str = z"/>
      <argument value="acts_container_z: double = -800*mm"/>
    </plugin>
    <plugin name="DD4hep_ParametersPlugin">
      <argument value="PixelBarrel"/>
      <argument value="acts_container: bool = true"/>
      <argument value="acts_container_type: int = 3"/>
      <argument value="acts_container_bvalues_n: int = 3"/>
      <argument value="acts_container_bvalues_0: double = 20*mm"/>
      <argument value="acts_container_bvalues_1: double = 150*mm"/>
      <argument value="acts_container_bvalues_2: double = 600*mm"/>
      <argument value="acts_container_binning: str = r"/>
    </plugin>
    <plugin name="DD4hep_ParametersPlugin">
      <argument value="PixelPositiveEndcap"/>
      <argument value="acts_container: bool = true"/>
      <argument value="acts_container_type: int = 3"/>
      <argument value="acts_container_bvalues_n: int = 3"/>
      <argument value="acts_container_bvalues_0: double = 20*mm"/>
      <argument value="acts_container_bvalues_1: double = 150*mm"/>
      <argument value="acts_container_bvalues_2: double = 200*mm"/>
      <argument value="acts_container_binning: str = z"/>
      <argument value="acts_container_z: double = 800*mm"/>
    </plugin>
  </plugins>
  )""";

const char* indent_4_xml = "    ";
const char* indent_8_xml = "        ";
const char* indent_12_xml = "            ";

void generateXML(const std::filesystem::path& xmlPath) {
  CylindricalTrackingGeometry::DetectorStore dStore;

  // Nec surfaces
  double necZ = -800.;
  auto necR0Surfaces = cGeometry.surfacesRing(dStore, 6.4, 12.4, 18., 0.125, 0.,
                                              42., necZ, 2., 22u);

  auto necR1Surfaces = cGeometry.surfacesRing(dStore, 12.4, 20.4, 30., 0.125,
                                              0., 80., necZ, 2., 22u);

  std::vector<std::vector<Surface*>> necSurfaces = {necR0Surfaces,
                                                    necR1Surfaces};

  // Barrel surfaces
  std::vector<std::array<double, 2u>> innerOuter = {
      {25., 35.}, {65., 75.}, {110., 120.}};
  auto b0Surfaces = cGeometry.surfacesCylinder(dStore, 8.4, 36., 0.15, 0.14,
                                               31., 3., 2., {16, 14});

  auto b1Surfaces = cGeometry.surfacesCylinder(dStore, 8.4, 36., 0.15, 0.14,
                                               71., 3., 2., {32, 14});

  auto b2Surfaces = cGeometry.surfacesCylinder(dStore, 8.4, 36., 0.15, 0.14,
                                               116., 3., 2., {52, 14});

  std::vector<std::vector<Surface*>> barrelSurfaces = {b0Surfaces, b1Surfaces,
                                                       b2Surfaces};

  // Nec surfaces
  double pecZ = 800.;
  auto pecR0Surfaces = cGeometry.surfacesRing(dStore, 6.4, 12.4, 18., 0.125, 0.,
                                              42., pecZ, 2., 22u);

  auto pecR1Surfaces = cGeometry.surfacesRing(dStore, 12.4, 20.4, 30., 0.125,
                                              0., 80., pecZ, 2., 22u);

  std::vector<std::vector<Surface*>> pecSurfaces = {pecR0Surfaces,
                                                    pecR1Surfaces};

  // Create an XML from it
  std::ofstream cxml;
  cxml.open(xmlPath);
  cxml << head_xml;
  cxml << segmentation_xml;

  // Add the beam pipe first
  cxml << beampipe_head_xml << '\n';
  cxml << indent_12_xml << "<acts_passive_surface>" << '\n';
  cxml << indent_12_xml
       << "<tubs rmin=\"19*mm\" rmax=\"20*mm\" dz=\"800*mm\" "
          "material=\"Air\"/>"
       << '\n';
  cxml << indent_12_xml << "</acts_passive_surface>" << '\n';
  cxml << indent_12_xml << "</layer> " << '\n';
  cxml << indent_8_xml << "</layers>" << '\n';
  cxml << indent_8_xml << "</detector>" << '\n';

  // Declare the assembly
  const char* pixel_assembly_xml =
      R"""(<detector id="4" name="Pixel" type="DD4hep_SubdetectorAssembly" vis="invisible">
            <shape name="PixelVolume" type="Tube" rmin="20*mm" rmax="150*mm" dz="1000*mm" material="Air"/>
            <composite name="PixelNegativeEndcap"/>
            <composite name="PixelBarrel"/>
            <composite name="PixelPositiveEndcap"/>
        </detector>)""";

  cxml << pixel_assembly_xml << '\n';
  cxml << "</detectors>" << '\n';

  // The Pixel detector

  // Nec - 1 Layer only
  cxml << "<detectors>" << '\n';

  cxml << nec_head_xml << '\n';
  cxml << indent_8_xml << "<layers>" << '\n';
  cxml << indent_8_xml << "<acts_container/> " << '\n';
  cxml << indent_4_xml << "<layer name=\"NegEndcapLayer_0\" id=\"0\">" << '\n';
  cxml << indent_4_xml << "<acts_volume dz=\"20*mm\" cz=\"-800.*mm\"/>" << '\n';
  cxml << indent_8_xml << "<modules>" << '\n';
  for (const auto& ring : necSurfaces) {
    for (const auto& s : ring) {
      cxml << indent_12_xml
           << DD4hepTestsHelper::surfaceToXML(tContext, *s,
                                              Transform3::Identity())
           << "\n";
    }
  }
  cxml << indent_8_xml << "</modules>" << '\n';
  cxml << indent_8_xml << "</layer> " << '\n';
  cxml << indent_8_xml << "</layers>" << '\n';
  cxml << indent_8_xml << "</detector>" << '\n';

  // Barrel
  cxml << barrel_head_xml << '\n';
  cxml << indent_8_xml << "<layers>" << '\n';
  cxml << indent_8_xml << "<acts_container/> " << '\n';
  for (const auto [ib, bs] : enumerate(barrelSurfaces)) {
    cxml << indent_4_xml << "<layer name=\"PixelBarrel_" << ib << "\" id=\""
         << ib << "\">" << '\n';
    cxml << indent_4_xml << "<acts_volume rmin=\"" << innerOuter[ib][0u]
         << "*mm\" rmax=\"" << innerOuter[ib][1u] << "*mm\"/>";
    cxml << indent_8_xml << "<modules>" << '\n';
    for (const auto& s : bs) {
      cxml << indent_12_xml
           << DD4hepTestsHelper::surfaceToXML(tContext, *s,
                                              Transform3::Identity())
           << "\n";
    }
    cxml << indent_8_xml << "</modules>" << '\n';
    cxml << indent_8_xml << "</layer>" << '\n';
  }
  cxml << indent_8_xml << "</layers>" << '\n';
  cxml << indent_8_xml << "</detector>" << '\n';

  // Pec - 1 layer only
  cxml << pec_head_xml << '\n';
  cxml << indent_8_xml << "<layers>" << '\n';
  cxml << indent_8_xml << "<acts_container/> " << '\n';
  cxml << indent_4_xml << "<layer name=\"PosEndcapLayer_0\" id=\"0\">" << '\n';
  cxml << indent_4_xml << "<acts_volume dz=\"20*mm\" cz=\"800.*mm\"/>" << '\n';
  cxml << indent_8_xml << "<modules>" << '\n';
  for (const auto& ring : pecSurfaces) {
    for (const auto& s : ring) {
      cxml << indent_12_xml
           << DD4hepTestsHelper::surfaceToXML(tContext, *s,
                                              Transform3::Identity())
           << "\n";
    }
  }
  cxml << indent_8_xml << "</modules>" << '\n';
  cxml << indent_8_xml << "</layer> " << '\n';
  cxml << indent_8_xml << "</layers>" << '\n';
  cxml << indent_8_xml << "</detector>" << '\n';
  cxml << indent_8_xml << "</detectors>" << '\n';

  cxml << plugin_xml << '\n';
  cxml << end_xml << '\n';
  cxml.close();
}

using namespace dd4hep;
using namespace UnitLiterals;

using enum AxisBoundaryType;
using enum AxisDirection;
using enum CylinderVolumeBounds::Face;
using enum SurfaceArrayNavigationPolicy::LayerType;
using AttachmentStrategy = VolumeAttachmentStrategy;
using ResizeStrategy = VolumeResizeStrategy;

// --------- Helper functions --------------
// Print Element Tree
void print_elements(const DetElement& el, unsigned int level = 0) {
  for (const auto& [name, child] : el.children()) {
    for (unsigned int i = 0; i < level; ++i) {
      std::cout << "\t";
    }
    std::cout << "-> " << name << std::endl;
    print_elements(child, level + 1);
  }
}

// Find the first element with a given name
const DetElement* find_element(const DetElement& el, const std::string& el_name,
                               unsigned int search_depth = 0) {
  static unsigned int level = 0;
  static bool abort_search = false;
  for (const auto& [name, child] : el.children()) {
    if (el_name == name) {
      abort_search = true;
      return &child;
    }
    if (++level <= search_depth) {
      auto el_child = find_element(child, el_name, search_depth);
      if (abort_search) {
        return el_child;
      }
    }
  }
  return nullptr;
}

// Helper function to convert shared_ptr vector to const ptr vector
std::vector<const Surface*> makeConstPtrVector(
    const std::vector<std::shared_ptr<Surface>>& surfs) {
  std::vector<const Surface*> constPtrs;
  constPtrs.reserve(surfs.size());
  for (const auto& surf : surfs) {
    constPtrs.push_back(surf.get());
  }
  return constPtrs;
}

// Helper struct to keep ProtoLayer and its associated surfaces together
struct LayerData {
  ProtoLayer protoLayer;
  std::vector<std::shared_ptr<Surface>> surfaces;

  LayerData(const GeometryContext& gctx,
            std::vector<std::shared_ptr<Surface>> surfs)
      : protoLayer(gctx, makeConstPtrVector(surfs)),
        surfaces(std::move(surfs)) {}
};

// Helper function to merge layers that overlap in z
std::vector<LayerData> mergeLayers(const GeometryContext& gctx,
                                   std::vector<LayerData> layers) {
  using enum AxisDirection;

  std::vector<LayerData> mergedLayers;
  if (layers.empty()) {
    return mergedLayers;
  }

  mergedLayers.push_back(std::move(layers.front()));

  for (std::size_t i = 1; i < layers.size(); i++) {
    auto& current = layers[i];
    auto& prev = mergedLayers.back();

    // Check if they overlap in z
    bool overlap =
        (current.protoLayer.min(AxisZ) <= prev.protoLayer.max(AxisZ) &&
         current.protoLayer.max(AxisZ) >= prev.protoLayer.min(AxisZ));

    if (overlap) {
      // Merge surfaces
      std::vector<std::shared_ptr<Surface>> mergedSurfaces;
      mergedSurfaces.reserve(current.surfaces.size() + prev.surfaces.size());
      mergedSurfaces.insert(mergedSurfaces.end(), current.surfaces.begin(),
                            current.surfaces.end());
      mergedSurfaces.insert(mergedSurfaces.end(), prev.surfaces.begin(),
                            prev.surfaces.end());

      mergedLayers.pop_back();
      mergedLayers.emplace_back(gctx, std::move(mergedSurfaces));
      auto& merged = mergedLayers.back();
      merged.protoLayer.envelope[AxisR] = current.protoLayer.envelope[AxisR];
      merged.protoLayer.envelope[AxisZ] = current.protoLayer.envelope[AxisZ];
    } else {
      mergedLayers.push_back(std::move(current));
    }
  }

  return mergedLayers;
}

// --------- Helper functions --------------

BOOST_AUTO_TEST_SUITE(DD4hepPlugin)

BOOST_AUTO_TEST_CASE(DD4hepCylidricalDetectorExplicit) {
  TemporaryDirectory tempDir{};
  auto xmlPath = tempDir.path() / "CylindricalDetector.xml";
  generateXML(xmlPath);

  auto lcdd = &(dd4hep::Detector::getInstance());
  lcdd->fromCompact(xmlPath.string());
  lcdd->volumeManager();
  lcdd->apply("DD4hepVolumeManager", 0, nullptr);

  constexpr std::size_t s_beamPipeVolumeId =
      1;  // where do volume IDs come from?
  constexpr std::size_t s_pixelVolumeId = 10;

  DetElement world = lcdd->world();

  // std::cout << "DD4Hep Elements: " << std::endl;
  // print_elements(world);

  auto worldSolidDim = lcdd->worldVolume().solid().dimensions();  // better way?

  Experimental::Blueprint::Config cfg;
  // Is the following correct?
  cfg.envelope[AxisX] = {worldSolidDim[0], worldSolidDim[0]};
  cfg.envelope[AxisY] = {worldSolidDim[1], worldSolidDim[1]};
  cfg.envelope[AxisZ] = {worldSolidDim[2], worldSolidDim[2]};
  // cfg.envelope[AxisR] = {0_mm, worldSolidDim[1]};

  auto blueprint = std::make_unique<Experimental::Blueprint>(cfg);
  auto& cylinder = blueprint->addCylinderContainer("Detector", AxisR);

  // ------- Add Beam Pipe to Blueprint -------
  auto level1_elements = world.children();
  DetElement bpipe_top = level1_elements["BeamPipe"];
  auto bpipe_top_position = bpipe_top.placement().position();

  Transform3 beamPipeTransform;
  beamPipeTransform.setIdentity();
  beamPipeTransform = Translation3(
      bpipe_top_position.x(), bpipe_top_position.y(), bpipe_top_position.z());

  auto beamPipeDim = bpipe_top.solid().dimensions();
  // Do I use the values retrieved from DD4Hep correctly?
  double beamPipeRMax = beamPipeDim[1] * 1_cm;
  double beamPipeHalfZ = beamPipeDim[2] * 1_cm;

  std::vector<std::shared_ptr<DD4hepDetectorElement>> detectorElements;

  cylinder.withGeometryIdentifier([&](auto& geoId) {
    geoId.setAllVolumeIdsTo(s_beamPipeVolumeId);
    geoId.addMaterial("BeamPipe_Material", [&](auto& mat) {
      mat.configureFace(
          OuterCylinder,
          {AxisRPhi, Bound,
           20},  // Where do these 20-s come from? (here and on the next line)
          {AxisZ, Bound, 20});
      mat.addStaticVolume(beamPipeTransform,
                          std::make_shared<CylinderVolumeBounds>(
                              0, beamPipeRMax, beamPipeHalfZ),
                          "BeamPipe");
    });
  });
  // ------- Add Beam Pipe to Blueprint -------

  // ------- Add Pixel to Blueprint -------
  auto pixelElement = find_element(world, "Pixel");
  auto pixelBarrelElement = find_element(*pixelElement, "PixelBarrel");

  cylinder.addMaterial("Pixel_Material", [&](auto& mat) {
    mat.configureFace(OuterCylinder, {AxisRPhi, Bound, 20}, {AxisZ, Bound, 20});
    auto& pixelContainer = mat.addCylinderContainer("Pixel", AxisZ);

    // Add barrel container
    auto& barrelGeoId = pixelContainer.withGeometryIdentifier();
    barrelGeoId.setAllVolumeIdsTo(s_pixelVolumeId)
        .incrementLayerIds(1)
        .sortBy([](auto& a, auto& b) {
          auto& boundsA =
              dynamic_cast<const CylinderVolumeBounds&>(a.volumeBounds());
          auto& boundsB =
              dynamic_cast<const CylinderVolumeBounds&>(b.volumeBounds());
          using enum CylinderVolumeBounds::BoundValues;
          double aMidR = (boundsA.get(eMinR) + boundsA.get(eMaxR)) / 2.0;
          double bMidR = (boundsB.get(eMinR) + boundsB.get(eMaxR)) / 2.0;
          return aMidR < bMidR;
        });

    auto& barrel = barrelGeoId.addCylinderContainer("Pixel_Barrel", AxisR);
    barrel.setAttachmentStrategy(AttachmentStrategy::Gap);
    barrel.setResizeStrategy(ResizeStrategy::Gap);

    std::map<int, std::vector<std::shared_ptr<Surface>>> layers{};
    int layerId{0};
    for (const auto& [nameLayer, layer] : pixelBarrelElement->children()) {
      for (const auto& [nameModule, module] : layer.children()) {
        std::string detAxis =
            getParamOr<std::string>("axis_definitions", module, "XYZ");
        auto dd4hepDetEl = std::make_shared<DD4hepDetectorElement>(
            module, detAxis, 1_cm, false, nullptr);
        detectorElements.push_back(dd4hepDetEl);
        layers[layerId].push_back(dd4hepDetEl->surface().getSharedPtr());
      }
      layerId++;
    }

    for (const auto& [ilayer, surfaces] : layers) {
      Experimental::BlueprintNode* lparent = nullptr;

      // Outermost layer can't have material, because it will get merged
      // with the outer cylinder shell of the endcap cylinders
      //
      //                          Material here is fine   |
      //                                                  |
      //                                                  |
      //                Will get merged. Therefore none   |
      //                   of them can have material      |
      //                                                  |
      //                               |                  |
      //                               |                  |
      //                               |                  |
      //       +-----------------------+-+----------------+---------+
      //       |                         |                |         |
      //       |                         |                |         |
      //       v                         v                |         v
      // +----------+-------------------------------------+---+----------+
      // |          |                Barrel L2            |   |          |
      // |          +-------------------------------------v---+          |
      // |          |                   Gap                   |          |
      // |          +-----------------------------------------+          |
      // |  Neg EC  |                Barrel L1                |  Pos EC  |
      // |          +-----------------------------------------+          |
      // |          |                   Gap                   |          |
      // |          +-----------------------------------------+          |
      // |          |                Barrel L0                |          |
      // +----------+-----------------------------------------+----------+
      if (ilayer < static_cast<int>(layers.size() - 1)) {
        auto& lmat = barrel.addMaterial(
            std::format("Pixel_Barrel_L{}_Material", ilayer));
        lmat.configureFace(OuterCylinder, {AxisRPhi, Bound, 40},
                           {AxisZ, Bound, 20});
        lparent = &lmat;
      } else {
        lparent = &barrel;
      }

      // Add layer with surfaces
      auto& layer = lparent->addLayer(std::format("Pixel_Barrel_L{}", ilayer));

      // Set navigation policy for efficient surface lookup
      layer.setNavigationPolicyFactory(
          NavigationPolicyFactory{}
              .add<SurfaceArrayNavigationPolicy>(
                  SurfaceArrayNavigationPolicy::Config{.layerType = Cylinder,
                                                       .bins = {30, 10}})
              .add<TryAllNavigationPolicy>(
                  TryAllNavigationPolicy::Config{.sensitives = false})
              .asUniquePtr());

      layer.setSurfaces(surfaces);
      layer.setEnvelope(ExtentEnvelope{{
          .z = {5_mm, 5_mm},  // ???
          .r = {2_mm, 2_mm},  // ???
      }});
    }

    // Add endcap containers
    for (int ecid : {-1, 1}) {
      const std::string_view s = ecid == 1 ? "p" : "n";
      auto& ecGeoId = pixelContainer.withGeometryIdentifier();
      ecGeoId.setAllVolumeIdsTo(s_pixelVolumeId + ecid).incrementLayerIds(1);
      auto& ec =
          ecGeoId.addCylinderContainer(std::format("Pixel_{}EC", s), AxisZ);
      ec.setAttachmentStrategy(AttachmentStrategy::Gap);
      ec.setResizeStrategy(ResizeStrategy::Expand);

      std::map<int, std::vector<std::shared_ptr<Surface>>> initialLayers{};
      const DetElement* pixelEndcapElement =
          ecid == 1 ? find_element(*pixelElement, "PixelPositiveEndcap")
                    : find_element(*pixelElement, "PixelNegativeEndcap");
      layerId = 0;
      for (const auto& [nameLayer, layer] : pixelEndcapElement->children()) {
        for (const auto& [nameModule, module] : layer.children()) {
          std::string detAxis =
              getParamOr<std::string>("axis_definitions", module, "XYZ");
          auto dd4hepDetEl = std::make_shared<DD4hepDetectorElement>(
              module, detAxis, 1_cm, false, nullptr);
          detectorElements.push_back(dd4hepDetEl);
          initialLayers[layerId].push_back(
              dd4hepDetEl->surface().getSharedPtr());
        }
        layerId++;
      }

      // Create proto layers from surfaces
      std::vector<LayerData> protoLayers;
      protoLayers.reserve(initialLayers.size());
      for (const auto& [key, surfaces] : initialLayers) {
        auto& layer = protoLayers.emplace_back(GeometryContext(), surfaces);
        layer.protoLayer.envelope[AxisR] = {2_mm, 2_mm};
        layer.protoLayer.envelope[AxisZ] = {1_mm, 1_mm};
      }
      // Sort by z position
      std::ranges::sort(protoLayers,
                        [](const LayerData& a, const LayerData& b) {
                          return std::abs(a.protoLayer.medium(AxisZ)) <
                                 std::abs(b.protoLayer.medium(AxisZ));
                        });

      std::vector<LayerData> mergedLayers =
          mergeLayers(GeometryContext(), protoLayers);

      // Create layers from proto layers
      for (const auto& [key, pl] : enumerate(mergedLayers)) {
        pl.protoLayer.medium(AxisZ);
        auto layerName = std::format("Pixel_{}EC_L{}", s, key);
        auto addLayer = [&layerName, &pl](auto& parent) {
          // Add layer with surfaces
          auto& layer = parent.addLayer(layerName);

          layer.setNavigationPolicyFactory(
              NavigationPolicyFactory{}
                  .add<SurfaceArrayNavigationPolicy>(
                      SurfaceArrayNavigationPolicy::Config{.layerType = Disc,
                                                           .bins = {30, 30}})
                  .add<TryAllNavigationPolicy>(
                      TryAllNavigationPolicy::Config{.sensitives = false})
                  .asUniquePtr());

          layer.setSurfaces(pl.surfaces);
          layer.setEnvelope(ExtentEnvelope{{
              .z = {1_mm, 1_mm},
              .r = {2_mm, 2_mm},
          }});
        };

        if (key < mergedLayers.size() - 1) {
          ec.addMaterial(layerName + "_Material", [&](auto& lmat) {
            lmat.configureFace(ecid < 0 ? NegativeDisc : PositiveDisc,
                               {AxisR, Bound, 40}, {AxisPhi, Bound, 40});
            addLayer(lmat);
          });
        } else {
          addLayer(ec);
        }
      }
    }
  });
  // ------- Add Pixel to Blueprint -------

  std::cout << tempDir.path() << std::endl;
  std::ofstream dotOut{tempDir.path() / "CylindricalDetector.dot"};
  blueprint->graphviz(dotOut);

  // Final step
  GeometryContext gctx;
  auto logger = getDefaultLogger("Geo", Logging::VERBOSE);
  auto trackingGeometry = blueprint->construct({}, gctx, *logger);

  BOOST_REQUIRE_NE(&world, nullptr);
  BOOST_REQUIRE_NE(trackingGeometry.get(), nullptr);

  // Kill that instance before going into the next test
  lcdd->destroyInstance();
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
