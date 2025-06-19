// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorSurfaceFactory.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <fstream>
#include <iostream>

#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hepTestsHelper.hpp"
#include "XML/Utilities.h"
#include "XMLFragments.hpp"

namespace {
Acts::GeometryContext tContext;
Acts::Test::CylindricalTrackingGeometry cGeometry =
    Acts::Test::CylindricalTrackingGeometry(tContext);
}  // namespace

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

const std::string indent_4_xml(4, ' ');
const std::string indent_8_xml(8, ' ');
const std::string indent_12_xml(12, ' ');

namespace {

Acts::Test::CylindricalTrackingGeometry::DetectorStore generateXML() {
  Acts::Test::CylindricalTrackingGeometry::DetectorStore dStore;

  // Nec surfaces
  double necZ = -800.;
  auto necR0Surfaces = cGeometry.surfacesRing(dStore, 6.4, 12.4, 18., 0.125, 0.,
                                              42., necZ, 2., 22u);

  auto necR1Surfaces = cGeometry.surfacesRing(dStore, 12.4, 20.4, 30., 0.125,
                                              0., 80., necZ, 2., 22u);

  std::vector<std::vector<const Acts::Surface*>> necSurfaces = {necR0Surfaces,
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

  std::vector<std::vector<const Acts::Surface*>> barrelSurfaces = {
      b0Surfaces, b1Surfaces, b2Surfaces};

  // Nec surfaces
  double pecZ = 800.;
  auto pecR0Surfaces = cGeometry.surfacesRing(dStore, 6.4, 12.4, 18., 0.125, 0.,
                                              42., pecZ, 2., 22u);

  auto pecR1Surfaces = cGeometry.surfacesRing(dStore, 12.4, 20.4, 30., 0.125,
                                              0., 80., pecZ, 2., 22u);

  std::vector<std::vector<const Acts::Surface*>> pecSurfaces = {pecR0Surfaces,
                                                                pecR1Surfaces};

  // Create an XML from it
  std::ofstream cxml;
  cxml.open("CylindricalDetectorTestSF.xml");
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
                                              Acts::Transform3::Identity())
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
  for (const auto [ib, bs] : Acts::enumerate(barrelSurfaces)) {
    cxml << indent_4_xml << "<layer name=\"PixelBarrel_" << ib << "\" id=\""
         << ib << "\">" << '\n';
    cxml << indent_4_xml << "<acts_volume rmin=\"" << innerOuter[ib][0u]
         << "*mm\" rmax=\"" << innerOuter[ib][1u] << "*mm\"/>";
    cxml << indent_8_xml << "<modules>" << '\n';
    for (const auto& s : bs) {
      cxml << indent_12_xml
           << DD4hepTestsHelper::surfaceToXML(tContext, *s,
                                              Acts::Transform3::Identity())
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
                                              Acts::Transform3::Identity())
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

  return dStore;
}

}  // namespace

auto store = generateXML();

BOOST_AUTO_TEST_SUITE(DD4hepPlugin)

BOOST_AUTO_TEST_CASE(ConvertSensitivesDefault) {
  auto lcdd = &(dd4hep::Detector::getInstance());
  lcdd->fromCompact("CylindricalDetectorTestSF.xml");
  lcdd->volumeManager();
  lcdd->apply("DD4hepVolumeManager", 0, nullptr);

  auto world = lcdd->world();

  // Test starts here - with nonimal detector construction
  Acts::DD4hepDetectorSurfaceFactory::Config sFactoryConfig;
  auto surfaceFactory = Acts::DD4hepDetectorSurfaceFactory(
      sFactoryConfig, Acts::getDefaultLogger("DD4hepDetectorSurfaceFactory",
                                             Acts::Logging::VERBOSE));

  Acts::DD4hepDetectorSurfaceFactory::Cache sFactoryCache;
  Acts::DD4hepDetectorSurfaceFactory::Options sFactoryOptions;

  surfaceFactory.construct(sFactoryCache, tContext, world, sFactoryOptions);
  // Check the number of surfaces
  BOOST_CHECK_EQUAL(sFactoryCache.sensitiveSurfaces.size(), 1488u);

  // Kill that instance before going into the next test
  lcdd->destroyInstance();
}

BOOST_AUTO_TEST_CASE(ConvertSensitivesextended) {
  auto lcdd = &(dd4hep::Detector::getInstance());
  lcdd->fromCompact("CylindricalDetectorTestSF.xml");
  lcdd->volumeManager();
  lcdd->apply("DD4hepVolumeManager", 0, nullptr);

  auto world = lcdd->world();

  // A typical extension would be overriding the `tranform(const
  // GeometryContext&)` method in order change how the detector element is
  // handled in alignment, for simplicity here we show a simple extension that
  // overrides the  thickness
  class ExtendedDetectorElement : public Acts::DD4hepDetectorElement {
   public:
    using Acts::DD4hepDetectorElement::DD4hepDetectorElement;

    double thickness() const final {
      // Return a fixed thickness for testing purposes
      return 42. * Acts::UnitConstants::mm;
    }
  };

  auto extendedFactory =
      [](const dd4hep::DetElement& detElem, const std::string& axes,
         double scalor, bool isDisc,
         const std::shared_ptr<const Acts::ISurfaceMaterial>& material)
      -> std::shared_ptr<Acts::DD4hepDetectorElement> {
    return std::make_shared<ExtendedDetectorElement>(detElem, axes, scalor,
                                                     isDisc, material);
  };

  // Test starts here - with nonimal detector construction
  Acts::DD4hepDetectorSurfaceFactory::Config sFactoryConfig;
  sFactoryConfig.detectorElementFactory = extendedFactory;
  auto surfaceFactory = Acts::DD4hepDetectorSurfaceFactory(
      sFactoryConfig, Acts::getDefaultLogger("DD4hepDetectorSurfaceFactory",
                                             Acts::Logging::VERBOSE));

  Acts::DD4hepDetectorSurfaceFactory::Cache sFactoryCache;
  Acts::DD4hepDetectorSurfaceFactory::Options sFactoryOptions;

  surfaceFactory.construct(sFactoryCache, tContext, world, sFactoryOptions);
  // Check the number of surfaces
  BOOST_CHECK_EQUAL(sFactoryCache.sensitiveSurfaces.size(), 1488u);
  for (const auto& [detElem, surface] : sFactoryCache.sensitiveSurfaces) {
    // Check that the extended detector element is used
    BOOST_CHECK_NE(dynamic_cast<const ExtendedDetectorElement*>(detElem.get()),
                   nullptr);
    // Check that the thickness is 42 mm
    CHECK_CLOSE_ABS(detElem->thickness(), 42. * Acts::UnitConstants::mm, 1e-10);
  }

  // Kill that instance before going into the next test
  lcdd->destroyInstance();
}

BOOST_AUTO_TEST_SUITE_END()
