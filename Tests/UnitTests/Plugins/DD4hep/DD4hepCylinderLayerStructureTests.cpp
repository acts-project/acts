// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorSurfaceFactory.hpp"
#include "Acts/Plugins/DD4hep/DD4hepLayerStructure.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <tuple>

#include <DD4hep/DetElement.h>
#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/Detector.h>
#include <XML/Utilities.h>
#include <XMLFragments.hpp>

#include "DD4hepTestsHelper.hpp"

Acts::GeometryContext tContext;
Acts::Test::CylindricalTrackingGeometry cGeometry =
    Acts::Test::CylindricalTrackingGeometry(tContext);

// ******* TO BE REMOVED **********

#include "Acts/Plugins/ActSVG/IndexedSurfacesSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/SvgUtils.hpp"

Acts::Svg::IndexedSurfacesConverter::Options generateDrawOptions() {
  // The converter options
  Acts::Svg::IndexedSurfacesConverter::Options isOptions;
  // Sensitive surface stuyle
  Acts::Svg::Style sensitiveStyle;
  sensitiveStyle.fillColor = {51, 153, 255};
  sensitiveStyle.fillOpacity = 0.9;
  sensitiveStyle.highlightColor = {255, 153, 51};
  sensitiveStyle.highlights = {"onmouseover", "onmouseout"};
  sensitiveStyle.strokeWidth = 0.5;
  sensitiveStyle.strokeColor = {0, 0, 0};
  sensitiveStyle.nSegments = 72u;
  std::pair<Acts::GeometryIdentifier, Acts::Svg::Style> allSensitives = {
      Acts::GeometryIdentifier(0u), sensitiveStyle};

  // Hierarchy map of styles
  Acts::GeometryHierarchyMap<Acts::Svg::Style> surfaceStyles({allSensitives});
  isOptions.surfaceStyles = surfaceStyles;

  // The grid style
  Acts::Svg::GridConverter::Options gridOptions;
  Acts::Svg::Style gridStyle;
  gridStyle.fillOpacity = 0.;
  gridStyle.strokeColor = {0, 0, 255};
  gridStyle.strokeWidth = 1.;
  gridStyle.highlightStrokeWidth = 3;
  gridStyle.highlightStrokeColor = {255, 0, 0};
  gridOptions.style = gridStyle;

  isOptions.gridOptions = gridOptions;
  return isOptions;
};

auto drawOptions = generateDrawOptions();

// ****************

const char* cylinder_layer_head_xml =
    R""""(
    <detectors>
        <detector id="1" name="CylinderLayer" type="CylinderLayer" readout="PixelReadout">
            <type_flags type="DetType_TRACKER + DetType_BARREL"/>
)"""";

const char* cylinder_layer_tail_xml =
    R""""(
        </detector>
    </detectors>
)"""";

const char* indent_12_xml = "            ";

BOOST_AUTO_TEST_SUITE(DD4hepPlugin)

BOOST_AUTO_TEST_CASE(DD4hepPluginCylinderLayerStructure) {
  // First create some test surfaces
  Acts::Test::CylindricalTrackingGeometry::DetectorStore dStore;
  auto cSurfaces = cGeometry.surfacesCylinder(dStore, 8.4, 36., 0.15, 0.145,
                                              116., 3., 2., {52, 14});

  std::vector<std::array<unsigned int, 4u> > zphiBinning = {
      {14u, 52u, 1u, 1u}, {28u, 104u, 0u, 0u}};

  for (auto [nz, nphi, ez, ephi] : zphiBinning) {
    // Create an XML from it
    std::ofstream cxml;

    std::string fNameBase = "CylinderLayer_nz";
    fNameBase += std::to_string(nz);
    fNameBase += "_nphi";
    fNameBase += std::to_string(nphi);

    cxml.open(fNameBase + ".xml");
    cxml << head_xml;
    cxml << segmentation_xml;
    cxml << cylinder_layer_head_xml;

    cxml << indent_12_xml
         << "<envelope rmin=\"180*mm\" rmax=\"240*mm\" dz=\"800*cm\" "
            "cz=\"0*mm\" "
            "material=\"Air\"/>"
         << '\n';

    cxml << indent_12_xml << "<surface_binning";
    cxml << " ztype=\"equidistant\"";
    cxml << " phitype=\"equidistant\"";
    cxml << " nz=\"" << nz << "\" zmin=\"-500\" zmax=\"500\"";
    cxml << " zexpansion= \"" << ez << "\"";
    cxml << " nphi=\"" << nphi << "\"  phimin=\"-3.1415\" phimax=\"3.1415\"";
    cxml << " phiexpansion= \"" << ephi << "\"/>";
    cxml << "<modules>" << '\n';

    for (const auto& s : cSurfaces) {
      cxml << indent_12_xml
           << DD4hepTestsHelper::surfaceToXML(tContext, *s,
                                              Acts::Transform3::Identity())
           << "\n";
    }

    cxml << "</modules>" << '\n';
    cxml << cylinder_layer_tail_xml;
    cxml << end_xml;
    cxml.close();

    auto lcdd = &(dd4hep::Detector::getInstance());
    lcdd->fromCompact(fNameBase + ".xml");
    lcdd->volumeManager();
    lcdd->apply("DD4hepVolumeManager", 0, nullptr);

    auto world = lcdd->world();

    // Now the test starts ...
    auto sFactory = std::make_shared<Acts::DD4hepDetectorSurfaceFactory>(
        Acts::getDefaultLogger("DD4hepDetectorSurfaceFactory",
                               Acts::Logging::DEBUG));

    Acts::Experimental::DD4hepLayerStructure barrelStructure(
        std::move(sFactory),
        Acts::getDefaultLogger("DD4hepLayerStructure", Acts::Logging::VERBOSE));

    Acts::DD4hepDetectorElement::Store dd4hepStore;

    Acts::Experimental::DD4hepLayerStructure::Options lsOptions;
    lsOptions.name = "BarrelLayer";
    lsOptions.logLevel = Acts::Logging::INFO;

    auto barrelInternalsBuilder =
        barrelStructure.builder(dd4hepStore, world, lsOptions);

    // Build the internal volume structure
    auto [surfaces, volumes, surfacesUpdator, volumeUpdator] =
        barrelInternalsBuilder->construct(tContext);

    // All surfaces are filled
    BOOST_CHECK(surfaces.size() == 14u * 52u);
    // No volumes are added
    BOOST_CHECK(volumes.empty());
    // The surface updator is connected
    BOOST_CHECK(surfacesUpdator.connected());
    // The volume updator is connected
    BOOST_CHECK(volumeUpdator.connected());

    // *** TO BE REMOVED ****
    // Draw the thing

    // The displaying
    auto pIndexeCylinder = Acts::Svg::IndexedSurfacesConverter::convert(
        tContext, surfaces, surfacesUpdator, drawOptions);
    auto pIndexCylinderView = Acts::Svg::View::zphi(pIndexeCylinder, fNameBase);

    Acts::Svg::toFile({pIndexCylinderView}, fNameBase + ".svg");
    // **********************

    // Kill that instance before going into the next test
    lcdd->destroyInstance();
  }
}

BOOST_AUTO_TEST_SUITE_END()
