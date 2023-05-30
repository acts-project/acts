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

const char* disc_layer_head_xml =
    R""""(
    <detectors>
        <detector id="1" name="DiscLayer" type="DiscLayer" readout="PixelReadout">
            <type_flags type="DetType_TRACKER + DetType_EDNCAP"/>
)"""";

const char* disc_layer_tail_xml =
    R""""(
        </detector>
    </detectors>
)"""";

const char* indent_8_xml = "        ";
const char* indent_12_xml = "            ";

BOOST_AUTO_TEST_SUITE(DD4hepPlugin)

// This test creates a DD4hep::compact.xml file to describe various disc
// setups with different binnings.
//
// It then reads this xml, creates a DD4hep geometry and converts the
// layer structure.
//
// This example is without bin expansion, bin expansion is showcased in
// the cylinder layer example.
//
// All combiniations of equidistant/variable are tested in this test.
BOOST_AUTO_TEST_CASE(DD4hepPluginDiscLayerStructure) {
  // binR type, outerR, hlafY, binPhi, bphi-bins
  using DiscSpecs = std::tuple<std::string, Acts::ActsScalar, Acts::ActsScalar,
                               std::string, std::string>;
  std::vector<DiscSpecs> discSpecifications = {
      {"equidistant", 80., 20., "equidistant", "noop"},
      {"variable", 100., 40., "equidistant", "noop"},
      {"equidistant", 80., 20., "variable",
       "-3.1415,-2.0,-0.5,0.7,1.45,3.1415"},
      {"variable", 100., 40., "variable", "-3.1415,-2.0,-0.5,0.7,1.45,3.1415"}};

  for (auto [nR, oR, oHy, nPhi, bPhi] : discSpecifications) {
    // Detector store
    Acts::Test::CylindricalTrackingGeometry::DetectorStore dStore;
    auto rSurfacesR0 = cGeometry.surfacesRing(dStore, 6.4, 11.4, 18., 0.125, 0.,
                                              42., 95., 2., 22u);

    auto rSurfacesR1 = cGeometry.surfacesRing(dStore, 12.4, 18.4, oHy, 0.125,
                                              0., oR, 100., 2., 22u);

    std::vector<std::vector<const Acts::Surface*>> rSurfaces = {rSurfacesR0,
                                                                rSurfacesR1};

    std::string fNameBase = "DiscLayer_" + nR + "_" + nPhi;

    std::ofstream dxml;
    dxml.open(fNameBase + ".xml");
    dxml << head_xml;
    dxml << segmentation_xml;
    dxml << disc_layer_head_xml;

    dxml
        << indent_12_xml
        << "<envelope rmin=\"20*mm\" rmax=\"200*mm\" dz=\"5*cm\" cz=\"100*mm\" "
           "material=\"Air\"/>"
        << '\n';
    dxml << indent_12_xml << "<surface_binning";
    if (nR == "equidistant" and nPhi == "equidistant") {
      dxml << " nr=\"2\" rmin=\"25\" rmax=\"100\"";
      dxml << " nphi=\"22\"  phimin=\"-3.1415\" phimax=\"3.1415\"/>";
    } else if (nR == "variable" and nPhi == "equidistant") {
      dxml << " rboundaries=\"2.5*cm,60*mm,140*mm\"";
      dxml << " nphi=\"22\"  phimin=\"-3.1415\" phimax=\"3.1415\"/>";
    } else if (nR == "equidistant") {
      dxml << " nr=\"2\" rmin=\"25\" rmax=\"100\"";
      dxml << " phiboundaries=\"" + bPhi + "\"/>";
    } else {
      dxml << " rboundaries=\"25*mm,60*mm,140*mm\"";
      dxml << " phiboundaries=\"" + bPhi + "\"/>";
    }
    dxml << '\n';

    Acts::Transform3 layerTransform = Acts::Transform3::Identity();
    layerTransform.pretranslate(Acts::Vector3{0., 0., 100.});
    Acts::Transform3 invLayerT = layerTransform.inverse();

    dxml << indent_8_xml << "<modules>" << '\n';

    for (const auto& ss : rSurfaces) {
      for (const auto& s : ss) {
        dxml << indent_12_xml
             << DD4hepTestsHelper::surfaceToXML(tContext, *s, invLayerT)
             << "\n";
      }
    }

    dxml << indent_8_xml << "</modules>" << '\n';
    dxml << disc_layer_tail_xml;
    dxml << end_xml;
    dxml.close();

    auto lcdd = &(dd4hep::Detector::getInstance());
    lcdd->fromCompact(fNameBase + ".xml");
    lcdd->volumeManager();
    lcdd->apply("DD4hepVolumeManager", 0, nullptr);

    auto world = lcdd->world();

    // Now the test starts ...
    auto sFactory = std::make_shared<Acts::DD4hepDetectorSurfaceFactory>(
        Acts::getDefaultLogger("DD4hepDetectorSurfaceFactory",
                               Acts::Logging::DEBUG));

    Acts::Experimental::DD4hepLayerStructure endcapStructure(
        std::move(sFactory),
        Acts::getDefaultLogger("DD4hepLayerStructure", Acts::Logging::DEBUG));

    Acts::DD4hepDetectorElement::Store dd4hepStore;

    Acts::Experimental::DD4hepLayerStructure::Options lsOptions;
    lsOptions.name = "EndcapLayer";
    lsOptions.logLevel = Acts::Logging::INFO;

    auto endcapInternalsBuilder =
        endcapStructure.builder(dd4hepStore, world, lsOptions);

    // Build the internal volume structure
    auto [surfaces, volumes, surfacesUpdator, volumeUpdator] =
        endcapInternalsBuilder->construct(tContext);

    // All surfaces are filled
    BOOST_CHECK(surfaces.size() == 44u);
    // No volumes are added
    BOOST_CHECK(volumes.empty());
    // The surface updator is connected
    BOOST_CHECK(surfacesUpdator.connected());
    // The volume updator is connected
    BOOST_CHECK(volumeUpdator.connected());

    // Kill that instance before going into the next test
    lcdd->destroyInstance();
  }
}

BOOST_AUTO_TEST_SUITE_END()
