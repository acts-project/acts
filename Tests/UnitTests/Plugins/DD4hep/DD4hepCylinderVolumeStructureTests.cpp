// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/DD4hep/DD4hepVolumeStructure.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <fstream>
#include <string>

#include <DD4hep/DetElement.h>
#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/Detector.h>
#include <XML/Utilities.h>
#include <XMLFragments.hpp>

#include "DD4hepTestsHelper.hpp"

Acts::GeometryContext tContext;

const char* beampipe_head_xml =
    R""""(
    <detectors>
        <detector id="0" name="BeamPipe" type="BarrelDetector">
            <type_flags type="DetType_TRACKER + DetType_BEAMPIPE"/>
)"""";

const char* tail_xml =
    R""""(
        </detector>
    </detectors>
)"""";

const char* indent_12_xml = "            ";

BOOST_AUTO_TEST_SUITE(DD4hepPlugin)

// This tests creates a cylindrical structure
BOOST_AUTO_TEST_CASE(DD4hepPluginCylindricalVolumeStructure) {
  // Create an XML from it
  std::ofstream cxml;

  std::string fNameBase = "CylindricalVolume";

  cxml.open(fNameBase + ".xml");
  cxml << head_xml;
  cxml << beampipe_head_xml;
  // DD4hep description
  cxml << indent_12_xml
       << "<envelope rmin=\"0*mm\" rmax=\"30*mm\" dz=\"2000*mm\" "
          "cz=\"0*mm\" "
          "material=\"Air\"/>"
       << '\n';
  // External description
  cxml << indent_12_xml << "<layers>" << '\n';
  cxml << indent_12_xml << "<layer name=\"CylinderLayer\" id=\"0\">" << '\n';
  cxml << indent_12_xml
       << "<acts_volume name=\"CylinderLayerVolume\" sequence=\"1\">" << '\n';
  cxml << indent_12_xml
       << "    <tubs rmin=\"0*mm\" rmax=\"30.0*mm\" dz=\"2000*mm\" />" << '\n';
  cxml << indent_12_xml << "</acts_volume>" << '\n';

  // Internal description
  cxml << indent_12_xml << "<passive_surface>" << '\n';
  cxml << indent_12_xml
       << "    <tubs rmin=\"25*mm\" rmax=\"25.8*mm\" dz=\"1800*mm\" "
          "material=\"Air\"/>"
       << '\n';
  cxml << indent_12_xml << "</passive_surface>" << '\n';
  cxml << indent_12_xml << "</layer>" << '\n';
  cxml << indent_12_xml << "</layers>" << '\n';

  cxml << tail_xml;
  cxml << end_xml;
  cxml.close();

  auto lcdd = &(dd4hep::Detector::getInstance());
  lcdd->fromCompact(fNameBase + ".xml");
  lcdd->volumeManager();
  lcdd->apply("DD4hepVolumeManager", 0, nullptr);

  auto world = lcdd->world();

  // Create the structure & the builder
  Acts::Experimental::DD4hepVolumeStructure volumeStructure(
      Acts::getDefaultLogger("DD4hepVolumeStructure", Acts::Logging::VERBOSE));

  Acts::Experimental::DD4hepVolumeStructure::Options vOptions;
  vOptions.logLevel = Acts::Logging::VERBOSE;
  vOptions.name = "CylinderLayerVolume";

  auto volumeBuilder = volumeStructure.builder(world, vOptions);

  auto [transform, bounds, portalGenerator] =
      volumeBuilder->construct(tContext);

  BOOST_CHECK(transform.isApprox(Acts::Transform3::Identity()));
  BOOST_CHECK(bounds != nullptr);
  BOOST_CHECK(bounds->type() == Acts::VolumeBounds::BoundsType::eCylinder);
  BOOST_CHECK(bounds->values().size() == 7u);
  BOOST_CHECK(bounds->values()[0u] == 0.);
  BOOST_CHECK(bounds->values()[1u] == 30.);
  BOOST_CHECK(bounds->values()[2u] == 1000.);

  // Kill that instance before going into the next test
  lcdd->destroyInstance();
}

// This tests creates a cylindrical structure
BOOST_AUTO_TEST_CASE(DD4hepPluginShiftedCylindricalVolumeStructure) {
  // Create an XML from it
  std::ofstream cxml;

  std::string fNameBase = "ShiftedCylindricalVolume";

  cxml.open(fNameBase + ".xml");
  cxml << head_xml;
  cxml << beampipe_head_xml;
  // DD4hep description
  cxml << indent_12_xml
       << "<envelope rmin=\"0*mm\" rmax=\"30*mm\" dz=\"2000*mm\" "
          "cz=\"0*mm\" "
          "material=\"Air\"/>"
       << '\n';
  // External description
  cxml << indent_12_xml << "<layers>" << '\n';
  cxml << indent_12_xml << "<layer name=\"ShiftedCylinderLayer\" id=\"0\">"
       << '\n';
  cxml << indent_12_xml
       << "<acts_volume name=\"ShiftedCylinderLayerVolume\" sequence=\"1\">"
       << '\n';
  cxml << indent_12_xml
       << "    <tubs rmin=\"0*mm\" rmax=\"30.0*mm\" dz=\"2000*mm\" "
          "x=\"10.*mm\" y=\"1.5*cm\" z=\"-1300*mm\"/>"
       << '\n';
  cxml << indent_12_xml << "</acts_volume>" << '\n';

  // Internal description
  cxml << indent_12_xml << "<passive_surface>" << '\n';
  cxml << indent_12_xml
       << "    <tubs rmin=\"25*mm\" rmax=\"25.8*mm\" dz=\"1800*mm\" "
          "material=\"Air\"/>"
       << '\n';
  cxml << indent_12_xml << "</passive_surface>" << '\n';
  cxml << indent_12_xml << "</layer>" << '\n';
  cxml << indent_12_xml << "</layers>" << '\n';

  cxml << tail_xml;
  cxml << end_xml;
  cxml.close();

  auto lcdd = &(dd4hep::Detector::getInstance());
  lcdd->fromCompact(fNameBase + ".xml");
  lcdd->volumeManager();
  lcdd->apply("DD4hepVolumeManager", 0, nullptr);

  auto world = lcdd->world();

  // Create the structure & the builder
  Acts::Experimental::DD4hepVolumeStructure volumeStructure(
      Acts::getDefaultLogger("DD4hepVolumeStructure", Acts::Logging::VERBOSE));

  Acts::Experimental::DD4hepVolumeStructure::Options vOptions;
  vOptions.logLevel = Acts::Logging::VERBOSE;
  vOptions.name = "ShiftedCylinderLayerVolume";

  auto volumeBuilder = volumeStructure.builder(world, vOptions);

  auto [transform, bounds, portalGenerator] =
      volumeBuilder->construct(tContext);

  BOOST_CHECK(bounds != nullptr);
  BOOST_CHECK(bounds->type() == Acts::VolumeBounds::BoundsType::eCylinder);
  BOOST_CHECK(bounds->values().size() == 7u);
  BOOST_CHECK(bounds->values()[0u] == 0.);
  BOOST_CHECK(bounds->values()[1u] == 30.);
  BOOST_CHECK(bounds->values()[2u] == 1000.);

  // It is shifted - so let us check this
  auto shift = transform.translation();
  CHECK_CLOSE_ABS(shift.x(), 10., 1e-10);
  CHECK_CLOSE_ABS(shift.y(), 15., 1e-10);
  CHECK_CLOSE_ABS(shift.z(), -1300., 1e-10);

  // Kill that instance before going into the next test
  lcdd->destroyInstance();
}

BOOST_AUTO_TEST_SUITE_END()
