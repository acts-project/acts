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

const char* container_head_xml =
    R""""(
    <detectors>
)"""";

const char* beampipe_xml =
    R""""(
        <detector id="0" name="BeamPipe" type="BarrelDetector">
            <type_flags type="DetType_TRACKER + DetType_BEAMPIPE"/>
            <envelope rmin="0*mm" rmax="3*cm" dz="200*cm" cz="0*mm" material="Air"/>
            <proto_container name="ContainerInR" binning="r" order="0"/>
            <layers>
                <layer id="0" name="inner" rmin="0*cm" rmax="30*mm" dz="200*cm" material="Air" vis="None">
                    <proto_volume>            
                        <tubs rmin="0*mm" rmax="30.0*mm" dz="2000*mm"/>
                    </proto_volume>
                    <passive_surface>
                        <tubs rmin="25*mm" rmax="25.8*mm" dz="180*cm" material="Air"/>
                    </passive_surface>
                </layer>
            </layers>            
        </detector>
    )"""";

const char* inner_cylinder_xml =
    R""""(
        <detector id="1" name="CylinderLayer" type="CylinderLayer">
            <type_flags type="DetType_TRACKER + DetType_BEAMPIPE"/>
            <envelope rmin="0*mm" rmax="3*cm" dz="200*cm" cz="0*mm" material="Air"/>
            <proto_container name="ContainerInR" binning="r" order="0"/>
            <proto_volume>            
              <tubs rmin="0*mm" rmax="30.0*mm" dz="2000*mm"/>
            </proto_volume>
            <passive_surface>
               <tubs rmin="25*mm" rmax="25.8*mm" dz="180*cm" material="Air"/>
            </passive_surface>
            <layers>
                <layer name="inner"/>
                <layer name="outer"/>
            </layers>
        </detector>
    )"""";

const char* tail_xml =
    R""""(
    </detectors>
)"""";

BOOST_AUTO_TEST_SUITE(DD4hepPlugin);

// This tests creates a few volumes and puts them into a container
BOOST_AUTO_TEST_CASE(ContainerInR) {
  // Create an XML from it
  ofstream cxml;
  cxml.open("ContainerInR.xml");

  cxml << head_xml;
  cxml << container_head_xml;
  cxml << beampipe_xml;
  cxml << tail_xml;
  cxml << end_xml;

  cxml.close();

  // Create the DD4hep geometry
  auto lcdd = &(dd4hep::Detector::getInstance());
  lcdd->fromCompact("ContainerInR.xml");
  lcdd->volumeManager();
  lcdd->apply("DD4hepVolumeManager", 0, nullptr);
}

BOOST_AUTO_TEST_SUITE_END();
