// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Utilities/Units.hpp>
#include <cstdlib>
#include <utility>

#include "ACTFW/DD4hepDetector/DD4hepGeometryService.hpp"
#include "ACTFW/Utilities/Options.hpp"

namespace po = boost::program_options;

namespace au = Acts::units;

namespace FW {

namespace Options {

void sortFCChhDetElements(std::vector<dd4hep::DetElement>& det) {
  std::vector<dd4hep::DetElement> tracker;
  std::vector<dd4hep::DetElement> eCal;
  std::vector<dd4hep::DetElement> hCal;
  std::vector<dd4hep::DetElement> muon;
  for (auto& detElement : det) {
    std::string detName = detElement.name();
    if (detName.find("Muon") != std::string::npos)
      muon.push_back(detElement);
    else if (detName.find("ECal") != std::string::npos)
      eCal.push_back(detElement);
    else if (detName.find("HCal") != std::string::npos)
      hCal.push_back(detElement);
    else
      tracker.push_back(detElement);
  }
  sort(muon.begin(), muon.end(),
       [](const dd4hep::DetElement& a, const dd4hep::DetElement& b) {
         return (a.id() < b.id());
       });
  sort(eCal.begin(), eCal.end(),
       [](const dd4hep::DetElement& a, const dd4hep::DetElement& b) {
         return (a.id() < b.id());
       });
  sort(hCal.begin(), hCal.end(),
       [](const dd4hep::DetElement& a, const dd4hep::DetElement& b) {
         return (a.id() < b.id());
       });
  sort(tracker.begin(), tracker.end(),
       [](const dd4hep::DetElement& a, const dd4hep::DetElement& b) {
         return (a.id() < b.id());
       });
  det.clear();
  det = tracker;

  det.insert(det.end(), eCal.begin(), eCal.end());
  det.insert(det.end(), hCal.begin(), hCal.end());
  det.insert(det.end(), muon.begin(), muon.end());
}

/// the particle gun options, the are prefixes with gp
template <typename aopt_t>
void addDD4hepOptions(aopt_t& opt) {
  opt.add_options()(
      "dd4hep-input",
      po::value<read_strings>()->multitoken()->default_value(
          {"file:Detectors/DD4hepDetector/compact/OpenDataDetector/"
           "OpenDataDetector.xml"}),
      "The locations of the input DD4hep files, use 'file:foo.xml'. In case "
      "you want to read in multiple files, just seperate the strings by "
      "space.")("dd4hep-envelopeR",
                po::value<double>()->default_value(1. * Acts::units::_mm),
                "The envelop cover in R for DD4hep volumes.")(
      "dd4hep-envelopeR",
      po::value<double>()->default_value(1. * Acts::units::_mm),
      "The tolerance added to the geometrical extension in r of the "
      "layers contained to build the volume envelope around in mm.")(
      "dd4hep-envelopeZ",
      po::value<double>()->default_value(1. * Acts::units::_mm),
      "The tolerance added to the geometrical extension in z of the "
      "layers contained to build the volume envelope around in mm.")(
      "dd4hep-layerThickness", po::value<double>()->default_value(10e-10),
      "In case no surfaces (to be contained by the layer) are handed over, "
      "the layer thickness will be set to this value.")(
      "dd4hep-buildFCChh", po::value<bool>()->default_value(true),
      "If you are not building the FCChh detector please set this flag to "
      "false.")("dd4hep-loglevel", po::value<size_t>()->default_value(2),
                "The output log level of the geometry building. Please set the "
                "wished "
                "number (0 = VERBOSE, 1 = "
                "DEBUG, 2 = INFO, 3 = WARNING, 4 = ERROR, 5 = FATAL).");
}

/// read the particle gun options and return a Config file
template <typename amap_t>
FW::DD4hep::DD4hepGeometryService::Config readDD4hepConfig(const amap_t& vm) {
  FW::DD4hep::DD4hepGeometryService::Config gsConfig;
  gsConfig.logLevel =
      Acts::Logging::Level(vm["dd4hep-loglevel"].template as<size_t>());
  gsConfig.xmlFileNames = vm["dd4hep-input"].template as<read_strings>();
  gsConfig.bTypePhi = Acts::equidistant;
  gsConfig.bTypeR = Acts::arbitrary;
  gsConfig.bTypeZ = Acts::equidistant;
  gsConfig.envelopeR = vm["dd4hep-envelopeR"].template as<double>();
  gsConfig.envelopeZ = vm["dd4hep-envelopeZ"].template as<double>();
  gsConfig.defaultLayerThickness =
      vm["dd4hep-layerThickness"].template as<double>();
  if (vm["dd4hep-buildFCChh"].template as<bool>()) {
    gsConfig.sortDetectors = sortFCChhDetElements;
  }
  return gsConfig;
}
}  // namespace Options
}  // namespace FW
