// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"
#include "ActsPlugins/DD4hep/ConvertDD4hepDetector.hpp"

#include <algorithm>
#include <memory>
#include <stdexcept>

#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>
#include <DD4hep/Handle.h>
#include <DD4hep/Volumes.h>
#include <Parsers/Printout.h>
#include <TError.h>

using namespace Acts;
using namespace ActsPlugins;

namespace ActsExamples {

DD4hepDetectorBase::DD4hepDetectorBase(const Config& cfg)
    : Detector(Acts::getDefaultLogger("DD4hepDetector", cfg.logLevel)) {
  if (cfg.xmlFileNames.empty()) {
    throw std::invalid_argument("Missing DD4hep XML filenames");
  }

  m_nominalGeometryContext = GeometryContext::dangerouslyDefaultConstruct();

  m_detector = buildDD4hepGeometry(cfg);

  auto logger = getDefaultLogger("DD4hepConversion", cfg.logLevel);
}

dd4hep::Detector& DD4hepDetectorBase::dd4hepDetector() {
  return *m_detector;
}

const dd4hep::Detector& DD4hepDetectorBase::dd4hepDetector() const {
  return *m_detector;
}

std::shared_ptr<DD4hepFieldAdapter> DD4hepDetectorBase::field() const {
  throw_assert(m_detector != nullptr, "Detector not initialized");
  return std::make_shared<DD4hepFieldAdapter>(m_detector->field());
}

TGeoNode& DD4hepDetectorBase::tgeoGeometry() {
  return *m_detector->world().placement().ptr();
}

std::unique_ptr<dd4hep::Detector> DD4hepDetectorBase::buildDD4hepGeometry(
    const Config& cfg) const {
  const int old_gErrorIgnoreLevel = gErrorIgnoreLevel;

  switch (cfg.dd4hepLogLevel) {
    case Logging::Level::VERBOSE:
      dd4hep::setPrintLevel(dd4hep::PrintLevel::VERBOSE);
      break;
    case Logging::Level::DEBUG:
      dd4hep::setPrintLevel(dd4hep::PrintLevel::DEBUG);
      break;
    case Logging::Level::INFO:
      dd4hep::setPrintLevel(dd4hep::PrintLevel::INFO);
      break;
    case Logging::Level::WARNING:
      dd4hep::setPrintLevel(dd4hep::PrintLevel::WARNING);
      gErrorIgnoreLevel = kWarning;
      break;
    case Logging::Level::ERROR:
      dd4hep::setPrintLevel(dd4hep::PrintLevel::ERROR);
      gErrorIgnoreLevel = kError;
      break;
    case Logging::Level::FATAL:
      dd4hep::setPrintLevel(dd4hep::PrintLevel::FATAL);
      gErrorIgnoreLevel = kFatal;
      break;
    case Logging::Level::MAX:
      dd4hep::setPrintLevel(dd4hep::PrintLevel::ALWAYS);
      break;
  }
  // completely silence std::cout as DD4HEP is using it for logging
  if (cfg.dd4hepLogLevel >= Logging::Level::WARNING) {
    std::cout.setstate(std::ios_base::failbit);
  }

  std::unique_ptr<dd4hep::Detector> detector =
      dd4hep::Detector::make_unique(cfg.name);
  for (const auto& file : cfg.xmlFileNames) {
    detector->fromCompact(file);
  }
  detector->volumeManager();
  detector->apply("DD4hepVolumeManager", 0, nullptr);

  // restore the logging
  gErrorIgnoreLevel = old_gErrorIgnoreLevel;
  std::cout.clear();

  return detector;
}

DD4hepDetector::DD4hepDetector(const Config& cfg)
    : DD4hepDetectorBase{cfg}, m_cfg{cfg} {
  m_trackingGeometry = convertDD4hepDetector(
      m_detector->world(), logger(), m_cfg.bTypePhi, m_cfg.bTypeR, m_cfg.bTypeZ,
      m_cfg.envelopeR, m_cfg.envelopeZ, m_cfg.defaultLayerThickness,
      m_cfg.sortDetectors, m_nominalGeometryContext, m_cfg.materialDecorator,
      m_cfg.geometryIdentifierHook, m_cfg.detectorElementFactory);
}

auto DD4hepDetector::config() const -> const Config& {
  return m_cfg;
}

}  // namespace ActsExamples

void ActsExamples::sortFCChhDetElements(std::vector<dd4hep::DetElement>& det) {
  std::vector<dd4hep::DetElement> tracker;
  std::vector<dd4hep::DetElement> eCal;
  std::vector<dd4hep::DetElement> hCal;
  std::vector<dd4hep::DetElement> muon;
  for (const auto& detElement : det) {
    std::string detName = detElement.name();
    if (detName.find("Muon") != std::string::npos) {
      muon.push_back(detElement);
    } else if (detName.find("ECal") != std::string::npos) {
      eCal.push_back(detElement);
    } else if (detName.find("HCal") != std::string::npos) {
      hCal.push_back(detElement);
    } else {
      tracker.push_back(detElement);
    }
  }
  std::ranges::sort(
      muon, [](const dd4hep::DetElement& a, const dd4hep::DetElement& b) {
        return (a.id() < b.id());
      });
  std::ranges::sort(
      eCal, [](const dd4hep::DetElement& a, const dd4hep::DetElement& b) {
        return (a.id() < b.id());
      });
  std::ranges::sort(
      hCal, [](const dd4hep::DetElement& a, const dd4hep::DetElement& b) {
        return (a.id() < b.id());
      });
  std::ranges::sort(
      tracker, [](const dd4hep::DetElement& a, const dd4hep::DetElement& b) {
        return (a.id() < b.id());
      });
  det.clear();
  det = tracker;

  det.insert(det.end(), eCal.begin(), eCal.end());
  det.insert(det.end(), hCal.begin(), hCal.end());
  det.insert(det.end(), muon.begin(), muon.end());
}
