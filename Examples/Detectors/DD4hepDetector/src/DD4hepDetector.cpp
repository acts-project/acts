// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"

#include "Acts/Plugins/DD4hep/ConvertDD4hepDetector.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/DetectorCommons/DetectorBase.hpp"

#include <algorithm>
#include <memory>
#include <stdexcept>

#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>
#include <DD4hep/Handle.h>
#include <DD4hep/Volumes.h>
#include <Parsers/Printout.h>
#include <TError.h>

namespace ActsExamples {

DD4hepDetector::DD4hepDetector(
    Acts::GeometryContext geometryContext,
    std::vector<std::shared_ptr<const Acts::DetectorElementBase>> detectorStore,
    std::shared_ptr<const Acts::TrackingGeometry> gen1Geometry,
    std::shared_ptr<Acts::Experimental::Detector> gen2Geometry,
    std::vector<std::shared_ptr<IContextDecorator>> contextDecorators,
    std::shared_ptr<dd4hep::Detector> detector)
    : PreConstructedDetector(std::move(geometryContext),
                             std::move(detectorStore), std::move(gen1Geometry),
                             std::move(gen2Geometry),
                             std::move(contextDecorators)),
      m_detector(std::move(detector)) {}

dd4hep::Detector& DD4hepDetector::dd4hepDetector() {
  return *m_detector;
}

dd4hep::DetElement DD4hepDetector::dd4hepGeometry() {
  return m_detector->world();
}

TGeoNode& DD4hepDetector::tgeoGeometry() {
  return *dd4hepGeometry().placement().ptr();
}

DD4hepDetectorFactory::DD4hepDetectorFactory(const Config& cfg)
    : DetectorFactoryBase(
          Acts::getDefaultLogger("DD4hepDetector", cfg.logLevel)),
      m_cfg(cfg) {
  if (m_cfg.xmlFileNames.empty()) {
    throw std::invalid_argument("Missing DD4hep XML filenames");
  }
}

DD4hepDetectorFactory::DD4hepDetectorFactory(DD4hepDetectorFactory&&) = default;

DD4hepDetectorFactory::~DD4hepDetectorFactory() = default;

DD4hepDetectorFactory& DD4hepDetectorFactory::operator=(
    DD4hepDetectorFactory&&) = default;

std::unique_ptr<dd4hep::Detector> DD4hepDetectorFactory::buildDD4hepGeometry()
    const {
  const int old_gErrorIgnoreLevel = gErrorIgnoreLevel;
  switch (m_cfg.dd4hepLogLevel) {
    case Acts::Logging::Level::VERBOSE:
      dd4hep::setPrintLevel(dd4hep::PrintLevel::VERBOSE);
      break;
    case Acts::Logging::Level::DEBUG:
      dd4hep::setPrintLevel(dd4hep::PrintLevel::DEBUG);
      break;
    case Acts::Logging::Level::INFO:
      dd4hep::setPrintLevel(dd4hep::PrintLevel::INFO);
      break;
    case Acts::Logging::Level::WARNING:
      dd4hep::setPrintLevel(dd4hep::PrintLevel::WARNING);
      gErrorIgnoreLevel = kWarning;
      break;
    case Acts::Logging::Level::ERROR:
      dd4hep::setPrintLevel(dd4hep::PrintLevel::ERROR);
      gErrorIgnoreLevel = kError;
      break;
    case Acts::Logging::Level::FATAL:
      dd4hep::setPrintLevel(dd4hep::PrintLevel::FATAL);
      gErrorIgnoreLevel = kFatal;
      break;
    case Acts::Logging::Level::MAX:
      dd4hep::setPrintLevel(dd4hep::PrintLevel::ALWAYS);
      break;
  }
  // completely silence std::cout as DD4HEP is using it for logging
  if (m_cfg.dd4hepLogLevel >= Acts::Logging::Level::WARNING) {
    std::cout.setstate(std::ios_base::failbit);
  }

  std::unique_ptr<dd4hep::Detector> detector =
      dd4hep::Detector::make_unique(m_cfg.name);
  for (auto& file : m_cfg.xmlFileNames) {
    detector->fromCompact(file.c_str());
  }
  detector->volumeManager();
  detector->apply("DD4hepVolumeManager", 0, nullptr);

  // restore the logging
  gErrorIgnoreLevel = old_gErrorIgnoreLevel;
  std::cout.clear();

  return detector;
}

std::shared_ptr<DetectorBase> DD4hepDetectorFactory::buildDetector() const {
  Acts::GeometryContext geometryContext;
  std::vector<std::shared_ptr<const Acts::DetectorElementBase>> detectorStore;
  std::shared_ptr<const Acts::TrackingGeometry> gen1Geometry;
  std::shared_ptr<Acts::Experimental::Detector> gen2Geometry;
  std::vector<std::shared_ptr<IContextDecorator>> contextDecorators;
  std::unique_ptr<dd4hep::Detector> detector;

  geometryContext = Acts::GeometryContext();

  detector = buildDD4hepGeometry();

  auto logger = Acts::getDefaultLogger("DD4hepConversion", m_cfg.logLevel);
  gen1Geometry = Acts::convertDD4hepDetector(
      detector->world(), *logger, m_cfg.bTypePhi, m_cfg.bTypeR, m_cfg.bTypeZ,
      m_cfg.envelopeR, m_cfg.envelopeZ, m_cfg.defaultLayerThickness,
      m_cfg.sortDetectors, geometryContext, m_cfg.materialDecorator,
      m_cfg.geometryIdentifierHook);

  return std::make_shared<DD4hepDetector>(
      geometryContext, detectorStore, gen1Geometry, gen2Geometry,
      contextDecorators, std::move(detector));
}

}  // namespace ActsExamples

void ActsExamples::sortFCChhDetElements(std::vector<dd4hep::DetElement>& det) {
  std::vector<dd4hep::DetElement> tracker;
  std::vector<dd4hep::DetElement> eCal;
  std::vector<dd4hep::DetElement> hCal;
  std::vector<dd4hep::DetElement> muon;
  for (auto& detElement : det) {
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
