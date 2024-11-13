// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootAthenaDumpGeoIdCollector.hpp"

#include "ActsExamples/ITkHelpers/ITkDetectorElement.hpp"

#include <fstream>
#include <ranges>

namespace ActsExamples {

RootAthenaDumpGeoIdCollector::RootAthenaDumpGeoIdCollector(
    const Config& config, Acts::Logging::Level level)
    : IReader(),
      m_cfg(config),
      m_logger(Acts::getDefaultLogger(name(), level)) {
  if (m_cfg.inputfile.empty()) {
    throw std::invalid_argument("Empty input file list");
  }
  if (m_cfg.treename.empty()) {
    throw std::invalid_argument("Missing tree name");
  }
  if (!m_cfg.trackingGeometry) {
    throw std::invalid_argument("Need tracking geometry!");
  }
  if (!m_cfg.geometryIdMap) {
    throw std::invalid_argument("Need map to write into!");
  }

  m_inputchain = std::make_shared<TChain>(m_cfg.treename.c_str());
  m_inputchain->SetBranchStatus("*", false);

  auto setBranchAddress = [&](const char* name, auto ptr) {
    m_inputchain->SetBranchStatus(name, true);
    m_inputchain->SetBranchAddress(name, ptr);
  };

  setBranchAddress("run_number", &run_number);
  setBranchAddress("event_number", &event_number);
  setBranchAddress("nCL", &nCL);
  setBranchAddress("CLhardware", &CLhardware);
  setBranchAddress("CLbarrel_endcap", CLbarrel_endcap);
  setBranchAddress("CLlayer_disk", CLlayer_disk);
  setBranchAddress("CLeta_module", CLeta_module);
  setBranchAddress("CLphi_module", CLphi_module);
  setBranchAddress("CLside", CLside);
  setBranchAddress("CLmoduleID", CLmoduleID);
  for (const auto& f : m_cfg.inputfile) {
    m_inputchain->Add(f.c_str());
    ACTS_DEBUG("Adding file '" << f << "' to tree " << m_cfg.treename);
  }
  m_events = m_inputchain->GetEntries();
  ACTS_DEBUG("Found " << m_events << " to read");

  m_inputchain->GetEntry(0);

  m_cfg.trackingGeometry->visitSurfaces([&](const Acts::Surface* surface) {
    auto detEl = dynamic_cast<const ActsExamples::ITkDetectorElement*>(
        surface->associatedDetectorElement());
    if (detEl == nullptr) {
      throw std::runtime_error("Could not convert to ITkDetectorElement!");
    }

    m_detectorElementMap.emplace(detEl->identifier().value(), detEl);
  });

  ACTS_INFO("Collected " << m_detectorElementMap.size()
                         << " detector elements for mapping");
  // throw std::runtime_error("stop here");
}

ProcessCode RootAthenaDumpGeoIdCollector::read(const AlgorithmContext& ctx) {
  ACTS_DEBUG("Reading event " << ctx.eventNumber);
  auto entry = ctx.eventNumber;
  if (entry >= m_events) {
    ACTS_ERROR("event out of bounds");
    return ProcessCode::ABORT;
  }

  std::lock_guard<std::mutex> lock(m_read_mutex);
  m_inputchain->GetEntry(entry);

  auto& athenaToActsGeoId = m_cfg.geometryIdMap->left;

  const auto prev = athenaToActsGeoId.size();
  ACTS_DEBUG("Read data from " << nCL << " measurements");

  std::vector<const ActsExamples::ITkDetectorElement*> matched;
  std::vector<ActsExamples::ITkIdentifier> missed;

  for (int im = 0; im < nCL; im++) {
    const auto athenaRepresentation = CLmoduleID[im];
    const auto hardware = CLhardware->at(im) == "PIXEL" ? 0 : 1;

    if (athenaToActsGeoId.find(athenaRepresentation) ==
        athenaToActsGeoId.end()) {
      // Here we use the same identifier that is used during geometry
      // construction to match the Acts::GeometryIdentifier and the athena
      // module id
      ActsExamples::ITkIdentifier itkId(hardware, CLbarrel_endcap[im],
                                        CLlayer_disk[im], CLeta_module[im],
                                        CLphi_module[im], CLside[im]);

      if (!m_detectorElementMap.contains(itkId.value())) {
        ACTS_WARNING("Missing sensitive surface for "
                     << itkId << ", cannot map geometry ids");
        missed.push_back(itkId);
        continue;
      }

      matched.push_back(m_detectorElementMap.at(itkId.value()));
      const auto gid =
          m_detectorElementMap.at(itkId.value())->surface().geometryId();

      ACTS_VERBOSE("Insert " << athenaRepresentation << " -> " << gid);
      athenaToActsGeoId.insert({athenaRepresentation, gid});
    }
  }

  {
    std::ranges::sort(matched);
    auto ret = std::ranges::unique(matched);
    matched.erase(ret.begin(), ret.end());
  }

  {
    std::ranges::sort(missed, {}, &ActsExamples::ITkIdentifier::value);
    auto ret =
        std::ranges::unique(missed, {}, &ActsExamples::ITkIdentifier::value);
    missed.erase(ret.begin(), ret.end());
  }

  std::ofstream missedFile("missed_surfaces.csv");
  missedFile << "hardware,bec,lw,em,pm,side\n";
  for (const auto& m : missed) {
    missedFile << m.hardware() << "," << m.barrelEndcap() << ","
               << m.layerWheel() << "," << m.etaModule() << "," << m.phiModule()
               << "," << m.side() << "\n";
  }

  std::ofstream matchedFile("matched_surfaces.csv");
  matchedFile << "acts_geoid,hardware,bec,lw,em,pm,side\n";
  for (auto d : matched) {
    const auto& m = d->identifier();
    matchedFile << d->surface().geometryId().value() << ",";
    matchedFile << m.hardware() << "," << m.barrelEndcap() << ","
                << m.layerWheel() << "," << m.etaModule() << ","
                << m.phiModule() << "," << m.side() << "\n";
  }

  ACTS_DEBUG("Added " << athenaToActsGeoId.size() - prev << " entries");

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
