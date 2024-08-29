// This file is part of the Acts project.
//
// Copyright (C) 2022-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/ModuleMapCpp.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorElementITk.hpp"
#include "ActsExamples/Io/Root/RootAthenaDumpGeoIdCollecter.hpp"

namespace ActsExamples {

RootAthenaDumpGeoIdCollecter::RootAthenaDumpGeoIdCollecter(
    const Config& config, Acts::Logging::Level level)
    : IReader(),
      m_cfg(config),
      m_logger(Acts::getDefaultLogger(name(), level)) {
  if (m_cfg.inputfile.empty()) {
    throw std::invalid_argument("Missing input filename");
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

  m_inputchain->Add(m_cfg.inputfile.c_str());
  ACTS_DEBUG("Adding file " << m_cfg.inputfile << " to tree" << m_cfg.treename);

  m_events = m_inputchain->GetEntries();

  m_cfg.trackingGeometry->visitSurfaces([&](const Acts::Surface* surface) {
    auto detEl = dynamic_cast<const Acts::GeoModelDetectorElementITk*>(
        surface->associatedDetectorElement());
    if (detEl == nullptr) {
      throw std::runtime_error(
          "Could not convert to GeoModelDetectorElementITk!");
    }

    m_detectorElementMap.emplace(detEl->identifier().value(), detEl);
  });
}

ProcessCode RootAthenaDumpGeoIdCollecter::read(const AlgorithmContext& ctx) {
  ACTS_DEBUG("Reading event " << ctx.eventNumber);
  auto entry = ctx.eventNumber;
  if (entry >= m_events) {
    ACTS_ERROR("event out of bounds");
    return ProcessCode::ABORT;
  }

  std::lock_guard<std::mutex> lock(m_read_mutex);

  auto& athenaToActsGeoId = m_cfg.geometryIdMap->left;

  for (int im = 0; im < nCL; im++) {
    const auto athenaRepresentation = CLmoduleID[im];
    const auto hardware = CLhardware->at(im) == "PIXEL" ? 0 : 1;

    if (athenaToActsGeoId.find(athenaRepresentation) ==
        athenaToActsGeoId.end()) {
      // Here we use the same identifier that is used during geometry
      // construction to match the Acts::GeometryIdentifier and the athena
      // module id
      Acts::ITkIdentifier itkId(hardware, CLbarrel_endcap[im], CLlayer_disk[im],
                                CLeta_module[im], CLphi_module[im], CLside[im]);

      if (!m_detectorElementMap.contains(itkId.value())) {
        ACTS_ERROR("Missing sensitive surface, cannot map geometry ids");
        return ProcessCode::ABORT;
      }

      const auto gid =
          m_detectorElementMap.at(itkId.value())->surface().geometryId();

      ACTS_VERBOSE("Insert " << athenaRepresentation << " -> " << gid);
      athenaToActsGeoId.insert({athenaRepresentation, gid});
    }
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
