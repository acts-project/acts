// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Json/JsonMaterialWriter.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"

#include <fstream>
#include <ios>
#include <iostream>
#include <stdexcept>

ActsExamples::JsonMaterialWriter::JsonMaterialWriter(
    const ActsExamples::JsonMaterialWriter::Config& cfg)
    : m_cfg(cfg) {
  // Validate the configuration
  if (m_cfg.converterCfg.name.empty()) {
    throw std::invalid_argument("Missing service name");
  }
}

ActsExamples::JsonMaterialWriter::~JsonMaterialWriter() {}

void ActsExamples::JsonMaterialWriter::write(
    const Acts::DetectorMaterialMaps& detMaterial) {
  // Evoke the converter
  Acts::MaterialMapJsonConverter jmConverter(m_cfg.converterCfg);
  auto jOut = jmConverter.materialMapsToJson(detMaterial);
  // And write the file(s)
  if (ACTS_CHECK_BIT(m_cfg.writeFormat, ActsExamples::JsonFormat::Json)) {
    std::ofstream ofj(m_cfg.fileName + ".json");
    ofj << std::setw(4) << jOut << std::endl;
  }
  if (ACTS_CHECK_BIT(m_cfg.writeFormat, ActsExamples::JsonFormat::Cbor)) {
    std::vector<uint8_t> cborOut = nlohmann::json::to_cbor(jOut);
    std::ofstream ofj(m_cfg.fileName + ".cbor",
                      std::ios::out | std::ios::binary);
    ofj.write((char*)cborOut.data(), cborOut.size());
  }
}

void ActsExamples::JsonMaterialWriter::write(
    const Acts::TrackingGeometry& tGeometry) {
  // Evoke the converter
  Acts::MaterialMapJsonConverter jmConverter(m_cfg.converterCfg);
  auto jOut = jmConverter.trackingGeometryToJson(tGeometry);
  // And write the file(s)
  if (ACTS_CHECK_BIT(m_cfg.writeFormat, ActsExamples::JsonFormat::Json)) {
    std::ofstream ofj(m_cfg.fileName + ".json");
    ofj << std::setw(4) << jOut << std::endl;
  }
  if (ACTS_CHECK_BIT(m_cfg.writeFormat, ActsExamples::JsonFormat::Cbor)) {
    std::vector<uint8_t> cborOut = nlohmann::json::to_cbor(jOut);
    std::ofstream ofj(m_cfg.fileName + ".cbor",
                      std::ios::out | std::ios::binary);
    ofj.write((char*)cborOut.data(), cborOut.size());
  }
}
