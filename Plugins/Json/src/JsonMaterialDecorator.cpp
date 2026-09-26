// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/JsonMaterialDecorator.hpp"

#include "Acts/Geometry/TrackingVolume.hpp"

#include <fstream>
#include <stdexcept>

namespace Acts {

JsonMaterialDecorator::JsonMaterialDecorator(
    const MaterialMapJsonConverter::Config& rConfig,
    const std::string& jFileName, Acts::Logging::Level level)
    : m_readerConfig(rConfig),
      m_logger{getDefaultLogger("JsonMaterialDecorator", level)} {
  // the material reader
  Acts::MaterialMapJsonConverter jmConverter(rConfig, level);

  ACTS_VERBOSE("Reading JSON material description from: " << jFileName);
  std::ifstream ifj(jFileName.c_str());
  if (!ifj.good()) {
    throw std::runtime_error{"Unable to open input JSON material file: " +
                             jFileName};
  }
  nlohmann::json jin;

  if (jFileName.find(".cbor") != std::string::npos) {
    std::vector<std::uint8_t> iCbor((std::istreambuf_iterator<char>(ifj)),
                                    std::istreambuf_iterator<char>());
    jin = nlohmann::json::from_cbor(iCbor);
  } else {
    ifj >> jin;
  }

  m_materialMaps = jmConverter.jsonToMaterialMaps(jin);
  ACTS_VERBOSE("JSON material description read complete");
}

void JsonMaterialDecorator::decorate(Surface& surface) const {
  m_materialMaps.apply(surface);
}

void JsonMaterialDecorator::decorate(TrackingVolume& volume) const {
  if (!m_materialMaps.keyedSurfaces.empty() && volume.portals().empty()) {
    throw std::invalid_argument(
        "Cannot apply a keyed material map to Gen1 geometry: stable material "
        "keys require Gen3 material designators. Use a Gen1 map indexed by "
        "geometry ID or apply this map to the matching Gen3 geometry.");
  }
  m_materialMaps.apply(volume);
}

}  // namespace Acts
