// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/JsonMaterialDecorator.hpp"

#include "Acts/Geometry/TrackingVolume.hpp"
#include "ActsPlugins/Json/TrackingGeometryMaterialJsonConverter.hpp"
#include "ActsPlugins/Json/detail/JsonIo.hpp"

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
  const auto jin = detail::readJsonFile(jFileName);
  if (jin.contains("header")) {
    m_materialMaps = TrackingGeometryMaterialJsonConverter{}.fromJson(jin);
  } else if (jin.contains("Surfaces") && jin.contains("Volumes")) {
    m_materialMaps = jmConverter.jsonToMaterialMaps(jin);
  } else {
    throw std::invalid_argument("Unrecognized material map format: " +
                                jFileName);
  }
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
