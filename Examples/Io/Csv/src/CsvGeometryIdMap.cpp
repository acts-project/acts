// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvGeometryIdMap.hpp"

#include "ActsExamples/Io/Csv/CsvInputOutput.hpp"

#include <stdexcept>

#include "CsvOutputData.hpp"

namespace ActsExamples {

std::shared_ptr<GeometryIdMapActsAthena> readCsvGeometryIdMap(
    const std::string& path) {
  auto map = std::make_shared<GeometryIdMapActsAthena>();

  BoostDescribeCsvReader<GeometryIdMapData> reader(path);
  GeometryIdMapData data;
  while (reader.read(data)) {
    const bool inserted =
        map->insert({data.athena_id, Acts::GeometryIdentifier(data.acts_id)})
            .second;
    if (!inserted) {
      throw std::runtime_error("Duplicate geometry identifier in " + path);
    }
  }

  return map;
}

}  // namespace ActsExamples
