// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/Conversion/CellMapConversion.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/Cluster.hpp"

#include <cstdint>
#include <cstdlib>
#include <map>
#include <vector>

#include "traccc/edm/cell.hpp"

namespace ActsExamples::Traccc::Common::Conversion {

std::map<Acts::GeometryIdentifier, std::vector<traccc::cell>> tracccCellsMap(
    const std::map<Acts::GeometryIdentifier, std::vector<Cluster::Cell>>& map) {
  std::map<Acts::GeometryIdentifier, std::vector<traccc::cell>> tracccCellMap;
  for (const auto& [geometryID, cells] : map) {
    std::vector<traccc::cell> tracccCells;
    for (const auto& cell : cells) {
      traccc::cell::link_type moduleLink = 0;

      tracccCells.emplace_back(static_cast<unsigned int>(cell.bin[0]),
                               static_cast<unsigned int>(cell.bin[1]),
                               static_cast<float>(cell.activation), 0.f,
                               moduleLink);
    }
    std::ranges::sort(tracccCells,
                      [](const traccc::cell& lhs, const traccc::cell& rhs) {
                        if (lhs.module_link != rhs.module_link) {
                          return lhs.module_link < rhs.module_link;
                        } else if (lhs.channel1 != rhs.channel1) {
                          return (lhs.channel1 < rhs.channel1);
                        } else {
                          return (lhs.channel0 < rhs.channel0);
                        }
                      });
    tracccCellMap.insert({geometryID, std::move(tracccCells)});
  }
  return tracccCellMap;
}

}  // namespace ActsExamples::Traccc::Common::Conversion
