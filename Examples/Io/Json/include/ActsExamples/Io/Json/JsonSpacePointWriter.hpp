// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @date 2016-05-23 Initial version
/// @date 2017-08-07 Rewrite with new interfaces

#pragma once

#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <fstream>

namespace ActsExamples {

/// Write out a space point collection in JSON format.
///
/// This writes one file per event into the configured output directory. By
/// default it writes to the current working directory. Files are named
/// using the following schema
///
///     event000000001-spacepoints.json
///     event000000002-spacepoints.json
///
template <class T>
class JsonSpacePointWriter : public WriterT<GeometryIdMultimap<T>> {
 public:
  struct Config {
    std::string collection;      ///< which collection to write
    std::string outputDir;       ///< where to place output files
    size_t outputPrecision = 6;  ///< floating point precision
  };

  JsonSpacePointWriter(const Config& cfg,
                       Acts::Logging::Level level = Acts::Logging::INFO);

 protected:
  ActsExamples::ProcessCode writeT(
      const ActsExamples::AlgorithmContext& context,
      const GeometryIdMultimap<T>& spacePoints) final override;

 private:
  // since class itself is templated, base class template must be fixed
  using Base = WriterT<GeometryIdMultimap<T>>;

  Config m_cfg;
};

}  // namespace ActsExamples

template <class T>
ActsExamples::JsonSpacePointWriter<T>::JsonSpacePointWriter(
    const ActsExamples::JsonSpacePointWriter<T>::Config& cfg,
    Acts::Logging::Level level)
    : Base(cfg.collection, "JsonSpacePointWriter", level), m_cfg(cfg) {
  if (m_cfg.collection.empty()) {
    throw std::invalid_argument("Missing input collection");
  }
}

template <class T>
ActsExamples::ProcessCode ActsExamples::JsonSpacePointWriter<T>::writeT(
    const ActsExamples::AlgorithmContext& context,
    const GeometryIdMultimap<T>& spacePoints) {
  // open per-event file
  std::string path = perEventFilepath(m_cfg.outputDir, "spacepoints.json",
                                      context.eventNumber);
  std::ofstream os(path, std::ofstream::out | std::ofstream::trunc);
  if (!os) {
    throw std::ios_base::failure("Could not open '" + path + "' to write");
  }

  os << std::setprecision(m_cfg.outputPrecision);
  os << "{\n";

  bool firstVolume = true;
  for (auto& volumeData : spacePoints) {
    geo_id_value volumeID = volumeData.first;

    if (!firstVolume)
      os << ",\n";
    os << "  \"SpacePoints_" << volumeID << "\" : [\n";

    bool firstPoint = true;
    for (auto& layerData : volumeData.second) {
      for (auto& moduleData : layerData.second) {
        for (auto& data : moduleData.second) {
          // set the comma correctly
          if (!firstPoint)
            os << ",\n";
          // write the space point
          os << "    [" << data.x() << ", " << data.y() << ", " << data.z()
             << "]";
          firstPoint = false;
        }
      }
    }
    os << "]";
    firstVolume = false;
  }
  os << "\n}\n";

  return ProcessCode::SUCCESS;
}
