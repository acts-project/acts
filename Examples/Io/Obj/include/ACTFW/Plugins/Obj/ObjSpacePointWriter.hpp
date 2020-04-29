// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <fstream>

#include "ACTFW/EventData/GeometryContainers.hpp"
#include "ACTFW/Framework/WriterT.hpp"
#include "ACTFW/Utilities/Paths.hpp"

namespace FW {
namespace Obj {

/// Write out a space point collection in OBJ format.
///
/// This writes one file per event into the configured output directory. By
/// default it writes to the current working directory. Files are named
/// using the following schema
///
///     event000000001-spacepoints.obj
///     event000000002-spacepoints.obj
///
/// One write call per thread and hence thread safe.
template <typename T>
class ObjSpacePointWriter : public WriterT<GeometryIdMultimap<T>> {
 public:
  struct Config {
    std::string collection;      ///< which collection to write
    std::string outputDir;       ///< where to place output files
    double outputScalor = 1.0;   ///< scale output values
    size_t outputPrecision = 6;  ///< floating point precision
  };

  ObjSpacePointWriter(const Config& cfg,
                      Acts::Logging::Level level = Acts::Logging::INFO);

 protected:
  ProcessCode writeT(const AlgorithmContext& context,
                     const GeometryIdMultimap<T>& spacePoints);

 private:
  // since class iitself is templated, base class template must be fixed
  using Base = WriterT<GeometryIdMultimap<T>>;

  Config m_cfg;
};

}  // namespace Obj
}  // namespace FW

template <typename T>
inline FW::Obj::ObjSpacePointWriter<T>::ObjSpacePointWriter(
    const ObjSpacePointWriter<T>::Config& cfg, Acts::Logging::Level level)
    : Base(cfg.collection, "ObjSpacePointWriter", level), m_cfg(cfg) {
  if (m_cfg.collection.empty()) {
    throw std::invalid_argument("Missing input collection");
  }
}

template <typename T>
inline FW::ProcessCode FW::Obj::ObjSpacePointWriter<T>::writeT(
    const FW::AlgorithmContext& context,
    const FW::GeometryIdMultimap<T>& spacePoints) {
  // open per-event file
  std::string path = FW::perEventFilepath(m_cfg.outputDir, "spacepoints.obj",
                                          context.eventNumber);
  std::ofstream os(path, std::ofstream::out | std::ofstream::trunc);
  if (!os) {
    throw std::ios_base::failure("Could not open '" + path + "' to write");
  }

  os << std::setprecision(m_cfg.outputPrecision);
  // count the vertex
  size_t vertex = 0;
  // loop and fill the space point data
  for (auto& volumeData : spacePoints) {
    for (auto& layerData : volumeData.second) {
      for (auto& moduleData : layerData.second) {
        for (auto& data : moduleData.second) {
          // write the space point
          os << "v " << m_cfg.outputScalor * data.x() << " "
             << m_cfg.outputScalor * data.y() << " "
             << m_cfg.outputScalor * data.z() << '\n';
          os << "p " << ++vertex << '\n';
        }
      }
    }
  }
  return ProcessCode::SUCCESS;
}
