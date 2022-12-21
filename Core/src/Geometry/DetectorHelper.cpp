// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/DetectorHelper.hpp"
#include "Acts/Geometry/CylindricalDetectorHelper.hpp"

#include <exception>

Acts::Experimental::ProtoContainer Acts::Experimental::connectDetectorVolumes(
    const GeometryContext& gctx, BinningValue bValue,
    std::vector<std::shared_ptr<Experimental::DetectorVolume>>& volumes,
    const std::vector<unsigned int>& selectedOnly,
    Acts::Logging::Level logLevel) {
  // Simple dispatcher to make client code a bit more readable
  switch (bValue) {
    case binR: {
      return connectDetectorVolumesInR(gctx, volumes, selectedOnly, logLevel);
    };
    case binPhi: {
      return connectDetectorVolumesInPhi(gctx, volumes, selectedOnly, logLevel);
    };
    case binZ: {
      return connectDetectorVolumesInZ(gctx, volumes, selectedOnly, logLevel);
    };
    default: {
      throw std::invalid_argument("DetectorHelper: connector mode not implemented.");
    }; 
  }
  return ProtoContainer{};
}

Acts::Experimental::ProtoContainer Acts::Experimental::connectContainers(
    const GeometryContext& gctx, BinningValue bValue,
    const std::vector<ProtoContainer>& containers,
    const std::vector<unsigned int>& selectedOnly,
    Acts::Logging::Level logLevel) {
  // Simple dispatcher to make client code a bit more readable
  switch (bValue) {
    case binR: {
      return connectContainersInR(gctx, containers, selectedOnly, logLevel);
    };
    case binZ: {
      return connectContainersInZ(gctx, containers, selectedOnly, logLevel);
    };
    default: {
      throw std::invalid_argument("DetectorHelper: connector mode not implemented.");
    }; 
  }
  return ProtoContainer{};
}