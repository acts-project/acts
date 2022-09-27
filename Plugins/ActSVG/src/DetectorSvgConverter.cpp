// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ActSVG/DetectorSvgConverter.hpp"

#include "Acts/Experimental/Detector.hpp"
#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/Portal.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Plugins/ActSVG/DetectorVolumeSvgConverter.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <map>

Acts::Svg::ProtoDetector Acts::Svg::DetectorConverter::convert(
    const GeometryContext& gctx, const Experimental::Detector& detector,
    const DetectorConverter::Options& detectorOptions) {
  // The local logger
  ACTS_LOCAL_LOGGER(
      getDefaultLogger("DetectorVolumeConverter", detectorOptions.logLevel));

  ProtoDetector pDetector;

  const auto volumes = detector.volumes();

  // Create the index map vor volumes
  std::map<const Experimental::DetectorVolume*, unsigned int> volumeIndices;
  for (auto [iv, v] : enumerate(volumes)) {
    volumeIndices[v] = iv;
  }

  // Create the index map for portals
  unsigned int ip = 0;
  std::map<const Experimental::Portal*, unsigned int> portalIndices;
  for (const auto& v : volumes) {
    for (const auto& p : v->portals()) {
      portalIndices[p] = ip++;
    }
  }

  DetectorVolumeConverter::Options volumeOptions =
      detectorOptions.volumeOptions;
  volumeOptions.portalIndices = portalIndices;
  volumeOptions.portalOptions.volumeIndices = volumeIndices;

  // Loop over the list of volumes, remember internal ones
  for (auto [iv, v] : enumerate(volumes)) {
    auto pVolume = convert(gctx, *v, volumeOptions);
    pVolume._index = iv;
    pDetector._volumes.push_back(pVolume);
  }

  return pDetector;
}
