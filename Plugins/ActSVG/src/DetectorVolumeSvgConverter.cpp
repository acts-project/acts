// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ActSVG/DetectorVolumeSvgConverter.hpp"

#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/Portal.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <set>

Acts::Svg::ProtoVolume Acts::Svg::DetectorVolumeConverter::convert(
    const GeometryContext& gctx, const Experimental::DetectorVolume& dVolume,
    const DetectorVolumeConverter::Options& volumeOptions) {
  // The local logger
  ACTS_LOCAL_LOGGER(
      getDefaultLogger("DetectorVolumeSvgConverter", volumeOptions.logLevel));

  ProtoVolume pVolume;

  auto volumeTransform = dVolume.transform(gctx);

  // The detector volume is of cylindrical shape
  const auto& boundValues = dVolume.volumeBounds().values();
  if (dVolume.volumeBounds().type() == Acts::VolumeBounds::eCylinder) {
    // we keep 6 for the moment
    for (unsigned int ib = 0; ib < 2u; ++ib) {
      pVolume._bound_values.push_back(
          static_cast<actsvg::scalar>(boundValues[ib]));
    }
    pVolume._bound_values.push_back(
        static_cast<actsvg::scalar>(volumeTransform.translation().z()));
    for (unsigned int ib = 2u; ib < 5u; ++ib) {
      pVolume._bound_values.push_back(
          static_cast<actsvg::scalar>(boundValues[ib]));
    }
  }

  // Adding the portals to the lot
  for (auto [ip, p] : enumerate(dVolume.portals())) {
    auto pPortal = PortalConverter::convert(gctx, *p, volumeOptions.portalOptions);
    auto pPortalCandidate = volumeOptions.portalIndices.find(p);
    if (pPortalCandidate != volumeOptions.portalIndices.end()){
      pPortal._name = "portal_" + std::to_string(pPortalCandidate->second);
    } else {
      pPortal._name = dVolume.name() + "_portal_" + std::to_string(ip);
    }
    pVolume._portals.push_back(pPortal);
  }

  return pVolume;
}
