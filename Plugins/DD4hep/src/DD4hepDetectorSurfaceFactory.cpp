// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/DD4hepDetectorSurfaceFactory.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Plugins/DD4hep/DD4hepBinningHelpers.hpp"
#include "Acts/Plugins/DD4hep/DD4hepConversionHelpers.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
#include "Acts/Plugins/TGeo/TGeoPrimitivesHelper.hpp"
#include "Acts/Plugins/TGeo/TGeoSurfaceConverter.hpp"

#include "DD4hep/DetElement.h"

using namespace Acts::detail;

Acts::DD4hepDetectorSurfaceFactory::DD4hepDetectorSurfaceFactory(
    std::unique_ptr<const Logger> mlogger)
    : m_logger(std::move(mlogger)) {
  ACTS_DEBUG("UnitLength conversion factor (DD4hep -> Acts): " << unitLength);
}

void Acts::DD4hepDetectorSurfaceFactory::construct(
    Cache& cache, const dd4hep::DetElement& dd4hepElement) {
  ACTS_DEBUG("Constructing DD4hepDetectorElements - tree level call from  "
             << dd4hepElement.name() << ".");
  recursiveConstruct(cache, dd4hepElement, 1);
  ACTS_DEBUG("Recursive search did yield: "
             << cache.sensitiveSurfaces.size() << " sensitive surfaces, "
             << cache.passiveSurfaces.size() << " passive surfaces, and "
             << cache.passiveSurfaceProxies.size()
             << " passive surface proxies.");
}

void Acts::DD4hepDetectorSurfaceFactory::recursiveConstruct(
    Cache& cache, const dd4hep::DetElement& dd4hepElement, int level) {
  ACTS_VERBOSE("Conversion call at level " << level << " for element "
                                           << dd4hepElement.name());

  // Check if any surface binnning can be detected
  int sBinning = getParamOr<int>("surface_binning_dim", dd4hepElement, 0);
  if (sBinning > 0) {
    cache.binnings = convertBinning(dd4hepElement, "surface_binning");
  }

  // Deal with passive surface if detected
  // Passive surfaces can be given through:
  // - [ direct ] translate of the DD4hepElement
  // - [ proxy ] placement through proxy parameters
  // - [ inner, outer, negative, positive, representing ]
  //    proxies to be determined at layer structure building
  bool pSurface = getParamOr<bool>("passive_surface", dd4hepElement, false);
  if (pSurface) {
    ACTS_VERBOSE("Passive surface(s) detected.");
    int pSurfaceCount =
        getParamOr<int>("passive_surface_count", dd4hepElement, 1);
    for (int ips = 0; ips < pSurfaceCount; ips++) {
      std::string pSurfaceType = getParamOr<std::string>(
          "passive_surface_n" + std::to_string(ips) + "_type", dd4hepElement,
          "direct");
      // The passive surface is a proxy, the actual dimension will be determined
      // during layer structure building
      if (std::find(allowedPassiveProxies.begin(), allowedPassiveProxies.end(),
                    pSurfaceType) != allowedPassiveProxies.end()) {
        cache.passiveSurfaceProxies.push_back(
            DD4hepPassiveSurfaceProxy{pSurfaceType});
      } else if (pSurfaceType == "direct") {
        cache.passiveSurfaces.push_back(constructPassiveElement(dd4hepElement));
      }
    }
  }

  const dd4hep::DetElement::Children& children = dd4hepElement.children();
  if (!children.empty()) {
    ACTS_VERBOSE(children.size() << " child(ren) detected.");
    for (auto& child : children) {
      dd4hep::DetElement childDetElement = child.second;
      ACTS_VERBOSE("Processing child " << childDetElement.name());
      if (childDetElement.volume().isSensitive()) {
        ACTS_VERBOSE("Sensitive surface detected.");
        cache.sensitiveSurfaces.push_back(
            constructSensitiveElement(childDetElement));
      }
      recursiveConstruct(cache, childDetElement, level + 1);
    }
  } else {
    ACTS_VERBOSE("No children detected.");
  }
}

Acts::DD4hepDetectorSurfaceFactory::DD4hepSensitiveSurface
Acts::DD4hepDetectorSurfaceFactory::constructSensitiveElement(
    const dd4hep::DetElement& dd4hepElement) const {
  // Extract the axis definition
  std::string detAxis =
      getParamOr<std::string>("axis_definitions", dd4hepElement, "XYZ");
  // Create the corresponding detector element
  auto dd4hepDetElement = std::make_shared<Acts::DD4hepDetectorElement>(
      dd4hepElement, detAxis, unitLength, false, nullptr, nullptr);
  // return the surface
  return {dd4hepDetElement, dd4hepDetElement->surface().getSharedPtr()};
}

Acts::DD4hepDetectorSurfaceFactory::DD4hepPassiveSurface
Acts::DD4hepDetectorSurfaceFactory::constructPassiveElement(
    const dd4hep::DetElement& dd4hepElement) const {
  // Underlying TGeo node, shape & transform
  const auto& tgeoNode = *(dd4hepElement.placement().ptr());
  auto tgeoShape = tgeoNode.GetVolume()->GetShape();
  const auto tgeoTransform = dd4hepElement.nominal().worldTransformation();
  // Extract the axis definition
  std::string detAxis =
      getParamOr<std::string>("axis_definitions", dd4hepElement, "XYZ");
  // Return a passive surface
  return TGeoSurfaceConverter::toSurface(*tgeoShape, tgeoTransform, detAxis,
                                         unitLength);
}
