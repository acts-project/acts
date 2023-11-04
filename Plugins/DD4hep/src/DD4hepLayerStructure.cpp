// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/DD4hepLayerStructure.hpp"

Acts::Experimental::DD4hepLayerStructure::DD4hepLayerStructure(
    std::shared_ptr<DD4hepDetectorSurfaceFactory> surfaceFactory,
    std::unique_ptr<const Logger> logger)
    : m_surfaceFactory(std::move(surfaceFactory)), m_logger(std::move(logger)) {
  if (m_surfaceFactory == nullptr) {
    throw std::invalid_argument(
        "DD4hepLayerStructure: no surface factory provided");
  }
}

std::shared_ptr<Acts::Experimental::LayerStructureBuilder>
Acts::Experimental::DD4hepLayerStructure::builder(
    DD4hepDetectorElement::Store& dd4hepStore,
    const dd4hep::DetElement& dd4hepElement, const Options& options) const {
  // Check for misconfiguration with double naming
  if (dd4hepStore.find(options.name) != dd4hepStore.end()) {
    std::string reMessage = "DD4hepLayerStructure: structure with name '";
    reMessage += options.name;
    reMessage += "' already registered in DetectorElementStore";
    throw std::runtime_error(reMessage.c_str());
  }

  // This object is going to be filled with the created surfaces
  DD4hepDetectorSurfaceFactory::Cache fCache;
  m_surfaceFactory->construct(fCache, dd4hepElement, options.conversionOptions);

  ACTS_DEBUG("Conversion from DD4Hep : " << fCache.sensitiveSurfaces.size()
                                         << " sensitive surfaces");

  ACTS_DEBUG("Conversion from DD4Hep : " << fCache.passiveSurfaces.size()
                                         << " passive surfaces");

  // Check if binning was provided or detected
  if (fCache.binnings.empty() and
      (fCache.sensitiveSurfaces.size() + fCache.passiveSurfaces.size()) > 0u) {
    ACTS_VERBOSE(
        "Surface binning neither provided nor found, navigation will be "
        "'tryAll' (could result in slow navigation).");
  }

  // Surfaces are prepared for creating the builder
  LayerStructureBuilder::Config lsbConfig;
  lsbConfig.auxiliary = "*** DD4hep driven builder for: ";
  lsbConfig.auxiliary += options.name;
  // Translate binings and supports
  lsbConfig.binnings = fCache.binnings;
  lsbConfig.supports = fCache.supports;

  std::vector<std::shared_ptr<Surface>> lSurfaces;
  lSurfaces.reserve(fCache.sensitiveSurfaces.size() +
                    fCache.passiveSurfaces.size());

  std::vector<std::shared_ptr<DD4hepDetectorElement>> cElements;
  cElements.reserve(fCache.sensitiveSurfaces.size());

  // Fill them in to the surface provider struct and detector store
  for (auto [de, ds] : fCache.sensitiveSurfaces) {
    lSurfaces.push_back(ds);
    cElements.push_back(de);
  }
  dd4hepStore[options.name] = cElements;

  // Passive surfaces to be added
  for (auto [ps, toAll] : fCache.passiveSurfaces) {
    // Passive surface is not declared to be added to all navigation bins
    if (not toAll) {
      lSurfaces.push_back(ps);
    } else {
      // Passive surface is indeed declared to be added to all navigaiton bins
      Experimental::ProtoSupport pSupport;
      pSupport.surface = ps;
      pSupport.assignToAll = true;
      lsbConfig.supports.push_back(pSupport);
    }
  }

  lsbConfig.surfacesProvider =
      std::make_shared<Experimental::LayerStructureBuilder::SurfacesHolder>(
          lSurfaces);

  ACTS_DEBUG("Configured with " << lsbConfig.binnings.size() << " binnings.");
  ACTS_DEBUG("Configured to build " << lsbConfig.supports.size()
                                    << " supports.");

  // Return the structure builder
  return std::make_shared<LayerStructureBuilder>(
      lsbConfig, getDefaultLogger(options.name, options.logLevel));
}
