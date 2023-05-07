// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/DD4hepLayerStructure.hpp"

#include "Acts/Plugins/DD4hep/DD4hepDetectorSurfaceFactory.hpp"

Acts::Experimental::DD4hepLayerStructure::DD4hepLayerStructure(
    std::shared_ptr<DD4hepDetectorSurfaceFactory> surfaceFactory)
    : m_surfaceFactory(std::move(surfaceFactory)) {
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
  DD4hepDetectorSurfaceFactory::Options fOptions;
  fOptions.unitLength = options.unitLength;
  m_surfaceFactory->construct(fCache, dd4hepElement, fOptions);

  // The constructed surfaces and detector elements
  DD4hepLayerStructure::Surfaces cStructure;
  cStructure.surfaces.reserve(fCache.sensitiveSurfaces.size() +
                              fCache.passiveSurfaces.size());

  std::vector<std::shared_ptr<DD4hepDetectorElement>> cElements;
  cElements.reserve(fCache.sensitiveSurfaces.size());

  // Fill them in to the surface provider struct and detector store
  for (auto [de, ds] : fCache.sensitiveSurfaces) {
    cStructure.surfaces.push_back(ds);
    cElements.push_back(de);
  }
  dd4hepStore[options.name] = cElements;

  // Passive surfaces
  cStructure.surfaces.insert(cStructure.surfaces.end(),
                             fCache.passiveSurfaces.begin(),
                             fCache.passiveSurfaces.end());

  // Surfaces are prepared for creating the builder
  LayerStructureBuilder::Config lsbConfig;
  lsbConfig.auxilliary = "*** DD4hep driven builder for: ";
  lsbConfig.auxilliary += options.name;
  lsbConfig.surfaces = cStructure;

  // Translate the binnings
  if (!options.binnings.empty()) {
    lsbConfig.binnings = options.binnings;
  }
  lsbConfig.supports = options.supports;

  return std::make_shared<LayerStructureBuilder>(
      lsbConfig, getDefaultLogger(options.name, options.logLevel));
}
