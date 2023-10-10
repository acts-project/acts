// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/DD4hepVolumeStructure.hpp"

#include "Acts/Plugins/DD4hep/DD4hepConversionHelpers.hpp"

#include <DD4hep/DetElement.h>

Acts::Experimental::DD4hepVolumeStructure::DD4hepVolumeStructure(
    std::unique_ptr<const Logger> mlogger)
    : m_logger(std::move(mlogger)) {
  ACTS_DEBUG("UnitLength conversion factor (DD4hep -> Acts): " << unitLength);
}

std::shared_ptr<Acts::Experimental::VolumeStructureBuilder>
Acts::Experimental::DD4hepVolumeStructure::builder(
    const dd4hep::DetElement& dd4hepElement, const Options& options) const {
  // The configuration
  VolumeStructureBuilder::Config vsbConfig;
  if (recursiveParse(vsbConfig, dd4hepElement)) {
    ACTS_DEBUG("ProtoVolume description successfully parsed.");
  }
  // Return the structure builder
  return std::make_shared<VolumeStructureBuilder>(
      vsbConfig, getDefaultLogger(options.name, options.logLevel));
}

bool Acts::Experimental::DD4hepVolumeStructure::recursiveParse(
    VolumeStructureBuilder::Config& vsbConfig,
    const dd4hep::DetElement& dd4hepElement) const {
  // Deal with a proto volume if detected
  bool actsVolume = getParamOr<bool>("acts_volume", dd4hepElement, false);
  if (actsVolume) {
    auto bValueInt =
        getParamOr<int>("acts_volume_type", dd4hepElement,
                        static_cast<int>(VolumeBounds::BoundsType::eOther));
    auto bValueType = static_cast<VolumeBounds::BoundsType>(bValueInt);
    if (bValueType < VolumeBounds::BoundsType::eOther) {
      auto bValues = extractSeries<ActsScalar>(
          dd4hepElement, "acts_volume_bvalues", unitLength);
      // Set the parameters for the config
      vsbConfig.boundsType = bValueType;
      vsbConfig.boundValues = bValues;
    }
    vsbConfig.transform =
        extractTransform(dd4hepElement, "acts_volume_pos", unitLength);
    return true;
  } else {
    ACTS_VERBOSE("No ProtoVolume description detected.");
  }

  // Parse through children
  const dd4hep::DetElement::Children& children = dd4hepElement.children();
  if (!children.empty()) {
    ACTS_VERBOSE(children.size() << " children detected.");
    for (auto& child : children) {
      dd4hep::DetElement childDetElement = child.second;
      if (recursiveParse(vsbConfig, childDetElement)) {
        return true;
      }
    }
  }

  return false;
}
