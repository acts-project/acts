// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/ProtoDetector.hpp"

#include "Acts/Utilities/Enumerate.hpp"

#include <exception>
#include <sstream>

void Acts::ProtoVolume::extendUp(Acts::ProtoVolume& ptVolume) {
  ptVolume.extent.extend(extent);

  for (auto& cv : constituentVolumes) {
    ptVolume.extent.extend(cv.extent);
    cv.extendUp(ptVolume);
  }
}

void Acts::ProtoVolume::propagateMinDown(BinningValue bValue) {
  for (auto& cv : constituentVolumes) {
    cv.extent.set(bValue, extent.min(bValue), cv.extent.max(bValue));
    cv.propagateMinDown(bValue);
  }
}

void Acts::ProtoVolume::propagateMaxDown(BinningValue bValue) {
  for (auto& cv : constituentVolumes) {
    cv.extent.set(bValue, cv.extent.min(bValue), extent.max(bValue));
    cv.propagateMaxDown(bValue);
  }
}

void Acts::ProtoVolume::constrainDown(const Acts::ProtoVolume& ptVolume) {
  extent.addConstrain(ptVolume.extent);
  for (auto& cv : constituentVolumes) {
    cv.extent.addConstrain(extent);
  }
}

void Acts::ProtoVolume::harmonize(bool legacy) {
  std::vector<BinningValue> otherConstrains;

  // Deal with the constituents
  if (not constituentVolumes.empty()) {
    if (constituentBinning.empty()) {
      std::string errorMsg = std::string("ProtoVolume '") + name +
                             std::string("' with constituents, but no binning");
      throw std::runtime_error(errorMsg);
    }

    // For legacy volumes, check if layers are present
    bool layersPresent = layerContainer;
    for (const auto& cv : constituentVolumes) {
      layersPresent =
          layersPresent or cv.layerType != Surface::SurfaceType::Other;
      if (layersPresent) {
        break;
      }
    }

    // If layers are present, it can't be a container in the legacy style
    layerContainer = layersPresent;
    auto binValue = constituentBinning[0].binvalue;
    // Set the first last
    auto& fVolume = *constituentVolumes.begin();
    auto& lVolume = constituentVolumes.back();

    std::vector<float> borders = {};

    // The volumes should be harmonized in all other constraining values
    for (auto obValue : s_binningValues) {
      if (obValue != binValue and extent.constrains(obValue)) {
        otherConstrains.push_back(obValue);
      }
    }

    // Legacy  conversion - layers are kept untouched
    if (not layerContainer) {
      // Set the outer boundaries
      fVolume.extent.set(binValue, extent.min(binValue),
                         fVolume.extent.max(binValue));
      lVolume.extent.set(binValue, lVolume.extent.min(binValue),
                         extent.max(binValue));
      // Align the containers
      borders.push_back(static_cast<float>(fVolume.extent.min(binValue)));
      for (unsigned int iv = 1; iv < constituentVolumes.size(); ++iv) {
        auto& lv = constituentVolumes[iv - 1u];
        ActsScalar zero = lv.extent.min(binValue);
        ActsScalar low = lv.extent.max(binValue);

        auto& hv = constituentVolumes[iv];
        ActsScalar high = hv.extent.min(binValue);
        ActsScalar mid = 0.5 * (low + high);
        ActsScalar max = hv.extent.max(binValue);
        lv.extent.set(binValue, zero, mid);
        hv.extent.set(binValue, mid, max);
        borders.push_back(mid);
      }
      borders.push_back(constituentVolumes.back().extent.max(binValue));

    } else if (layerContainer and not legacy) {
      // Count the gaps
      std::size_t gaps = 0;
      std::vector<float> boundaries = {};
      // New container vector
      std::vector<ProtoVolume> updatedConstituents;
      ActsScalar containerMin = extent.min(binValue);
      if (fVolume.extent.min(binValue) > containerMin) {
        ProtoVolume gap;
        gap.name = name + "-gap-" + std::to_string(gaps++);
        gap.extent.set(binValue, containerMin, fVolume.extent.min(binValue));
        updatedConstituents.push_back(gap);
        borders.push_back(static_cast<float>(containerMin));
      }
      // Fill the gaps
      for (unsigned int iv = 1; iv < constituentVolumes.size(); ++iv) {
        auto& lv = constituentVolumes[iv - 1u];
        // This volume is one to save
        updatedConstituents.push_back(lv);
        borders.push_back(static_cast<float>(lv.extent.min(binValue)));
        // check if a gap to the next is needed
        ActsScalar low = lv.extent.max(binValue);
        auto& hv = constituentVolumes[iv];
        ActsScalar high = hv.extent.min(binValue);
        if (high > low) {
          ProtoVolume gap;
          gap.name = name + "-gap-" + std::to_string(gaps++);
          gap.extent.set(binValue, low, high);
          updatedConstituents.push_back(gap);
          borders.push_back(static_cast<float>(low));
        }
      }
      ActsScalar constituentsMax = lVolume.extent.max(binValue);
      updatedConstituents.push_back(lVolume);
      borders.push_back(static_cast<float>(constituentsMax));
      // Check the container min/max setting
      ActsScalar containerMax = extent.max(binValue);
      if (constituentsMax < containerMax) {
        ProtoVolume gap;
        gap.name = name + "-gap-" + std::to_string(gaps++);
        gap.extent.set(binValue, constituentsMax, containerMax);
        updatedConstituents.push_back(gap);
        borders.push_back(static_cast<float>(containerMax));
      }
      constituentVolumes = updatedConstituents;
    } else if (legacy and layerContainer) {
      borders = {0., 1.};
    }
    constituentBinning = {
        BinningData(constituentBinning[0].option, binValue, borders)};
  }

  // Harmonize downwards
  for (auto& cv : constituentVolumes) {
    cv.extent.extend(extent, otherConstrains);
    cv.harmonize(legacy);
  }
}

bool Acts::ProtoVolume::operator==(const Acts::ProtoVolume& pv) const {
  // Simple checks
  if (name != pv.name or extent != pv.extent or
      layerContainer != pv.layerContainer or layerType != pv.layerType) {
    return false;
  }
  if (layerSurfaceBinning != pv.layerSurfaceBinning) {
    return false;
  }
  if (constituentVolumes != pv.constituentVolumes or
      constituentBinning != pv.constituentBinning) {
    return false;
  }
  return true;
}

std::string Acts::ProtoVolume::toString(const std::string& indent) const {
  std::string subIndent("  ");
  std::stringstream ss;
  ss << indent << "> volume: " << name << '\n';
  ss << indent << "  extent: ";
  ss << extent.toString(indent) << '\n';
  if (not constituentVolumes.empty()) {
    ss << indent << "  container of " << constituentVolumes.size()
       << " constituents. " << '\n';
    ss << indent << "  constituent binning:" << '\n';
    for (const auto& cb : constituentBinning) {
      ss << cb.toString(indent) << '\n';
    }
    ss << indent << "  constituents are:" << '\n';
    for (auto cv : constituentVolumes) {
      ss << cv.toString(indent + subIndent) << '\n';
    }
  }
  return ss.str();
}

std::string Acts::ProtoDetector::toString(const std::string& indent) const {
  std::string subIndent("  ");
  std::stringstream ss;
  ss << indent << "> detector: " << name << '\n';
  ss << worldVolume.toString(indent + subIndent) << '\n';
  return ss.str();
}