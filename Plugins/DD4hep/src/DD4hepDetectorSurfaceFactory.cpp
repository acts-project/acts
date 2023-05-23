// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/DD4hepDetectorSurfaceFactory.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
#include "Acts/Plugins/DD4hep/detail/DD4hepConversionHelpers.hpp"
#include "Acts/Plugins/TGeo/TGeoPrimitivesHelper.hpp"

#include "DD4hep/DetElement.h"

Acts::DD4hepDetectorSurfaceFactory::DD4hepDetectorSurfaceFactory(
    std::unique_ptr<const Logger> logger)
    : m_logger(std::move(logger)) {}

void Acts::DD4hepDetectorSurfaceFactory::construct(
    Cache& cache, const dd4hep::DetElement& dd4hepElement,
    const Options& option) {
  ACTS_DEBUG("Constructing DD4hepDetectorElements - tree level call.");
  recursiveConstruct(cache, dd4hepElement, option, 1);
}

void Acts::DD4hepDetectorSurfaceFactory::recursiveConstruct(
    Cache& cache, const dd4hep::DetElement& dd4hepElement,
    const Options& option, int level) {
  using Binning = Experimental::LayerStructureBuilder::Binning;

  ACTS_VERBOSE("Conversion call at level " << level << " for element "
                                           << dd4hepElement.name());
  bool detAxis = getParamOr<bool>("surface_binning", dd4hepElement, false);
  if (detAxis) {
    // Check the bins to be be converted, cascade throught he options
    using SurfaceBinning =
        std::tuple<Acts::BinningValue, Acts::BinningOption, std::string>;
    // The possible binnings
    std::vector<SurfaceBinning> possibleBinnings = {{binX, open, "x"},
                                                    {binY, open, "y"},
                                                    {binZ, open, "z"},
                                                    {binR, open, "r"},
                                                    {binPhi, closed, "phi"}};
    // Check the possible bindings
    for (auto [bVal, bOpt, bStr] : possibleBinnings) {
      // Get the number of bins
      int nB =
          getParamOr<int>("surface_binning_" + bStr + "_n", dd4hepElement, 0);
      // Get the expansion
      int nE =
          getParamOr<int>("surface_binning_" + bStr + "_exp", dd4hepElement, 0);
      // Get the binning type
      std::string binningType = getParamOr<std::string>(
          "surface_binning_" + bStr + "_type", dd4hepElement, "equidistant");
      // Equidistant case
      if (nB > 0 and binningType == "equidistant") {
        ACTS_DEBUG("Binning " << bStr << " is equidistant.");
        ActsScalar minVal = getParamOr<ActsScalar>(
            "surface_binning_" + bStr + "_min", dd4hepElement, 0.);
        ActsScalar maxVal = getParamOr<ActsScalar>(
            "surface_binning_" + bStr + "_max", dd4hepElement, 1.);
        // Force to [-pi, pi)]
        if (bVal == binPhi) {
          minVal = -M_PI;
          maxVal = M_PI;
        }
        cache.binnings.push_back(
            Binning{BinningData(bOpt, bVal, nB, minVal, maxVal),
                    static_cast<std::size_t>(nE)});
      } else if (nB > 0) {
        ACTS_DEBUG("Binning " << bStr << " is variable.");

        // Reconstruct the boundary vector
        std::vector<float> boundaries = {};
        boundaries.reserve(nB + 1);
        for (auto ib = 0; ib <= nB; ++ib) {
          float val = getParamOr<ActsScalar>(
              "surface_binning_" + bStr + "_b" + std::to_string(ib),
              dd4hepElement, 0.);
          boundaries.push_back(val);
        }
        // Force to [-pi, pi)
        if (bVal == binPhi) {
          boundaries.front() = -M_PI;
          boundaries.back() = M_PI;
        }
        cache.binnings.push_back(Binning{BinningData(bOpt, bVal, boundaries),
                                         static_cast<std::size_t>(nE)});
      }
    }
    ACTS_VERBOSE("Surface binning in " << cache.binnings.size()
                                       << " dimensions.");
  }

  const dd4hep::DetElement::Children& children = dd4hepElement.children();
  if (!children.empty()) {
    ACTS_VERBOSE(children.size() << " children detected.");
    for (auto& child : children) {
      dd4hep::DetElement childDetElement = child.second;
      if (childDetElement.volume().isSensitive()) {
        cache.sensitiveSurfaces.push_back(
            constructSensitiveElement(childDetElement, option));
      }
      recursiveConstruct(cache, childDetElement, option, level + 1);
    }
  }
}

Acts::DD4hepDetectorSurfaceFactory::DD4hepSensitiveSurface
Acts::DD4hepDetectorSurfaceFactory::constructSensitiveElement(
    const dd4hep::DetElement& dd4hepElement, const Options& options) const {
  std::string detAxis =
      getParamOr<std::string>("axis_definitions", dd4hepElement, "XYZ");
  // Create the corresponding detector element
  auto dd4hepDetElement = std::make_shared<Acts::DD4hepDetectorElement>(
      dd4hepElement, detAxis, options.unitLength, false, nullptr, nullptr);
  // return the surface
  return {dd4hepDetElement, dd4hepDetElement->surface().getSharedPtr()};
}

Acts::DD4hepDetectorSurfaceFactory::DD4hepPassiveSurface
Acts::DD4hepDetectorSurfaceFactory::constructPassiveElement(
    const dd4hep::DetElement& /*dd4hepElemen*/,
    const Options& /*options*/) const {
  return {};
}
