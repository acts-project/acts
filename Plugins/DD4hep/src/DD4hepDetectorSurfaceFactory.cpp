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
#include "Acts/Plugins/TGeo/TGeoSurfaceConverter.hpp"

#include "DD4hep/DetElement.h"

Acts::DD4hepDetectorSurfaceFactory::DD4hepDetectorSurfaceFactory(
    std::unique_ptr<const Logger> mlogger)
    : m_logger(std::move(mlogger)) {
  ACTS_DEBUG("UnitLength conversion factor (DD4hep -> Acts): " << unitLength);
}

void Acts::DD4hepDetectorSurfaceFactory::construct(
    Cache& cache, const dd4hep::DetElement& dd4hepElement) {
  ACTS_DEBUG("Constructing DD4hepDetectorElements - tree level call.");
  recursiveConstruct(cache, dd4hepElement, 1);
}

void Acts::DD4hepDetectorSurfaceFactory::recursiveConstruct(
    Cache& cache, const dd4hep::DetElement& dd4hepElement, int level) {
  using Binning = Experimental::LayerStructureBuilder::Binning;

  ACTS_VERBOSE("Conversion call at level " << level << " for element "
                                           << dd4hepElement.name());
  // Deal with surface binning if detected
  bool sBinning = getParamOr<bool>("surface_binning", dd4hepElement, false);
  if (sBinning) {
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
  // Deal with passive surface if detected
  bool pSurface = getParamOr<bool>("passive_surface", dd4hepElement, false);
  if (pSurface) {
    cache.passiveSurfaces.push_back(constructPassiveElement(dd4hepElement));
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
