// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/DD4hepDetectorSurfaceFactory.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Plugins/DD4hep/DD4hepBinningHelpers.hpp"
#include "Acts/Plugins/DD4hep/DD4hepConversionHelpers.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
#include "Acts/Plugins/TGeo/TGeoMaterialConverter.hpp"
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
    Cache& cache, const dd4hep::DetElement& dd4hepElement,
    const Options& options) {
  ACTS_DEBUG("Configured to convert "
             << (options.convertSensitive ? "sensitive components and " : "")
             << (options.convertPassive
                     ? "passive surfaces."
                     : (!options.convertSensitive
                            ? "nothing (this is likely a configuration error)."
                            : "")));
  ACTS_DEBUG("Constructing DD4hepDetectorElements - tree level call from  "
             << dd4hepElement.name() << ".");
  recursiveConstruct(cache, dd4hepElement, options, 1);
  ACTS_DEBUG("Recursive search did yield: "
             << cache.sensitiveSurfaces.size() << " sensitive surface(s), "
             << cache.passiveSurfaces.size() << " passive surface(s)");
}

void Acts::DD4hepDetectorSurfaceFactory::recursiveConstruct(
    Cache& cache, const dd4hep::DetElement& dd4hepElement,
    const Options& options, int level) {
  ACTS_VERBOSE("Conversion call at level " << level << " for element "
                                           << dd4hepElement.name());

  // Check if any surface binnning can be detected
  int sBinning = getParamOr<int>("acts_surface_binning_dim", dd4hepElement, 0);
  if (sBinning > 0) {
    cache.binnings = convertBinning(dd4hepElement, "acts_surface_binning");
  }

  // Deal with passive surface if detected
  bool pSurface =
      getParamOr<bool>("acts_passive_surface", dd4hepElement, false);
  if (pSurface && options.convertPassive) {
    ACTS_VERBOSE("Passive surface(s) detected.");
    cache.passiveSurfaces.push_back(
        constructPassiveComponents(dd4hepElement, options));
  }

  const dd4hep::DetElement::Children& children = dd4hepElement.children();
  if (!children.empty()) {
    ACTS_VERBOSE(children.size() << " child(ren) detected.");
    for (auto& child : children) {
      dd4hep::DetElement childDetElement = child.second;
      ACTS_VERBOSE("Processing child " << childDetElement.name());
      if (childDetElement.volume().isSensitive() && options.convertSensitive) {
        ACTS_VERBOSE("Sensitive surface detected.");
        cache.sensitiveSurfaces.push_back(
            constructSensitiveComponents(childDetElement, options));
      }
      recursiveConstruct(cache, childDetElement, options, level + 1);
    }
  } else {
    ACTS_VERBOSE("No children detected.");
  }
}

Acts::DD4hepDetectorSurfaceFactory::DD4hepSensitiveSurface
Acts::DD4hepDetectorSurfaceFactory::constructSensitiveComponents(
    const dd4hep::DetElement& dd4hepElement, const Options& options) const {
  // Extract the axis definition
  std::string detAxis =
      getParamOr<std::string>("axis_definitions", dd4hepElement, "XYZ");
  std::shared_ptr<const Acts::ISurfaceMaterial> surfaceMaterial = nullptr;

  // Create the corresponding detector element
  auto dd4hepDetElement = std::make_shared<Acts::DD4hepDetectorElement>(
      dd4hepElement, detAxis, unitLength, false, nullptr);
  auto sSurface = dd4hepDetElement->surface().getSharedPtr();
  attachSurfaceMaterial(dd4hepElement, *sSurface.get(),
                        dd4hepDetElement->thickness(), options);
  // return the surface
  return {dd4hepDetElement, sSurface};
}

Acts::DD4hepDetectorSurfaceFactory::DD4hepPassiveSurface
Acts::DD4hepDetectorSurfaceFactory::constructPassiveComponents(
    const dd4hep::DetElement& dd4hepElement, const Options& options) const {
  // Underlying TGeo node, shape & transform
  const auto& tgeoNode = *(dd4hepElement.placement().ptr());
  auto tgeoShape = tgeoNode.GetVolume()->GetShape();
  const auto tgeoTransform = dd4hepElement.nominal().worldTransformation();
  // Extract the axis definition
  auto detAxis =
      getParamOr<std::string>("axis_definitions", dd4hepElement, "XYZ");
  bool assignToAll = getParamOr<bool>("assign_to_all", dd4hepElement, true);

  auto [pSurface, thickness] =
      TGeoSurfaceConverter::toSurface(*tgeoShape, tgeoTransform, detAxis);
  attachSurfaceMaterial(dd4hepElement, *pSurface.get(), thickness, options);
  // Return a passive surface
  return {pSurface, assignToAll};
}

void Acts::DD4hepDetectorSurfaceFactory::attachSurfaceMaterial(
    const dd4hep::DetElement& dd4hepElement, Acts::Surface& surface,
    ActsScalar thickness, const Options& options) const {
  if (options.convertMaterial) {
    // Extract the material
    const auto& tgeoNode = *(dd4hepElement.placement().ptr());
    auto tgeoMaterial = tgeoNode.GetMedium()->GetMaterial();
    // Convert the material
    TGeoMaterialConverter::Options materialOptions;
    materialOptions.unitLengthScalor = unitLength;
    auto materialSlab = TGeoMaterialConverter::materialSlab(
        *tgeoMaterial, thickness, options.surfaceMaterialThickness,
        materialOptions);
    auto surfaceMaterial =
        std::make_shared<HomogeneousSurfaceMaterial>(materialSlab);
    // Assign the material to the surface
    surface.assignSurfaceMaterial(std::move(surfaceMaterial));
  }
}
