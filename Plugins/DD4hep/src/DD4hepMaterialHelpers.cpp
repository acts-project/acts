// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/DD4hep/DD4hepMaterialHelpers.hpp"

#include "Acts/Geometry/ApproachDescriptor.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AxisSpec.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "ActsPlugins/DD4hep/DD4hepConversionHelpers.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <ostream>
#include <stdexcept>

#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

using namespace Acts;

std::shared_ptr<ProtoGridSurfaceMaterial> ActsPlugins::createProtoMaterial(
    const dd4hep::rec::VariantParameters& params, const std::string& valueTag,
    const std::vector<std::pair<const std::string, BinningOption> >& binning,
    const Logger& logger) {
  if (binning.size() != 2) {
    throw std::invalid_argument("DD4hep surface material requires two axes");
  }
  const bool cylinder = std::ranges::any_of(binning, [](const auto& axis) {
    return axisDirectionFromName(axis.first) == AxisDirection::AxisZ;
  });
  auto axisSpec = [&](const auto& axis) {
    auto direction = axisDirectionFromName(axis.first);
    if (cylinder && direction == AxisDirection::AxisPhi) {
      direction = AxisDirection::AxisRPhi;
    }
    const int bins = params.get<int>(valueTag + "_" + axis.first);
    ACTS_VERBOSE("  - material binning for " << axis.first << " on " << valueTag
                                             << ": " << bins);
    // An omitted legacy dimension is homogeneous, represented by one bin.
    // Keep ranges deferred; material mapping resolves them from the surface.
    return AxisSpec::Equidistant(
        static_cast<std::size_t>(std::max(1, bins)), std::nullopt, std::nullopt,
        axis.second == closed ? AxisBoundaryType::Closed
                              : AxisBoundaryType::Bound,
        direction);
  };
  return std::make_shared<ProtoGridSurfaceMaterial>(
      MultiAxisSpec2D({axisSpec(binning[0]), axisSpec(binning[1])}));
}

void ActsPlugins::addLayerProtoMaterial(
    const dd4hep::rec::VariantParameters& params, Layer& layer,
    const std::vector<std::pair<const std::string, BinningOption> >& binning,
    const Logger& logger) {
  ACTS_VERBOSE("addLayerProtoMaterial");
  // Start with the representing surface
  std::vector<std::string> materialOptions = {"layer_material_representing"};
  std::vector<const Surface*> materialSurfaces = {
      &(layer.surfaceRepresentation())};
  // Now fill (optionally) with the approach surfaces
  auto aDescriptor = layer.approachDescriptor();
  if (aDescriptor != nullptr && aDescriptor->containedSurfaces().size() >= 2) {
    // Add the inner and outer approach surface
    const std::vector<const Surface*>& aSurfaces =
        aDescriptor->containedSurfaces();
    materialOptions.push_back("layer_material_inner");
    materialSurfaces.push_back(aSurfaces[0]);
    materialOptions.push_back("layer_material_outer");
    materialSurfaces.push_back(aSurfaces[1]);
  }

  // Now loop over it and create the ProtoMaterial
  for (unsigned int is = 0; is < materialOptions.size(); ++is) {
    // if (actsExtension.hasValue(materialOptions[is])) {
    ACTS_VERBOSE(" - checking material for: " << materialOptions[is]);
    if (params.contains(materialOptions[is])) {
      ACTS_VERBOSE(" - have material");
      // Create the material and assign it
      auto psMaterial =
          createProtoMaterial(params, materialOptions[is], binning, logger);
      // const_cast (ugly - to be changed after internal geometry stored
      // non-const)
      Surface* surface = const_cast<Surface*>(materialSurfaces[is]);
      surface->assignSurfaceMaterial(psMaterial);
    }
  }
}

void ActsPlugins::addCylinderLayerProtoMaterial(dd4hep::DetElement detElement,
                                                Layer& cylinderLayer,
                                                const Logger& logger) {
  ACTS_VERBOSE(
      "Translating DD4hep material into Acts material for CylinderLayer : "
      << detElement.name());
  if (hasParams(detElement)) {
    ACTS_VERBOSE(" params: " << getParams(detElement));
  } else {
    ACTS_VERBOSE(" NO params");
  }
  if (getParamOr<bool>("layer_material", detElement, false)) {
    addLayerProtoMaterial(getParams(detElement), cylinderLayer,
                          {{"binPhi", closed}, {"binZ", open}}, logger);
  }
}

void ActsPlugins::addDiscLayerProtoMaterial(dd4hep::DetElement detElement,
                                            Layer& discLayer,
                                            const Logger& logger) {
  ACTS_VERBOSE("Translating DD4hep material into Acts material for DiscLayer : "
               << detElement.name());

  if (hasParams(detElement)) {
    ACTS_VERBOSE(" params: " << getParams(detElement));
  } else {
    ACTS_VERBOSE(" NO params");
  }
  if (getParamOr<bool>("layer_material", detElement, false)) {
    addLayerProtoMaterial(getParams(detElement), discLayer,
                          {{"binPhi", closed}, {"binR", open}}, logger);
  }
}
