// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"
#include "Acts/Plugins/TGeo/TGeoItkModuleSplitter.hpp"


Acts::TGeoItkModuleSplitter::TGeoItkModuleSplitter(
    const TGeoItkModuleSplitter::Config& cfg,
    std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {}


/// If applicable, returns a split detector element 
inline std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>> Acts::TGeoItkModuleSplitter::split(
  const GeometryContext& /*gctx*/,
  std::shared_ptr<const Acts::TGeoDetectorElement> detElement) const {

  // Is the current node covered by this splitter?
  const TGeoNode& node = detElement->tgeoNode();
  auto sensorName = std::string(node.GetName());
  if(sensorName.find("Brl") != std::string::npos) {
    // First parameter for every node is the splitting multiplicity
    return Acts::TGeoItkModuleSplitter::splitBarrelModule(detElement, m_cfg.paramMap.at("Brl"));
  }

  ACTS_DEBUG("Node not found in splitting config.");
  return {detElement};
}

/// If applicable, returns a split detector element 
inline std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>> Acts::TGeoItkModuleSplitter::splitBarrelModule(
  std::shared_ptr<const Acts::TGeoDetectorElement> detElement, 
  unsigned int nSegments) const {

  // Retrive the surface
  auto identifier = detElement->identifier();
  const Surface& surface = detElement->surface();
  const SurfaceBounds& bounds = surface.bounds();
  if (bounds.type() != SurfaceBounds::eRectangle or nSegments <= 1u) {
    ACTS_DEBUG("Invalid splitting config for node");
    return {detElement};
  }

  std::cout << nSegments << std::endl;
  // Get the geometric information
  double thickness = detElement->thickness();
  const Transform3& transform = surface.transform({});

  // Determine the new bounds (Y should be along global z)
  const std::vector<double> boundsValues = bounds.values();
  double lengthY = boundsValues[RectangleBounds::eMaxY] 
                  - boundsValues[RectangleBounds::eMinY];
  double lengthX = (boundsValues[RectangleBounds::eMaxX] 
                  - boundsValues[RectangleBounds::eMinX]) / nSegments;
  // Add the translation for every new surface to the transform
  std::vector<std::shared_ptr<const TGeoDetectorElement>> detElements = {};

  auto localTranslation = Vector3(-nSegments*lengthX, 0., 0.);
  const auto step = Vector3(lengthX, 0., 0.);
  for (size_t i = 0; i < nSegments; i++) {
    auto rectBounds = std::make_shared<RectangleBounds>(lengthX, lengthY);
    auto elemTransform = Transform3(transform).translate(localTranslation);
    localTranslation += step;

    auto element = std::make_shared<const TGeoDetectorElement>(identifier, detElement->tgeoNode(), elemTransform, rectBounds, thickness);
    detElements.push_back(element);
  }
  return detElements;
}

/// If applicable, returns a split detector element 
/*inline std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>> Acts::TGeoItkModuleSplitter::splitDiskModule(
    std::shared_ptr<const Acts::TGeoDetectorElement> detElement, std::vector<double>& params) const {
    return {detElement};
  }*/

