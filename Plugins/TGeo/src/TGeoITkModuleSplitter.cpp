// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"
#include "Acts/Plugins/TGeo/TGeoITkModuleSplitter.hpp"


Acts::TGeoITkModuleSplitter::TGeoITkModuleSplitter(
    const TGeoITkModuleSplitter::Config& cfg,
    std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {}

/// If applicable, returns a split detector element 
inline std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>> Acts::TGeoITkModuleSplitter::split(
  const GeometryContext& gctx,
  std::shared_ptr<const Acts::TGeoDetectorElement> detElement) const {

  // Is the current node covered by this splitter?
  const TGeoNode& node = detElement->tgeoNode();
  auto sensorName = std::string(node.GetName());
  
  if(sensorName.find("BRL") != std::string::npos) {
    if(sensorName.find("MS") != std::string::npos) {
      return Acts::TGeoITkModuleSplitter::splitBarrelModule(gctx, detElement,m_cfg.barrelMap.at("MS"));
    }
    if(sensorName.find("SS") != std::string::npos) {
      return Acts::TGeoITkModuleSplitter::splitBarrelModule(gctx, detElement,m_cfg.barrelMap.at("SS"));
    }
  }
  if(sensorName.find("EC") != std::string::npos) {
    if(sensorName.find("Sensor0") != std::string::npos) {
      return Acts::TGeoITkModuleSplitter::splitDiscModule(gctx, detElement,
      m_cfg.discMap.at("EC0"));
    }
    if(sensorName.find("Sensor1") != std::string::npos) {
      return Acts::TGeoITkModuleSplitter::splitDiscModule(gctx, detElement,
      m_cfg.discMap.at("EC1"));
    }
    if(sensorName.find("Sensor2") != std::string::npos) {
      return Acts::TGeoITkModuleSplitter::splitDiscModule(gctx, detElement,
      m_cfg.discMap.at("EC2"));
    }
    if(sensorName.find("Sensor3") != std::string::npos) {
      return Acts::TGeoITkModuleSplitter::splitDiscModule(gctx, detElement,
      m_cfg.discMap.at("EC3"));
    }
    if(sensorName.find("Sensor4") != std::string::npos) {
      return Acts::TGeoITkModuleSplitter::splitDiscModule(gctx, detElement,
      m_cfg.discMap.at("EC4"));
    }
    if(sensorName.find("Sensor5") != std::string::npos) {
      return Acts::TGeoITkModuleSplitter::splitDiscModule(gctx, detElement,
      m_cfg.discMap.at("EC5"));
    }
  }

  return {detElement};
}

/// If applicable, returns a split detector element 
inline std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>> Acts::TGeoITkModuleSplitter::splitBarrelModule(
  const GeometryContext& gctx,
  std::shared_ptr<const Acts::TGeoDetectorElement> detElement,
  unsigned int nSegments) const {

  // Retrive the surface
  auto identifier = detElement->identifier();
  const Surface& surface = detElement->surface();
  const SurfaceBounds& bounds = surface.bounds();
  if (bounds.type() != SurfaceBounds::eRectangle or nSegments <= 1u) {
    ACTS_WARNING("Invalid splitting config for barrel node: " 
                 + std::string(detElement->tgeoNode().GetName())
                 + "! Node will not be slpit.");
    return {detElement};
  }

  // Output container for the submodules
  std::vector<std::shared_ptr<const TGeoDetectorElement>> detElements = {};
  detElements.reserve(nSegments);

  // Get the geometric information
  double thickness = detElement->thickness();
  const Transform3& transform = surface.transform(gctx);
  // Determine the new bounds
  const std::vector<double> boundsValues = bounds.values();
  double lengthX = boundsValues[RectangleBounds::eMaxX] 
                  - boundsValues[RectangleBounds::eMinX];
  double lengthY = (boundsValues[RectangleBounds::eMaxY] 
                  - boundsValues[RectangleBounds::eMinY]) / nSegments;
  auto rectBounds = std::make_shared<RectangleBounds>(0.5*lengthX,
                                                      0.5*lengthY);
  // Translation for every subelement
  auto localTranslation = Vector2(0., -0.5*lengthY * (nSegments - 1));
  const auto step = Vector2(0., lengthY);

  for (size_t i = 0; i < nSegments; i++) {
    Vector3 globalTranslation = surface.localToGlobal(gctx, localTranslation,
                                                      {})
                                - transform.translation();
    auto elemTransform = Transform3(transform).pretranslate(globalTranslation);
    auto element = std::make_shared<const TGeoDetectorElement>(identifier, detElement->tgeoNode(), elemTransform, rectBounds, thickness);
    detElements.push_back(std::move(element));

    localTranslation += step;
  }
  return detElements;
}

/// If applicable, returns a split detector element 
inline std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>> Acts::TGeoITkModuleSplitter::splitDiscModule(
    const GeometryContext& gctx,
    std::shared_ptr<const Acts::TGeoDetectorElement> detElement,
    const std::vector<SplitRange>& splitRanges) const {

  // Retrive the surface
  auto identifier = detElement->identifier();
  const Surface& surface = detElement->surface();
  const SurfaceBounds& bounds = surface.bounds();

  if (bounds.type() != SurfaceBounds::eAnnulus or splitRanges.empty()) {
    ACTS_WARNING("Invalid splitting config for disk node: " 
                 + std::string(detElement->tgeoNode().GetName())
                 + "! Node will not be slpit.");
    return {detElement};
  }

  auto nSegments = splitRanges.size();

  // Output container for the submodules
  std::vector<std::shared_ptr<const TGeoDetectorElement>> detElements = {};
  detElements.reserve(nSegments);

  // Get the geometric information
  double thickness = detElement->thickness();
  const Transform3& transform = surface.transform(gctx);
  const std::vector<double> boundsValues = bounds.values();
  std::array<double, AnnulusBounds::eSize> values;
  std::copy_n(boundsValues.begin(), AnnulusBounds::eSize, values.begin());

  for (size_t i = 0; i < nSegments; i++) {
    values[AnnulusBounds::eMinR] = splitRanges[i].first;
    values[AnnulusBounds::eMaxR] = splitRanges[i].second;
    auto annulusBounds = std::make_shared<AnnulusBounds>(values);

    auto element = std::make_shared<const TGeoDetectorElement>(identifier, detElement->tgeoNode(), transform, annulusBounds, thickness);
    detElements.push_back(std::move(element));
  }
  return detElements;
}

