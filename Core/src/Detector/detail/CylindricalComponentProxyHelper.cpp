// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/detail/CylindricalComponentProxyHelper.hpp"

#include "Acts/Detector/DetectorVolumeBuilder.hpp"
#include "Acts/Detector/VolumeStructureBuilder.hpp"
#include "Acts/Utilities/BinningData.hpp"

void Acts::Experimental::detail::CylindricalComponentProxyHelper::sortChildren(
    ComponentBuilderProxy::ContainerProxy& container) {
  // Sort accordingly - in z
  if (container.binning == std::vector<BinningValue>{binZ}) {
    std::sort(container.children.begin(), container.children.end(),
              [](const auto& a, const auto& b) {
                return a->transform.translation().z() <
                       b->transform.translation().z();
              });
  } else if (container.binning == std::vector<BinningValue>{binR}) {
    // Sort - in r
    std::sort(container.children.begin(), container.children.end(),
              [](const auto& a, const auto& b) {
                return 0.5 * (a->boundaryValues[0] + a->boundaryValues[1]) <
                       0.5 * (b->boundaryValues[0] + b->boundaryValues[1]);
              });
  }
}

std::array<Acts::Vector3, 2u>
Acts::Experimental::detail::CylindricalComponentProxyHelper::extractVolumeEndsZ(
    const ComponentBuilderProxy& proxy) {
  Vector3 zAxis = proxy.transform.rotation().col(2);
  Vector3 negZ =
      proxy.transform.translation() - proxy.boundaryValues[2] * zAxis;
  Vector3 posZ =
      proxy.transform.translation() + proxy.boundaryValues[2] * zAxis;
  return {negZ, posZ};
};

std::shared_ptr<Acts::Experimental::ComponentBuilderProxy>
Acts::Experimental::detail::CylindricalComponentProxyHelper::createGapProxyZ(
    const std::string& gapName, ActsScalar innerR, ActsScalar outerR,
    const Acts::Transform3& transform, const Vector3& z0, const Vector3& z1,
    Acts::Logging::Level logLevel) {
  ActsScalar gap = (z1 - z0).norm();
  Transform3 gapTransform = Transform3::Identity();
  gapTransform.pretranslate(0.5 * (z0 + z1));
  gapTransform.rotate(transform.rotation());
  // Create the gap proxy & builder
  VolumeStructureBuilder::Config gapExtConfig{
      VolumeBounds::BoundsType::eCylinder,
      gapTransform,
      {innerR, outerR, 0.5 * gap},
      std::nullopt,
      "*** acts auto-generated gap shape ***"};
  auto externalsBuilder = std::make_shared<VolumeStructureBuilder>(
      gapExtConfig, getDefaultLogger(gapName + "_shape", logLevel));
  DetectorVolumeBuilder::Config gapConfig;
  gapConfig.externalsBuilder = externalsBuilder;
  auto componentBuilder = std::make_shared<DetectorVolumeBuilder>(
      gapConfig, getDefaultLogger(gapName, logLevel));
  return ComponentBuilderProxy::createRootProxy(
      gapName, gapTransform, VolumeBounds::eCylinder,
      {innerR, outerR, 0.5 * gap}, {}, componentBuilder);
};

std::shared_ptr<Acts::Experimental::ComponentBuilderProxy>
Acts::Experimental::detail::CylindricalComponentProxyHelper::createGapProxyR(
    const std::string& gapName, ActsScalar innerR, ActsScalar outerR,
    ActsScalar halfZ, const Acts::Transform3& transform,
    Acts::Logging::Level logLevel) {
  // Create the gap proxy & builder
  VolumeStructureBuilder::Config gapExtConfig{
      VolumeBounds::BoundsType::eCylinder,
      transform,
      {innerR, outerR, halfZ},
      std::nullopt,
      "*** acts auto-generated gap shape ***"};
  auto externalsBuilder = std::make_shared<VolumeStructureBuilder>(
      gapExtConfig, getDefaultLogger(gapName + "_shape", logLevel));
  DetectorVolumeBuilder::Config gapConfig;
  gapConfig.externalsBuilder = externalsBuilder;
  auto componentBuilder = std::make_shared<DetectorVolumeBuilder>(
      gapConfig, getDefaultLogger(gapName, logLevel));
  return ComponentBuilderProxy::createRootProxy(
      gapName, transform, VolumeBounds::eCylinder, {innerR, outerR, halfZ}, {},
      componentBuilder);
};

void Acts::Experimental::detail::CylindricalComponentProxyHelper::addGapProxies(
    ComponentBuilderProxy& proxy, Acts::Logging::Level logLevel) {
  ComponentBuilderProxy::ContainerProxy& container =
      std::get<ComponentBuilderProxy::ContainerProxy>(proxy.holder);

  // Sort them first
  sortChildren(container);
  // If  container values are present, gap filling can be attempted
  if (not proxy.boundaryValues.empty()) {
    std::vector<std::shared_ptr<ComponentBuilderProxy>> gapProxies;
    ActsScalar innerR = proxy.boundaryValues[0];
    ActsScalar outerR = proxy.boundaryValues[1];
    // Binning in Z
    if (container.binning == std::vector<BinningValue>{binZ}) {
      auto [negC, posC] = extractVolumeEndsZ(proxy);
      // Assume sorted along the local z axis
      unsigned int igap = 0;
      for (const auto& child : container.children) {
        auto [neg, pos] = extractVolumeEndsZ(*child);
        if (not neg.isApprox(negC)) {
          gapProxies.push_back(createGapProxyZ(
              proxy.name + "_gap_" + std::to_string(igap), innerR, outerR,
              proxy.transform, negC, neg, logLevel));
          ++igap;
        }
        // Set to new current negative value
        negC = pos;
      }
      // Check if a last one needs to be filled
      if (not negC.isApprox(posC)) {
        gapProxies.push_back(
            createGapProxyZ(proxy.name + "_gap_" + std::to_string(igap), innerR,
                            outerR, proxy.transform, negC, posC, logLevel));
      }
    } else if (container.binning == std::vector<BinningValue>{binR}) {
      // Binning in R
      // Assume sorted along the local z axis
      unsigned int igap = 0;
      ActsScalar lastR = innerR;
      ActsScalar halfZ = proxy.boundaryValues[2];
      for (const auto& child : container.children) {
        ActsScalar cInnerR = child->boundaryValues[0];
        ActsScalar cOuterR = child->boundaryValues[1];
        if (cInnerR > lastR) {
          gapProxies.push_back(createGapProxyR(
              proxy.name + "_gap_" + std::to_string(igap), lastR, cInnerR,
              halfZ, proxy.transform, logLevel));
          ++igap;
        }
        // Set to new outer radius
        lastR = cOuterR;
      }
      // Check if a last one needs to be filled
      if (lastR < outerR) {
        gapProxies.push_back(
            createGapProxyR(proxy.name + "_gap_" + std::to_string(igap), lastR,
                            outerR, halfZ, proxy.transform, logLevel));
      }
    }
    // Insert gaps
    container.children.insert(container.children.end(), gapProxies.begin(),
                              gapProxies.end());
    // Re-sort after adding the gap proxies
    sortChildren(container);
  }
}
