// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Implementation of BlueprintBuilder template methods.
/// Include this header only in TUs that perform explicit template instantiation
/// (e.g. geometry plugin .cpp files). Do not include from the main
/// BlueprintBuilder.hpp.

#pragma once

#include "Acts/Geometry/BlueprintBuilder.hpp"
#include "Acts/Utilities/FunctionComposition.hpp"

namespace Acts::Experimental {

namespace detail {

template <typename ElementT>
struct LayerBuildInputs {
  std::vector<ElementT> layerElements;
  std::optional<std::string> deducedContainerName;
};

template <detail::BlueprintBackend BackendT>
LayerBuildInputs<typename BackendT::Element> resolveLayerBuildInputs(
    const BlueprintBuilder<BackendT>& builder,
    const std::optional<typename BackendT::Element>& container,
    const std::optional<std::vector<typename BackendT::Element>>&
        explicitLayerElements,
    const std::optional<std::regex>& filter) {
  using Element = typename BackendT::Element;
  LayerBuildInputs<Element> inputs;

  if (container.has_value()) {
    inputs.deducedContainerName = builder.backend().nameOf(container.value());
  }

  if (explicitLayerElements.has_value()) {
    inputs.layerElements = *explicitLayerElements;

    if (filter.has_value()) {
      std::vector<Element> filteredLayerElements;
      filteredLayerElements.reserve(inputs.layerElements.size());
      for (const auto& layerElement : inputs.layerElements) {
        const std::string layerElementName =
            builder.backend().nameOf(layerElement);
        if (std::regex_match(layerElementName, filter.value())) {
          filteredLayerElements.push_back(layerElement);
        }
      }
      inputs.layerElements = std::move(filteredLayerElements);
    }
    return inputs;
  }

  if (!container.has_value()) {
    throw std::runtime_error("Container not set in ElementLayerAssembler");
  }
  if (!filter.has_value()) {
    throw std::runtime_error("Pattern not set in ElementLayerAssembler");
  }

  inputs.layerElements =
      builder.findDetElementByNamePattern(container.value(), filter.value());
  return inputs;
}

// Both ElementLayerAssembler and SensorLayerAssembler use the same underlying
// LayerCustomizer type (a std::function with identical signature), so a single
// template deducing CustomizerT covers both.
template <detail::BlueprintBackend BackendT, typename CustomizerT>
std::shared_ptr<Acts::Experimental::LayerBlueprintNode> finalizeLayer(
    const std::optional<typename BackendT::Element>& layerElement,
    std::shared_ptr<Acts::Experimental::LayerBlueprintNode> layer,
    const std::optional<Acts::ExtentEnvelope>& envelope,
    const CustomizerT& onLayer) {
  if (envelope.has_value()) {
    layer->setEnvelope(envelope.value());
  }
  if (onLayer) {
    layer = onLayer(layerElement, std::move(layer));
  }
  return layer;
}

}  // namespace detail

// ElementLayerAssembler
template <detail::BlueprintBackend BackendT>
ElementLayerAssembler<BackendT>::ElementLayerAssembler(const Builder& builder)
    : m_builder{&builder} {}

template <detail::BlueprintBackend BackendT>
ElementLayerAssembler<BackendT>&& ElementLayerAssembler<BackendT>::setLayerType(
    LayerType layerType) && {
  m_layerType = layerType;
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
ElementLayerAssembler<BackendT>&& ElementLayerAssembler<BackendT>::endcap() && {
  return std::move(*this).setLayerType(LayerType::Disc);
}

template <detail::BlueprintBackend BackendT>
ElementLayerAssembler<BackendT>&& ElementLayerAssembler<BackendT>::barrel() && {
  return std::move(*this).setLayerType(LayerType::Cylinder);
}

template <detail::BlueprintBackend BackendT>
ElementLayerAssembler<BackendT>&& ElementLayerAssembler<BackendT>::planar() && {
  return std::move(*this).setLayerType(LayerType::Plane);
}

template <detail::BlueprintBackend BackendT>
ElementLayerAssembler<BackendT>&&
ElementLayerAssembler<BackendT>::setLayerFilter(const std::string& pattern) && {
  return std::move(*this).setLayerFilter(std::regex{pattern});
}

template <detail::BlueprintBackend BackendT>
ElementLayerAssembler<BackendT>&&
ElementLayerAssembler<BackendT>::setLayerFilter(const std::regex& pattern) && {
  m_filter = pattern;
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
ElementLayerAssembler<BackendT>&& ElementLayerAssembler<BackendT>::setContainer(
    const Element& container) && {
  m_container = container;
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
ElementLayerAssembler<BackendT>&&
ElementLayerAssembler<BackendT>::setContainerName(
    std::string containerName) && {
  m_containerName = std::move(containerName);
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
ElementLayerAssembler<BackendT>&&
ElementLayerAssembler<BackendT>::setLayerNameSuffix(
    const std::optional<std::string>& layerNameSuffix) && {
  m_layerSpec.layerName = layerNameSuffix;
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
ElementLayerAssembler<BackendT>&&
ElementLayerAssembler<BackendT>::setLayerElements(
    std::vector<Element> layerElements) && {
  m_layerElements = std::move(layerElements);
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
ElementLayerAssembler<BackendT>&& ElementLayerAssembler<BackendT>::setContainer(
    const std::string& name) && {
  m_container = m_builder->findDetElementByName(name);
  if (!m_container.has_value()) {
    throw std::runtime_error("Could not find DetElement with name " + name +
                             " in ElementLayerAssembler");
  }
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
ElementLayerAssembler<BackendT>&& ElementLayerAssembler<BackendT>::setEnvelope(
    const Acts::ExtentEnvelope& envelope) && {
  m_envelope = envelope;
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
ElementLayerAssembler<BackendT>&& ElementLayerAssembler<BackendT>::setEmptyOk(
    bool emptyOk) && {
  m_emptyOk = emptyOk;
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
ElementLayerAssembler<BackendT>&&
ElementLayerAssembler<BackendT>::setAttachmentStrategy(
    std::optional<Acts::VolumeAttachmentStrategy> strategy) && {
  m_attachmentStrategy = strategy;
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
void ElementLayerAssembler<BackendT>::addTo(
    Acts::Experimental::BlueprintNode& node) const&& {
  node.addChild(build());
}

template <detail::BlueprintBackend BackendT>
std::shared_ptr<Acts::Experimental::ContainerBlueprintNode>
ElementLayerAssembler<BackendT>::build() const {
  const auto& logger = m_builder->logger();

  if (!m_layerType.has_value()) {
    throw std::runtime_error("Layer type not set in ElementLayerAssembler");
  }

  if constexpr (detail::HasAxisDefinition<BackendT>) {
    if (!m_layerSpec.axes.has_value()) {
      throw std::runtime_error("Axes not set in ElementLayerAssembler");
    }
  }

  if (!m_filter.has_value() && !m_layerElements.has_value()) {
    throw std::runtime_error(
        "Neither filter nor layer elements set in ElementLayerAssembler");
  }

  // Resolve the concrete layer-element set and optional deduced container name.
  auto inputs = detail::resolveLayerBuildInputs<BackendT>(
      *m_builder, m_container, m_layerElements, m_filter);

  // Resolve the final output container name (manual override wins).
  std::string containerName;
  if (m_containerName.has_value()) {
    containerName = m_containerName.value();
  } else if (inputs.deducedContainerName.has_value()) {
    containerName = inputs.deducedContainerName.value();
  } else {
    throw std::runtime_error(
        "Container name is not set in ElementLayerAssembler. Provide "
        "setContainerName() or setContainer().");
  }

  if (inputs.layerElements.empty()) {
    ACTS_LOG(m_emptyOk ? Acts::Logging::INFO : Acts::Logging::ERROR,
             "No layers found in container " << containerName
                                             << " matching pattern");
    if (!m_emptyOk) {
      throw std::runtime_error(std::format(
          "No layers found in container {} matching pattern", containerName));
    }
  }

  std::shared_ptr<Acts::Experimental::ContainerBlueprintNode> node;
  if (m_layerType != LayerType::Plane) {
    const Acts::AxisDirection axisDir = m_layerType == LayerType::Cylinder
                                            ? Acts::AxisDirection::AxisR
                                            : Acts::AxisDirection::AxisZ;
    node = std::make_shared<Acts::Experimental::CylinderContainerBlueprintNode>(
        containerName, axisDir);
  } else {
    node = std::make_shared<Acts::Experimental::CuboidContainerBlueprintNode>(
        containerName, Acts::AxisDirection::AxisZ);
  }

  if (m_attachmentStrategy.has_value()) {
    node->setAttachmentStrategy(m_attachmentStrategy.value());
  }

  for (const auto& layerElement : inputs.layerElements) {
    LayerSpec resolvedLayerSpec = m_layerSpec;
    std::string fullLayerName = m_builder->getPathToElementName(layerElement);
    if (resolvedLayerSpec.layerName.has_value() &&
        !resolvedLayerSpec.layerName->empty()) {
      fullLayerName += "|" + *resolvedLayerSpec.layerName;
    }
    resolvedLayerSpec.layerName = std::move(fullLayerName);
    auto layer = m_builder->makeLayer(layerElement, resolvedLayerSpec);
    layer->setLayerType(m_layerType.value());
    node->addChild(detail::finalizeLayer<BackendT>(
        std::optional<Element>{layerElement}, std::move(layer), m_envelope,
        m_onLayer));
  }

  return node;
}

// SensorLayerAssembler
template <detail::BlueprintBackend BackendT>
SensorLayerAssembler<BackendT>::SensorLayerAssembler(const Builder& builder)
    : m_builder{&builder} {}

template <detail::BlueprintBackend BackendT>
SensorLayerAssembler<BackendT>&& SensorLayerAssembler<BackendT>::setLayerType(
    LayerType layerType) && {
  m_layerType = layerType;
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
SensorLayerAssembler<BackendT>&& SensorLayerAssembler<BackendT>::endcap() && {
  return std::move(*this).setLayerType(LayerType::Disc);
}

template <detail::BlueprintBackend BackendT>
SensorLayerAssembler<BackendT>&& SensorLayerAssembler<BackendT>::barrel() && {
  return std::move(*this).setLayerType(LayerType::Cylinder);
}

template <detail::BlueprintBackend BackendT>
SensorLayerAssembler<BackendT>&& SensorLayerAssembler<BackendT>::planar() && {
  return std::move(*this).setLayerType(LayerType::Plane);
}

template <detail::BlueprintBackend BackendT>
SensorLayerAssembler<BackendT>&& SensorLayerAssembler<BackendT>::setSensors(
    std::vector<Element> sensors) && {
  m_sensors = std::move(sensors);
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
SensorLayerAssembler<BackendT>&& SensorLayerAssembler<BackendT>::groupBy(
    typename SensorLayerAssembler<BackendT>::LayerGrouper grouper) && {
  m_groupBy = std::move(grouper);
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
SensorLayerAssembler<BackendT>&&
SensorLayerAssembler<BackendT>::setContainerName(std::string containerName) && {
  m_containerName = std::move(containerName);
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
SensorLayerAssembler<BackendT>&& SensorLayerAssembler<BackendT>::setEnvelope(
    const Acts::ExtentEnvelope& envelope) && {
  m_envelope = envelope;
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
SensorLayerAssembler<BackendT>&&
SensorLayerAssembler<BackendT>::setAttachmentStrategy(
    std::optional<Acts::VolumeAttachmentStrategy> strategy) && {
  m_attachmentStrategy = strategy;
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
void SensorLayerAssembler<BackendT>::addTo(
    Acts::Experimental::BlueprintNode& node) const&& {
  node.addChild(build());
}

template <detail::BlueprintBackend BackendT>
std::shared_ptr<Acts::Experimental::ContainerBlueprintNode>
SensorLayerAssembler<BackendT>::build() const {
  using enum Acts::AxisDirection;

  if (!m_layerType.has_value()) {
    throw std::runtime_error("Layer type not set in SensorLayerAssembler");
  }

  if constexpr (detail::HasAxisDefinition<BackendT>) {
    if (!m_layerSpec.axes.has_value()) {
      throw std::runtime_error(
          std::format("Axes not set in SensorLayerAssembler (backend: {})",
                      BackendT::kIdentifier));
    }
  }

  if (!m_sensors.has_value()) {
    throw std::runtime_error("Sensors not set in SensorLayerAssembler");
  }
  if (!m_containerName.has_value()) {
    throw std::runtime_error(
        "Container name not set in SensorLayerAssembler. "
        "Call setContainerName().");
  }
  if (!m_groupBy) {
    throw std::runtime_error(
        "SensorLayerAssembler requires groupBy(). For a single layer without "
        "grouping, use BlueprintBuilder::layerFromSensors() instead.");
  }

  std::shared_ptr<Acts::Experimental::ContainerBlueprintNode> node;
  if (m_layerType != LayerType::Plane) {
    const AxisDirection axisDir =
        m_layerType == LayerType::Cylinder ? AxisR : AxisZ;
    node = std::make_shared<Acts::Experimental::CylinderContainerBlueprintNode>(
        m_containerName.value(), axisDir);
  } else {
    node = std::make_shared<Acts::Experimental::CuboidContainerBlueprintNode>(
        m_containerName.value(), AxisZ);
  }

  if (m_attachmentStrategy.has_value()) {
    node->setAttachmentStrategy(m_attachmentStrategy.value());
  }

  struct GroupData {
    std::string key;
    std::vector<Element> sensors;
  };

  // Keep first-seen group order deterministic.
  std::vector<GroupData> groups;
  for (const auto& sensor : *m_sensors) {
    const std::string key = m_groupBy(sensor);
    auto it = std::ranges::find_if(
        groups, [&](const GroupData& g) { return g.key == key; });
    if (it == groups.end()) {
      groups.push_back(GroupData{.key = key, .sensors = {}});
      it = std::prev(groups.end());
    }
    it->sensors.push_back(sensor);
  }

  for (const auto& group : groups) {
    if (group.key.empty()) {
      throw std::runtime_error(
          "groupBy() key must be non-empty for all sensors in "
          "SensorLayerAssembler");
    }
    LayerSpec layerSpec = m_layerSpec;
    layerSpec.layerName = group.key;
    auto layer = m_builder->makeLayer(std::span<const Element>{group.sensors},
                                      layerSpec);
    layer->setLayerType(m_layerType.value());
    node->addChild(detail::finalizeLayer<BackendT>(
        std::nullopt, std::move(layer), m_envelope, m_onLayer));
  }

  return node;
}

// SensorLayer
template <detail::BlueprintBackend BackendT>
SensorLayer<BackendT>::SensorLayer(const Builder& builder)
    : m_builder{&builder} {}

template <detail::BlueprintBackend BackendT>
SensorLayer<BackendT>&& SensorLayer<BackendT>::setLayerType(
    LayerType layerType) && {
  m_layerType = layerType;
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
SensorLayer<BackendT>&& SensorLayer<BackendT>::endcap() && {
  return std::move(*this).setLayerType(LayerType::Disc);
}

template <detail::BlueprintBackend BackendT>
SensorLayer<BackendT>&& SensorLayer<BackendT>::barrel() && {
  return std::move(*this).setLayerType(LayerType::Cylinder);
}

template <detail::BlueprintBackend BackendT>
SensorLayer<BackendT>&& SensorLayer<BackendT>::planar() && {
  return std::move(*this).setLayerType(LayerType::Plane);
}

template <detail::BlueprintBackend BackendT>
SensorLayer<BackendT>&& SensorLayer<BackendT>::setSensors(
    std::vector<Element> sensors) && {
  m_sensors = std::move(sensors);
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
SensorLayer<BackendT>&& SensorLayer<BackendT>::setLayerName(
    std::string name) && {
  m_layerName = std::move(name);
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
SensorLayer<BackendT>&& SensorLayer<BackendT>::setEnvelope(
    const Acts::ExtentEnvelope& envelope) && {
  m_envelope = envelope;
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
void SensorLayer<BackendT>::addTo(
    Acts::Experimental::BlueprintNode& node) const&& {
  node.addChild(build());
}

template <detail::BlueprintBackend BackendT>
std::shared_ptr<Acts::Experimental::LayerBlueprintNode>
SensorLayer<BackendT>::build() const {
  if (!m_layerType.has_value()) {
    throw std::runtime_error("Layer type not set in SensorLayer");
  }

  if constexpr (detail::HasAxisDefinition<BackendT>) {
    if (!m_layerSpec.axes.has_value()) {
      throw std::runtime_error("Axes not set in SensorLayer");
    }
  }

  if (!m_sensors.has_value()) {
    throw std::runtime_error("Sensors not set in SensorLayer");
  }
  if (!m_layerName.has_value() || m_layerName->empty()) {
    throw std::runtime_error(
        "Layer name not set in SensorLayer. Call setLayerName().");
  }

  LayerSpec layerSpec = m_layerSpec;
  layerSpec.layerName = m_layerName;
  auto layer =
      m_builder->makeLayer(std::span<const Element>{*m_sensors}, layerSpec);
  layer->setLayerType(m_layerType.value());
  return detail::finalizeLayer<BackendT>(std::nullopt, std::move(layer),
                                         m_envelope, m_onLayer);
}

// BarrelEndcapAssembler
template <detail::BlueprintBackend BackendT>
BarrelEndcapAssembler<BackendT>::BarrelEndcapAssembler(const Builder& builder)
    : m_builder{&builder} {}

template <detail::BlueprintBackend BackendT>
void BarrelEndcapAssembler<BackendT>::addTo(
    Acts::Experimental::BlueprintNode& node) const&&
  requires(detail::HasBarrelEndcapClassifier<BackendT>)
{
  node.addChild(build());
}

template <detail::BlueprintBackend BackendT>
BarrelEndcapAssembler<BackendT>&& BarrelEndcapAssembler<BackendT>::setAssembly(
    const Element& assembly) && {
  m_assembly = assembly;
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
BarrelEndcapAssembler<BackendT>&&
BarrelEndcapAssembler<BackendT>::setSensorAxes(
    typename BarrelEndcapAssembler<BackendT>::AxisDefinition barrel,
    typename BarrelEndcapAssembler<BackendT>::AxisDefinition endcap) &&
  requires(detail::HasAxisDefinition<BackendT>)
{
  m_barrelAxes = std::move(barrel);
  m_endcapAxes = std::move(endcap);
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
BarrelEndcapAssembler<BackendT>&&
BarrelEndcapAssembler<BackendT>::setEndcapAxes(
    typename BarrelEndcapAssembler<BackendT>::AxisDefinition axes) &&
  requires(detail::HasAxisDefinition<BackendT>)
{
  m_endcapAxes = std::move(axes);
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
BarrelEndcapAssembler<BackendT>&&
BarrelEndcapAssembler<BackendT>::setLayerFilter(const std::regex& pattern) && {
  m_layerFilter = pattern;
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
std::shared_ptr<Acts::Experimental::CylinderContainerBlueprintNode>
BarrelEndcapAssembler<BackendT>::build() const
  requires(detail::HasBarrelEndcapClassifier<BackendT>)
{
  using enum Acts::AxisDirection;

  const auto& logger = m_builder->logger();

  if (!m_assembly.has_value()) {
    throw std::runtime_error(
        "Assembly detector element not set in BarrelEndcapAssembler");
  }

  if constexpr (detail::HasAxisDefinition<BackendT>) {
    if (!m_barrelAxes.has_value()) {
      throw std::runtime_error("Barrel axes not set in BarrelEndcapAssembler");
    }

    if (!m_endcapAxes.has_value()) {
      throw std::runtime_error("Endcap axes not set in BarrelEndcapAssembler");
    }
  }

  if (!m_layerFilter.has_value()) {
    throw std::runtime_error("Layer pattern not set in BarrelEndcapAssembler");
  }

  const auto& assembly = m_assembly.value();
  const std::string assemblyName = m_builder->backend().nameOf(assembly);

  ACTS_INFO("Converting barrel-endcap assembly from element: " << assemblyName);
  auto barrels = m_builder->findBarrelElements(assembly);

  ACTS_DEBUG("Have " << barrels.size() << " barrel elements in assembly "
                     << assemblyName);
  if (barrels.size() > 1) {
    ACTS_ERROR("Expected exactly zero or one barrel in assembly "
               << assemblyName << ", found " << barrels.size());
    throw std::runtime_error(std::format(
        "Expected exactly zero or one barrel in assembly {}", assemblyName));
  }

  auto endcaps = m_builder->findEndcapElements(assembly);

  ACTS_DEBUG("Have " << endcaps.size() << " endcap elements in assembly "
                     << assemblyName);

  auto node =
      std::make_shared<Acts::Experimental::CylinderContainerBlueprintNode>(
          assemblyName, Acts::AxisDirection::AxisZ);

  auto maybeAddAxes = [](const auto& axes) {
    return [&axes]<typename T>(T&& assembler) {
      if constexpr (detail::HasAxisDefinition<BackendT>) {
        return std::forward<T>(assembler).setSensorAxes(axes.value());
      } else {
        return std::forward<T>(assembler);
      }
    };
  };

  auto build = []<typename T>(T&& assembler) {
    return std::forward<T>(assembler).build();
  };

  auto addTo =
      std::bind_front(&Acts::Experimental::BlueprintNode::addChild, node.get());

  for (const auto& barrel : barrels) {
    auto compose = Acts::compose(addTo, std::bind_front(m_onContainer, barrel),
                                 build, maybeAddAxes(m_barrelAxes));

    compose(m_builder->layers()
                .barrel()
                .setLayerFilter(m_layerFilter.value())
                .setContainer(barrel)
                .onLayer(m_onLayer));
  }

  for (const auto& endcap : endcaps) {
    auto compose = Acts::compose(addTo, std::bind_front(m_onContainer, endcap),
                                 build, maybeAddAxes(m_endcapAxes));

    compose(m_builder->layers()
                .endcap()
                .setLayerFilter(m_layerFilter.value())
                .setContainer(endcap)
                .onLayer(m_onLayer));
  }

  return node;
}

// BlueprintBuilder
template <detail::BlueprintBackend BackendT>
BlueprintBuilder<BackendT>::BlueprintBuilder(
    const typename Backend::Config& cfg,
    std::unique_ptr<const Acts::Logger> logger_)
    : m_logger(logger_ ? std::move(logger_)
                       : Acts::getDefaultLogger("BlueprintBuilder",
                                                Acts::Logging::INFO)),
      m_backend(cfg, *m_logger) {}

template <detail::BlueprintBackend BackendT>
std::shared_ptr<Acts::Experimental::LayerBlueprintNode>
BlueprintBuilder<BackendT>::makeLayer(const Element& layerElement,
                                      const LayerSpec& layerSpec) const {
  auto sensitives = resolveSensitives(layerElement);
  return makeLayer(layerElement, sensitives, layerSpec);
}

template <detail::BlueprintBackend BackendT>
typename BlueprintBuilder<BackendT>::ElementLayerAssembler
BlueprintBuilder<BackendT>::layers() const {
  return ElementLayerAssembler(*this);
}

template <detail::BlueprintBackend BackendT>
typename BlueprintBuilder<BackendT>::SensorLayerAssembler
BlueprintBuilder<BackendT>::layersFromSensors() const {
  return SensorLayerAssembler(*this);
}

template <detail::BlueprintBackend BackendT>
typename BlueprintBuilder<BackendT>::SensorLayer
BlueprintBuilder<BackendT>::layerFromSensors() const {
  return SensorLayer(*this);
}

template <detail::BlueprintBackend BackendT>
typename BlueprintBuilder<BackendT>::BarrelEndcapAssembler
BlueprintBuilder<BackendT>::barrelEndcap() const
  requires(detail::HasBarrelEndcapClassifier<BackendT>)
{
  return BarrelEndcapAssembler(*this);
}

template <detail::BlueprintBackend BackendT>
const Acts::Logger& BlueprintBuilder<BackendT>::logger() const {
  return *m_logger;
}

template <detail::BlueprintBackend BackendT>
const typename BlueprintBuilder<BackendT>::Backend&
BlueprintBuilder<BackendT>::backend() const {
  return m_backend;
}

template <detail::BlueprintBackend BackendT>
std::optional<typename BlueprintBuilder<BackendT>::Element>
BlueprintBuilder<BackendT>::findDetElementByName(
    const Element& parent, const std::string& name) const {
  if (m_backend.nameOf(parent) == name) {
    return parent;
  }

  for (const auto& child : m_backend.children(parent)) {
    auto result = findDetElementByName(child, name);
    if (result.has_value()) {
      return result;
    }
  }

  return std::nullopt;
}

template <detail::BlueprintBackend BackendT>
std::optional<typename BlueprintBuilder<BackendT>::Element>
BlueprintBuilder<BackendT>::findDetElementByName(
    const std::string& name) const {
  return findDetElementByName(m_backend.world(), name);
}

template <detail::BlueprintBackend BackendT>
std::string BlueprintBuilder<BackendT>::getPathToElementName(
    const Element& elem, std::string_view separator) const {
  std::vector<std::string> names;
  names.emplace_back(m_backend.nameOf(elem));

  const auto world = m_backend.world();
  auto current = elem;
  while (current != world) {
    current = m_backend.parent(current);
    if (current == world) {
      break;
    }
    names.emplace_back(m_backend.nameOf(current));
  }

  std::ranges::reverse(names);
  std::string path;
  for (std::size_t i = 0; i < names.size(); ++i) {
    if (i > 0) {
      path += separator;
    }
    path += names[i];
  }
  return path;
}

template <detail::BlueprintBackend BackendT>
std::vector<typename BlueprintBuilder<BackendT>::Element>
BlueprintBuilder<BackendT>::findDetElementByNamePattern(
    const Element& parent, const std::regex& pattern) const {
  std::vector<Element> matches;

  std::function<void(const Element&)> visit = [&](const Element& elem) {
    if (const std::string elemName = m_backend.nameOf(elem);
        std::regex_match(elemName, pattern)) {
      matches.push_back(elem);
    }
    for (const auto& child : m_backend.children(elem)) {
      visit(child);
    }
  };
  visit(parent);

  return matches;
}

template <detail::BlueprintBackend BackendT>
std::vector<typename BlueprintBuilder<BackendT>::Element>
BlueprintBuilder<BackendT>::findBarrelElements(const Element& assembly) const
  requires(detail::HasBarrelEndcapClassifier<BackendT>)
{
  std::vector<Element> barrels;

  std::function<void(const Element&)> visit = [&](const Element& elem) {
    if (m_backend.isTracker(elem) && m_backend.isBarrel(elem)) {
      barrels.push_back(elem);
    }
    for (const auto& child : m_backend.children(elem)) {
      visit(child);
    }
  };
  visit(assembly);
  return barrels;
}

template <detail::BlueprintBackend BackendT>
std::vector<typename BlueprintBuilder<BackendT>::Element>
BlueprintBuilder<BackendT>::findEndcapElements(const Element& assembly) const
  requires(detail::HasBarrelEndcapClassifier<BackendT>)
{
  std::vector<Element> endcaps;

  std::function<void(const Element&)> visit = [&](const Element& elem) {
    if (m_backend.isTracker(elem) && m_backend.isEndcap(elem)) {
      endcaps.push_back(elem);
    }
    for (const auto& child : m_backend.children(elem)) {
      visit(child);
    }
  };
  visit(assembly);
  return endcaps;
}

template <detail::BlueprintBackend BackendT>
std::shared_ptr<Acts::Experimental::LayerBlueprintNode>
BlueprintBuilder<BackendT>::makeLayer(const Element& parent,
                                      std::span<const Element> sensitives,
                                      const LayerSpec& layerSpec) const {
  const std::string nodeName =
      layerSpec.layerName.value_or(m_backend.nameOf(parent));
  auto node =
      std::make_shared<Acts::Experimental::LayerBlueprintNode>(nodeName);
  node->setSurfaces(m_backend.makeSurfaces(sensitives, layerSpec));

  if constexpr (detail::HasLayerTransformLookup<BackendT>) {
    if (const auto transform =
            m_backend.lookupLayerTransform(parent, layerSpec);
        transform.has_value()) {
      node->setTransform(transform.value());
    }
  }

  return node;
}

template <detail::BlueprintBackend BackendT>
std::shared_ptr<Acts::Experimental::LayerBlueprintNode>
BlueprintBuilder<BackendT>::makeLayer(std::span<const Element> sensitives,
                                      const LayerSpec& layerSpec) const {
  if (!layerSpec.layerName.has_value() || layerSpec.layerName->empty()) {
    throw std::runtime_error(
        "BlueprintBuilder::makeLayer(sensitives, layerSpec): "
        "layerSpec.layerName must be set");
  }

  auto node = std::make_shared<Acts::Experimental::LayerBlueprintNode>(
      layerSpec.layerName.value());
  node->setSurfaces(m_backend.makeSurfaces(sensitives, layerSpec));
  return node;
}

template <detail::BlueprintBackend BackendT>
std::vector<typename BlueprintBuilder<BackendT>::Element>
BlueprintBuilder<BackendT>::resolveSensitives(const Element& detElement) const {
  std::vector<Element> sensitives;

  std::function<void(const Element&)> visit = [&](const Element& elem) {
    if (m_backend.isSensitive(elem)) {
      sensitives.push_back(elem);
    }
    for (const auto& child : m_backend.children(elem)) {
      visit(child);
    }
  };
  visit(detElement);
  return sensitives;
}

}  // namespace Acts::Experimental
