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

namespace Acts::Experimental {

namespace detail {

template <typename F, typename G, typename... Rest>
auto compose(F f, G g, Rest... rest) {
  auto fg = [=](auto&& x) -> decltype(auto) {
    return std::invoke(f, g(std::forward<decltype(x)>(x)));
  };
  if constexpr (sizeof...(rest) == 0) {
    return fg;
  } else {
    return compose(fg, rest...);
  }
}

}  // namespace detail

// LayerAssembler
template <detail::BlueprintBackend BackendT>
LayerAssembler<BackendT>::LayerAssembler(const Builder& builder)
    : m_builder{&builder} {}

template <detail::BlueprintBackend BackendT>
LayerAssembler<BackendT>&& LayerAssembler<BackendT>::setLayerType(
    LayerType layerType) && {
  m_layerType = layerType;
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
LayerAssembler<BackendT>&& LayerAssembler<BackendT>::endcap() && {
  return std::move(*this).setLayerType(LayerType::Disc);
}

template <detail::BlueprintBackend BackendT>
LayerAssembler<BackendT>&& LayerAssembler<BackendT>::barrel() && {
  return std::move(*this).setLayerType(LayerType::Cylinder);
}

template <detail::BlueprintBackend BackendT>
LayerAssembler<BackendT>&& LayerAssembler<BackendT>::setFilter(
    const std::string& pattern) && {
  return std::move(*this).setFilter(std::regex{pattern});
}

template <detail::BlueprintBackend BackendT>
LayerAssembler<BackendT>&& LayerAssembler<BackendT>::setFilter(
    const std::regex& pattern) && {
  m_filter = pattern;
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
LayerAssembler<BackendT>&& LayerAssembler<BackendT>::setContainer(
    const Element& container) && {
  m_container = container;
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
LayerAssembler<BackendT>&& LayerAssembler<BackendT>::setContainer(
    const std::string& name) && {
  m_container = m_builder->findDetElementByName(name);
  if (!m_container.has_value()) {
    throw std::runtime_error("Could not find DetElement with name " + name +
                             " in LayerAssembler");
  }
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
LayerAssembler<BackendT>&& LayerAssembler<BackendT>::setEnvelope(
    const Acts::ExtentEnvelope& envelope) && {
  m_envelope = envelope;
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
LayerAssembler<BackendT>&& LayerAssembler<BackendT>::setEmptyOk(
    bool emptyOk) && {
  m_emptyOk = emptyOk;
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
LayerAssembler<BackendT>&& LayerAssembler<BackendT>::setAttachmentStrategy(
    std::optional<Acts::VolumeAttachmentStrategy> strategy) && {
  m_attachmentStrategy = strategy;
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
LayerAssembler<BackendT>&& LayerAssembler<BackendT>::onLayer(
    typename LayerAssembler<BackendT>::LayerCustomizer customizer) && {
  m_onLayer = std::move(customizer);
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
LayerAssembler<BackendT>&& LayerAssembler<BackendT>::onLayer(
    typename LayerAssembler<BackendT>::ReplacingLayerCustomizer customizer) && {
  m_onLayer = [customizer = std::move(customizer)](
                  const typename LayerAssembler<BackendT>::Element& elem,
                  std::shared_ptr<Acts::Experimental::LayerBlueprintNode>
                      layer) mutable {
    customizer(elem, *layer);
    return layer;
  };
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
void LayerAssembler<BackendT>::addTo(
    Acts::Experimental::BlueprintNode& node) const&& {
  node.addChild(build());
}

template <detail::BlueprintBackend BackendT>
std::shared_ptr<Acts::Experimental::CylinderContainerBlueprintNode>
LayerAssembler<BackendT>::build() const {
  const auto& logger = m_builder->logger();

  if (!m_layerType.has_value()) {
    throw std::runtime_error("Layer type not set in LayerAssembler");
  }

  if constexpr (detail::HasAxisDefinition<BackendT>) {
    if (!m_layerSpec.axes.has_value()) {
      throw std::runtime_error("Axes not set in LayerAssembler");
    }
  }

  if (!m_filter.has_value()) {
    throw std::runtime_error("Pattern not set in LayerAssembler");
  }

  if (!m_container.has_value()) {
    throw std::runtime_error("Container not set in LayerAssembler");
  }

  Acts::AxisDirection axisDir = m_layerType == LayerType::Cylinder
                                    ? Acts::AxisDirection::AxisR
                                    : Acts::AxisDirection::AxisZ;

  const auto& container = m_container.value();
  const std::string containerName{m_builder->backend().nameOf(container)};
  auto node =
      std::make_shared<Acts::Experimental::CylinderContainerBlueprintNode>(
          containerName, axisDir);

  if (m_attachmentStrategy.has_value()) {
    node->setAttachmentStrategy(m_attachmentStrategy.value());
  }

  const auto& layerElements =
      m_builder->findDetElementByNamePattern(container, m_filter.value());

  if (layerElements.empty()) {
    ACTS_LOG(m_emptyOk ? Acts::Logging::INFO : Acts::Logging::ERROR,
             "No layers found in container " << containerName
                                             << " matching pattern");
    if (!m_emptyOk) {
      throw std::runtime_error(std::format(
          "No layers found in container {} matching pattern", containerName));
    }
  }
  for (const auto& element : layerElements) {
    auto layer = m_builder->makeLayer(element, m_layerSpec);
    if (m_envelope.has_value()) {
      layer->setEnvelope(m_envelope.value());
    }

    if (m_onLayer) {
      layer = m_onLayer(element, std::move(layer));
    }

    node->addChild(layer);
  }

  return node;
}

// BarrelEndcapAssembler
template <detail::BlueprintBackend BackendT>
BarrelEndcapAssembler<BackendT>::BarrelEndcapAssembler(const Builder& builder)
    : m_builder{&builder} {}

template <detail::BlueprintBackend BackendT>
void BarrelEndcapAssembler<BackendT>::addTo(
    Acts::Experimental::BlueprintNode& node) const&& {
  node.addChild(build());
}

template <detail::BlueprintBackend BackendT>
BarrelEndcapAssembler<BackendT>&& BarrelEndcapAssembler<BackendT>::setAssembly(
    const Element& assembly) && {
  m_assembly = assembly;
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
BarrelEndcapAssembler<BackendT>&& BarrelEndcapAssembler<BackendT>::onLayer(
    typename BarrelEndcapAssembler<BackendT>::LayerAssembler::LayerCustomizer
        customizer) && {
  m_onLayer = std::move(customizer);
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
BarrelEndcapAssembler<BackendT>&& BarrelEndcapAssembler<BackendT>::onLayer(
    typename BarrelEndcapAssembler<
        BackendT>::LayerAssembler::ReplacingLayerCustomizer customizer) && {
  m_onLayer = [customizer = std::move(customizer)](
                  const typename BarrelEndcapAssembler<BackendT>::Element& elem,
                  std::shared_ptr<Acts::Experimental::LayerBlueprintNode>
                      layer) mutable {
    customizer(elem, *layer);
    return layer;
  };
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
BarrelEndcapAssembler<BackendT>&& BarrelEndcapAssembler<BackendT>::onContainer(
    typename BarrelEndcapAssembler<BackendT>::ContainerCustomizer
        customizer) && {
  m_onContainer = std::move(customizer);
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
BarrelEndcapAssembler<BackendT>&& BarrelEndcapAssembler<BackendT>::onContainer(
    typename BarrelEndcapAssembler<BackendT>::ReplacingContainerCustomizer
        customizer) && {
  m_onContainer =
      [customizer = std::move(customizer)](
          const typename BarrelEndcapAssembler<BackendT>::Element& elem,
          std::shared_ptr<Acts::Experimental::CylinderContainerBlueprintNode>
              node) mutable {
        customizer(elem, *node);
        return node;
      };
  return std::move(*this);
}

template <detail::BlueprintBackend BackendT>
BarrelEndcapAssembler<BackendT>&& BarrelEndcapAssembler<BackendT>::setAxes(
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
BarrelEndcapAssembler<BackendT>::build() const {
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
  const std::string assemblyName{m_builder->backend().nameOf(assembly)};

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

  std::shared_ptr node =
      std::make_shared<Acts::Experimental::CylinderContainerBlueprintNode>(
          assemblyName, Acts::AxisDirection::AxisZ);

  auto maybeAddAxes = [](const auto& axes) {
    return [&axes](auto&& assembler) {
      if constexpr (detail::HasAxisDefinition<BackendT>) {
        return std::move(assembler).setAxes(axes.value());
      } else {
        return std::move(assembler);
      }
    };
  };

  auto build = [](auto&& assembler) { return assembler.build(); };

  auto addTo =
      std::bind_front(&Acts::Experimental::BlueprintNode::addChild, node.get());

  for (const auto& barrel : barrels) {
    auto compose =
        detail::compose(addTo, std::bind_front(m_onContainer, barrel), build,
                        maybeAddAxes(m_barrelAxes));

    compose(m_builder->layers()
                .barrel()
                .setFilter(m_layerFilter.value())
                .setContainer(barrel)
                .onLayer(m_onLayer));
  }

  for (const auto& endcap : endcaps) {
    auto compose =
        detail::compose(addTo, std::bind_front(m_onContainer, endcap), build,
                        maybeAddAxes(m_endcapAxes));

    compose(m_builder->layers()
                .endcap()
                .setFilter(m_layerFilter.value())
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
BlueprintBuilder<BackendT>::makeLayer(const Element& detElement,
                                      const LayerSpec& layerSpec) const {
  auto sensitives = resolveSensitives(detElement);
  return makeLayer(detElement, sensitives, layerSpec);
}

template <detail::BlueprintBackend BackendT>
typename BlueprintBuilder<BackendT>::LayerAssembler
BlueprintBuilder<BackendT>::layers() const {
  return LayerAssembler(*this);
}

template <detail::BlueprintBackend BackendT>
typename BlueprintBuilder<BackendT>::BarrelEndcapAssembler
BlueprintBuilder<BackendT>::barrelEndcap() const {
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
    const Element& elem, const std::string& separator) const {
  std::vector<std::string> names;
  names.emplace_back(std::string{m_backend.nameOf(elem)});

  const auto world = m_backend.world();
  auto current = elem;
  while (current != world) {
    current = m_backend.parent(current);
    if (current == world) {
      break;
    }
    names.emplace_back(std::string{m_backend.nameOf(current)});
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
    const std::string elemName{m_backend.nameOf(elem)};
    if (std::regex_match(elemName, pattern)) {
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
BlueprintBuilder<BackendT>::findBarrelElements(const Element& assembly) const {
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
BlueprintBuilder<BackendT>::findEndcapElements(const Element& assembly) const {
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
BlueprintBuilder<BackendT>::makeLayer(
    const Element& parent, std::span<const Element> sensitives,
    const LayerSpec& layerSpec) const {
  std::string fullLayerName = getPathToElementName(parent);
  LayerSpec resolvedLayerSpec = layerSpec;
  if (resolvedLayerSpec.layerName.has_value()) {
    fullLayerName += "|" + *resolvedLayerSpec.layerName;
  }
  resolvedLayerSpec.layerName = fullLayerName;

  return m_backend.makeLayer(parent, sensitives, resolvedLayerSpec);
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
