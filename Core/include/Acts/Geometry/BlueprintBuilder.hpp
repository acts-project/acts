// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/LayerBlueprintNode.hpp"
#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Geometry/VolumeResizeStrategy.hpp"
#include "Acts/Navigation/CylinderNavigationPolicy.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <concepts>
#include <format>
#include <functional>
#include <memory>
#include <optional>
#include <ranges>
#include <regex>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace Acts::Experimental {

namespace detail {

template <typename BackendT>
concept HasLayerSpec = requires { typename BackendT::LayerSpec; };

template <typename BackendT>
concept HasAxisDefinition =
    HasLayerSpec<BackendT> && requires { typename BackendT::AxisDefinition; } &&
    requires(typename BackendT::LayerSpec layerSpec,
             typename BackendT::AxisDefinition axes,
             std::optional<typename BackendT::AxisDefinition> layerAxes) {
      layerSpec.axes = std::move(axes);
      { layerSpec.axes.has_value() } -> std::convertible_to<bool>;
      layerSpec.layerAxes = std::move(axes);
      layerSpec.layerAxes = std::move(layerAxes);
    };

using LayerNodePtr = std::shared_ptr<Acts::Experimental::LayerBlueprintNode>;

template <typename BackendT>
concept HasLayerFactory =
    HasLayerSpec<BackendT> &&
    requires(const BackendT& backend, const typename BackendT::Element& parent,
             std::span<const typename BackendT::Element> sensitives,
             const typename BackendT::LayerSpec& layerSpec) {
      {
        backend.makeLayer(parent, sensitives, layerSpec)
      } -> std::same_as<LayerNodePtr>;
    };

template <typename BackendT>
concept HasLayerNameMember =
    HasLayerSpec<BackendT> && requires(typename BackendT::LayerSpec layerSpec,
                                       std::optional<std::string> layerName) {
      layerSpec.layerName = std::move(layerName);
      { layerSpec.layerName.has_value() } -> std::convertible_to<bool>;
    };

template <typename F, typename ElementT>
concept LayerCustomizerCallable =
    (std::invocable<F&, const ElementT&, LayerNodePtr> &&
     std::same_as<std::invoke_result_t<F&, const ElementT&, LayerNodePtr>,
                  LayerNodePtr>) ||
    (std::invocable<F&, const ElementT&,
                    Acts::Experimental::LayerBlueprintNode&> &&
     std::same_as<std::invoke_result_t<F&, const ElementT&,
                                       Acts::Experimental::LayerBlueprintNode&>,
                  void>);

template <typename...>
inline constexpr bool AlwaysFalse = false;

template <typename BackendT>
concept BlueprintBackend =
    HasLayerFactory<BackendT> && HasLayerNameMember<BackendT> &&
    requires(const typename BackendT::Config& cfg, const Acts::Logger& logger,
             const BackendT& backend,
             const typename BackendT::Element& element) {
      BackendT{cfg, logger};
      {
        backend.makeBeampipe()
      }
      -> std::same_as<std::shared_ptr<Acts::Experimental::StaticBlueprintNode>>;
      { backend.world() } -> std::same_as<typename BackendT::Element>;
      { backend.nameOf(element) } -> std::convertible_to<std::string_view>;
      {
        backend.children(element)
      } -> std::same_as<std::vector<typename BackendT::Element>>;
      { backend.parent(element) } -> std::same_as<typename BackendT::Element>;
      { backend.isSensitive(element) } -> std::same_as<bool>;
      { backend.isBarrel(element) } -> std::same_as<bool>;
      { backend.isEndcap(element) } -> std::same_as<bool>;
      { backend.isTracker(element) } -> std::same_as<bool>;
      { backend.isWorld(element) } -> std::same_as<bool>;
    };

}  // namespace detail

template <detail::BlueprintBackend BackendT>
class BlueprintBuilder;

template <detail::BlueprintBackend BackendT>
class LayerAssembler {
 public:
  using Builder = BlueprintBuilder<BackendT>;
  using LayerType = Acts::Experimental::LayerBlueprintNode::LayerType;
  using Element = typename BackendT::Element;
  using LayerSpec = typename BackendT::LayerSpec;
  using LayerCustomizer =
      std::function<std::shared_ptr<Acts::Experimental::LayerBlueprintNode>(
          const Element&,
          std::shared_ptr<Acts::Experimental::LayerBlueprintNode>)>;

  explicit LayerAssembler(const Builder& builder) : m_builder{&builder} {}

  LayerAssembler&& setLayerType(LayerType layerType) && {
    m_layerType = layerType;
    return std::move(*this);
  }

  LayerAssembler&& endcap() && {
    return std::move(*this).setLayerType(LayerType::Disc);
  }
  LayerAssembler&& barrel() && {
    return std::move(*this).setLayerType(LayerType::Cylinder);
  }

  template <typename B = BackendT>
  LayerAssembler&& setAxes(typename B::AxisDefinition axes) &&
    requires(detail::HasAxisDefinition<B>)
  {
    m_layerSpec.axes = std::move(axes);
    return std::move(*this);
  }

  template <typename B = BackendT>
  LayerAssembler&& setLayerAxes(typename B::AxisDefinition layerAxes) &&
    requires(detail::HasAxisDefinition<B>)
  {
    m_layerSpec.layerAxes = std::move(layerAxes);
    return std::move(*this);
  }

  LayerAssembler&& setFilter(const std::string& pattern) && {
    return std::move(*this).setFilter(std::regex{pattern});
  }

  LayerAssembler&& setFilter(const std::regex& pattern) && {
    m_filter = pattern;
    return std::move(*this);
  }

  LayerAssembler&& setContainer(const Element& container) && {
    m_container = container;
    return std::move(*this);
  }

  LayerAssembler&& setContainer(const std::string& name) && {
    m_container = m_builder->findDetElementByName(name);
    if (!m_container.has_value()) {
      throw std::runtime_error("Could not find DetElement with name " + name +
                               " in LayerAssembler");
    }
    return std::move(*this);
  }

  LayerAssembler&& setEnvelope(const Acts::ExtentEnvelope& envelope) && {
    m_envelope = envelope;
    return std::move(*this);
  }

  LayerAssembler&& setEmptyOk(bool emptyOk) && {
    m_emptyOk = emptyOk;
    return std::move(*this);
  }

  template <detail::LayerCustomizerCallable<Element> F>
  LayerAssembler&& onLayer(F&& customizer) && {
    using LayerNode = Acts::Experimental::LayerBlueprintNode;
    using LayerNodePtr = std::shared_ptr<LayerNode>;
    using Fn = std::decay_t<F>;

    if constexpr (std::invocable<Fn&, const Element&, LayerNodePtr>) {
      using ResultT = std::invoke_result_t<Fn&, const Element&, LayerNodePtr>;
      static_assert(std::same_as<ResultT, LayerNodePtr>,
                    "Layer customizer with pointer signature must return "
                    "std::shared_ptr<LayerBlueprintNode>");
      m_onLayer = std::forward<F>(customizer);
    } else if constexpr (std::invocable<Fn&, const Element&, LayerNode&>) {
      using ResultT = std::invoke_result_t<Fn&, const Element&, LayerNode&>;
      static_assert(std::same_as<ResultT, void>,
                    "Layer customizer with ref signature must return void");
      m_onLayer = [customizer = std::forward<F>(customizer)](
                      const Element& elem, LayerNodePtr layer) mutable {
        customizer(elem, *layer);
        return layer;
      };
    }
    return std::move(*this);
  }

  LayerAssembler&& setAttachmentStrategy(
      std::optional<Acts::VolumeAttachmentStrategy> strategy) && {
    m_attachmentStrategy = strategy;
    return std::move(*this);
  }

  std::shared_ptr<Acts::Experimental::CylinderContainerBlueprintNode> build()
      const {
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

  void addTo(Acts::Experimental::BlueprintNode& node) const&& {
    node.addChild(build());
  }

 private:
  const Builder* m_builder;
  std::optional<LayerType> m_layerType;
  LayerSpec m_layerSpec{};
  std::optional<std::regex> m_filter;
  std::optional<Element> m_container;
  std::optional<Acts::ExtentEnvelope> m_envelope;
  std::optional<Acts::VolumeAttachmentStrategy> m_attachmentStrategy;
  bool m_emptyOk = false;

  LayerCustomizer m_onLayer;
};

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

template <detail::BlueprintBackend BackendT>
class BarrelEndcapAssembler {
 public:
  using Builder = BlueprintBuilder<BackendT>;
  using Element = typename BackendT::Element;
  using AxisDefinition =
      std::conditional_t<detail::HasAxisDefinition<BackendT>,
                         typename BackendT::AxisDefinition, std::monostate>;
  using LayerAssembler = Acts::Experimental::LayerAssembler<BackendT>;
  using ContainerCustomizer = std::function<
      std::shared_ptr<Acts::Experimental::CylinderContainerBlueprintNode>(
          const Element&,
          std::shared_ptr<Acts::Experimental::CylinderContainerBlueprintNode>)>;

  explicit BarrelEndcapAssembler(const Builder& builder)
      : m_builder{&builder} {}

  std::shared_ptr<Acts::Experimental::CylinderContainerBlueprintNode> build()
      const {
    using enum Acts::AxisDirection;

    const auto& logger = m_builder->logger();

    if (!m_assembly.has_value()) {
      throw std::runtime_error(
          "Assembly detector element not set in BarrelEndcapAssembler");
    }

    if constexpr (detail::HasAxisDefinition<BackendT>) {
      if (!m_barrelAxes.has_value()) {
        throw std::runtime_error(
            "Barrel axes not set in BarrelEndcapAssembler");
      }

      if (!m_endcapAxes.has_value()) {
        throw std::runtime_error(
            "Endcap axes not set in BarrelEndcapAssembler");
      }
    }

    if (!m_layerFilter.has_value()) {
      throw std::runtime_error(
          "Layer pattern not set in BarrelEndcapAssembler");
    }

    const auto& assembly = m_assembly.value();
    const std::string assemblyName{m_builder->backend().nameOf(assembly)};

    ACTS_INFO(
        "Converting barrel-endcap assembly from element: " << assemblyName);
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

    // Composable helper to add axes to a layer assembler
    auto maybeAddAxes = [](const auto& axes) {
      return [&axes](auto&& assembler) {
        if constexpr (detail::HasAxisDefinition<BackendT>) {
          return std::move(assembler).setAxes(axes.value());
        } else {
          return std::move(assembler);
        }
      };
    };

    // Composable helper to build a layer assembler
    auto build = [](auto&& assembler) { return assembler.build(); };

    // Composable helper to add a node to a parent node
    auto addTo = std::bind_front(&Acts::Experimental::BlueprintNode::addChild,
                                 node.get());

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

  void addTo(Acts::Experimental::BlueprintNode& node) const&& {
    node.addChild(build());
  }

  template <typename F>
  BarrelEndcapAssembler&& onLayer(F&& customizer) &&
    requires(detail::LayerCustomizerCallable<std::decay_t<F>, Element>)
  {
    using LayerNode = Acts::Experimental::LayerBlueprintNode;
    using LayerNodePtr = std::shared_ptr<LayerNode>;
    using Fn = std::decay_t<F>;

    if constexpr (std::invocable<Fn&, const Element&, LayerNodePtr>) {
      using ResultT = std::invoke_result_t<Fn&, const Element&, LayerNodePtr>;
      static_assert(std::same_as<ResultT, LayerNodePtr>,
                    "Layer customizer with pointer signature must return "
                    "std::shared_ptr<LayerBlueprintNode>");
      m_onLayer = std::forward<F>(customizer);
    } else if constexpr (std::invocable<Fn&, const Element&, LayerNode&>) {
      using ResultT = std::invoke_result_t<Fn&, const Element&, LayerNode&>;
      static_assert(std::same_as<ResultT, void>,
                    "Layer customizer with ref signature must return void");
      m_onLayer = [customizer = std::forward<F>(customizer)](
                      const Element& elem, LayerNodePtr layer) mutable {
        customizer(elem, *layer);
        return layer;
      };
    }
    return std::move(*this);
  }

  template <typename F>
  BarrelEndcapAssembler&& onContainer(F&& customizer) && {
    using ContainerNode = Acts::Experimental::CylinderContainerBlueprintNode;
    using ContainerNodePtr = std::shared_ptr<ContainerNode>;
    using Fn = std::decay_t<F>;

    if constexpr (std::invocable<Fn&, const Element&, ContainerNodePtr>) {
      using ResultT =
          std::invoke_result_t<Fn&, const Element&, ContainerNodePtr>;
      static_assert(std::same_as<ResultT, ContainerNodePtr>,
                    "Container customizer with pointer signature must return "
                    "std::shared_ptr<CylinderContainerBlueprintNode>");
      m_onContainer = std::forward<F>(customizer);
    } else if constexpr (std::invocable<Fn&, const Element&, ContainerNode&>) {
      using ResultT = std::invoke_result_t<Fn&, const Element&, ContainerNode&>;
      static_assert(std::same_as<ResultT, void>,
                    "Container customizer with ref signature must return void");
      m_onContainer = [customizer = std::forward<F>(customizer)](
                          const Element& elem, ContainerNodePtr node) mutable {
        customizer(elem, *node);
        return node;
      };
    } else {
      static_assert(
          detail::AlwaysFalse<Fn>,
          "Container customizer must be invocable as either "
          "(const Element&, "
          "std::shared_ptr<CylinderContainerBlueprintNode>) -> "
          "std::shared_ptr<CylinderContainerBlueprintNode> or "
          "(const Element&, CylinderContainerBlueprintNode&) -> void");
    }
    return std::move(*this);
  }

  BarrelEndcapAssembler&& setAssembly(const Element& assembly) && {
    m_assembly = assembly;
    return std::move(*this);
  }

  BarrelEndcapAssembler&& setAxes(AxisDefinition barrel,
                                  AxisDefinition endcap) &&
    requires(detail::HasAxisDefinition<BackendT>)
  {
    m_barrelAxes = std::move(barrel);
    m_endcapAxes = std::move(endcap);
    return std::move(*this);
  }

  BarrelEndcapAssembler&& setEndcapAxes(AxisDefinition axes) &&
    requires(detail::HasAxisDefinition<BackendT>)
  {
    m_endcapAxes = std::move(axes);
    return std::move(*this);
  }

  BarrelEndcapAssembler&& setLayerFilter(const std::regex& pattern) && {
    m_layerFilter = pattern;
    return std::move(*this);
  }

 private:
  typename LayerAssembler::LayerCustomizer m_onLayer;
  ContainerCustomizer m_onContainer =
      [](const Element&,
         std::shared_ptr<Acts::Experimental::CylinderContainerBlueprintNode>
             node) { return node; };

  std::optional<Element> m_assembly;
  std::optional<AxisDefinition> m_barrelAxes;
  std::optional<AxisDefinition> m_endcapAxes;
  std::optional<std::regex> m_layerFilter;
  const Builder* m_builder;
};

template <detail::BlueprintBackend BackendT>
class BlueprintBuilder {
 public:
  using Backend = BackendT;
  using Element = typename Backend::Element;
  using LayerSpec = typename Backend::LayerSpec;
  using AxisDefinition =
      std::conditional_t<detail::HasAxisDefinition<Backend>,
                         typename Backend::AxisDefinition, std::monostate>;
  using LayerAssembler = Acts::Experimental::LayerAssembler<Backend>;
  using BarrelEndcapAssembler =
      Acts::Experimental::BarrelEndcapAssembler<Backend>;

  explicit BlueprintBuilder(const typename Backend::Config& cfg,
                            std::unique_ptr<const Acts::Logger> logger_ =
                                Acts::getDefaultLogger("BlueprintBuilder",
                                                       Acts::Logging::INFO))
      : m_logger{std::move(logger_)}, m_backend(cfg, *m_logger) {}

  std::shared_ptr<Acts::Experimental::LayerBlueprintNode> makeLayer(
      const Element& detElement, const LayerSpec& layerSpec) const {
    return makeLayerFromSensitives(
        detElement, resolveSensitives(m_backend, detElement), layerSpec);
  }

  std::shared_ptr<Acts::Experimental::StaticBlueprintNode> makeBeampipe()
      const {
    return m_backend.makeBeampipe();
  }

  LayerAssembler layers() const { return LayerAssembler(*this); }
  BarrelEndcapAssembler barrelEndcap() const {
    return BarrelEndcapAssembler(*this);
  }

  std::optional<Element> findDetElementByName(const Element& parent,
                                              const std::string& name) const {
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

  std::optional<Element> findDetElementByName(const std::string& name) const {
    return findDetElementByName(m_backend.world(), name);
  }

  std::string getPathToElementName(const Element& elem,
                                   const std::string& separator = "|") const {
    std::vector<std::string> names;
    names.emplace_back(std::string{m_backend.nameOf(elem)});

    auto current = elem;
    while (!m_backend.isWorld(current)) {
      current = m_backend.parent(current);
      if (m_backend.isWorld(current)) {
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

  std::vector<Element> findDetElementByNamePattern(
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

  std::vector<Element> findBarrelElements(const Element& assembly) const {
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

  std::vector<Element> findEndcapElements(const Element& assembly) const {
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

  const Acts::Logger& logger() const { return *m_logger; }
  const Backend& backend() const { return m_backend; }

 private:
  std::shared_ptr<Acts::Experimental::LayerBlueprintNode>
  makeLayerFromSensitives(const Element& parent,
                          std::span<const Element> sensitives,
                          const LayerSpec& layerSpec) const {
    std::string fullLayerName = getPathToElementName(parent);
    LayerSpec resolvedLayerSpec = layerSpec;
    if (resolvedLayerSpec.layerName.has_value()) {
      fullLayerName += "|" + *resolvedLayerSpec.layerName;
    }
    resolvedLayerSpec.layerName = fullLayerName;

    return m_backend.makeLayer(parent, sensitives, resolvedLayerSpec);
  }

  static std::vector<Element> resolveSensitives(const Backend& backend,
                                                const Element& detElement) {
    std::vector<Element> sensitives;

    std::function<void(const Element&)> visit = [&](const Element& elem) {
      if (backend.isSensitive(elem)) {
        sensitives.push_back(elem);
      }
      for (const auto& child : backend.children(elem)) {
        visit(child);
      }
    };
    visit(detElement);
    return sensitives;
  }

  std::unique_ptr<const Acts::Logger> m_logger;
  Backend m_backend;
};

}  // namespace Acts::Experimental
