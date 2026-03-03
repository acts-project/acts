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

template <typename BackendT>
concept BlueprintBackend =
    HasLayerFactory<BackendT> && HasLayerNameMember<BackendT> &&
    requires(const typename BackendT::Config& cfg, const Acts::Logger& logger,
             const BackendT& backend,
             const typename BackendT::Element& element) {
      BackendT{cfg, logger};
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
      requires requires(const typename BackendT::Element& a,
                        const typename BackendT::Element& b) {
        { a == b } -> std::convertible_to<bool>;
      };
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
  using ReplacingLayerCustomizer = std::function<void(
      const Element&, Acts::Experimental::LayerBlueprintNode&)>;

  explicit LayerAssembler(const Builder& builder);

  LayerAssembler&& setLayerType(LayerType layerType) &&;
  LayerAssembler&& endcap() &&;
  LayerAssembler&& barrel() &&;

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

  LayerAssembler&& setFilter(const std::string& pattern) &&;
  LayerAssembler&& setFilter(const std::regex& pattern) &&;

  LayerAssembler&& setContainer(const Element& container) &&;
  LayerAssembler&& setContainer(const std::string& name) &&;

  LayerAssembler&& setEnvelope(const Acts::ExtentEnvelope& envelope) &&;
  LayerAssembler&& setEmptyOk(bool emptyOk) &&;

  LayerAssembler&& onLayer(LayerCustomizer customizer) &&;
  LayerAssembler&& onLayer(ReplacingLayerCustomizer customizer) &&;

  LayerAssembler&& setAttachmentStrategy(
      std::optional<Acts::VolumeAttachmentStrategy> strategy) &&;

  std::shared_ptr<Acts::Experimental::CylinderContainerBlueprintNode> build()
      const;

  void addTo(Acts::Experimental::BlueprintNode& node) const&&;

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

// compose() and template method definitions live in
// detail/BlueprintBuilder_impl.hpp (included only for explicit instantiation)

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
  using ReplacingContainerCustomizer = std::function<void(
      const Element&, Acts::Experimental::CylinderContainerBlueprintNode&)>;

  explicit BarrelEndcapAssembler(const Builder& builder);

  std::shared_ptr<Acts::Experimental::CylinderContainerBlueprintNode> build()
      const;

  void addTo(Acts::Experimental::BlueprintNode& node) const&&;

  BarrelEndcapAssembler&& onLayer(
      typename LayerAssembler::LayerCustomizer customizer) &&;
  BarrelEndcapAssembler&& onLayer(
      typename LayerAssembler::ReplacingLayerCustomizer customizer) &&;

  BarrelEndcapAssembler&& onContainer(ContainerCustomizer customizer) &&;
  BarrelEndcapAssembler&& onContainer(
      ReplacingContainerCustomizer customizer) &&;

  BarrelEndcapAssembler&& setAssembly(const Element& assembly) &&;

  BarrelEndcapAssembler&& setAxes(AxisDefinition barrel,
                                  AxisDefinition endcap) &&
    requires(detail::HasAxisDefinition<BackendT>);

  BarrelEndcapAssembler&& setEndcapAxes(AxisDefinition axes) &&
    requires(detail::HasAxisDefinition<BackendT>);

  BarrelEndcapAssembler&& setLayerFilter(const std::regex& pattern) &&;

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
                                                       Acts::Logging::INFO));

  std::shared_ptr<Acts::Experimental::LayerBlueprintNode> makeLayer(
      const Element& detElement, const LayerSpec& layerSpec) const;

  LayerAssembler layers() const;
  BarrelEndcapAssembler barrelEndcap() const;

  std::optional<Element> findDetElementByName(const Element& parent,
                                              const std::string& name) const;
  std::optional<Element> findDetElementByName(const std::string& name) const;

  std::string getPathToElementName(const Element& elem,
                                   const std::string& separator = "|") const;

  std::vector<Element> findDetElementByNamePattern(
      const Element& parent, const std::regex& pattern) const;

  std::vector<Element> findBarrelElements(const Element& assembly) const;
  std::vector<Element> findEndcapElements(const Element& assembly) const;

  const Acts::Logger& logger() const;
  const Backend& backend() const;

 private:
  std::shared_ptr<Acts::Experimental::LayerBlueprintNode>
  makeLayerFromSensitives(const Element& parent,
                          std::span<const Element> sensitives,
                          const LayerSpec& layerSpec) const;

  static std::vector<Element> resolveSensitives(const Backend& backend,
                                                const Element& detElement);

  std::unique_ptr<const Acts::Logger> m_logger;
  Backend m_backend;
};

}  // namespace Acts::Experimental
