// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/LayerBlueprintNode.hpp"
#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Navigation/CylinderNavigationPolicy.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <concepts>
#include <functional>
#include <memory>
#include <optional>
#include <regex>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace Acts::Experimental {

namespace detail {

/// @brief Concept requiring @p BackendT to expose a nested `LayerSpec` type.
///
/// A backend satisfying this concept must define a `LayerSpec` struct that
/// carries the configuration needed to construct a single detector layer.
template <typename BackendT>
concept HasLayerSpec = requires { typename BackendT::LayerSpec; };

/// @brief Concept requiring a backend to expose axis-definition support in its
/// `LayerSpec`.
///
/// In addition to satisfying @ref HasLayerSpec, the backend must define an
/// `AxisDefinition` type and the corresponding `LayerSpec` must have
/// - `axes` : optional axes used to orient sensitive surfaces, and
/// - `layerAxes` : optional axes used to determine the layer transform from the
///   parent detector element shape.
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

using LayerNodePtr = std::shared_ptr<LayerBlueprintNode>;
using ContainerNodePtr = std::shared_ptr<ContainerBlueprintNode>;
using SurfacePtr = std::shared_ptr<Acts::Surface>;
using SurfaceVector = std::vector<SurfacePtr>;

/// @brief Callback type that can replace or wrap a @ref LayerBlueprintNode.
///
/// Receives the source layer element (or @c std::nullopt when no element
/// context exists) and the newly created layer node, and returns the
/// (possibly replaced) node to be added to the container.
template <typename ElementT>
using LayerCustomizer = std::function<std::shared_ptr<LayerBlueprintNode>(
    const std::optional<ElementT>&, std::shared_ptr<LayerBlueprintNode>)>;

/// @brief Callback type that can replace or wrap a
/// @ref CylinderContainerBlueprintNode.
template <typename ElementT>
using ContainerCustomizer =
    std::function<std::shared_ptr<ContainerBlueprintNode>(
        const ElementT&, std::shared_ptr<ContainerBlueprintNode>)>;

/// @brief Concept satisfied when @p CallableT can be called with an optional
/// element and a @ref LayerBlueprintNode shared pointer and returns a (possibly
/// different) @ref LayerBlueprintNode shared pointer.
///
/// Used to constrain the returning form of the `onLayer` callback accepted by
/// the assembler builders.
template <typename ElementT, typename CallableT>
concept LayerNodeReturningCallable =
    std::invocable<CallableT&, const std::optional<ElementT>&, LayerNodePtr> &&
    std::same_as<std::invoke_result_t<
                     CallableT&, const std::optional<ElementT>&, LayerNodePtr>,
                 LayerNodePtr>;

/// @brief Concept satisfied when @p CallableT can be called with an optional
/// element and a mutable @ref LayerBlueprintNode reference and returns `void`.
///
/// Used to constrain the in-place (mutating) form of the `onLayer` callback.
template <typename ElementT, typename CallableT>
concept LayerNodeReplacingCallable =
    std::invocable<CallableT&, const std::optional<ElementT>&,
                   LayerBlueprintNode&> &&
    std::same_as<
        std::invoke_result_t<CallableT&, const std::optional<ElementT>&,
                             LayerBlueprintNode&>,
        void>;

/// @brief Concept satisfied when @p CallableT can be called with an element and
/// a @ref ContainerBlueprintNode shared pointer and returns a (possibly
/// different) @ref ContainerBlueprintNode shared pointer.
///
/// Used to constrain the returning form of the `onContainer` callback.
template <typename ElementT, typename CallableT>
concept ContainerNodeReturningCallable =
    std::invocable<CallableT&, const ElementT&, ContainerNodePtr> &&
    std::same_as<
        std::invoke_result_t<CallableT&, const ElementT&, ContainerNodePtr>,
        ContainerNodePtr>;

/// @brief Concept satisfied when @p CallableT can be called with an element and
/// a mutable @ref ContainerBlueprintNode reference and returns `void`.
///
/// Used to constrain the in-place (mutating) form of the `onContainer`
/// callback.
template <typename ElementT, typename CallableT>
concept ContainerNodeReplacingCallable =
    std::invocable<CallableT&, const ElementT&, ContainerBlueprintNode&> &&
    std::same_as<std::invoke_result_t<CallableT&, const ElementT&,
                                      ContainerBlueprintNode&>,
                 void>;

/// @brief Concept requiring a backend to provide a surface-construction method.
///
/// The method must accept a span of sensitive child elements plus a `LayerSpec`
/// and return the corresponding ACTS surfaces.
template <typename BackendT>
concept HasSurfaceFactory =
    HasLayerSpec<BackendT> &&
    requires(const BackendT& backend,
             std::span<const typename BackendT::Element> sensitives,
             const typename BackendT::LayerSpec& layerSpec) {
      {
        backend.makeSurfaces(sensitives, layerSpec)
      } -> std::same_as<SurfaceVector>;
    };

/// @brief Optional backend capability to extract a layer transform from one
/// context element and a layer specification.
///
/// Backends that satisfy this concept may provide automatic layer-transform
/// extraction (for example from detector-element geometry). The interface layer
/// applies this only when a context element is available.
template <typename BackendT>
concept HasLayerTransformLookup =
    HasLayerSpec<BackendT> &&
    requires(const BackendT& backend, const typename BackendT::Element& elem,
             const typename BackendT::LayerSpec& layerSpec) {
      {
        backend.lookupLayerTransform(elem, layerSpec)
      } -> std::same_as<std::optional<Acts::Transform3>>;
    };

/// @brief Concept requiring `LayerSpec` to carry an optional `layerName` field.
///
/// The `layerName` member, when set, overrides the name derived from the
/// detector element hierarchy when constructing a layer node.
template <typename BackendT>
concept HasLayerNameMember =
    HasLayerSpec<BackendT> && requires(typename BackendT::LayerSpec layerSpec,
                                       std::optional<std::string> layerName) {
      layerSpec.layerName = std::move(layerName);
      { layerSpec.layerName.has_value() } -> std::convertible_to<bool>;
    };

/// @brief Optional backend capability for barrel/endcap assembly discovery.
template <typename BackendT>
concept HasBarrelEndcapClassifier = requires(
    const BackendT& backend, const typename BackendT::Element& element) {
  { backend.isBarrel(element) } -> std::same_as<bool>;
  { backend.isEndcap(element) } -> std::same_as<bool>;
  { backend.isTracker(element) } -> std::same_as<bool>;
};

/// @brief Concept that fully constrains a geometry backend usable with
/// @ref BlueprintBuilder.
///
/// A conforming backend must:
/// - satisfy @ref HasSurfaceFactory and @ref HasLayerNameMember,
/// - be constructible from a `Config` object and an `Acts::Logger` reference,
/// - expose the element-hierarchy query interface (`world`, `nameOf`,
///   `children`, `parent`),
/// - expose `isSensitive()`, and
/// - support equality comparison between `Element` instances, and
/// - define `static constexpr std::string_view kIdentifier`.
template <typename BackendT>
concept BlueprintBackend =
    HasSurfaceFactory<BackendT> && HasLayerNameMember<BackendT> &&
    requires(const typename BackendT::Config& cfg, const Acts::Logger& logger,
             const BackendT& backend,
             const typename BackendT::Element& element) {
      BackendT{cfg, logger};
      { BackendT::kIdentifier } -> std::convertible_to<std::string_view>;
      { backend.world() } -> std::same_as<typename BackendT::Element>;
      { backend.nameOf(element) } -> std::same_as<std::string>;
      {
        backend.children(element)
      } -> std::same_as<std::vector<typename BackendT::Element>>;
      { backend.parent(element) } -> std::same_as<typename BackendT::Element>;
      { backend.isSensitive(element) } -> std::same_as<bool>;
      requires requires(const typename BackendT::Element& a,
                        const typename BackendT::Element& b) {
        { a == b } -> std::convertible_to<bool>;
      };
    };

}  // namespace detail

template <detail::BlueprintBackend BackendT>
class BlueprintBuilder;

template <detail::BlueprintBackend BackendT>
class SensorLayerAssembler;

template <detail::BlueprintBackend BackendT>
class SensorLayer;

/// @brief Fluent builder that assembles a flat collection of cylindrical or
/// disc-like detector layers from layer-representative detector elements into a
/// @ref CylinderContainerBlueprintNode.
///
/// Each supplied element represents one layer: its hierarchy path is used as
/// the layer name and (optionally) its geometry drives the layer transform.
/// Obtained from @ref BlueprintBuilder::layers().
///
/// ```cpp
/// builder.layers()
///   .barrel()
///   .setSensorAxes(myAxes)
///   .setLayerFilter(layerPattern)
///   .setContainer(containerElement)
///   .addTo(parentNode);
/// ```
///
/// @tparam BackendT Geometry backend that provides detector elements, layer
///         specifications, hierarchy traversal, sensitive-element
///         classification, and surface construction.
template <detail::BlueprintBackend BackendT>
class ElementLayerAssembler {
 public:
  /// The associated @ref BlueprintBuilder type.
  using Builder = BlueprintBuilder<BackendT>;
  /// Distinguishes barrel (Cylinder) from endcap (Disc) layer geometry.
  using LayerType = Acts::Experimental::LayerBlueprintNode::LayerType;
  /// Backend detector element handle type.
  using Element = typename BackendT::Element;
  /// Backend layer-specification type.
  using LayerSpec = typename BackendT::LayerSpec;
  /// Axis definition type, or `std::monostate` when the backend does not
  /// support axis definitions.
  using AxisDefinition =
      std::conditional_t<detail::HasAxisDefinition<BackendT>,
                         typename BackendT::AxisDefinition, std::monostate>;
  /// Callback type that can replace or wrap a @ref LayerBlueprintNode.
  using LayerCustomizer = detail::LayerCustomizer<Element>;

  /// @brief Set the layer geometry type explicitly.
  /// @param layerType `LayerType::Cylinder` for barrel, `LayerType::Disc` for
  ///                  endcap.
  /// @return `*this` (rvalue).
  [[nodiscard]] ElementLayerAssembler&& setLayerType(LayerType layerType) &&;

  /// @brief Shorthand for `setLayerType(LayerType::Disc)`.
  /// @return `*this` (rvalue).
  [[nodiscard]] ElementLayerAssembler&& endcap() &&;

  /// @brief Shorthand for `setLayerType(LayerType::Cylinder)`.
  /// @return `*this` (rvalue).
  [[nodiscard]] ElementLayerAssembler&& barrel() &&;

  /// @brief Shorthand for `setLayerType(LayerType::Plane)`.
  /// @return `*this` (rvalue).
  [[nodiscard]] ElementLayerAssembler&& planar() &&;

  /// @brief Set the axis definition used to orient sensitive surfaces.
  ///
  /// Only available when the backend defines an @ref AxisDefinition type and
  /// stores optional surface-axis information in `LayerSpec`.
  /// @param axes Axis definition forwarded to `LayerSpec::axes`.
  /// @return `*this` (rvalue).
  template <typename B = BackendT>
      [[nodiscard]] ElementLayerAssembler&& setSensorAxes(
          typename B::AxisDefinition axes) &&
      requires(detail::HasAxisDefinition<B>) {
        m_layerSpec.axes = std::move(axes);
        return std::move(*this);
      }

      /// @brief Set the axis definition used to derive the layer transform from the
      /// parent element shape.
      ///
      /// Only available when the backend defines an @ref AxisDefinition type and
      /// stores optional layer-axis information in `LayerSpec`.
      /// When set, the layer transform is extracted automatically from the
      /// geometry of the enclosing detector element.
      /// @param layerAxes Axis definition forwarded to `LayerSpec::layerAxes`.
      /// @return `*this` (rvalue).
      template <typename B = BackendT>
      [[nodiscard]] ElementLayerAssembler&& setLayerAxes(
          typename B::AxisDefinition layerAxes) &&
      requires(detail::HasAxisDefinition<B>) {
        m_layerSpec.layerAxes = std::move(layerAxes);
        return std::move(*this);
      }

      /// @brief Set the regex filter used to select layer elements inside the
      /// container by name string.
      /// @param pattern Regular-expression string; converted to `std::regex`
      ///                internally.
      /// @return `*this` (rvalue).
      [[nodiscard]] ElementLayerAssembler&& setLayerFilter(
          const std::string& pattern) &&;

  /// @brief Set the regex filter used to select layer elements inside the
  /// container.
  /// @param pattern Pre-compiled regular expression matched against each
  ///                child element name.
  /// @return `*this` (rvalue).
  [[nodiscard]] ElementLayerAssembler&& setLayerFilter(
      const std::regex& pattern) &&;

  /// @brief Set the detector element that acts as the containing volume for the
  /// layer search.
  /// @param container Element whose subtree is searched for layers matching the
  ///                  filter.
  /// @return `*this` (rvalue).
  [[nodiscard]] ElementLayerAssembler&& setContainer(
      const Element& container) &&;

  /// @brief Override the output container node name.
  ///
  /// If set, this name is used verbatim for the produced
  /// @ref CylinderContainerBlueprintNode. Otherwise, the name is taken from the
  /// configured container element (when @ref setContainer is used).
  /// @param containerName Explicit container-node name to use.
  /// @return `*this` (rvalue).
  [[nodiscard]] ElementLayerAssembler&& setContainerName(
      std::string containerName) &&;

  /// @brief Set an explicit suffix for produced layer node names.
  ///
  /// The final node name is `"<anchor path>|<suffix>"`. Use `std::nullopt` to
  /// clear a previously configured suffix.
  /// @param layerNameSuffix Optional suffix appended to each produced layer
  ///                        name.
  /// @return `*this` (rvalue).
  [[nodiscard]] ElementLayerAssembler&& setLayerNameSuffix(
      const std::optional<std::string>& layerNameSuffix) &&;

  /// @brief Set an explicit list of layer-representative elements.
  ///
  /// When set, the assembler skips subtree discovery via `setContainer()` and
  /// uses these elements directly as layer representatives. `setLayerFilter()`
  /// still applies as an optional post-filter on this list. Set either
  /// @ref setContainerName or @ref setContainer to define the output container
  /// name.
  /// @param layerElements Layer-representative elements; one layer per element.
  /// @return `*this` (rvalue).
  [[nodiscard]] ElementLayerAssembler&& setLayerElements(
      std::vector<Element> layerElements) &&;

  /// @brief Set the container element by name, searching from the world root.
  ///
  /// @throws std::runtime_error if no element with @p name is found.
  /// @param name Name of the detector element to use as the container.
  /// @return `*this` (rvalue).
  [[nodiscard]] ElementLayerAssembler&& setContainer(
      const std::string& name) &&;

  /// @brief Set an envelope to be applied to every layer node produced.
  /// @param envelope Envelope margins added around each layer's extent.
  /// @return `*this` (rvalue).
  [[nodiscard]] ElementLayerAssembler&& setEnvelope(
      const Acts::ExtentEnvelope& envelope) &&;

  /// @brief Control whether an empty layer collection is an error.
  ///
  /// When @p emptyOk is `false` (the default) and no layers are found in the
  /// container, @ref build() throws. Setting it to `true` downgrades the
  /// failure to an informational log message.
  /// @param emptyOk If `true`, silently accept an empty result.
  /// @return `*this` (rvalue).
  [[nodiscard]] ElementLayerAssembler&& setEmptyOk(bool emptyOk) &&;

  /// @brief Register a callback invoked for each created layer node.
  ///
  /// The callback may either:
  /// - return a (possibly replaced/wrapped) @ref LayerBlueprintNode, or
  /// - mutate a @ref LayerBlueprintNode in-place and return `void`.
  ///
  /// In both cases, the first argument is the source layer element.
  /// @param customizer Callback applied to each created layer node.
  /// @return `*this` (rvalue).
  template <typename CustomizerT>
  [[nodiscard]] ElementLayerAssembler&& onLayer(CustomizerT customizer) &&
    requires(
        detail::LayerNodeReturningCallable<Element,
                                           std::decay_t<CustomizerT>> ||
        detail::LayerNodeReplacingCallable<Element, std::decay_t<CustomizerT>>)
  {
    if constexpr (detail::LayerNodeReturningCallable<
                      Element, std::decay_t<CustomizerT>>) {
      m_onLayer = std::move(customizer);
    } else {
      m_onLayer = [customizer = std::move(customizer)](
                      const std::optional<Element>& layerElement,
                      std::shared_ptr<LayerBlueprintNode> layer) mutable {
        customizer(layerElement, *layer);
        return layer;
      };
    }
    return std::move(*this);
  }

  /// @brief Override the attachment strategy for the container node.
  ///
  /// When unset the backend's default strategy is used.
  /// @param strategy Optional attachment strategy; pass `std::nullopt` to
  ///                 reset to the default.
  /// @return `*this` (rvalue).
  [[nodiscard]] ElementLayerAssembler&& setAttachmentStrategy(
      std::optional<Acts::VolumeAttachmentStrategy> strategy) &&;

  /// @brief Build and return the assembled container node.
  ///
  /// Each resolved layer element becomes exactly one @ref LayerBlueprintNode.
  /// The layer name is derived from the element's full path in the hierarchy
  /// (plus an optional suffix). The layer transform is deduced from the element
  /// when `setLayerAxes()` is configured.
  ///
  /// @throws std::runtime_error if the layer type has not been set, if the
  ///         backend requires axes and none were provided, if neither a layer
  ///         filter nor explicit layer elements have been provided, if no
  ///         container name is resolvable, or if the container yields no
  ///         matching elements and @p emptyOk is `false`.
  /// @return Shared pointer to the fully assembled container node.
  [[nodiscard]] std::shared_ptr<Acts::Experimental::ContainerBlueprintNode>
  build() const;

  /// @brief Build the container node and attach it as a child of @p node.
  ///
  /// Equivalent to `node.addChild(build())`.
  /// @param node Blueprint node that will receive the built container as a
  ///             child.
  void addTo(Acts::Experimental::BlueprintNode& node) const&&;

 private:
  friend class BlueprintBuilder<BackendT>;

  /// @brief Construct an @ref ElementLayerAssembler bound to @p builder.
  /// @param builder The owning @ref BlueprintBuilder; must outlive this object.
  explicit ElementLayerAssembler(const Builder& builder);

  const Builder* m_builder = nullptr;
  std::optional<LayerType> m_layerType;
  LayerSpec m_layerSpec{};
  std::optional<std::regex> m_filter;
  std::optional<Element> m_container;
  std::optional<std::string> m_containerName;
  std::optional<std::vector<Element>> m_layerElements;
  std::optional<Acts::ExtentEnvelope> m_envelope;
  std::optional<Acts::VolumeAttachmentStrategy> m_attachmentStrategy;
  bool m_emptyOk = false;
  LayerCustomizer m_onLayer;
};

/// @brief Fluent builder that assembles multiple cylindrical or disc-like
/// detector layers directly from sensor elements into a
/// @ref CylinderContainerBlueprintNode.
///
/// Unlike @ref ElementLayerAssembler, no layer-representative element is
/// assumed to exist. Names and transforms are never deduced from the element
/// hierarchy; they come solely from @ref groupBy keys.
/// Obtained from @ref BlueprintBuilder::layersFromSensors().
///
/// A `groupBy` function is required: sensors mapped to the same key are merged
/// into one layer, and the key becomes the layer name.
///
/// For a single layer (no grouping), use @ref SensorLayer via
/// @ref BlueprintBuilder::layerFromSensors() instead.
///
/// ```cpp
/// builder.layersFromSensors()
///   .barrel()
///   .setSensorAxes(myAxes)
///   .setSensors(sensorElements)
///   .groupBy(keyExtractor)
///   .setContainerName("MyBarrel")
///   .addTo(parentNode);
/// ```
///
/// @tparam BackendT Geometry backend that provides detector elements, layer
///         specifications, hierarchy traversal, sensitive-element
///         classification, and surface construction.
template <detail::BlueprintBackend BackendT>
class SensorLayerAssembler {
 public:
  /// The associated @ref BlueprintBuilder type.
  using Builder = BlueprintBuilder<BackendT>;
  /// Distinguishes barrel (Cylinder) from endcap (Disc) layer geometry.
  using LayerType = LayerBlueprintNode::LayerType;
  /// Backend detector element handle type.
  using Element = typename BackendT::Element;
  /// Backend layer-specification type.
  using LayerSpec = typename BackendT::LayerSpec;
  /// Axis definition type, or `std::monostate` when the backend does not
  /// support axis definitions.
  using AxisDefinition =
      std::conditional_t<detail::HasAxisDefinition<BackendT>,
                         typename BackendT::AxisDefinition, std::monostate>;
  /// Callback type that can replace or wrap a @ref LayerBlueprintNode.
  using LayerCustomizer = detail::LayerCustomizer<Element>;
  /// Callable that maps a sensor element to a string group key.
  using LayerGrouper = std::function<std::string(const Element&)>;

  /// @brief Set the layer geometry type explicitly.
  /// @param layerType `LayerType::Cylinder` for barrel, `LayerType::Disc` for
  ///                  endcap.
  /// @return `*this` (rvalue).
  [[nodiscard]] SensorLayerAssembler&& setLayerType(LayerType layerType) &&;

  /// @brief Shorthand for `setLayerType(LayerType::Disc)`.
  /// @return `*this` (rvalue).
  [[nodiscard]] SensorLayerAssembler&& endcap() &&;

  /// @brief Shorthand for `setLayerType(LayerType::Cylinder)`.
  /// @return `*this` (rvalue).
  [[nodiscard]] SensorLayerAssembler&& barrel() &&;

  /// @brief Shorthand for `setLayerType(LayerType::Plane)`.
  /// @return `*this` (rvalue).
  [[nodiscard]] SensorLayerAssembler&& planar() &&;

  /// @brief Set the axis definition used to orient sensitive surfaces.
  ///
  /// Only available when the backend defines an @ref AxisDefinition type and
  /// stores optional surface-axis information in `LayerSpec`.
  /// @param axes Axis definition forwarded to `LayerSpec::axes`.
  /// @return `*this` (rvalue).
  template <typename B = BackendT>
      [[nodiscard]] SensorLayerAssembler&& setSensorAxes(
          typename B::AxisDefinition axes) &&
      requires(detail::HasAxisDefinition<B>) {
        m_layerSpec.axes = std::move(axes);
        return std::move(*this);
      }

      /// @brief Set the sensor elements to assemble into layers.
      /// @param sensors Sensor elements (leaf-level sensitives).
      /// @return `*this` (rvalue).
      [[nodiscard]] SensorLayerAssembler&& setSensors(
          std::vector<Element> sensors) &&;

  /// @brief Group sensors into layers by key (required).
  ///
  /// Sensors mapped to the same key are merged into one layer. The key becomes
  /// the layer name.
  /// @param grouper Callable `std::string(const Element& sensor)`.
  /// @return `*this` (rvalue).
  [[nodiscard]] SensorLayerAssembler&& groupBy(LayerGrouper grouper) &&;

  /// @brief Set the output container node name (required).
  /// @param containerName Name of the produced
  ///                      @ref CylinderContainerBlueprintNode.
  /// @return `*this` (rvalue).
  [[nodiscard]] SensorLayerAssembler&& setContainerName(
      std::string containerName) &&;

  /// @brief Set an envelope applied to every produced layer node.
  /// @param envelope Envelope margins added around each layer's extent.
  /// @return `*this` (rvalue).
  [[nodiscard]] SensorLayerAssembler&& setEnvelope(
      const Acts::ExtentEnvelope& envelope) &&;

  /// @brief Override the attachment strategy for the container node.
  ///
  /// When unset the backend's default strategy is used.
  /// @param strategy Optional attachment strategy; pass `std::nullopt` to
  ///                 reset to the default.
  /// @return `*this` (rvalue).
  [[nodiscard]] SensorLayerAssembler&& setAttachmentStrategy(
      std::optional<Acts::VolumeAttachmentStrategy> strategy) &&;

  /// @brief Register a callback invoked for each created layer node.
  ///
  /// The callback may either return a (possibly replaced/wrapped) layer node,
  /// or mutate a layer node in-place and return `void`.
  /// @param customizer Callback applied to each created layer node.
  /// @return `*this` (rvalue).
  template <typename CustomizerT>
  [[nodiscard]] SensorLayerAssembler&& onLayer(CustomizerT customizer) &&
    requires(
        detail::LayerNodeReturningCallable<Element,
                                           std::decay_t<CustomizerT>> ||
        detail::LayerNodeReplacingCallable<Element, std::decay_t<CustomizerT>>)
  {
    if constexpr (detail::LayerNodeReturningCallable<
                      Element, std::decay_t<CustomizerT>>) {
      m_onLayer = std::move(customizer);
    } else {
      m_onLayer = [customizer = std::move(customizer)](
                      const std::optional<Element>& elem,
                      std::shared_ptr<LayerBlueprintNode> layer) mutable {
        customizer(elem, *layer);
        return layer;
      };
    }
    return std::move(*this);
  }

  /// @brief Build and return the assembled container node.
  ///
  /// @throws std::runtime_error if the layer type is not set, if the backend
  ///         requires axes and none were provided, if sensors are not set,
  ///         if
  ///         the container name is not set, or if @ref groupBy has not been
  ///         configured.
  /// @return Shared pointer to the assembled container node.
  [[nodiscard]] std::shared_ptr<ContainerBlueprintNode> build() const;

  /// @brief Build the container node and attach it as a child of @p node.
  ///
  /// Equivalent to `node.addChild(build())`.
  /// @param node Blueprint node that will receive the built container as a
  ///             child.
  void addTo(BlueprintNode& node) const&&;

 private:
  friend class BlueprintBuilder<BackendT>;

  /// @brief Construct a @ref SensorLayerAssembler bound to @p builder.
  /// @param builder The owning @ref BlueprintBuilder; must outlive this object.
  explicit SensorLayerAssembler(const Builder& builder);

  const Builder* m_builder = nullptr;
  std::optional<LayerType> m_layerType;
  LayerSpec m_layerSpec{};
  std::optional<std::vector<Element>> m_sensors;
  LayerGrouper m_groupBy;
  std::optional<std::string> m_containerName;
  std::optional<Acts::ExtentEnvelope> m_envelope;
  std::optional<Acts::VolumeAttachmentStrategy> m_attachmentStrategy;
  LayerCustomizer m_onLayer;
};

/// @brief Fluent builder that assembles a single cylindrical or disc-like
/// detector layer directly from sensor elements, returning a
/// @ref LayerBlueprintNode (no container wrapper).
///
/// Unlike @ref SensorLayerAssembler, this builder produces exactly one layer.
/// The layer name must be provided explicitly via @ref setLayerName. No
/// grouping function is required or supported.
/// Obtained from @ref BlueprintBuilder::layerFromSensors().
///
/// ```cpp
/// builder.layerFromSensors()
///   .barrel()
///   .setSensorAxes(myAxes)
///   .setSensors(sensorElements)
///   .setLayerName("MyLayer")
///   .addTo(parentNode);
/// ```
///
/// @tparam BackendT Geometry backend that provides detector elements, layer
///         specifications, hierarchy traversal, sensitive-element
///         classification, and surface construction.
template <detail::BlueprintBackend BackendT>
class SensorLayer {
 public:
  /// The associated @ref BlueprintBuilder type.
  using Builder = BlueprintBuilder<BackendT>;
  /// Distinguishes barrel (Cylinder) from endcap (Disc) layer geometry.
  using LayerType = LayerBlueprintNode::LayerType;
  /// Backend detector element handle type.
  using Element = typename BackendT::Element;
  /// Backend layer-specification type.
  using LayerSpec = typename BackendT::LayerSpec;
  /// Axis definition type, or `std::monostate` when the backend does not
  /// support axis definitions.
  using AxisDefinition =
      std::conditional_t<detail::HasAxisDefinition<BackendT>,
                         typename BackendT::AxisDefinition, std::monostate>;
  /// Callback type that can replace or wrap a @ref LayerBlueprintNode.
  using LayerCustomizer = detail::LayerCustomizer<Element>;

  /// @brief Set the layer geometry type explicitly.
  /// @param layerType `LayerType::Cylinder` for barrel, `LayerType::Disc` for
  ///                  endcap.
  /// @return `*this` (rvalue).
  [[nodiscard]] SensorLayer&& setLayerType(LayerType layerType) &&;

  /// @brief Shorthand for `setLayerType(LayerType::Disc)`.
  /// @return `*this` (rvalue).
  [[nodiscard]] SensorLayer&& endcap() &&;

  /// @brief Shorthand for `setLayerType(LayerType::Cylinder)`.
  /// @return `*this` (rvalue).
  [[nodiscard]] SensorLayer&& barrel() &&;

  /// @brief Shorthand for `setLayerType(LayerType::Plane)`.
  /// @return `*this` (rvalue).
  [[nodiscard]] SensorLayer&& planar() &&;

  /// @brief Set the axis definition used to orient sensitive surfaces.
  ///
  /// Only available when the backend defines an @ref AxisDefinition type and
  /// stores optional surface-axis information in `LayerSpec`.
  /// @param axes Axis definition forwarded to `LayerSpec::axes`.
  /// @return `*this` (rvalue).
  template <typename B = BackendT>
      [[nodiscard]] SensorLayer&& setSensorAxes(
          typename B::AxisDefinition axes) &&
      requires(detail::HasAxisDefinition<B>) {
        m_layerSpec.axes = std::move(axes);
        return std::move(*this);
      }

      /// @brief Set the sensor elements to assemble into the layer.
      /// @param sensors Sensor elements (leaf-level sensitives).
      /// @return `*this` (rvalue).
      [[nodiscard]] SensorLayer&& setSensors(std::vector<Element> sensors) &&;

  /// @brief Set the name for the produced layer node (required).
  /// @param name Layer node name.
  /// @return `*this` (rvalue).
  [[nodiscard]] SensorLayer&& setLayerName(std::string name) &&;

  /// @brief Set an envelope applied to the produced layer node.
  /// @param envelope Envelope margins added around the layer's extent.
  /// @return `*this` (rvalue).
  [[nodiscard]] SensorLayer&& setEnvelope(
      const Acts::ExtentEnvelope& envelope) &&;

  /// @brief Register a callback invoked for the created layer node.
  ///
  /// The callback may either return a (possibly replaced/wrapped) layer node,
  /// or mutate a layer node in-place and return `void`.
  /// @param customizer Callback applied to the created layer node.
  /// @return `*this` (rvalue).
  template <typename CustomizerT>
  [[nodiscard]] SensorLayer&& onLayer(CustomizerT customizer) &&
    requires(
        detail::LayerNodeReturningCallable<Element,
                                           std::decay_t<CustomizerT>> ||
        detail::LayerNodeReplacingCallable<Element, std::decay_t<CustomizerT>>)
  {
    if constexpr (detail::LayerNodeReturningCallable<
                      Element, std::decay_t<CustomizerT>>) {
      m_onLayer = std::move(customizer);
    } else {
      m_onLayer = [customizer = std::move(customizer)](
                      const std::optional<Element>& elem,
                      std::shared_ptr<LayerBlueprintNode> layer) mutable {
        customizer(elem, *layer);
        return layer;
      };
    }
    return std::move(*this);
  }

  /// @brief Build and return the assembled layer node.
  ///
  /// @throws std::runtime_error if the layer type is not set, if the backend
  ///         requires axes and none were provided, if sensors are not set,
  ///         or
  ///         if @ref setLayerName has not been called.
  /// @return Shared pointer to the assembled @ref LayerBlueprintNode.
  [[nodiscard]] std::shared_ptr<LayerBlueprintNode> build() const;

  /// @brief Build the layer node and attach it as a child of @p node.
  ///
  /// Equivalent to `node.addChild(build())`.
  /// @param node Blueprint node that will receive the built layer as a child.
  void addTo(BlueprintNode& node) const&&;

 private:
  friend class BlueprintBuilder<BackendT>;

  /// @brief Construct a @ref SensorLayer bound to @p builder.
  /// @param builder The owning @ref BlueprintBuilder; must outlive this object.
  explicit SensorLayer(const Builder& builder);

  const Builder* m_builder = nullptr;
  std::optional<LayerType> m_layerType;
  LayerSpec m_layerSpec{};
  std::optional<std::vector<Element>> m_sensors;
  std::optional<std::string> m_layerName;
  std::optional<Acts::ExtentEnvelope> m_envelope;
  LayerCustomizer m_onLayer;
};

/// @brief Fluent builder that assembles a combined barrel + endcap subdetector
/// into a @ref CylinderContainerBlueprintNode arranged along the Z axis.
///
/// Instances are obtained from @ref BlueprintBuilder::barrelEndcap(). The
/// builder inspects the subtree of the provided assembly element for barrel and
/// endcap children (using the backend's `isBarrel` / `isEndcap` / `isTracker`
/// predicates, when available) and delegates individual layer assembly to
/// @ref ElementLayerAssembler internally.
///
/// Typical usage:
/// @code
/// builder.barrelEndcap()
///   .setAssembly(innerTrackerElement)
///   .setSensorAxes(barrelAxes, endcapAxes)
///   .setLayerFilter(layerPattern)
///   .addTo(rootNode);
/// @endcode
///
/// This builder requires backend predicates that classify elements as barrel,
/// endcap, and tracker components.
/// @tparam BackendT Geometry backend that provides detector elements, layer
///         specifications, hierarchy traversal, sensitive-element
///         classification, and surface construction.
template <detail::BlueprintBackend BackendT>
class BarrelEndcapAssembler {
 public:
  /// The associated @ref BlueprintBuilder type.
  using Builder = BlueprintBuilder<BackendT>;
  /// Backend detector element handle type.
  using Element = typename BackendT::Element;
  /// Axis definition type, or `std::monostate` when the backend does not
  /// support axis definitions.
  using AxisDefinition =
      std::conditional_t<detail::HasAxisDefinition<BackendT>,
                         typename BackendT::AxisDefinition, std::monostate>;
  /// The @ref ElementLayerAssembler specialisation for this backend.
  using ElementLayerAssembler =
      ::Acts::Experimental::ElementLayerAssembler<BackendT>;
  /// Callback type that can replace or wrap a
  /// @ref CylinderContainerBlueprintNode.
  using ContainerCustomizer = detail::ContainerCustomizer<Element>;

  /// @brief Construct a @ref BarrelEndcapAssembler bound to @p builder.
  /// @param builder The owning @ref BlueprintBuilder; must outlive this object.
  explicit BarrelEndcapAssembler(const Builder& builder);

  /// @brief Build and return the assembled barrel+endcap container node.
  ///
  /// Locates barrel and endcap sub-elements inside the assembly, creates one
  /// @ref ElementLayerAssembler -based barrel container and one or more endcap
  /// containers, then returns a Z-axis @ref CylinderContainerBlueprintNode
  /// holding them all.
  ///
  /// @throws std::runtime_error if the assembly element has not been set, if
  ///         axes are required by the backend but not provided, if the layer
  ///         filter has not been set, or if more than one barrel element is
  ///         found inside the assembly.
  /// @return Shared pointer to the assembled Z-axis container node.
  [[nodiscard]] std::shared_ptr<CylinderContainerBlueprintNode> build() const
    requires(detail::HasBarrelEndcapClassifier<BackendT>);

  /// @brief Build the container node and attach it as a child of @p node.
  ///
  /// Equivalent to `node.addChild(build())`.
  /// @param node Blueprint node that will receive the built container as a
  ///             child.
  void addTo(BlueprintNode& node) const&&
    requires(detail::HasBarrelEndcapClassifier<BackendT>);

  /// @brief Register a layer callback forwarded to each inner
  /// @ref ElementLayerAssembler.
  ///
  /// The callback may either return a (possibly replaced/wrapped) layer node,
  /// or mutate a layer node in-place and return `void`.
  /// @param customizer Callback applied to each created layer node.
  /// @return `*this` (rvalue).
  template <typename CustomizerT>
  [[nodiscard]] BarrelEndcapAssembler&& onLayer(CustomizerT customizer) &&
    requires(
        detail::LayerNodeReturningCallable<Element,
                                           std::decay_t<CustomizerT>> ||
        detail::LayerNodeReplacingCallable<Element, std::decay_t<CustomizerT>>)
  {
    if constexpr (detail::LayerNodeReturningCallable<
                      Element, std::decay_t<CustomizerT>>) {
      m_onLayer = std::move(customizer);
    } else {
      m_onLayer = [customizer = std::move(customizer)](
                      const std::optional<Element>& elem,
                      std::shared_ptr<LayerBlueprintNode> layer) mutable {
        customizer(elem, *layer);
        return layer;
      };
    }
    return std::move(*this);
  }

  /// @brief Register a callback invoked for each barrel or endcap container
  /// node.
  ///
  /// The callback may either return a (possibly replaced/wrapped) container
  /// node, or mutate a container node in-place and return `void`.
  /// @param customizer Callback applied to each created barrel or endcap
  ///        container node.
  /// @return `*this` (rvalue).
  template <typename CustomizerT>
  [[nodiscard]] BarrelEndcapAssembler&& onContainer(CustomizerT customizer) &&
    requires(detail::ContainerNodeReturningCallable<
                 Element, std::decay_t<CustomizerT>> ||
             detail::ContainerNodeReplacingCallable<Element,
                                                    std::decay_t<CustomizerT>>)
  {
    if constexpr (detail::ContainerNodeReturningCallable<
                      Element, std::decay_t<CustomizerT>>) {
      m_onContainer = std::move(customizer);
    } else {
      m_onContainer =
          [customizer = std::move(customizer)](
              const Element& elem,
              std::shared_ptr<ContainerBlueprintNode> node) mutable {
            customizer(elem, *node);
            return node;
          };
    }
    return std::move(*this);
  }

  /// @brief Set the top-level detector element whose subtree is searched for
  /// barrel and endcap elements.
  /// @param assembly Root element of the barrel+endcap sub-detector.
  /// @return `*this` (rvalue).
  [[nodiscard]] BarrelEndcapAssembler&& setAssembly(const Element& assembly) &&;

  /// @brief Set the axis definitions for both barrel and endcap layers at once.
  ///
  /// Only available when the backend defines an @ref AxisDefinition type and
  /// stores optional surface-axis information in `LayerSpec`.
  /// @param barrel Axis definition forwarded to barrel
  ///               @ref ElementLayerAssembler s.
  /// @param endcap Axis definition forwarded to endcap
  ///               @ref ElementLayerAssembler s.
  /// @return `*this` (rvalue).
  [[nodiscard]] BarrelEndcapAssembler&& setSensorAxes(AxisDefinition barrel,
                                                      AxisDefinition endcap) &&
    requires(detail::HasAxisDefinition<BackendT>);

  /// @brief Set the axis definition used for endcap layers only.
  ///
  /// Only available when the backend defines an @ref AxisDefinition type and
  /// stores optional surface-axis information in `LayerSpec`.
  /// @param axes Axis definition forwarded to endcap
  ///             @ref ElementLayerAssembler s.
  /// @return `*this` (rvalue).
  [[nodiscard]] BarrelEndcapAssembler&& setEndcapAxes(AxisDefinition axes) &&
    requires(detail::HasAxisDefinition<BackendT>);

  /// @brief Set the regex filter used to select individual layer elements
  /// within each barrel or endcap container.
  /// @param pattern Regular expression matched against child element names.
  /// @return `*this` (rvalue).
  [[nodiscard]] BarrelEndcapAssembler&& setLayerFilter(
      const std::regex& pattern) &&;

 private:
  typename ElementLayerAssembler::LayerCustomizer m_onLayer;
  ContainerCustomizer m_onContainer =
      [](const Element&, std::shared_ptr<ContainerBlueprintNode> node) {
        return node;
      };

  std::optional<Element> m_assembly;
  std::optional<AxisDefinition> m_barrelAxes;
  std::optional<AxisDefinition> m_endcapAxes;
  std::optional<std::regex> m_layerFilter;
  const Builder* m_builder = nullptr;
};

/// @brief High-level builder that converts a backend detector element hierarchy
/// into a blueprint node tree.
///
/// @ref BlueprintBuilder provides the entry points for blueprint construction:
/// - @ref BlueprintBuilder::layers() returns an @ref ElementLayerAssembler
///   for building a layer stack from layer-representative detector elements,
/// - @ref BlueprintBuilder::layersFromSensors() returns a
///   @ref SensorLayerAssembler for building multiple layers directly from
///   sensor elements (`groupBy` required),
/// - @ref BlueprintBuilder::layerFromSensors() returns a @ref SensorLayer for
///   building a single layer directly from sensor elements,
/// - @ref BlueprintBuilder::barrelEndcap() returns a
///   @ref BarrelEndcapAssembler for combined barrel+endcap sub-detectors.
///
/// It also exposes helpers for traversing and querying the detector element
/// hierarchy, which are used internally by the assembler classes.
///
/// The builder is parameterised on a backend type @p BackendT. The backend
/// must be constructible from its `Config` plus an `Acts::Logger`, expose
/// `Element` and `LayerSpec` types, traverse the detector hierarchy via
/// `world()` / `children()` / `parent()` / `nameOf()`, classify sensitive
/// elements with `isSensitive()`, and build ACTS surfaces from sensitive
/// elements and a layer specification. It encapsulates all
/// geometry-framework-specific knowledge (e.g. DD4hep `DetElement`
/// navigation, type-flag interpretation, surface conversion).
///
/// @tparam BackendT Geometry backend that provides detector elements, layer
///         specifications, hierarchy traversal, sensitive-element
///         classification, and surface construction.
template <detail::BlueprintBackend BackendT>
class BlueprintBuilder {
 public:
  /// The backend type.
  using Backend = BackendT;
  /// Backend detector element handle type.
  using Element = typename Backend::Element;
  /// Backend layer-specification type.
  using LayerSpec = typename Backend::LayerSpec;
  /// Axis definition type, or `std::monostate` when the backend does not
  /// support axis definitions.
  using AxisDefinition =
      std::conditional_t<detail::HasAxisDefinition<Backend>,
                         typename Backend::AxisDefinition, std::monostate>;
  /// The @ref ElementLayerAssembler specialisation for this backend.
  using ElementLayerAssembler =
      ::Acts::Experimental::ElementLayerAssembler<Backend>;
  /// The @ref SensorLayerAssembler specialisation for this backend.
  using SensorLayerAssembler =
      ::Acts::Experimental::SensorLayerAssembler<Backend>;
  /// The @ref SensorLayer specialisation for this backend.
  using SensorLayer = ::Acts::Experimental::SensorLayer<Backend>;
  /// The @ref BarrelEndcapAssembler specialisation for this backend.
  using BarrelEndcapAssembler =
      ::Acts::Experimental::BarrelEndcapAssembler<Backend>;

  /// @brief Construct a `BlueprintBuilder` from a backend configuration.
  ///
  /// @param cfg  Backend-specific configuration object passed directly to the
  ///             backend constructor.
  /// @param logger_ Optional logger; defaults to an `INFO`-level logger named
  ///                `"BlueprintBuilder"`.
  explicit BlueprintBuilder(const typename Backend::Config& cfg,
                            std::unique_ptr<const Acts::Logger> logger_ =
                                Acts::getDefaultLogger("BlueprintBuilder",
                                                       Acts::Logging::INFO));

  /// @brief Create a @ref LayerBlueprintNode from a single detector element.
  ///
  /// Recursively collects all sensitive descendants of @p layerElement and
  /// delegates to the two-argument overload.
  /// @param layerElement Layer element whose subtree is scanned for
  ///                   sensitive volumes.
  /// @param layerSpec  Specification controlling surface axes and optional
  ///                   name/transform overrides.
  /// @return Shared pointer to the constructed @ref LayerBlueprintNode.
  std::shared_ptr<LayerBlueprintNode> makeLayer(
      const Element& layerElement, const LayerSpec& layerSpec) const;

  /// @brief Create a @ref LayerBlueprintNode from an explicit list of sensitive
  /// elements.
  ///
  /// This is a low-level forwarding API: @p layerSpec is passed to the backend
  /// as-is.
  /// @param parent     Detector element that contextualises the layer
  ///                   (used for naming and transform extraction).
  /// @param sensitives Span of sensitive detector elements to be converted to
  ///                   surfaces and assigned to the node.
  /// @param layerSpec  Specification controlling surface axes and optional
  ///                   name/transform overrides.
  /// @return Shared pointer to the constructed @ref LayerBlueprintNode.
  std::shared_ptr<LayerBlueprintNode> makeLayer(
      const Element& parent, std::span<const Element> sensitives,
      const LayerSpec& layerSpec) const;

  /// @brief Create a @ref LayerBlueprintNode from sensitive elements without
  /// a parent/context element.
  ///
  /// This variant cannot perform parent-based transform extraction.
  /// `layerSpec.layerName` must be set and is used verbatim.
  /// @param sensitives Span of sensitive detector elements to be converted to
  ///                   surfaces and assigned to the node.
  /// @param layerSpec  Specification controlling surface axes and name.
  /// @throws std::runtime_error if `layerSpec.layerName` is not set.
  /// @return Shared pointer to the constructed @ref LayerBlueprintNode.
  std::shared_ptr<LayerBlueprintNode> makeLayer(
      std::span<const Element> sensitives, const LayerSpec& layerSpec) const;

  /// @brief Create an @ref ElementLayerAssembler bound to this builder.
  ///
  /// Use when layer-representative detector elements exist in the hierarchy.
  /// Each element becomes one layer; name and transform are deduced from it.
  /// @return A new @ref ElementLayerAssembler instance.
  [[nodiscard]] ElementLayerAssembler layers() const;

  /// @brief Create a @ref SensorLayerAssembler bound to this builder.
  ///
  /// Use when sensor elements are supplied directly and no layer-representative
  /// element exists. A @ref SensorLayerAssembler::groupBy function is required;
  /// sensors sharing the same key are merged into one layer per key. For a
  /// single layer without grouping, use @ref layerFromSensors() instead.
  /// @return A new @ref SensorLayerAssembler instance.
  [[nodiscard]] SensorLayerAssembler layersFromSensors() const;

  /// @brief Create a @ref SensorLayer bound to this builder.
  ///
  /// Use when all sensor elements belong to exactly one layer and no grouping
  /// is needed. The layer name must be set explicitly via
  /// @ref SensorLayer::setLayerName. The result is a single
  /// @ref LayerBlueprintNode (no container wrapper).
  /// @return A new @ref SensorLayer instance.
  [[nodiscard]] SensorLayer layerFromSensors() const;

  /// @brief Create a @ref BarrelEndcapAssembler bound to this builder.
  ///
  /// The returned assembler must be configured (assembly element, axes, layer
  /// filter) and then finalised via @ref BarrelEndcapAssembler::build() or
  /// @ref BarrelEndcapAssembler::addTo().
  /// @return A new @ref BarrelEndcapAssembler instance.
  [[nodiscard]] BarrelEndcapAssembler barrelEndcap() const
    requires(detail::HasBarrelEndcapClassifier<Backend>);

  /// @brief Search for a detector element by exact name within a subtree.
  ///
  /// Performs a depth-first search starting at @p parent.
  /// @param parent Starting element for the search.
  /// @param name   Exact element name to match.
  /// @return The first matching element, or `std::nullopt` if not found.
  std::optional<Element> findDetElementByName(const Element& parent,
                                              const std::string& name) const;

  /// @brief Search for a detector element by exact name starting from the world
  /// root.
  /// @param name Exact element name to match.
  /// @return The first matching element, or `std::nullopt` if not found.
  std::optional<Element> findDetElementByName(const std::string& name) const;

  /// @brief Build a separator-joined path string from the world root down to
  /// @p elem.
  ///
  /// The path consists of element names at each level of the hierarchy, from
  /// the immediate child of the world element down to @p elem, joined by
  /// @p separator.
  /// @param elem      Target element.
  /// @param separator String inserted between successive name components
  ///                  (default: `"|"`).
  /// @return The assembled path string.
  std::string getPathToElementName(const Element& elem,
                                   std::string_view separator = "|") const;

  /// @brief Collect all elements in a subtree whose names match a regex.
  ///
  /// Performs a depth-first traversal of the subtree rooted at @p parent and
  /// returns every element whose name fully matches @p pattern
  /// (`std::regex_match`).
  /// @param parent  Root of the subtree to search.
  /// @param pattern Regular expression matched against each element name.
  /// @return Vector of matching elements in depth-first order.
  std::vector<Element> findDetElementByNamePattern(
      const Element& parent, const std::regex& pattern) const;

  /// @brief Collect all barrel tracker elements within an assembly subtree.
  ///
  /// An element is included when the backend reports it as both a tracker
  /// element and a barrel element.
  /// @param assembly Root element of the assembly subtree to search.
  /// @return Vector of barrel elements in depth-first order.
  std::vector<Element> findBarrelElements(const Element& assembly) const
    requires(detail::HasBarrelEndcapClassifier<Backend>);

  /// @brief Collect all endcap tracker elements within an assembly subtree.
  ///
  /// An element is included when the backend reports it as both a tracker
  /// element and an endcap element.
  /// @param assembly Root element of the assembly subtree to search.
  /// @return Vector of endcap elements in depth-first order.
  std::vector<Element> findEndcapElements(const Element& assembly) const
    requires(detail::HasBarrelEndcapClassifier<Backend>);

  /// @brief Return the logger associated with this builder.
  /// @return Reference to the logger instance.
  const Acts::Logger& logger() const;

  /// @brief Return the backend associated with this builder.
  /// @return Const reference to the backend instance.
  const Backend& backend() const;

  /// @brief Recursively collect all sensitive descendant elements.
  ///
  /// Traverses the subtree rooted at @p detElement and returns every element
  /// for which the backend's `isSensitive()` predicate returns `true`.
  /// @param detElement Root of the subtree to scan.
  /// @return Vector of sensitive elements in depth-first order.
  std::vector<Element> resolveSensitives(const Element& detElement) const;

 private:
  std::unique_ptr<const Acts::Logger> m_logger;
  Backend m_backend;
};

}  // namespace Acts::Experimental
