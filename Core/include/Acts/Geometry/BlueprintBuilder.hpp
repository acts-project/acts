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

using LayerNodePtr = std::shared_ptr<Acts::Experimental::LayerBlueprintNode>;

/// @brief Concept requiring a backend to provide a `makeLayer` factory method.
///
/// The method must accept a parent element, a span of sensitive child elements,
/// and a `LayerSpec`, and return a `LayerBlueprintNode` shared pointer.
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

/// @brief Concept that fully constrains a geometry backend usable with
/// @ref BlueprintBuilder.
///
/// A conforming backend must:
/// - satisfy @ref HasLayerFactory and @ref HasLayerNameMember,
/// - be constructible from a `Config` object and an `Acts::Logger` reference,
/// - expose the element-hierarchy query interface (`world`, `nameOf`,
///   `children`, `parent`),
/// - expose element-classification predicates (`isSensitive`, `isBarrel`,
///   `isEndcap`, `isTracker`), and
/// - support equality comparison between `Element` instances.
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

/// @brief Fluent builder that assembles a flat collection of cylindrical or
/// disc-like detector layers into a @ref CylinderContainerBlueprintNode.
///
/// Instances are obtained from @ref BlueprintBuilder::layers() and configured
/// through a method-chaining API. All setters are rvalue-qualified and return
/// `LayerAssembler&&` so that calls can be chained without copies:
///
/// ```cpp
/// builder.layers()
///   .barrel()
///   .setAxes(myAxes)
///   .setFilter(layerPattern)
///   .setContainer(containerElement)
///   .addTo(parentNode);
/// ```
///
/// @tparam BackendT Geometry backend satisfying @ref detail::BlueprintBackend.
template <detail::BlueprintBackend BackendT>
class LayerAssembler {
 public:
  /// The associated @ref BlueprintBuilder type.
  using Builder = BlueprintBuilder<BackendT>;
  /// Distinguishes barrel (Cylinder) from endcap (Disc) layer geometry.
  using LayerType = Acts::Experimental::LayerBlueprintNode::LayerType;
  /// Backend detector element handle type.
  using Element = typename BackendT::Element;
  /// Backend layer-specification type.
  using LayerSpec = typename BackendT::LayerSpec;
  /// Callback type that can replace or wrap a @ref LayerBlueprintNode.
  ///
  /// Receives the source detector element and the newly created layer node,
  /// and returns the (possibly replaced) node to be added to the container.
  using LayerCustomizer =
      std::function<std::shared_ptr<Acts::Experimental::LayerBlueprintNode>(
          const Element&,
          std::shared_ptr<Acts::Experimental::LayerBlueprintNode>)>;
  /// Simplified callback type that mutates a @ref LayerBlueprintNode in-place.
  ///
  /// The return value is ignored; the original node is kept and added to the
  /// container after the callback returns.
  using ReplacingLayerCustomizer = std::function<void(
      const Element&, Acts::Experimental::LayerBlueprintNode&)>;

  /// @brief Construct a @ref LayerAssembler bound to @p builder.
  /// @param builder The owning @ref BlueprintBuilder; must outlive this object.
  explicit LayerAssembler(const Builder& builder);

  /// @brief Set the layer geometry type explicitly.
  /// @param layerType `LayerType::Cylinder` for barrel, `LayerType::Disc` for
  ///                  endcap.
  /// @return `*this` (rvalue).
  [[nodiscard]] LayerAssembler&& setLayerType(LayerType layerType) &&;

  /// @brief Shorthand for `setLayerType(LayerType::Disc)`.
  /// @return `*this` (rvalue).
  [[nodiscard]] LayerAssembler&& endcap() &&;

  /// @brief Shorthand for `setLayerType(LayerType::Cylinder)`.
  /// @return `*this` (rvalue).
  [[nodiscard]] LayerAssembler&& barrel() &&;

  /// @brief Set the axis definition used to orient sensitive surfaces.
  ///
  /// Only available when the backend satisfies @ref detail::HasAxisDefinition.
  /// @param axes Axis definition forwarded to `LayerSpec::axes`.
  /// @return `*this` (rvalue).
  template <typename B = BackendT>
  [[nodiscard]] LayerAssembler&& setAxes(typename B::AxisDefinition axes) &&
    requires(detail::HasAxisDefinition<B>)
  {
    m_layerSpec.axes = std::move(axes);
    return std::move(*this);
  }

  /// @brief Set the axis definition used to derive the layer transform from the
  /// parent element shape.
  ///
  /// Only available when the backend satisfies @ref detail::HasAxisDefinition.
  /// When set, the layer transform is extracted automatically from the geometry
  /// of the enclosing detector element.
  /// @param layerAxes Axis definition forwarded to `LayerSpec::layerAxes`.
  /// @return `*this` (rvalue).
  template <typename B = BackendT>
  [[nodiscard]] LayerAssembler&& setLayerAxes(
      typename B::AxisDefinition layerAxes) &&
    requires(detail::HasAxisDefinition<B>)
  {
    m_layerSpec.layerAxes = std::move(layerAxes);
    return std::move(*this);
  }

  /// @brief Set the regex filter used to select layer elements inside the
  /// container by name string.
  /// @param pattern Regular-expression string; converted to `std::regex`
  ///                internally.
  /// @return `*this` (rvalue).
  [[nodiscard]] LayerAssembler&& setFilter(const std::string& pattern) &&;

  /// @brief Set the regex filter used to select layer elements inside the
  /// container.
  /// @param pattern Pre-compiled regular expression matched against each
  ///                child element name.
  /// @return `*this` (rvalue).
  [[nodiscard]] LayerAssembler&& setFilter(const std::regex& pattern) &&;

  /// @brief Set the detector element that acts as the containing volume for the
  /// layer search.
  /// @param container Element whose subtree is searched for layers matching the
  ///                  filter.
  /// @return `*this` (rvalue).
  [[nodiscard]] LayerAssembler&& setContainer(const Element& container) &&;

  /// @brief Set the container element by name, searching from the world root.
  ///
  /// @throws std::runtime_error if no element with @p name is found.
  /// @param name Name of the detector element to use as the container.
  /// @return `*this` (rvalue).
  [[nodiscard]] LayerAssembler&& setContainer(const std::string& name) &&;

  /// @brief Set an envelope to be applied to every layer node produced.
  /// @param envelope Envelope margins added around each layer's extent.
  /// @return `*this` (rvalue).
  [[nodiscard]] LayerAssembler&& setEnvelope(
      const Acts::ExtentEnvelope& envelope) &&;

  /// @brief Control whether an empty layer collection is an error.
  ///
  /// When @p emptyOk is `false` (the default) and no layers are found in the
  /// container, @ref build() throws. Setting it to `true` downgrades the
  /// failure to an informational log message.
  /// @param emptyOk If `true`, silently accept an empty result.
  /// @return `*this` (rvalue).
  [[nodiscard]] LayerAssembler&& setEmptyOk(bool emptyOk) &&;

  /// @brief Register a customizer callback invoked for each created layer node.
  ///
  /// The callback receives the source detector element and the newly created
  /// @ref LayerBlueprintNode. It may return a different (or wrapped) node which
  /// will then be inserted into the container.
  /// @return `*this` (rvalue).
  [[nodiscard]] LayerAssembler&& onLayer(LayerCustomizer customizer) &&;

  /// @brief Register an in-place mutating callback invoked for each layer node.
  ///
  /// The callback receives the source detector element and a mutable reference
  /// to the @ref LayerBlueprintNode. The original node is kept regardless of
  /// what the callback does.
  /// @return `*this` (rvalue).
  [[nodiscard]] LayerAssembler&& onLayer(
      ReplacingLayerCustomizer customizer) &&;

  /// @brief Override the attachment strategy for the container node.
  ///
  /// When unset the backend's default strategy is used.
  /// @param strategy Optional attachment strategy; pass `std::nullopt` to
  ///                 reset to the default.
  /// @return `*this` (rvalue).
  [[nodiscard]] LayerAssembler&& setAttachmentStrategy(
      std::optional<Acts::VolumeAttachmentStrategy> strategy) &&;

  /// @brief Build and return the assembled container node.
  ///
  /// Searches the configured container element for children whose names match
  /// the filter, creates a @ref LayerBlueprintNode for each, applies optional
  /// envelope and customizer, then wraps them in a
  /// @ref CylinderContainerBlueprintNode oriented along the axis appropriate
  /// for the selected layer type.
  ///
  /// @throws std::runtime_error if the layer type, filter, or container has not
  ///         been set, or if the backend requires axes and none were provided,
  ///         or if the container yields no matching elements and @p emptyOk is
  ///         `false`.
  /// @return Shared pointer to the fully assembled container node.
  [[nodiscard]] std::shared_ptr<
      Acts::Experimental::CylinderContainerBlueprintNode>
  build() const;

  /// @brief Build the container node and attach it as a child of @p node.
  ///
  /// Equivalent to `node.addChild(build())`.
  /// @param node Blueprint node that will receive the built container as a
  ///             child.
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

/// @brief Fluent builder that assembles a combined barrel + endcap subdetector
/// into a @ref CylinderContainerBlueprintNode arranged along the Z axis.
///
/// Instances are obtained from @ref BlueprintBuilder::barrelEndcap(). The
/// builder inspects the subtree of the provided assembly element for barrel and
/// endcap children (using the backend's `isBarrel` / `isEndcap` / `isTracker`
/// predicates) and delegates individual layer assembly to @ref LayerAssembler
/// internally.
///
/// Typical usage:
/// @code
/// builder.barrelEndcap()
///   .setAssembly(innerTrackerElement)
///   .setAxes(barrelAxes, endcapAxes)
///   .setLayerFilter(layerPattern)
///   .addTo(rootNode);
/// @endcode
///
/// @tparam BackendT Geometry backend satisfying @ref detail::BlueprintBackend.
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
  /// The @ref LayerAssembler specialisation for this backend.
  using LayerAssembler = Acts::Experimental::LayerAssembler<BackendT>;
  /// Callback type that can replace or wrap a
  /// @ref CylinderContainerBlueprintNode.
  using ContainerCustomizer = std::function<
      std::shared_ptr<Acts::Experimental::CylinderContainerBlueprintNode>(
          const Element&,
          std::shared_ptr<Acts::Experimental::CylinderContainerBlueprintNode>)>;
  /// Simplified callback type that mutates a
  /// @ref CylinderContainerBlueprintNode in-place.
  using ReplacingContainerCustomizer = std::function<void(
      const Element&, Acts::Experimental::CylinderContainerBlueprintNode&)>;

  /// @brief Construct a @ref BarrelEndcapAssembler bound to @p builder.
  /// @param builder The owning @ref BlueprintBuilder; must outlive this object.
  explicit BarrelEndcapAssembler(const Builder& builder);

  /// @brief Build and return the assembled barrel+endcap container node.
  ///
  /// Locates barrel and endcap sub-elements inside the assembly, creates one
  /// @ref LayerAssembler -based barrel container and one or more endcap
  /// containers, then returns a Z-axis @ref CylinderContainerBlueprintNode
  /// holding them all.
  ///
  /// @throws std::runtime_error if the assembly element has not been set, if
  ///         axes are required by the backend but not provided, if the layer
  ///         filter has not been set, or if more than one barrel element is
  ///         found inside the assembly.
  /// @return Shared pointer to the assembled Z-axis container node.
  [[nodiscard]] std::shared_ptr<
      Acts::Experimental::CylinderContainerBlueprintNode>
  build() const;

  /// @brief Build the container node and attach it as a child of @p node.
  ///
  /// Equivalent to `node.addChild(build())`.
  /// @param node Blueprint node that will receive the built container as a
  ///             child.
  void addTo(Acts::Experimental::BlueprintNode& node) const&&;

  /// @brief Register a customizer callback forwarded to each inner
  /// @ref LayerAssembler.
  ///
  /// The callback may return a different or wrapped @ref LayerBlueprintNode.
  /// @param customizer See @ref LayerAssembler::onLayer(LayerCustomizer).
  /// @return `*this` (rvalue).
  [[nodiscard]] BarrelEndcapAssembler&& onLayer(
      typename LayerAssembler::LayerCustomizer customizer) &&;

  /// @brief Register an in-place mutating layer callback forwarded to each
  /// inner @ref LayerAssembler.
  /// @param customizer See
  ///        @ref LayerAssembler::onLayer(ReplacingLayerCustomizer).
  /// @return `*this` (rvalue).
  BarrelEndcapAssembler&& onLayer(
      typename LayerAssembler::ReplacingLayerCustomizer customizer) &&;

  /// @brief Register a customizer callback invoked for each barrel or endcap
  /// container node.
  ///
  /// The callback may return a different or wrapped
  /// @ref CylinderContainerBlueprintNode.
  /// @param customizer Customizer function for the container
  /// @return `*this` (rvalue).
  [[nodiscard]] BarrelEndcapAssembler&& onContainer(
      ContainerCustomizer customizer) &&;

  /// @brief Register an in-place mutating callback invoked for each barrel or
  /// endcap container node.
  /// @param customizer Customizer function for the container
  /// @return `*this` (rvalue).
  [[nodiscard]] BarrelEndcapAssembler&& onContainer(
      ReplacingContainerCustomizer customizer) &&;

  /// @brief Set the top-level detector element whose subtree is searched for
  /// barrel and endcap elements.
  /// @param assembly Root element of the barrel+endcap sub-detector.
  /// @return `*this` (rvalue).
  [[nodiscard]] BarrelEndcapAssembler&& setAssembly(const Element& assembly) &&;

  /// @brief Set the axis definitions for both barrel and endcap layers at once.
  ///
  /// Only available when the backend satisfies @ref detail::HasAxisDefinition.
  /// @param barrel Axis definition forwarded to barrel @ref LayerAssembler s.
  /// @param endcap Axis definition forwarded to endcap @ref LayerAssembler s.
  /// @return `*this` (rvalue).
  [[nodiscard]] BarrelEndcapAssembler&& setAxes(AxisDefinition barrel,
                                                AxisDefinition endcap) &&
    requires(detail::HasAxisDefinition<BackendT>);

  /// @brief Set the axis definition used for endcap layers only.
  ///
  /// Only available when the backend satisfies @ref detail::HasAxisDefinition.
  /// @param axes Axis definition forwarded to endcap @ref LayerAssembler s.
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

/// @brief High-level builder that converts a backend detector element hierarchy
/// into an ACTS blueprint node tree.
///
/// @ref BlueprintBuilder provides the entry points for blueprint construction:
/// - @ref BlueprintBuilder::layers() returns a @ref LayerAssembler
///   for building a flat layer stack,
/// - @ref BlueprintBuilder::barrelEndcap() returns a
///   @ref BarrelEndcapAssembler for combined barrel+endcap sub-detectors.
///
/// It also exposes helpers for traversing and querying the detector element
/// hierarchy, which are used internally by the assembler classes.
///
/// The builder is parameterised on a backend type @p BackendT which must
/// satisfy the @ref detail::BlueprintBackend concept. A backend encapsulates
/// all geometry-framework-specific knowledge (e.g. DD4hep `DetElement`
/// navigation, type-flag interpretation, surface conversion).
///
/// @tparam BackendT Geometry backend satisfying @ref detail::BlueprintBackend.
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
  /// The @ref LayerAssembler specialisation for this backend.
  using LayerAssembler = Acts::Experimental::LayerAssembler<Backend>;
  /// The @ref BarrelEndcapAssembler specialisation for this backend.
  using BarrelEndcapAssembler =
      Acts::Experimental::BarrelEndcapAssembler<Backend>;

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
  /// Recursively collects all sensitive descendants of @p detElement and
  /// delegates to the two-argument overload.
  /// @param detElement Parent detector element whose subtree is scanned for
  ///                   sensitive volumes.
  /// @param layerSpec  Specification controlling surface axes and optional
  ///                   name/transform overrides.
  /// @return Shared pointer to the constructed @ref LayerBlueprintNode.
  std::shared_ptr<Acts::Experimental::LayerBlueprintNode> makeLayer(
      const Element& detElement, const LayerSpec& layerSpec) const;

  /// @brief Create a @ref LayerBlueprintNode from an explicit list of sensitive
  /// elements.
  ///
  /// The layer node name is derived from the full path to @p parent in the
  /// element hierarchy (separator `"|"`), with an optional suffix from
  /// `layerSpec.layerName`.
  /// @param parent     Detector element that contextualises the layer
  ///                   (used for naming and transform extraction).
  /// @param sensitives Span of sensitive detector elements to be converted to
  ///                   surfaces and assigned to the node.
  /// @param layerSpec  Specification controlling surface axes and optional
  ///                   name/transform overrides.
  /// @return Shared pointer to the constructed @ref LayerBlueprintNode.
  std::shared_ptr<Acts::Experimental::LayerBlueprintNode> makeLayer(
      const Element& parent, std::span<const Element> sensitives,
      const LayerSpec& layerSpec) const;

  /// @brief Create a @ref LayerAssembler bound to this builder.
  ///
  /// The returned assembler must be configured (layer type, filter, container)
  /// and then finalised via @ref LayerAssembler::build() or
  /// @ref LayerAssembler::addTo().
  /// @return A new @ref LayerAssembler instance.
  [[nodiscard]] LayerAssembler layers() const;

  /// @brief Create a @ref BarrelEndcapAssembler bound to this builder.
  ///
  /// The returned assembler must be configured (assembly element, axes, layer
  /// filter) and then finalised via @ref BarrelEndcapAssembler::build() or
  /// @ref BarrelEndcapAssembler::addTo().
  /// @return A new @ref BarrelEndcapAssembler instance.
  [[nodiscard]] BarrelEndcapAssembler barrelEndcap() const;

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
                                   const std::string& separator = "|") const;

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
  std::vector<Element> findBarrelElements(const Element& assembly) const;

  /// @brief Collect all endcap tracker elements within an assembly subtree.
  ///
  /// An element is included when the backend reports it as both a tracker
  /// element and an endcap element.
  /// @param assembly Root element of the assembly subtree to search.
  /// @return Vector of endcap elements in depth-first order.
  std::vector<Element> findEndcapElements(const Element& assembly) const;

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
