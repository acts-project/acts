// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/BlueprintBuilder.hpp"
#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/LayerBlueprintNode.hpp"

#include <memory>
#include <optional>

using namespace Acts::Experimental;
using Acts::Experimental::detail::OnContainerMutatesContainer;
using Acts::Experimental::detail::OnContainerReturnsNode;
using Acts::Experimental::detail::OnLayerMutatesLayer;
using Acts::Experimental::detail::OnLayerReturnsNode;

namespace ActsTests {

// Regression coverage for the backward-compatible `onLayer`/`onContainer`
// customizer contract. The customizers were widened from returning the concrete
// node type (`LayerBlueprintNode` / `ContainerBlueprintNode`) to the base
// `BlueprintNode`. The accepting concepts use `std::convertible_to` (not
// `std::same_as`) so that pre-existing client callbacks returning the concrete
// derived pointer still satisfy the contract.

using Element = int;

// --- onLayer ---------------------------------------------------------------

// Old-style callback: returns the concrete LayerBlueprintNode pointer.
[[maybe_unused]] auto oldStyleLayer =
    [](const std::optional<Element>&, std::shared_ptr<LayerBlueprintNode> layer)
    -> std::shared_ptr<LayerBlueprintNode> { return layer; };
static_assert(OnLayerReturnsNode<Element, decltype(oldStyleLayer)>,
              "onLayer callback returning shared_ptr<LayerBlueprintNode> must "
              "remain accepted (backward compatibility)");

// New-style callback: returns the base BlueprintNode pointer (wrap/replace).
[[maybe_unused]] auto newStyleLayer =
    [](const std::optional<Element>&, std::shared_ptr<LayerBlueprintNode> layer)
    -> std::shared_ptr<BlueprintNode> { return layer; };
static_assert(OnLayerReturnsNode<Element, decltype(newStyleLayer)>,
              "onLayer callback returning shared_ptr<BlueprintNode> must be "
              "accepted");

// In-place mutating callback: returns void.
[[maybe_unused]] auto mutateLayer = [](const std::optional<Element>&,
                                       LayerBlueprintNode&) {};
static_assert(OnLayerMutatesLayer<Element, decltype(mutateLayer)>,
              "void-returning in-place onLayer callback must be accepted");
static_assert(!OnLayerReturnsNode<Element, decltype(mutateLayer)>,
              "void-returning callback must not match the node-returning form");

// --- onContainer -----------------------------------------------------------

// Old-style callback: returns the concrete ContainerBlueprintNode pointer.
[[maybe_unused]] auto oldStyleContainer =
    [](const Element&, std::shared_ptr<ContainerBlueprintNode> container)
    -> std::shared_ptr<ContainerBlueprintNode> { return container; };
static_assert(
    OnContainerReturnsNode<Element, decltype(oldStyleContainer)>,
    "onContainer callback returning shared_ptr<ContainerBlueprintNode> must "
    "remain accepted (backward compatibility)");

// New-style callback: returns the base BlueprintNode pointer (wrap/replace).
[[maybe_unused]] auto newStyleContainer =
    [](const Element&, std::shared_ptr<ContainerBlueprintNode> container)
    -> std::shared_ptr<BlueprintNode> { return container; };
static_assert(
    OnContainerReturnsNode<Element, decltype(newStyleContainer)>,
    "onContainer callback returning shared_ptr<BlueprintNode> must be "
    "accepted");

// In-place mutating callback: returns void.
[[maybe_unused]] auto mutateContainer = [](const Element&,
                                           ContainerBlueprintNode&) {};
static_assert(OnContainerMutatesContainer<Element, decltype(mutateContainer)>,
              "void-returning in-place onContainer callback must be accepted");
static_assert(!OnContainerReturnsNode<Element, decltype(mutateContainer)>,
              "void-returning callback must not match the node-returning form");

BOOST_AUTO_TEST_SUITE(GeometrySuite)

BOOST_AUTO_TEST_CASE(CustomizerReturnTypeBackwardCompatibility) {
  // The substantive checks are the static_asserts above; this case ensures the
  // translation unit is exercised by the test runner.
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
