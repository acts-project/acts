// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TypeTraits.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"

namespace Acts {

namespace Concepts {
namespace VertexFinder {

template <typename T>
using state_t = typename T::State;

METHOD_TRAIT(find_t, find);

// clang-format off
    template <typename S>
      struct VertexFinderConcept {
        constexpr static bool state_exists = exists<state_t, S>;
        static_assert(state_exists, "State type not found");
        
        constexpr static bool find_exists = has_method<const S, Result<std::vector<Vertex<typename S::InputTrack_t>>>,
         find_t, const std::vector<const typename S::InputTrack_t*>&, 
         const VertexingOptions<typename S::InputTrack_t>&, typename S::State&>;
        static_assert(find_exists, "find method not found");

        constexpr static bool value = require<state_exists, find_exists>;
      };
// clang-format on
}  // namespace VertexFinder
}  // namespace Concepts

template <typename finder>
constexpr bool VertexFinderConcept =
    Acts::Concepts ::VertexFinder::VertexFinderConcept<finder>::value;

}  // namespace Acts
