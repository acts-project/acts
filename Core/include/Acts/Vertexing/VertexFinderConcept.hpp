// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TypeTraits.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexFinderOptions.hpp"

namespace Acts {

namespace concept {
  namespace VertexFinder {

  METHOD_TRAIT(find_t, find);

  // clang-format off
    template <typename S>
      struct VertexFinderConcept {
        
        constexpr static bool find_exists = has_method<const S, Result<std::vector<Vertex<typename S::InputTrack>>>,
         find_t, const std::vector<typename S::InputTrack>&, 
         const VertexFinderOptions<typename S::InputTrack>&>;
        static_assert(find_exists, "find method not found");

        constexpr static bool value = require<find_exists>;
      };
  // clang-format on
  }  // namespace VertexFinder
}  // namespace concept

template <typename finder>
constexpr bool VertexFinderConcept =
    Acts::concept ::VertexFinder::VertexFinderConcept<finder>::value;

}  // namespace Acts