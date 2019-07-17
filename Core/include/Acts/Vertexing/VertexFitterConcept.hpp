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
#include "Acts/Vertexing/VertexFitterOptions.hpp"

namespace Acts {

namespace concept {
  namespace VertexFitter {

  METHOD_TRAIT(fit_t, fit);

  // clang-format off
    template <typename S>
      struct VertexFitterConcept {
        
        constexpr static bool fit_exists = has_method<const S, Result<Vertex<typename S::InputTrack>>,
         fit_t, const std::vector<typename S::InputTrack>&, 
         const VertexFitterOptions<typename S::InputTrack>&>;
        static_assert(fit_exists, "fit method not found");

        constexpr static bool value = require<fit_exists>;
      };
  // clang-format on
  }  // namespace VertexFitter
}  // namespace concept

template <typename fitter>
constexpr bool VertexFitterConcept =
    Acts::concept ::VertexFitter::VertexFitterConcept<fitter>::value;

}  // namespace Acts