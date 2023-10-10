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
namespace VertexFitter {

template <typename T>
using track_t = typename T::InputTrack_t;
template <typename T>
using propagator_t = typename T::Propagator_t;
template <typename T>
using linearizer_t = typename T::Linearizer_t;
template <typename T>
using state_t = typename T::State;

METHOD_TRAIT(fit_t, fit);

// clang-format off
    template <typename S>
      struct VertexFitterConcept {
        constexpr static bool fit_exists = has_method<const S, Result<Vertex<typename S::InputTrack_t>>,
         fit_t, 
         const std::vector<const typename S::InputTrack_t*>&, 
         const typename S::Linearizer_t&,
         const VertexingOptions<typename S::InputTrack_t>&,
         typename S::State&>;
        static_assert(fit_exists, "fit method not found");

        constexpr static bool track_exists = exists<track_t, S>;
        static_assert(track_exists, "Track type not found");
        constexpr static bool propagator_exists = exists<propagator_t, S>;
        static_assert(propagator_exists, "Propagator type not found");
        constexpr static bool linearizer_exists = exists<linearizer_t, S>;
        static_assert(linearizer_exists, "Linearizer type not found");
        constexpr static bool state_exists = exists<state_t, S>;
        static_assert(state_exists, "State type not found");

        constexpr static bool value = require<fit_exists,
                                              track_exists,
                                              propagator_exists,
                                              linearizer_exists,
                                              state_exists>;
      };
// clang-format on
}  // namespace VertexFitter
}  // namespace Concepts

template <typename fitter>
constexpr bool VertexFitterConcept =
    Acts::Concepts ::VertexFitter::VertexFitterConcept<fitter>::value;

}  // namespace Acts
