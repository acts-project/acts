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

  template <typename T>
  using track_t = typename T::InputTrack_t;
  template <typename T>
  using propagator_t = typename T::Propagator_t;
  template <typename T>
  using linearizer_t = typename T::Linearizer_t;
  template <typename T>
  using bfield_t = typename T::Bfield_t;

  METHOD_TRAIT(fit_t, fit);

  // clang-format off
    template <typename S>
      struct VertexFitterConcept {
        constexpr static bool fit_exists = has_method<const S, Result<Vertex<typename S::InputTrack_t>>,
         fit_t, const std::vector<typename S::InputTrack_t>&, 
         const VertexFitterOptions<typename S::InputTrack_t>&>;
        static_assert(fit_exists, "fit method not found");
        
        constexpr static bool track_exists = exists<track_t, S>;
        static_assert(track_exists, "Track type not found");
        constexpr static bool propagator_exists = exists<propagator_t, S>;
        static_assert(propagator_exists, "Propagator type not found");
        constexpr static bool linearizer_exists = exists<linearizer_t, S>;
        static_assert(linearizer_exists, "Linearizer type not found");
        constexpr static bool bfield_exists = exists<bfield_t, S>;
        static_assert(bfield_exists, "BField type not found");

        constexpr static bool value = require<fit_exists,
                                              track_exists,
                                              propagator_exists,
                                              linearizer_exists,
                                              bfield_exists>;
      };
  // clang-format on
  }  // namespace VertexFitter
}  // namespace concept

template <typename fitter>
constexpr bool VertexFitterConcept =
    Acts::concept ::VertexFitter::VertexFitterConcept<fitter>::value;

}  // namespace Acts