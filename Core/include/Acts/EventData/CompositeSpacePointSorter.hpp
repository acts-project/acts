// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/CompositeSpacePoint.hpp"
namespace Acts::Experimental {
/// @brief Define the concept of the space point measurement sorter. The sorter shall take a collection
///         of station space points and split them into straw && strip hits.
///         Hits
/// from each category are then subdivided further into the particular detector
/// layers.
///
///      A possible implementation of the CompositeSpacePointSorter needs to
///      have
/// the following attributes
///
///      using SpVec_t =  Standard container satisfiyng the
/// CompositeSpacePointContainer concept
///
///
///      const std::vector<SpVec_t>& strawHits();
///      const std::vector<SpVec_t>& stripHits();
///  Each SpVec_t contains all measurements from a particular detector layer
template <typename Splitter_t, typename SpacePointCont_t>
concept CompositeSpacePointSorter =
    CompositeSpacePointContainer<SpacePointCont_t> &&
    requires(Splitter_t sorter) {
      /// @brief Return the straw-hit space point sorted by straw layer
      {
        sorter.strawHits()
      } -> std::same_as<const std::vector<SpacePointCont_t>&>;
      /// @brief Return the strip-hit  space points sorted by detector layer
      {
        sorter.stripHits()
      } -> std::same_as<const std::vector<SpacePointCont_t>&>;
    };

}  // namespace Acts::Experimental
