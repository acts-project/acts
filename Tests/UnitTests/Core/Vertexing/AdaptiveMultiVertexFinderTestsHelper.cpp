// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFinder.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFitter.hpp"
#include "Acts/Vertexing/GaussianTrackDensity.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"

namespace Acts::Test {

namespace BoundAmvf {
using Propagator = Acts::Propagator<EigenStepper<>>;
using Linearizer = HelicalTrackLinearizer<Propagator>;
using Fitter = AdaptiveMultiVertexFitter<BoundTrackParameters, Linearizer>;
using SeedFinder =
    TrackDensityVertexFinder<Fitter,
                             GaussianTrackDensity<BoundTrackParameters>>;
using Finder = AdaptiveMultiVertexFinder<Fitter, SeedFinder>;
}  // namespace BoundAmvf

// template <typename input_track_t, typename finder_t>
// Result<std::vector<Vertex<input_track_t>>> find_impl(
// const finder_t& finder, const std::vector<const input_track_t*>& allTracks,
// const VertexingOptions<input_track_t>& vertexingOptions,
// typename finder_t::state_t& state) {
// return finder.find(allTracks, vertexingOptions, state);
// }

Result<std::vector<Vertex<BoundTrackParameters>>> find(
    const BoundAmvf::Finder& finder,
    const std::vector<const BoundTrackParameters*>& allTracks,
    const VertexingOptions<BoundTrackParameters>& vertexingOptions,
    BoundAmvf::Finder::State& state) {
  return finder.find(allTracks, vertexingOptions, state);
}

}  // namespace Acts::Test
