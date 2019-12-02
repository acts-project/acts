// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename vfitter_t>
auto Acts::TrackDensityVertexFinder<vfitter_t>::find(
    const std::vector<InputTrack_t>& trackVector,
    const VertexFinderOptions<InputTrack_t>& vFinderOptions) const
    -> Result<std::vector<Vertex<InputTrack_t>>> {

  return {};
}
