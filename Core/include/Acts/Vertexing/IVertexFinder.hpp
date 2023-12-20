// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Any.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"

#include <vector>

namespace Acts {

class Vertex;
struct InputTrack;
struct VertexingOptions;

class IVertexFinder {
 public:
  using State = Acts::AnyBase<128>;

  virtual Result<std::vector<Vertex>> find(
      const std::vector<InputTrack>& trackVector,
      const VertexingOptions& vertexingOptions, State& state) const = 0;

  virtual State makeState() const = 0;

  virtual bool hasTrivialState() const = 0;

  virtual void setTracksToRemove(
      State& state, const std::vector<InputTrack>& removedTracks) const = 0;
};
}  // namespace Acts
