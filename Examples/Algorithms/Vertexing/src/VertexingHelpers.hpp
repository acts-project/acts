// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Vertexing/Vertex.hpp"
#include "ActsExamples/EventData/ProtoVertex.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"

#include <vector>

namespace ActsExamples {

/// Create a pointers container for all track parameters in the input container.
///
/// @param trackParameters input examples track parameters container
/// @return track parameters pointer container referencing the input tracks
inline std::vector<const Acts::BoundTrackParameters*>
makeTrackParametersPointerContainer(
    const TrackParametersContainer& trackParameters) {
  std::vector<const Acts::BoundTrackParameters*> trackParametersPointers;
  trackParametersPointers.reserve(trackParameters.size());

  for (const auto& trackParam : trackParameters) {
    trackParametersPointers.push_back(&trackParam);
  }
  return trackParametersPointers;
}

inline auto makeParameterContainers(
    const ActsExamples::AlgorithmContext& ctx,
    const ReadDataHandle<std::vector<Acts::BoundTrackParameters>>&
        inputTrackParametersHandle,
    const ReadDataHandle<TrajectoriesContainer>& inputTrajectoriesHandle) {
  std::vector<Acts::BoundTrackParameters> inputTrackParameters;
  std::vector<const Acts::BoundTrackParameters*> inputTrackPointers;

  if (inputTrackParametersHandle.isInitialized()) {
    const auto& tmp = inputTrackParametersHandle(ctx);
    inputTrackParameters = tmp;
    inputTrackPointers = makeTrackParametersPointerContainer(tmp);
  } else {
    const auto& inputTrajectories = inputTrajectoriesHandle(ctx);

    for (const auto& trajectories : inputTrajectories) {
      for (auto tip : trajectories.tips()) {
        if (!trajectories.hasTrackParameters(tip)) {
          continue;
        }
        const auto& trackParam = trajectories.trackParameters(tip);
        inputTrackParameters.push_back(trackParam);
        inputTrackPointers.push_back(&trackParam);
      }
    }
  }

  return std::make_tuple(inputTrackParameters, inputTrackPointers);
}

/// Create proto vertices from reconstructed vertices.
///
/// @param trackParameters input track parameters container
/// @param vertices reconstructed vertices
/// @return proto vertices corresponding to the reconstructed vertices
///
/// Assumes that the original parameters pointers in the vertices point to
/// elements in the given input track parameters container. If that is not the
/// case the behaviour is undefined.
inline ProtoVertexContainer makeProtoVertices(
    const TrackParametersContainer& trackParameters,
    const std::vector<Acts::Vertex<Acts::BoundTrackParameters>>& vertices) {
  ProtoVertexContainer protoVertices;
  protoVertices.reserve(vertices.size());

  // use pointer arithmetic to compute the position/index of the original
  // parameters in the input track parameters container.
  const Acts::BoundTrackParameters* first = trackParameters.data();
  for (const auto& vertex : vertices) {
    ProtoVertex protoVertex;
    protoVertex.reserve(vertex.tracks().size());

    for (const auto& track : vertex.tracks()) {
      protoVertex.push_back(std::distance(first, track.originalParams));
    }
    protoVertices.push_back(std::move(protoVertex));
  }

  return protoVertices;
}

}  // namespace ActsExamples
