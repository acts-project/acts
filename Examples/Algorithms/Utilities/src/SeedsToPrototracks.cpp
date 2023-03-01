// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/SeedsToPrototracks.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/EventDataTransforms.hpp"

namespace ActsExamples {

ProcessCode SeedsToPrototracks::execute(const AlgorithmContext& ctx) const {
  const auto seeds = ctx.eventStore.get<SimSeedContainer>(m_cfg.inputSeeds);

  auto tracks = seedsToPrototracks(seeds);

  ctx.eventStore.add(m_cfg.outputProtoTracks, std::move(tracks));

  return ProcessCode::SUCCESS;
}
}  // namespace ActsExamples
