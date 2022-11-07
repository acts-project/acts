// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepMultiTrajectoryWriter.hpp"

#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <stdexcept>

namespace ActsExamples {

EDM4hepMultiTrajectoryWriter::EDM4hepMultiTrajectoryWriter(
    const EDM4hepMultiTrajectoryWriter::Config& config,
    Acts::Logging::Level level)
    : WriterT<TrajectoriesContainer>(config.inputTrajectories,
                                     "EDM4hepMultiTrajectoryWriter", level),
      m_cfg(config),
      m_writer(config.outputPath, &m_store) {
  if (m_cfg.inputTrajectories.empty()) {
    throw std::invalid_argument("Missing input trajectories collection");
  }

  m_trackCollection = &m_store.create<edm4hep::TrackCollection>("ActsTracks");
  m_writer.registerForWrite("ActsTracks");
}

ActsExamples::ProcessCode EDM4hepMultiTrajectoryWriter::endRun() {
  m_writer.finish();

  return ProcessCode::SUCCESS;
}

ProcessCode EDM4hepMultiTrajectoryWriter::writeT(
    const AlgorithmContext& ctx, const TrajectoriesContainer& trajectories) {
  const auto& hitParticlesMap =
      ctx.eventStore.get<IndexMultimap<ActsFatras::Barcode>>(
          m_cfg.inputMeasurementParticlesMap);

  for (const auto& from : trajectories) {
    for (const auto& trackTip : from.tips()) {
      auto to = m_trackCollection->create();
      EDM4hepUtil::writeTrajectory(from, to, trackTip, hitParticlesMap);
    }
  }

  m_writer.writeEvent();
  m_store.clearCollections();

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
