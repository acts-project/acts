// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/EventDataTransforms.hpp"

#include <mutex>
#include <string>
#include <vector>

#include <boost/histogram.hpp>

namespace ActsExamples {

class ProtoTrackEfficiencyPrinter final : public IAlgorithm {
 public:
  struct Config {
    std::string testProtoTracks;
    std::string refProtoTracks;
    // std::string spacePoints;
  };

  ProtoTrackEfficiencyPrinter(Config cfg, Acts::Logging::Level lvl)
      : IAlgorithm("ProtoTrackEfficencyPrinter", lvl),
        m_cfg(cfg),
        m_histogram(boost::histogram::make_histogram(
            boost::histogram::axis::regular<>(10, 0.0, 1.0))) {
    m_testProtoTracks.initialize(m_cfg.testProtoTracks);
    m_refProtoTracks.initialize(m_cfg.refProtoTracks);
    // m_inputSpacePoints.initialize(m_cfg.spacePoints);
  }

  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext &context) const override;

  ActsExamples::ProcessCode finalize() override;

  const Config &config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<ProtoTrackContainer> m_testProtoTracks{this,
                                                        "InputTestProtoTracks"};
  ReadDataHandle<ProtoTrackContainer> m_refProtoTracks{this,
                                                       "InputRefProtoTracks"};

  using Hist = decltype(boost::histogram::make_histogram(
      std::declval<boost::histogram::axis::regular<>>()));
  mutable Hist m_histogram;

  mutable std::mutex m_histogramMutex;
};

}  // namespace ActsExamples
