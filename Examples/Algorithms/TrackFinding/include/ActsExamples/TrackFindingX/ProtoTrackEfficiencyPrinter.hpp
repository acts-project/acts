// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/EventDataTransforms.hpp"
// #include "Acts/Utilities/KDTree.hpp"

#include <string>
#include <vector>

namespace ActsExamples {

class ProtoTrackEfficencyPrinter final : public IAlgorithm {
 public:
  struct Config {
    std::string testProtoTracks;
    std::string refProtoTracks;
    // std::string spacePoints;
  };

  ProtoTrackEfficencyPrinter(Config cfg, Acts::Logging::Level lvl)
      : IAlgorithm("ProtoTrackEfficencyPrinter", lvl), m_cfg(cfg) {
    m_testProtoTracks.initialize(m_cfg.testProtoTracks);
    m_refProtoTracks.initialize(m_cfg.refProtoTracks);
    // m_inputSpacePoints.initialize(m_cfg.spacePoints);
  }

  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext &context) const override {
    const auto testTracks = m_testProtoTracks(context);
    const auto refTracks = m_refProtoTracks(context);

    std::vector<double> effs(testTracks.size(), 0.0);

    for (auto [refTrack, eff] : Acts::zip(refTracks, effs)) {
      ProtoTrack intersection;
      for (const auto &testTrack : testTracks) {
        std::set_intersection(refTrack.begin(), refTrack.end(),
                              testTrack.begin(), testTrack.end(),
                              intersection.begin());
        eff = std::max(
            eff, static_cast<double>(intersection.size()) / refTrack.size());
        intersection.clear();
      }
    }

    std::sort(effs.begin(), effs.end());

    const static std::vector<double> thresholds = {0., .1, .2, .3, .4,
                                                   .5, .6, .7, .8, .9};

    auto it = effs.begin();
    std::vector<std::size_t> hist;
    for (double threshold : thresholds) {
      auto endIt = std::find_if(it, effs.end(),
                                [&](double eff) { return eff >= threshold; });
      hist.push_back(std::distance(it, endIt));
      it = endIt;
    }

    const auto max = *std::max_element(hist.begin(), hist.end());
    const int colMax = 40;

    ACTS_INFO("Prototrack efficiency histogram:");
    for (const auto &[v, th] : Acts::zip(hist, thresholds)) {
      auto rel = v / static_cast<float>(max));
      auto l = std::round(rel * colMax);
      std::string str(l, '#');
      ACTS_INFO(">=" << th << " | " << str);
    }

    return {};
  }

  const Config &config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<ProtoTrackContainer> m_testProtoTracks{this,
                                                        "InputTestProtoTracks"};
  ReadDataHandle<ProtoTrackContainer> m_refProtoTracks{this,
                                                       "InputRefProtoTracks"};
  // ReadDataHandle<SimSpacePointContainer> m_inputSpacePoints{
  // this, "InputSpacePoints"};
};

}  // namespace ActsExamples
