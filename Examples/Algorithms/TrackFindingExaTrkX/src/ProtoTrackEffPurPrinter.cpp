// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingExaTrkX/ProtoTrackEffPurPrinter.hpp"

#include <algorithm>
#include <vector>

#include <boost/histogram/ostream.hpp>

namespace bh = boost::histogram;

ActsExamples::ProtoTrackEffPurPrinter::ProtoTrackEffPurPrinter(
    Config cfg, Acts::Logging::Level lvl)
    : IAlgorithm("ProtoTrackEfficencyPrinter", lvl),
      m_cfg(cfg),
      m_effHistogram(bh::make_histogram(bh::axis::regular<>(10, 0.0, 1.0))),
      m_purHistogram(bh::make_histogram(bh::axis::regular<>(10, 0.0, 1.0))),
      m_countPerTrackHist(bh::make_histogram(bh::axis::integer<>(1, 5))) {
  m_testProtoTracks.initialize(m_cfg.testProtoTracks);
  m_refProtoTracks.initialize(m_cfg.refProtoTracks);
}

ActsExamples::ProcessCode ActsExamples::ProtoTrackEffPurPrinter::execute(
    const ActsExamples::AlgorithmContext &context) const {
  auto testTracks = m_testProtoTracks(context);
  auto truthTracks = m_refProtoTracks(context);

  ACTS_INFO("Receiving " << truthTracks.size() << " reference tracks");
  truthTracks.erase(std::remove_if(truthTracks.begin(), truthTracks.end(),
                                   [](const auto &t) { return t.size() < 3; }),
                    truthTracks.end());
  ACTS_INFO(" -> " << truthTracks.size() << " tracks with size >= 3");
  ACTS_INFO("Receiving " << testTracks.size() << " test tracks");

  // Build id-to-truth-track map
  // For now we assume that each space point only belongs to one truth track
  // (this might be wrong if we merge clusters)
  std::unordered_map<std::size_t, std::size_t> idToTruthTrack;

  for (auto i = 0ul; i < truthTracks.size(); ++i) {
    for (auto el : truthTracks[i]) {
      idToTruthTrack[el] = i;
    }
  }

  // Go through test tracks and search truth track ids
  constexpr static std::size_t invalid =
      std::numeric_limits<std::size_t>::max();

  std::vector<float> testTrackPurities(testTracks.size(), 0.f);
  std::vector<float> trueTrackEfficiencies(truthTracks.size(), 0.f);
  std::vector<std::size_t> trueTracksPerTestTrack;

  {
    // allocate once for memory optimization
    std::vector<std::size_t> truthTrackIds;
    std::vector<std::size_t> truthTrackUniqueIds;

    for (auto testId = 0ul; testId < testTracks.size(); ++testId) {
      const auto &testTrack = testTracks[testId];

      // Build vector of truth track ids
      for (const auto el : testTrack) {
        // invalid means, the points is not associated to any track (noise)
        auto tid = idToTruthTrack.count(el) > 0 ? idToTruthTrack[el] : invalid;
        truthTrackIds.push_back(tid);
      }

      // Find out which truth track ids are most fequent
      std::sort(truthTrackIds.begin(), truthTrackIds.end());
      std::unique_copy(truthTrackIds.begin(), truthTrackIds.end(),
                       std::back_inserter(truthTrackUniqueIds));

      std::sort(
          truthTrackUniqueIds.begin(), truthTrackUniqueIds.end(),
          [&](const auto &a, const auto &b) {
            auto ac = std::count(truthTrackIds.begin(), truthTrackIds.end(), a);
            auto bc = std::count(truthTrackIds.begin(), truthTrackIds.end(), b);
            return ac > bc;  // sort many-to-few
          });

      trueTracksPerTestTrack.push_back(truthTrackUniqueIds.size());

      // compute metrics
      if (truthTrackUniqueIds[0] == invalid &&
          truthTrackUniqueIds.size() == 1) {
        testTrackPurities[testId] = 0.f;
      } else {
        const auto truthId = truthTrackUniqueIds[0] == invalid
                                 ? truthTrackUniqueIds[1]
                                 : truthTrackUniqueIds[0];
        const auto nhits =
            std::count(truthTrackIds.begin(), truthTrackIds.end(), truthId);

        const auto &truthTrack = truthTracks[truthId];

        // Ensure eff,pur < 1.0, so the binning in the histogram is nice
        float eff = static_cast<float>(nhits) / truthTrack.size();
        eff = std::min(eff, 0.9999f);

        float pur = static_cast<float>(nhits) / testTrack.size();
        pur = std::min(pur, 0.9999f);

        trueTrackEfficiencies[truthId] =
            std::max(trueTrackEfficiencies[truthId], eff);
        testTrackPurities[testId] = pur;
      }

      // clear vectors
      truthTrackIds.clear();
      truthTrackUniqueIds.clear();
    }
  }

  std::lock_guard<std::mutex>{m_histogramMutex};
  m_effHistogram.fill(trueTrackEfficiencies);
  m_purHistogram.fill(testTrackPurities);
  m_countPerTrackHist.fill(trueTracksPerTestTrack);

  return {};
}

ActsExamples::ProcessCode ActsExamples::ProtoTrackEffPurPrinter::finalize() {
  ACTS_INFO("Truth track efficiency:\n" << m_effHistogram);
  ACTS_INFO("Test track purity:\n" << m_purHistogram);
  ACTS_INFO("Particles per test track:\n" << m_countPerTrackHist);
  return {};
}
