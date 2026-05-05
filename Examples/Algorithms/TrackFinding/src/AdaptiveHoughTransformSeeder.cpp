// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/AdaptiveHoughTransformSeeder.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/ScopedTimer.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <map>
#include <numeric>
#include <ostream>
#include <stdexcept>

namespace ActsExamples {

AdaptiveHoughTransformSeeder::AdaptiveHoughTransformSeeder(
    const Config &cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("AdaptiveHoughTransformSeeder", std::move(logger)),
      m_cfg(cfg) {
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument(
        "AdaptiveHoughTransformSeeder: Input space points collection name "
        "empty");
  }
  if (m_cfg.outputSeeds.empty()) {
    throw std::invalid_argument(
        "AdaptiveHoughTransformSeeder: Output seeds collection name empty");
  }

  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_outputSeeds.initialize(m_cfg.outputSeeds);
}

// When traversing parameters space in various directions it is useful
// to record the "path".
// Acts::Experimental::HoughAccumulatorSection has simple functionality allowing that in a
// form of a plain vector of floats. These numbers need to be accessed
// consistently, thus this indices.
namespace {
constexpr unsigned int phiSplitMinIndex = 0;
constexpr unsigned int phiSplitWidthIndex = 1;
constexpr unsigned int zSplitMinIndex = 2;
constexpr unsigned int zSplitWidthIndex = 3;
constexpr unsigned int cotThetaSplitMinIndex = 4;
constexpr unsigned int cotThetaSplitWidthIndex = 5;
}  // namespace

ProcessCode AdaptiveHoughTransformSeeder::execute(
    const AlgorithmContext &ctx) const {
  const SpacePointContainer &spacePoints = m_inputSpacePoints(ctx);

  // get inputs
  std::vector<PreprocessedMeasurement> measurements;
  preparePreprocessedMeasurements(spacePoints, measurements);

  // prepare initial stack
  std::vector<Acts::Experimental::HoughAccumulatorSection> stack1;
  fillStackPhiSplit(stack1, measurements);

  // split into regions in z_vertex cot theta, there is a lot of duplication and
  // it can be optimized better in the future
  for (auto &section : stack1) {
    section.updateDimensions(2.0f * m_cfg.zRange, 2.0f * m_cfg.cotThetaRange,
                             -m_cfg.zRange, -m_cfg.cotThetaRange);
  }
  std::vector<Acts::Experimental::HoughAccumulatorSection> stack2;
  {
    Acts::ScopedTimer st("splitInZCotTheta", logger());
    processStackZCotThetaSplit(stack1, stack2, measurements);
    if (m_cfg.doSecondPhase) {
      for (Acts::Experimental::HoughAccumulatorSection &section : stack2) {
        section.setHistory(zSplitMinIndex, section.xBegin());
        section.setHistory(zSplitWidthIndex, section.xSize());
        section.setHistory(cotThetaSplitMinIndex, section.yBegin());
        section.setHistory(cotThetaSplitWidthIndex, section.ySize());
      }
    }
  }
  ACTS_DEBUG("past StackZCotThetaSplit stack1 size "
             << stack1.size() << " stack2 size " << stack2.size());
  for (auto &section : stack2) {
    ACTS_DEBUG("Post z_vertex cot theta regions "
               << section.count() << " for region starting from "
               << section.xBegin() << " " << section.yBegin());
    // now need to change search space into phi - q/pT, section covers phi range
    // from initial splitting (therefore needed history)
    section.updateDimensions(
        section.history(phiSplitWidthIndex), 2.0f * m_cfg.qOverPtMin,
        section.history(phiSplitMinIndex), -m_cfg.qOverPtMin);
  }
  {
    Acts::ScopedTimer st("processQOverPtPhi", logger());
    processStackQOverPtPhi(stack2, stack1, measurements);
  }
  ACTS_DEBUG("past StackQOverPtPhi stack1 size "
             << stack1.size() << " stack2 size " << stack2.size());

  if (m_cfg.deduplicate) {
    Acts::ScopedTimer st("deduplication", logger());
    deduplicate(stack1);
    ACTS_DEBUG("past deduplication " << stack1.size());
  }

  // do scan in z_vertex - cot theta
  if (m_cfg.doSecondPhase) {
    Acts::ScopedTimer st("secondPhase", logger());

    for (auto &section : stack1) {
      section.updateDimensions(section.history(zSplitWidthIndex),
                               section.history(cotThetaSplitWidthIndex),
                               section.history(zSplitMinIndex),
                               section.history(cotThetaSplitMinIndex));
    }
    processStackZCotTheta(stack1, stack2, measurements);
    ACTS_DEBUG("past StackZCotTheta stack1 size "
               << stack1.size() << " stack2 size " << stack2.size());

    if (m_cfg.deduplicate) {
      Acts::ScopedTimer st2("deduplication", logger());
      deduplicate(stack2);
      ACTS_DEBUG("past second deduplication " << stack2.size());
    }
  }

  std::vector<Acts::Experimental::HoughAccumulatorSection> &solutions =
      m_cfg.doSecondPhase ? stack2 : stack1;
  Acts::ScopedTimer st("seedsMaking", logger());

  // post solutions
  SeedContainer seeds;
  seeds.assignSpacePointContainer(spacePoints);
  seeds.reserve(solutions.size());

  ACTS_VERBOSE("Number of solutions " << solutions.size());
  makeSeeds(seeds, solutions, measurements);

  m_outputSeeds(ctx, std::move(seeds));

  return ProcessCode::SUCCESS;
}

void AdaptiveHoughTransformSeeder::preparePreprocessedMeasurements(
    const SpacePointContainer &spacePoints,
    std::vector<PreprocessedMeasurement> &measurements) const {
  ACTS_DEBUG("Inserting " << spacePoints.size()
                          << " space points from collection\""
                          << m_cfg.inputSpacePoints << "\"");

  for (const ConstSpacePointProxy sp : spacePoints) {
    const float phi = std::atan2(sp.y(), sp.x());
    const float invr = 1.0f / sp.r();
    measurements.emplace_back(invr, phi, sp.z(), sp.index());
    // wrap phi by duplicating some seeds
    if (phi < -std::numbers::pi + m_cfg.phiWrap * std::numbers::pi) {
      measurements.emplace_back(invr, phi + 2 * std::numbers::pi, sp.z(),
                                sp.index());
    }
    if (phi > std::numbers::pi - m_cfg.phiWrap * std::numbers::pi) {
      measurements.emplace_back(invr, phi - 2 * std::numbers::pi, sp.z(),
                                sp.index());
    }
  }

  ACTS_DEBUG("Total number of " << measurements.size()
                                << " used to account for phi wrapping");
}

void AdaptiveHoughTransformSeeder::fillStackPhiSplit(
    std::vector<Acts::Experimental::HoughAccumulatorSection> &stack,
    const std::vector<PreprocessedMeasurement> &measurements) const {
  Acts::ScopedTimer st("splitInQuadrants", logger());
  const int nSplits = 8;
  const auto wedgeWidth = static_cast<float>(2.0 * std::numbers::pi / nSplits);

  for (int phiIndex = 0; phiIndex < nSplits; phiIndex++) {
    const auto startPhi = static_cast<float>(
        wedgeWidth * static_cast<float>(phiIndex) - std::numbers::pi);
    stack.emplace_back(1.2f * wedgeWidth, 2.0f * m_cfg.qOverPtMin, startPhi,
                       -m_cfg.qOverPtMin);
    stack.back().indices().resize(measurements.size());
    std::iota(std::begin(stack.back().indices()),
              std::end(stack.back().indices()), 0);
    stack.back().setHistory(phiSplitMinIndex, startPhi);
    stack.back().setHistory(phiSplitWidthIndex, wedgeWidth);

    updateSection(stack.back(), measurements, m_qOverPtPhiLineParams);
    ACTS_DEBUG("Initial split " << stack.back().count()
                                << " for region starting from " << startPhi);
    for (unsigned index : stack.back().indices()) {
      const PreprocessedMeasurement &m = measurements[index];
      ACTS_DEBUG("phi section "
                 << startPhi << " x: " << std::cos(m.phi) / m.invr
                 << " y: " << std::sin(m.phi) / m.invr << " z: " << m.z);
    }
  }
}

void AdaptiveHoughTransformSeeder::processStackQOverPtPhi(
    std::vector<Acts::Experimental::HoughAccumulatorSection> &input,
    std::vector<Acts::Experimental::HoughAccumulatorSection> &output,
    const std::vector<PreprocessedMeasurement> &measurements) const {
  struct Stats {
    double area{};
    int nSections{};
    int nLines{};
    int discardedByThresholdCut{};
    int discardedByCrossingCut{};
  };
  std::map<int, Stats> sStat;
  ExplorationOptions opt;
  opt.xMinBinSize = m_cfg.phiMinBinSize;
  opt.yMinBinSize = m_cfg.qOverPtMinBinSize;
  opt.lineFunctor = m_qOverPtPhiLineParams;
  opt.decisionFunctor = [&sStat, &cfg = m_cfg, &opt, this](
                            const Acts::Experimental::HoughAccumulatorSection &section,
                            const std::vector<PreprocessedMeasurement> &mes) {
    if (section.divisionLevel() <= 8) {
      return Decision::Drill;
    }

    if (section.count() < cfg.threshold) {
      sStat[section.divisionLevel()].discardedByThresholdCut += 1;
      return Decision::Discard;
    }
    if (section.count() < 3 * cfg.threshold &&
        !passIntersectionsCheck(section, mes, opt.lineFunctor,
                                cfg.threshold * (cfg.threshold - 1))) {
      sStat[section.divisionLevel()].discardedByCrossingCut += 1;
      return Decision::Discard;
    }

    if (section.count() >= cfg.threshold &&
        section.count() <= cfg.noiseThreshold &&
        section.xSize() <= cfg.phiMinBinSize &&
        section.ySize() <= cfg.qOverPtMinBinSize) {
      return Decision::Accept;
    }
    sStat[section.divisionLevel()].area += section.xSize() * section.ySize();
    sStat[section.divisionLevel()].nSections += 1;
    sStat[section.divisionLevel()].nLines += section.count();
    return Decision::Drill;
  };

  exploreHoughParametersSpace(input, measurements, opt, output);
  const unsigned nl =
      std::accumulate(output.begin(), output.end(), 0,
                      [](unsigned sum, const Acts::Experimental::HoughAccumulatorSection &s) {
                        return sum + s.count();
                      });
  ACTS_DEBUG("size " << sStat.size());
  for (const auto &[div, stats] : sStat) {
    ACTS_DEBUG("Area in used by div: "
               << div << " area: " << stats.area << " avglines: "
               << (stats.nSections > 0 ? stats.nLines / stats.nSections : 0)
               << " n sections: " << stats.nSections
               << " thresholdCut: " << stats.discardedByThresholdCut
               << " crossingCut: " << stats.discardedByCrossingCut);
  }
  ACTS_DEBUG("Phi - qOverPt scan finds "
             << output.size() << " solutions, included measurements " << nl);
}

void AdaptiveHoughTransformSeeder::processStackZCotTheta(
    std::vector<Acts::Experimental::HoughAccumulatorSection> &input,
    std::vector<Acts::Experimental::HoughAccumulatorSection> &output,
    const std::vector<PreprocessedMeasurement> &measurements) const {
  ExplorationOptions opt;
  opt.xMinBinSize = m_cfg.zMinBinSize;
  opt.yMinBinSize = m_cfg.cotThetaMinBinSize;
  opt.lineFunctor = m_zCotThetaLineParams;
  opt.decisionFunctor = [&cfg = m_cfg](
                            const Acts::Experimental::HoughAccumulatorSection &section,
                            const std::vector<PreprocessedMeasurement> &) {

    if (section.count() < cfg.threshold) {
      return Decision::Discard;
    }
    if (section.count() >= cfg.threshold &&
        section.xSize() <= cfg.zMinBinSize &&
        section.ySize() <= cfg.cotThetaMinBinSize) {
      return Decision::Accept;
    }
    return Decision::Drill;
  };

  exploreHoughParametersSpace(input, measurements, opt, output);
  const unsigned nl =
      std::accumulate(output.begin(), output.end(), 0,
                      [](unsigned sum, const Acts::Experimental::HoughAccumulatorSection &s) {
                        return sum + s.count();
                      });
  ACTS_DEBUG("Z - cotTheta scan finds " << output.size()
                                        << " included measurements " << nl);
}

void AdaptiveHoughTransformSeeder::processStackZCotThetaSplit(
    std::vector<Acts::Experimental::HoughAccumulatorSection> &input,
    std::vector<Acts::Experimental::HoughAccumulatorSection> &output,
    const std::vector<PreprocessedMeasurement> &measurements) const {
  ExplorationOptions opt;
  opt.xMinBinSize = 101.0f * Acts::UnitConstants::mm;
  opt.yMinBinSize = 10.1f;
  opt.lineFunctor = m_zCotThetaLineParams;
  opt.decisionFunctor = [&cfg = m_cfg, &opt](
                            const Acts::Experimental::HoughAccumulatorSection &section,
                            const std::vector<PreprocessedMeasurement> &) {

    if (section.count() < cfg.threshold) {
      return Decision::Discard;
    }
    if (section.count() >= cfg.threshold &&
        section.xSize() <= opt.xMinBinSize &&
        section.ySize() <= opt.yMinBinSize) {
      return Decision::Accept;
    }
    return Decision::Drill;
  };

  exploreHoughParametersSpace(input, measurements, opt, output);
  const unsigned nl =
      std::accumulate(output.begin(), output.end(), 0,
                      [](unsigned sum, const Acts::Experimental::HoughAccumulatorSection &s) {
                        return sum + s.count();
                      });
  ACTS_DEBUG("Z - cotTheta split scan finds "
             << output.size() << " included measurements " << nl);
}

void AdaptiveHoughTransformSeeder::makeSeeds(
    SeedContainer &seeds,
    const std::vector<Acts::Experimental::HoughAccumulatorSection> &solutions,
    const std::vector<PreprocessedMeasurement> &measurements) const {
  const SpacePointContainer &spacePoints = seeds.spacePointContainer();

  std::size_t seedIndex = 0;
  for (const Acts::Experimental::HoughAccumulatorSection &s : solutions) {
    std::vector<unsigned> sortedIndices = s.indices();
    std::ranges::sort(sortedIndices, [&measurements](unsigned i1, unsigned i2) {
      const auto &m1 = measurements[i1];
      const auto &m2 = measurements[i2];
      return m1.invr > m2.invr;
    });

    if (sortedIndices.size() < 3) {
      ACTS_VERBOSE("this solution has less than 3 SP, ignoring");
      continue;
    }

    unsigned spIndex = 0;
    std::array<std::optional<ConstSpacePointProxy>, 3> sp = {
        std::nullopt, std::nullopt, std::nullopt};
    for (unsigned sidx : sortedIndices) {
      const PreprocessedMeasurement &m = measurements[sidx];
      sp[spIndex] = spacePoints.at(m.sp);
      if (spIndex == 0 || std::abs(sp[spIndex]->r() - sp[spIndex - 1]->r()) >
                              5. * Acts::UnitConstants::mm) {
        spIndex++;
      }
      if (spIndex >= sp.size()) {
        break;
      }
    }

    auto cotThetaEstimate = static_cast<float>((sp[2]->z() - sp[0]->z()) /
                                               (sp[2]->r() - sp[0]->r()));
    auto cotThetaEstimate01 = static_cast<float>((sp[1]->z() - sp[0]->z()) /
                                                 (sp[1]->r() - sp[0]->r()));
    auto cotThetaEstimate12 = static_cast<float>((sp[2]->z() - sp[1]->z()) /
                                                 (sp[2]->r() - sp[1]->r()));
    if (std::abs(cotThetaEstimate01 - cotThetaEstimate12) > 0.1f) {
      ACTS_VERBOSE("Large difference in cotTheta " << cotThetaEstimate01 << " "
                                                   << cotThetaEstimate12);
      continue;
    }

    auto seed = seeds.createSeed();
    seed.assignSpacePointIndices(std::array<SpacePointIndex, 3>{
        sp[0]->index(), sp[1]->index(), sp[2]->index()});

    // for the time the quality is fixed
    // in the future we can use section count for instance
    seed.quality() = 1.0;

    const auto z =
        static_cast<float>(sp[1]->z() - sp[1]->r() * cotThetaEstimate);
    seed.vertexZ() = z;

    ACTS_VERBOSE(seedIndex << ": solution x: " << s.xBegin() << " " << s.xSize()
                           << " y: " << s.yBegin() << " " << s.ySize()
                           << " nlines: " << s.count() << " vertex z: " << z
                           << " r0,1,2: " << sp[0]->r() << " " << sp[1]->r()
                           << " " << sp[2]->r() << " cot,01,12 "
                           << cotThetaEstimate << " " << cotThetaEstimate01
                           << " " << cotThetaEstimate12);

    seedIndex++;
  }
}

bool AdaptiveHoughTransformSeeder::passIntersectionsCheck(
    const Acts::Experimental::HoughAccumulatorSection &section,
    const std::vector<PreprocessedMeasurement> &measurements,
    const LineFunctor &lineFunctor, unsigned threshold) const {
  unsigned inside = 0;
  for (std::size_t idx1 = 0; idx1 < section.count(); ++idx1) {
    const auto &m1 = measurements[section.indices()[idx1]];
    std::function<float(float)> line1 =
        std::bind_front(lineFunctor, std::cref(m1));
    for (std::size_t idx2 = idx1 + 1; idx2 < section.count(); ++idx2) {
      const auto &m2 = measurements[section.indices()[idx2]];
      std::function<float(float)> line2 =
          std::bind_front(lineFunctor, std::cref(m2));
      if (section.isCrossingInside(line1, line2)) {
        inside++;
        if (inside >= threshold) {
          return true;
        }
      }
    }
  }
  ACTS_VERBOSE("Number of crossings inside of section " << inside);
  return inside >= threshold;
}

void AdaptiveHoughTransformSeeder::deduplicate(
    std::vector<Acts::Experimental::HoughAccumulatorSection> &input) const {
  std::vector<const Acts::Experimental::HoughAccumulatorSection *> op;
  op.reserve(input.size());
  std::transform(input.begin(), input.end(), std::back_inserter(op),
                 [](const Acts::Experimental::HoughAccumulatorSection &s) { return &s; });

  auto binaryPredSort = [](const Acts::Experimental::HoughAccumulatorSection *a,
                           const Acts::Experimental::HoughAccumulatorSection *b) {
    return a->indices() < b->indices();
  };
  auto binaryPredUnique = [](const Acts::Experimental::HoughAccumulatorSection *a,
                             const Acts::Experimental::HoughAccumulatorSection *b) {
    return a->indices() == b->indices();
  };

  std::ranges::sort(op, binaryPredSort);
  auto [rbegin, rend] = std::ranges::unique(op, binaryPredUnique);
  std::vector<Acts::Experimental::HoughAccumulatorSection> temp;
  for (auto sPtr = std::begin(op); sPtr != rbegin; ++sPtr) {
    temp.push_back(**sPtr);
  }
  input = std::move(temp);
}

}  // namespace ActsExamples
