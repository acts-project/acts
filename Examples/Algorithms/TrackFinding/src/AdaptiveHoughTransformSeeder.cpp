// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/AdaptiveHoughTransformSeeder.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/ScopedTimer.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <map>
#include <numeric>
#include <ostream>
#include <stdexcept>

namespace ActsExamples {

// Helper class describing one section of the accumulator space
AccumulatorSection::AccumulatorSection(float xw, float yw, float xBegin,
                                       float yBegin, int div,
                                       const std::vector<unsigned> &indices,
                                       const std::vector<float> &history)
    : m_xSize(xw),
      m_ySize(yw),
      m_xBegin(xBegin),
      m_yBegin(yBegin),
      m_divisionLevel(div),
      m_indices(indices),
      m_history(history) {}

void AccumulatorSection::updateDimensions(float xw, float yw, float xBegin,
                                          float yBegin) {
  m_xSize = xw;
  m_ySize = yw;
  m_xBegin = xBegin;
  m_yBegin = yBegin;
}

void AccumulatorSection::expand(float xs, float ys) {
  m_xBegin = m_xBegin + 0.5f * m_xSize - m_xSize * xs * 0.5f;
  m_yBegin = m_yBegin + 0.5f * m_ySize - m_ySize * ys * 0.5f;
  m_xSize *= xs;
  m_ySize *= ys;
}

AccumulatorSection AccumulatorSection::bottomLeft(float xFraction,
                                                  float yFraction) const {
  return AccumulatorSection(m_xSize * xFraction, m_ySize * yFraction, m_xBegin,
                            m_yBegin, m_divisionLevel + 1, m_indices,
                            m_history);
}
AccumulatorSection AccumulatorSection::topLeft(float xFraction,
                                               float yFraction) const {
  return AccumulatorSection(m_xSize * xFraction, m_ySize * yFraction, m_xBegin,
                            m_yBegin + m_ySize - m_ySize * yFraction,
                            m_divisionLevel + 1, m_indices, m_history);
}
AccumulatorSection AccumulatorSection::topRight(float xFraction,
                                                float yFraction) const {
  return AccumulatorSection(m_xSize * xFraction, m_ySize * yFraction,
                            m_xBegin + m_xSize - m_xSize * xFraction,
                            m_yBegin + m_ySize - m_ySize * yFraction,
                            m_divisionLevel + 1, m_indices, m_history);
}
AccumulatorSection AccumulatorSection::bottomRight(float xFraction,
                                                   float yFraction) const {
  return AccumulatorSection(m_xSize * xFraction, m_ySize * yFraction,
                            m_xBegin + m_xSize - m_xSize * xFraction, m_yBegin,
                            m_divisionLevel + 1, m_indices, m_history);
}
AccumulatorSection AccumulatorSection::bottom(float yFraction) const {
  return bottomLeft(1.0, yFraction);
}
AccumulatorSection AccumulatorSection::top(float yFraction) const {
  return topLeft(1.0, yFraction);
}
AccumulatorSection AccumulatorSection::left(float xFraction) const {
  return bottomLeft(xFraction, 1.0);
}
AccumulatorSection AccumulatorSection::right(float xFraction) const {
  return bottomRight(xFraction, 1.0);
}

AdaptiveHoughTransformSeeder::AdaptiveHoughTransformSeeder(
    const Config &cfg, Acts::Logging::Level lvl)
    : IAlgorithm("AdaptiveHoughTransformSeeder", lvl),
      m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("AdaptiveHoughTransformSeeder", lvl)) {
  for (const auto &spName : config().inputSpacePoints) {
    const auto &handle = m_inputSpacePoints.emplace_back(
        std::make_unique<ReadDataHandle<SimSpacePointContainer>>(
            this,
            "InputSpacePoints#" + std::to_string(m_inputSpacePoints.size())));
    handle->initialize(spName);
  }
  if (config().outputSeeds.empty()) {
    throw std::invalid_argument(
        "AdaptiveHoughTransformSeeder: Output seeds collection name empty");
  }
  m_outputSeeds.initialize(config().outputSeeds);
}

// When traversing parameters space in various directions it is useful
// to record the "path".
// AccumulatorSection has simple functionality allowing that in a form of a
// plain vector of floats. These numbers need to be accessed consistently, thus
// this indices.
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
  // get inputs
  std::vector<PreprocessedMeasurement> measurements;
  preparePreprocessedMeasurements(ctx, measurements);

  // prepare initial stack
  std::vector<AccumulatorSection> stack1;
  fillStackPhiSplit(stack1, measurements);

  // split into regions in z_vertex cot theta, there is a lot of duplication and
  // it can be optimized better in the future
  for (auto &section : stack1) {
    section.updateDimensions(2.0f * config().zRange,
                             2.0f * config().cotThetaRange, -config().zRange,
                             -config().cotThetaRange);
  }
  std::vector<AccumulatorSection> stack2;
  {
    Acts::ScopedTimer st("splitInZCotTheta", logger());
    processStackZCotThetaSplit(stack1, stack2, measurements);
    if (config().doSecondPhase) {
      for (AccumulatorSection &section : stack2) {
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
        section.history(phiSplitWidthIndex), 2.0f * config().qOverPtMin,
        section.history(phiSplitMinIndex), -config().qOverPtMin);
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
  if (config().doSecondPhase) {
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

  std::vector<AccumulatorSection> &solutions =
      config().doSecondPhase ? stack2 : stack1;
  Acts::ScopedTimer st("seedsMaking", logger());

  // post solutions
  SimSeedContainer seeds;
  seeds.reserve(solutions.size());

  ACTS_VERBOSE("Number of solutions " << solutions.size());
  makeSeeds(seeds, solutions, measurements);

  m_outputSeeds(ctx, std::move(seeds));

  return ProcessCode::SUCCESS;
}

void AdaptiveHoughTransformSeeder::preparePreprocessedMeasurements(
    const AlgorithmContext &ctx,
    std::vector<PreprocessedMeasurement> &measurements) const {
  for (const auto &isp : m_inputSpacePoints) {
    const auto &spContainer = (*isp)(ctx);
    ACTS_DEBUG("Inserting " << spContainer.size()
                            << " space points from collection\"" << isp->key()
                            << "\"");

    for (const SimSpacePoint &sp : spContainer) {
      const double phi = std::atan2(sp.y(), sp.x());
      const double invr = 1.0 / sp.r();
      measurements.emplace_back(
          invr, phi, sp.z(),
          Acts::SourceLink(static_cast<const SimSpacePoint *>(&sp)));
      // wrap phi by duplicating some seeds
      if (phi < -std::numbers::pi + config().phiWrap * std::numbers::pi) {
        measurements.emplace_back(
            invr, phi + 2 * std::numbers::pi, sp.z(),
            Acts::SourceLink(static_cast<const SimSpacePoint *>(&sp)));
      }
      if (phi > std::numbers::pi - config().phiWrap * std::numbers::pi) {
        measurements.emplace_back(
            invr, phi - 2 * std::numbers::pi, sp.z(),
            Acts::SourceLink(static_cast<const SimSpacePoint *>(&sp)));
      }
    }
  }
  ACTS_DEBUG("Total number of " << measurements.size()
                                << " used to account for phi wrapping");
}

void AdaptiveHoughTransformSeeder::fillStackPhiSplit(
    std::vector<AccumulatorSection> &stack,
    const std::vector<PreprocessedMeasurement> &measurements) const {
  Acts::ScopedTimer st("splitInQuadrants", logger());
  const int nSplits = 8;
  const auto wedgeWidth = static_cast<float>(2.0 * std::numbers::pi / nSplits);

  for (int phiIndex = 0; phiIndex < nSplits; phiIndex++) {
    const auto startPhi = static_cast<float>(
        wedgeWidth * static_cast<float>(phiIndex) - std::numbers::pi);
    stack.emplace_back(1.2f * wedgeWidth, 2.0f * config().qOverPtMin, startPhi,
                       -config().qOverPtMin);
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
    std::vector<AccumulatorSection> &input,
    std::vector<AccumulatorSection> &output,
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
  opt.xMinBinSize = config().phiMinBinSize;
  opt.yMinBinSize = config().qOverPtMinBinSize;
  opt.lineParamFunctor = m_qOverPtPhiLineParams;
  opt.decisionFunctor = [&sStat, &cfg = m_cfg, &opt, this](
                            const AccumulatorSection &section,
                            const std::vector<PreprocessedMeasurement> &mes) {
    using enum ExplorationOptions<PreprocessedMeasurement>::Decision;
    if (section.divisionLevel() <= 8) {
      return Drill;
    }

    if (section.count() < cfg.threshold) {
      sStat[section.divisionLevel()].discardedByThresholdCut += 1;
      return Discard;
    }
    if (section.count() < 3 * cfg.threshold &&
        !passIntersectionsCheck(section, mes, opt.lineParamFunctor,
                                cfg.threshold * (cfg.threshold - 1))) {
      sStat[section.divisionLevel()].discardedByCrossingCut += 1;
      return Discard;
    }

    if (section.count() >= cfg.threshold &&
        section.count() <= cfg.noiseThreshold &&
        section.xSize() <= cfg.phiMinBinSize &&
        section.ySize() <= cfg.qOverPtMinBinSize) {
      return Accept;
    }
    sStat[section.divisionLevel()].area += section.xSize() * section.ySize();
    sStat[section.divisionLevel()].nSections += 1;
    sStat[section.divisionLevel()].nLines += section.count();
    return Drill;
  };

  exploreParametersSpace(input, measurements, opt, output);
  const unsigned nl =
      std::accumulate(output.begin(), output.end(), 0,
                      [](unsigned sum, const AccumulatorSection &s) {
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
    std::vector<AccumulatorSection> &input,
    std::vector<AccumulatorSection> &output,
    const std::vector<PreprocessedMeasurement> &measurements) const {
  ExplorationOptions opt;
  opt.xMinBinSize = config().zMinBinSize;
  opt.yMinBinSize = config().cotThetaMinBinSize;
  opt.lineParamFunctor = m_zCotThetaLineParams;
  opt.decisionFunctor = [&cfg = m_cfg](
                            const AccumulatorSection &section,
                            const std::vector<PreprocessedMeasurement> &) {
    using enum ExplorationOptions<PreprocessedMeasurement>::Decision;

    if (section.count() < cfg.threshold) {
      return Discard;
    }
    if (section.count() >= cfg.threshold &&
        section.xSize() <= cfg.zMinBinSize &&
        section.ySize() <= cfg.cotThetaMinBinSize) {
      return Accept;
    }
    return Drill;
  };

  exploreParametersSpace(input, measurements, opt, output);
  const unsigned nl =
      std::accumulate(output.begin(), output.end(), 0,
                      [](unsigned sum, const AccumulatorSection &s) {
                        return sum + s.count();
                      });
  ACTS_DEBUG("Z - cotTheta scan finds " << output.size()
                                        << " included measurements " << nl);
}

void AdaptiveHoughTransformSeeder::processStackZCotThetaSplit(
    std::vector<AccumulatorSection> &input,
    std::vector<AccumulatorSection> &output,
    const std::vector<PreprocessedMeasurement> &measurements) const {
  ExplorationOptions opt;
  opt.xMinBinSize = 101.0f * Acts::UnitConstants::mm;
  opt.yMinBinSize = 10.1f;
  opt.lineParamFunctor = m_zCotThetaLineParams;
  opt.decisionFunctor = [&cfg = m_cfg, &opt](
                            const AccumulatorSection &section,
                            const std::vector<PreprocessedMeasurement> &) {
    using enum ExplorationOptions<PreprocessedMeasurement>::Decision;
    if (section.count() < cfg.threshold) {
      return Discard;
    }
    if (section.count() >= cfg.threshold &&
        section.xSize() <= opt.xMinBinSize &&
        section.ySize() <= opt.yMinBinSize) {
      return Accept;
    }
    return Drill;
  };

  exploreParametersSpace(input, measurements, opt, output);
  const unsigned nl =
      std::accumulate(output.begin(), output.end(), 0,
                      [](unsigned sum, const AccumulatorSection &s) {
                        return sum + s.count();
                      });
  ACTS_DEBUG("Z - cotTheta split scan finds "
             << output.size() << " included measurements " << nl);
}

void AdaptiveHoughTransformSeeder::makeSeeds(
    SimSeedContainer &seeds, const std::vector<AccumulatorSection> &solutions,
    const std::vector<PreprocessedMeasurement> &measurements) const {
  std::size_t seedIndex = 0;
  for (const AccumulatorSection &s : solutions) {
    unsigned spIndex = 0;
    std::array<const SimSpacePoint *, 3> sp = {nullptr, nullptr, nullptr};
    std::vector<unsigned> sortedIndices = s.indices();
    std::ranges::sort(sortedIndices, [&measurements](unsigned i1, unsigned i2) {
      const auto &m1 = measurements[i1];
      const auto &m2 = measurements[i2];
      return m1.invr > m2.invr;
    });

    for (unsigned sidx : sortedIndices) {
      const PreprocessedMeasurement &m = measurements[sidx];
      sp[spIndex] = m.link.get<const SimSpacePoint *>();
      if (spIndex == 0 || std::abs(sp[spIndex]->r() - sp[spIndex - 1]->r()) >
                              5. * Acts::UnitConstants::mm) {
        spIndex++;
      }
      if (spIndex >= sp.size()) {
        break;
      }
    }
    if (sp[2] == nullptr) {
      ACTS_VERBOSE("this solution has less than 3 SP, ignoring");
      continue;
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
    seeds.emplace_back(*sp[0], *sp[1], *sp[2]);
    auto z = static_cast<float>(sp[1]->z() - sp[1]->r() * cotThetaEstimate);
    seeds.back().setVertexZ(z);

    ACTS_VERBOSE(seedIndex << ": solution x: " << s.xBegin() << " " << s.xSize()
                           << " y: " << s.yBegin() << " " << s.ySize()
                           << " nlines: " << s.count() << " vertex z: " << z
                           << " r0,1,2: " << sp[0]->r() << " " << sp[1]->r()
                           << " " << sp[2]->r() << " cot,01,12 "
                           << cotThetaEstimate << " " << cotThetaEstimate01
                           << " " << cotThetaEstimate12);

    // for the time the quality is fixed
    // in the future we can use section count for instance
    seeds.back().setQuality(1.0);
    seedIndex++;
  }
}

bool AdaptiveHoughTransformSeeder::passIntersectionsCheck(
    const AccumulatorSection &section,
    const std::vector<PreprocessedMeasurement> &measurements,
    const LineParamFunctor &lineFunctor, unsigned threshold) const {
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
    std::vector<AccumulatorSection> &input) const {
  std::vector<const AccumulatorSection *> op;
  op.reserve(input.size());
  std::transform(input.begin(), input.end(), std::back_inserter(op),
                 [](const AccumulatorSection &s) { return &s; });

  auto binaryPredSort = [](const AccumulatorSection *a,
                           const AccumulatorSection *b) {
    return a->indices() < b->indices();
  };
  auto binaryPredUnique = [](const AccumulatorSection *a,
                             const AccumulatorSection *b) {
    return a->indices() == b->indices();
  };

  std::ranges::sort(op, binaryPredSort);
  auto [rbegin, rend] = std::ranges::unique(op, binaryPredUnique);
  std::vector<AccumulatorSection> temp;
  for (auto sPtr = std::begin(op); sPtr != rbegin; ++sPtr) {
    temp.push_back(**sPtr);
  }
  input = std::move(temp);
}

}  // namespace ActsExamples
