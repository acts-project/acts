// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/AdaptiveHoughTransformSeeder.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/ScopedTimer.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Utilities/GroupBy.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>
#include <ostream>
#include <stack>
#include <stdexcept>
#include <variant>

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
  m_xBegin = m_xBegin + 0.5 * m_xSize - m_xSize * xs * 0.5;
  m_yBegin = m_yBegin + 0.5 * m_ySize - m_ySize * ys * 0.5;
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

float AccumulatorSection::distCC(float a, float b) const {
  const float y = fma(a, (m_xBegin + m_xSize), b);
  const float yEnd = m_yBegin + m_ySize;
  return y <= yEnd ? (m_xSize + (yEnd - y)) : ((yEnd - b) / a - m_xBegin);
}
float AccumulatorSection::distACC(float a, float b) const {
  const float y = fma(a, m_xBegin, b);
  return y <= m_yBegin ? (m_ySize + (m_yBegin - b) / a - m_xBegin)
                       : (m_yBegin + m_ySize - y);
}

AdaptiveHoughTransformSeeder::AdaptiveHoughTransformSeeder(
    ActsExamples::AdaptiveHoughTransformSeeder::Config cfg,
    Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("AdaptiveHoughTransformSeeder", lvl),
      m_cfg(std::move(cfg)),
      m_logger(Acts::getDefaultLogger("AdaptiveHoughTransformSeeder", lvl)) {
  for (const auto &spName : config().inputSpacePoints) {
    auto &handle = m_inputSpacePoints.emplace_back(
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

ProcessCode AdaptiveHoughTransformSeeder::execute(
    const AlgorithmContext &ctx) const {
  // get inputs
  std::vector<PreprocessedMeasurement> measurements;
  for (const auto &isp : m_inputSpacePoints) {
    const auto &spContainer = (*isp)(ctx);
    ACTS_DEBUG("Inserting " << spContainer.size()
                            << " space points from collection\"" << isp->key()
                            << "\"");
    for (const SimSpacePoint &sp : spContainer) {
      const double phi = std::atan2(sp.y(), sp.x());
      const double invr = 1.0 / sp.r();
      measurements.push_back(PreprocessedMeasurement(
          invr, phi, sp.z(),
          Acts::SourceLink(static_cast<const SimSpacePoint *>(&sp))));
      // wrap phi by duplicating some seeds
      if (phi < -std::numbers::pi + config().phiWrap * std::numbers::pi) {
        measurements.push_back(PreprocessedMeasurement(
            invr, phi + 2 * std::numbers::pi, sp.z(),
            Acts::SourceLink(static_cast<const SimSpacePoint *>(&sp))));
      }
      if (phi > std::numbers::pi - config().phiWrap * std::numbers::pi) {
        measurements.push_back(PreprocessedMeasurement(
            invr, phi - 2 * std::numbers::pi, sp.z(),
            Acts::SourceLink(static_cast<const SimSpacePoint *>(&sp))));
      }
    }
  }
  ACTS_DEBUG("Total number of " << measurements.size()
                                << " used to account for phi wrapping");

  const unsigned phiSplitMinIndex = 0;
  const unsigned phiSplitWidthIndex = 1;
  const unsigned zSplitMinIndex = 2;
  const unsigned zSplitWidthIndex = 3;
  const unsigned cotThetaSplitMinIndex = 4;
  const unsigned cotThetaSplitWidthIndex = 5;

  // prepare initial stack
  // will become more complicated with slicing
  std::deque<AccumulatorSection> stack1;
  {
    Acts::ScopedTimer st("splitInQuadrants", logger());
    const int nSplits = 8;
    float wedgeWidth = 2.0 * std::numbers::pi / nSplits;

    for (int phiIndex = 0; phiIndex < nSplits; phiIndex++) {
      const float startPhi = wedgeWidth * phiIndex - std::numbers::pi;
      stack1.push_back(AccumulatorSection(1.05 * wedgeWidth,
                                          2.0 * config().qOverPtMin, startPhi,
                                          -config().qOverPtMin));
      stack1.back().indices().resize(measurements.size());
      std::iota(std::begin(stack1.back().indices()),
                std::end(stack1.back().indices()), 0);
      stack1.back().setHistory(phiSplitMinIndex, startPhi);
      stack1.back().setHistory(phiSplitWidthIndex, wedgeWidth);

      updateSection(stack1.back(), measurements, m_qOverPtPhiLineParams);
      ACTS_DEBUG("Initial split " << stack1.back().count()
                                  << " for region starting from " << startPhi);
      for (unsigned index : stack1.back().indices()) {
        const PreprocessedMeasurement &m = measurements[index];
        ACTS_DEBUG("phi section "
                   << startPhi << " x: " << std::cos(m.phi) / m.invr
                   << " y: " << std::sin(m.phi) / m.invr << " z: " << m.z);
      }
    }
  }

  // split into regions in z_vertex cot theta, there is a lot of duplication and
  // it can be optimized better in the future
  for (auto &section : stack1) {
    section.updateDimensions(2.0 * config().zRange, 2. * config().cotThetaRange,
                             -config().zRange, -config().cotThetaRange);
  }
  std::deque<AccumulatorSection> stack2;
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
        section.history(phiSplitWidthIndex), 2.0 * config().qOverPtMin,
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

  std::deque<AccumulatorSection> &solutions =
      config().doSecondPhase ? stack2 : stack1;
  Acts::ScopedTimer st("seedsMaking", logger());

  // post solutions
  SimSeedContainer seeds;
  seeds.reserve(solutions.size());

  ProtoTrackContainer protoTracks;
  protoTracks.reserve(solutions.size());
  ACTS_VERBOSE("Number of solutions " << solutions.size());
  std::size_t seedIndex = 0;
  for (const AccumulatorSection &s : solutions) {
    unsigned spIndex = 0;
    std::array<const SimSpacePoint *, 3> sp = {nullptr, nullptr, nullptr};
    std::vector<unsigned> sortedIndices = s.indices();
    std::sort(sortedIndices.begin(), sortedIndices.end(),
              [&measurements](unsigned i1, unsigned i2) {
                auto &m1 = measurements[i1];
                auto &m2 = measurements[i2];
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

    float cotThetaEstimate =
        (sp[2]->z() - sp[0]->z()) / (sp[2]->r() - sp[0]->r());
    float cotThetaEstimate01 =
        (sp[1]->z() - sp[0]->z()) / (sp[1]->r() - sp[0]->r());
    float cotThetaEstimate12 =
        (sp[2]->z() - sp[1]->z()) / (sp[2]->r() - sp[1]->r());
    if (std::abs(cotThetaEstimate01 - cotThetaEstimate12) > 0.1) {
      ACTS_VERBOSE("Large difference in cotTheta " << cotThetaEstimate01 << " "
                                                   << cotThetaEstimate12);
      continue;
    }
    seeds.emplace_back(*sp[0], *sp[1], *sp[2]);
    float z = sp[1]->z() - sp[1]->r() * cotThetaEstimate;
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

  m_outputSeeds(ctx, std::move(seeds));

  return ActsExamples::ProcessCode::SUCCESS;
}

void AdaptiveHoughTransformSeeder::processStackQOverPtPhi(
    std::deque<AccumulatorSection> &input,
    std::deque<AccumulatorSection> &output,
    const std::vector<PreprocessedMeasurement> &measurements) const {
  struct Stats {
    double area{};
    int nSections{};
    int nLines{};
    int discardedByThresholdCut{};
    int discardedByCrossingCut{};
  };

  std::map<int, Stats> sStat;
  AHTExplorationOptions opt;
  opt.xMinBinSize = config().phiMinBinSize;
  opt.yMinBinSize = config().qOverPtMinBinSize;
  opt.lineParamFunctor = m_qOverPtPhiLineParams;
  opt.decisionFunctor = [&sStat, this](
                            const AccumulatorSection &section,
                            const std::vector<PreprocessedMeasurement> &mes) {
    if (section.divisionLevel() <= 8) {
      return AHTExplorationOptions<PreprocessedMeasurement>::Drill;
    }

    if (section.count() < this->m_cfg.threshold) {
      sStat[section.divisionLevel()].discardedByThresholdCut += 1;
      return AHTExplorationOptions<PreprocessedMeasurement>::Discard;
    }
    if (section.count() < 3 * this->m_cfg.threshold) {
      if (!passIntersectionsCheck(
              section, mes, m_qOverPtPhiLineParams,
              this->m_cfg.threshold * (this->m_cfg.threshold - 1))) {
        sStat[section.divisionLevel()].discardedByCrossingCut += 1;
        return AHTExplorationOptions<PreprocessedMeasurement>::Discard;
      }
    }
    // if (section.count() == m_cfg.threshold) {

    // }

    if (section.count() >= this->m_cfg.threshold &&
        section.count() <= this->m_cfg.noiseThreshold &&
        section.xSize() <= this->m_cfg.phiMinBinSize &&
        section.ySize() <= this->m_cfg.qOverPtMinBinSize) {
      return AHTExplorationOptions<PreprocessedMeasurement>::Accept;
    }
    sStat[section.divisionLevel()].area += section.xSize() * section.ySize();
    sStat[section.divisionLevel()].nSections += 1;
    sStat[section.divisionLevel()].nLines += section.count();
    return AHTExplorationOptions<PreprocessedMeasurement>::Drill;
  };

  exploreParametersSpace(input, measurements, opt, output);
  const unsigned nl =
      std::accumulate(output.begin(), output.end(), 0,
                      [](unsigned sum, const AccumulatorSection &s) {
                        return sum + s.count();
                      });
  ACTS_DEBUG("size " << sStat.size());
  for (const auto &data : sStat) {
    ACTS_DEBUG("Area in used by div: "
               << data.first << " area: " << data.second.area << " avglines: "
               << (data.second.nSections > 0
                       ? data.second.nLines / data.second.nSections
                       : 0)
               << " n sections: " << data.second.nSections
               << " thresholdCut: " << data.second.discardedByThresholdCut
               << " crossingCut: " << data.second.discardedByCrossingCut);
  }
  ACTS_DEBUG("Phi - qOverPt scan finds "
             << output.size() << " solutions, included measurements " << nl);
}

void AdaptiveHoughTransformSeeder::processStackZCotTheta(
    std::deque<AccumulatorSection> &input,
    std::deque<AccumulatorSection> &output,
    const std::vector<PreprocessedMeasurement> &measurements) const {
  AHTExplorationOptions opt;
  opt.xMinBinSize = config().zMinBinSize;
  opt.yMinBinSize = config().cotThetaMinBinSize;
  opt.lineParamFunctor = m_zCotThetaLineParams;
  opt.decisionFunctor = [&m_cfg = m_cfg](
                            const AccumulatorSection &section,
                            const std::vector<PreprocessedMeasurement> &) {
    if (section.count() < m_cfg.threshold) {
      return AHTExplorationOptions<PreprocessedMeasurement>::Discard;
    }
    if (section.count() >= m_cfg.threshold &&
        section.xSize() <= m_cfg.zMinBinSize &&
        section.ySize() <= m_cfg.cotThetaMinBinSize) {
      return AHTExplorationOptions<PreprocessedMeasurement>::Accept;
    }
    return AHTExplorationOptions<PreprocessedMeasurement>::Drill;
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
    std::deque<AccumulatorSection> &input,
    std::deque<AccumulatorSection> &output,
    const std::vector<PreprocessedMeasurement> &measurements) const {
  AHTExplorationOptions opt;
  opt.xMinBinSize = 101 * Acts::UnitConstants::mm;
  opt.yMinBinSize = 10.1;
  opt.lineParamFunctor = m_zCotThetaLineParams;
  opt.decisionFunctor = [&m_cfg = m_cfg, &opt](
                            const AccumulatorSection &section,
                            const std::vector<PreprocessedMeasurement> &) {
    if (section.count() < m_cfg.threshold) {
      return AHTExplorationOptions<PreprocessedMeasurement>::Discard;
    }
    if (section.count() >= m_cfg.threshold &&
        section.xSize() <= opt.xMinBinSize &&
        section.ySize() <= opt.yMinBinSize) {
      return AHTExplorationOptions<PreprocessedMeasurement>::Accept;
    }
    return AHTExplorationOptions<PreprocessedMeasurement>::Drill;
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

bool AdaptiveHoughTransformSeeder::passIntersectionsCheck(
    const AccumulatorSection &section,
    const std::vector<PreprocessedMeasurement> &measurements,
    const LineParamFunctor &lineFunctor, unsigned threshold) const {
  using namespace std::placeholders;
  unsigned inside = 0;
  for (std::size_t idx1 = 0; idx1 < section.count(); ++idx1) {
    const auto &m1 = measurements[section.indices()[idx1]];
    std::function<float(float)> line1 = std::bind(lineFunctor, m1, _1);
    for (std::size_t idx2 = idx1 + 1; idx2 < section.count(); ++idx2) {
      const auto &m2 = measurements[section.indices()[idx2]];
      std::function<float(float)> line2 = std::bind(lineFunctor, m2, _1);
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
    std::deque<AccumulatorSection> &input) const {
  std::vector<const AccumulatorSection *> op;
  op.reserve(input.size());
  for (AccumulatorSection &s : input) {
    op.push_back(&s);
  }

  auto binaryPredSort = [](const AccumulatorSection *a,
                           const AccumulatorSection *b) {
    return a->indices() < b->indices();
  };
  auto binaryPredUnique = [](const AccumulatorSection *a,
                             const AccumulatorSection *b) {
    return a->indices() == b->indices();
  };

  std::sort(std::begin(op), std::end(op), binaryPredSort);
  auto endOfUniqueSections =
      std::unique(std::begin(op), std::end(op), binaryPredUnique);
  std::deque<AccumulatorSection> temp;
  for (auto sPtr = std::begin(op); sPtr != endOfUniqueSections; ++sPtr) {
    temp.push_back(**sPtr);
  }
  input = std::move(temp);
}

}  // namespace ActsExamples
