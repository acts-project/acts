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
AccumulatorSection::AccumulatorSection(double xw, double yw, double xBegin,
                                       double yBegin, int div,
                                       const std::vector<unsigned> &indices)
    : m_xSize(xw),
      m_ySize(yw),
      m_xBegin(xBegin),
      m_yBegin(yBegin),
      m_divisionLevel(div),
      m_indices(indices) {}

AccumulatorSection AccumulatorSection::bottomLeft(float xFraction,
                                                  float yFraction) const {
  return AccumulatorSection(m_xSize * xFraction, m_ySize * yFraction, m_xBegin,
                            m_yBegin, m_divisionLevel + 1, m_indices);
}
AccumulatorSection AccumulatorSection::topLeft(float xFraction,
                                               float yFraction) const {
  return AccumulatorSection(m_xSize * xFraction, m_ySize * yFraction, m_xBegin,
                            m_yBegin + m_ySize - m_ySize * yFraction,
                            m_divisionLevel + 1, m_indices);
}
AccumulatorSection AccumulatorSection::topRight(float xFraction,
                                                float yFraction) const {
  return AccumulatorSection(m_xSize * xFraction, m_ySize * yFraction,
                            m_xBegin + m_xSize - m_xSize * xFraction,
                            m_yBegin + m_ySize - m_ySize * yFraction,
                            m_divisionLevel + 1, m_indices);
}
AccumulatorSection AccumulatorSection::bottomRight(float xFraction,
                                                   float yFraction) const {
  return AccumulatorSection(m_xSize * xFraction, m_ySize * yFraction,
                            m_xBegin + m_xSize - m_xSize * xFraction, m_yBegin,
                            m_divisionLevel + 1, m_indices);
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
    ACTS_DEBUG("Inserting " << spContainer.size() << " space points from "
                            << isp->key());
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
  ACTS_DEBUG("Collected " << measurements.size() << " space points");

  // prepare initial stack
  // will become more complicated with slicing
  std::stack<AccumulatorSection> sectionsStack(
      {AccumulatorSection(2. * std::numbers::pi, 2.0 * config().qOverPtMin,
                          -std::numbers::pi, -config().qOverPtMin)});

  // all the measurements should be considered at start
  // therefore indices of all of them are stored
  // additional loop over regions will wrap this
  sectionsStack.top().indices().resize(measurements.size());
  std::iota(std::begin(sectionsStack.top().indices()),
            std::end(sectionsStack.top().indices()), 0);

  std::vector<AccumulatorSection> solutions;
  processStackHeadQOverPtPhi(sectionsStack, solutions, measurements);

  if (config().doSecondPhase) {
    for (auto &s : solutions) {
      sectionsStack.push(AccumulatorSection(
          2. * config().zRange, 2. * config().cotThetaRange, -config().zRange,
          -config().cotThetaRange, s.divisionLevel(), s.indices()));
    }
    solutions.clear();
    processStackHeadZCotTheta(sectionsStack, solutions, measurements);
  }

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
      if (spIndex >= sp.size())
        break;
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

void AdaptiveHoughTransformSeeder::processStackHeadQOverPtPhi(
    std::stack<AccumulatorSection> &sections,
    std::vector<AccumulatorSection> &solutions,
    const std::vector<PreprocessedMeasurement> &measurements) const {
  AHTExplorationOptions opt;
  opt.xMinBinSize = config().phiMinBinSize;
  opt.yMinBinSize = config().qOverPtMinBinSize;
  opt.lineParamFunctors = m_qOverPtPhiLineParams;
  opt.decisionFunctor = [&m_cfg = m_cfg, &opt, this](
                            const AccumulatorSection &section,
                            const std::vector<PreprocessedMeasurement> &mes) {
    if (section.count() < m_cfg.threshold) {
      return AHTExplorationOptions<PreprocessedMeasurement>::Discard;
    }
    if (!passIntersectionsCheck(section, mes, m_qOverPtPhiLineParams,
                                m_cfg.threshold + 1)) {
      return AHTExplorationOptions<PreprocessedMeasurement>::Discard;
    }

    // count lines skipping first pixel layers, only when we are deep enough
    if (section.count() < 2 * m_cfg.threshold) {
      std::size_t countNoPix = 0;
      for (unsigned i : section.indices()) {
        const PreprocessedMeasurement &m = mes[i];
        if (m.invr < 1. / 300.)
          countNoPix++;
      }
      // this threshold could be made configurable as well
      if (countNoPix < m_cfg.threshold - 2) {
        return AHTExplorationOptions<PreprocessedMeasurement>::Discard;
      }
    }
    if (section.count() >= m_cfg.threshold &&
        section.xSize() <= m_cfg.phiMinBinSize &&
        section.ySize() <= m_cfg.qOverPtMinBinSize) {
      return AHTExplorationOptions<PreprocessedMeasurement>::Accept;
    }
    return AHTExplorationOptions<PreprocessedMeasurement>::Drill;
  };

  exploreParametersSpace(sections, measurements, opt, solutions);
  const unsigned nl =
      std::accumulate(solutions.begin(), solutions.end(), 0,
                      [](unsigned sum, const AccumulatorSection &s) {
                        return sum + s.count();
                      });
  ACTS_DEBUG("Phi - qOverPt scan finds " << solutions.size()
                                         << " included measurements " << nl);
}

void AdaptiveHoughTransformSeeder::processStackHeadZCotTheta(
    std::stack<AccumulatorSection> &sections,
    std::vector<AccumulatorSection> &solutions,
    const std::vector<PreprocessedMeasurement> &measurements) const {
  AHTExplorationOptions opt;
  opt.xMinBinSize = config().zMinBinSize;
  opt.yMinBinSize = config().cotThetaMinBinSize;
  opt.lineParamFunctors = m_zCotThetaLineParams;
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

  exploreParametersSpace(sections, measurements, opt, solutions);
  const unsigned nl =
      std::accumulate(solutions.begin(), solutions.end(), 0,
                      [](unsigned sum, const AccumulatorSection &s) {
                        return sum + s.count();
                      });
  ACTS_DEBUG("Z - cotTheta scan finds " << solutions.size()
                                        << " included measurements " << nl);
}

bool AdaptiveHoughTransformSeeder::passIntersectionsCheck(
    const AccumulatorSection &section,
    const std::vector<PreprocessedMeasurement> &measurements,
    const LineParamFunctors &lineParams, unsigned threshold) const {
  unsigned inside = 0;
  for (std::size_t idx1 = 0; idx1 < section.count(); ++idx1) {
    const auto &m1 = measurements[section.indices()[idx1]];
    const double a1 = lineParams.first(m1);
    const double b1 = lineParams.second(m1);
    for (std::size_t idx2 = idx1 + 1; idx2 < section.count(); ++idx2) {
      const auto &m2 = measurements[section.indices()[idx2]];
      const double a2 = lineParams.first(m2);
      const double b2 = lineParams.second(m2);
      if (section.isCrossingInside(a1, b1, a2, b2)) {
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

void AdaptiveHoughTransformSeeder::addSolution(
    AccumulatorSection &&s, std::vector<AccumulatorSection> &output) const {
  if (config().deduplicate) {
    std::size_t i = output.size() > 10 ? output.size() - 10 : 0;
    for (; i < std::size(output); ++i) {
      // compare first 4 measurements
      const std::size_t firstMeasurementsToCompare = 4;
      if (std::equal(output[i].indices().begin(),
                     output[i].indices().begin() + firstMeasurementsToCompare,
                     s.indices().begin(),
                     s.indices().begin() + firstMeasurementsToCompare)) {
        ACTS_VERBOSE(
            "There is a nearby duplicate solution. Skipping this one.");
        return;
      }
    }
  }
  ACTS_VERBOSE("solution section "
               << s.count() << " section " << s.xBegin() << " - "
               << s.xBegin() + s.xSize() << " " << s.yBegin() << " - "
               << s.yBegin() + s.ySize() << " nlines: " << s.count()
               << " div: " << s.divisionLevel());

  output.push_back(std::move(s));
}

}  // namespace ActsExamples
