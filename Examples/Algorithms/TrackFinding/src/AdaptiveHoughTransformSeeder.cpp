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
#include <fstream>

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
// Acts::Experimental::HoughAccumulatorSection has simple functionality allowing
// that in a form of a plain vector of floats. These numbers need to be accessed
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

  ACTS_VERBOSE("Starting maxDivision value: " << m_cfg.maxDivision);
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
        section.setDivisionLevel(0);
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

  // additional stack for smaller structer as HoughAccumulatorSection is no longer needed
  std::vector<SeedCandidate> stack3;

  // do scan in z_vertex - cot theta
  if (m_cfg.doSecondPhase) {
    Acts::ScopedTimer st("secondPhase", logger());

    for (auto &section : stack1) {
      section.updateDimensions(section.history(zSplitWidthIndex),
                               section.history(cotThetaSplitWidthIndex),
                               section.history(zSplitMinIndex),
                               section.history(cotThetaSplitMinIndex));
    }
    processStackZCotTheta(stack1, stack3, measurements);
    ACTS_DEBUG("past StackZCotTheta stack1 size "
               << stack1.size() << " stack3 size " << stack3.size());
  }

  // vector for final solutions with or without second phase
  std::vector<SeedCandidate> finalSolutions;

  if (m_cfg.doSecondPhase){
    finalSolutions = std::move(stack3);
  }
  else{
    finalSolutions.reserve(stack1.size());
    for (auto& section : stack1){
      finalSolutions.push_back({std::move(section.indices())});
    }
  }

  Acts::ScopedTimer st("seedsMaking", logger());

  // post solutions
  SeedContainer seeds;
  seeds.assignSpacePointContainer(spacePoints);
  seeds.reserve(finalSolutions.size());

  ACTS_VERBOSE("Number of solutions " << finalSolutions.size());
  makeSeeds(seeds, finalSolutions, measurements);

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
    double drilledArea{};
    double acceptedArea{};
    int nDrilledSections_standard{};
    int nDrilledSections_extra{};
    int nAcceptedLines{};
    int nAcceptedSections{};
    int discardedByThresholdCut{};
    int discardedByCrossingCut{};
  };
  std::map<int, Stats> sStat;
  ExplorationOptions opt;
  opt.xMinBinSize = m_cfg.phiMinBinSize;
  opt.yMinBinSize = m_cfg.qOverPtMinBinSize;
  opt.lineFunctor = m_qOverPtPhiLineParams;
  opt.decisionFunctor =
      [&sStat, &cfg = m_cfg, &opt,this](
          const Acts::Experimental::HoughAccumulatorSection &section,
          const std::vector<PreprocessedMeasurement> &mes) {
        using enum HoughAccumulatorSection::Decision;

        if (section.count() < cfg.threshold) {
          sStat[section.divisionLevel()].discardedByThresholdCut += 1;
          return Discard;
        }

        if (section.count() < (m_cfg.noiseThreshold) &&
            !passIntersectionsCheck(section, mes, opt.lineFunctor,
                                    cfg.threshold * (cfg.threshold - 1)/2)) {
          sStat[section.divisionLevel()].discardedByCrossingCut ++;
          return Discard;
        }

        if (section.divisionLevel() <= cfg.maxDivision-1) {
          sStat[section.divisionLevel()].nDrilledSections_standard++;
          sStat[section.divisionLevel()].drilledArea += section.xSize() * section.ySize();
          return Drill;
        }

        if (section.count() >= cfg.threshold &&
            section.xSize() <= cfg.phiMinBinSize &&
            section.ySize() <= cfg.qOverPtMinBinSize) {
          sStat[section.divisionLevel()].nAcceptedSections++;
          sStat[section.divisionLevel()].nAcceptedLines+=section.count(); 
          sStat[section.divisionLevel()].acceptedArea += section.xSize() * section.ySize();
          return Accept;
        }
        sStat[section.divisionLevel()].drilledArea +=
            section.xSize() * section.ySize();
        sStat[section.divisionLevel()].nDrilledSections_extra ++;
        return Drill;
      };

  exploreHoughParametersSpace(input, measurements, opt, output);
  const unsigned nl = std::accumulate(
      output.begin(), output.end(), 0,
      [](unsigned sum, const Acts::Experimental::HoughAccumulatorSection &s) {
        return sum + s.count();
      });
  ACTS_DEBUG("Phi - qOverPt scan finds "
             << output.size() << " solutions, included measurements " << nl);
  ACTS_DEBUG("size " << sStat.size());
  ACTS_DEBUG("===========================================================");
  ACTS_DEBUG(" [Q/pT - Phi] process stats ");
  ACTS_DEBUG("===========================================================");
  for (const auto &[div, stats] : sStat) {
    ACTS_DEBUG("Lvl: " << div 
              << " | Cut(Thresh): " << stats.discardedByThresholdCut  
              << " | Cut(Cross): " << stats.discardedByCrossingCut 
              << " | Area(Drill): " << stats.drilledArea
              << " | nSections(Drill_standard): " << stats.nDrilledSections_standard
              << " | nSections(Drill_extra): " << stats.nDrilledSections_extra
              << " | Area(Accept): " << stats.acceptedArea
              << " | nSections(Accept): " << stats.nAcceptedSections
              << " | nLines(Accept): " << stats.nAcceptedLines);
  }
  ACTS_DEBUG("===========================================================");
}

void AdaptiveHoughTransformSeeder::processStackZCotTheta(
    std::vector<Acts::Experimental::HoughAccumulatorSection> &input,
    std::vector<SeedCandidate> &output,
    const std::vector<PreprocessedMeasurement> &measurements) const {
    
    Acts::ScopedTimer st("processStackZCotTheta", logger());

    // set created to root out the duplication
    std::unordered_set<std::vector<unsigned>, VectorHash> unique_clusters;

    const float max_y_spread = m_cfg.cotThetaMinBinSize * 2.5f;
    const unsigned req_lines = m_cfg.threshold;
    struct LinePoint {
        float y;
        unsigned idx;
        bool operator<(const LinePoint& other) const {
            return y < other.y;
        }
    };

    for (const auto &section : input) {
        const unsigned cnt = section.count();
        
        if (cnt < req_lines) {
            continue;
        }

        const float xLeft = section.xBegin();
        const float xRight = xLeft + section.xSize();

        // 7 x point option
        const float xM1 = xLeft+(xRight-xLeft)/6*1;
        const float xM2 = xLeft+(xRight-xLeft)/6*2;
        const float xM3 = xLeft+(xRight-xLeft)/6*3;
        const float xM4 = xLeft+(xRight-xLeft)/6*4;
        const float xM5 = xLeft+(xRight-xLeft)/6*5;
        const float check_points[7] = {xLeft,xM1, xM2,xM3, xM4, xM5, xRight};

        std::vector<LinePoint> yVals(cnt);

        // loop calculating yValue for each line in selected x points
        // then it checks if section contains clusters which we take as solutions
        for (float x_check : check_points) {
            for (unsigned i = 0; i < cnt; ++i) {
                unsigned meas_idx = section.indices()[i];
                const auto &m = measurements[meas_idx];
                
                yVals[i].y = -m.invr * x_check + m.z * m.invr;
                yVals[i].idx = meas_idx;
            }

            std::sort(yVals.begin(), yVals.end());

            for (unsigned i = 0; i <= cnt - req_lines; ++i) {
              if (yVals[i + req_lines - 1].y - yVals[i].y <= max_y_spread) {
                std::vector<unsigned> clustered_indices;
                clustered_indices.reserve(cnt);
                float base_y = yVals[i].y;
                for(unsigned j = i; j < cnt; ++j) {
                    if (yVals[j].y - base_y <= max_y_spread) {
                        clustered_indices.push_back(yVals[j].idx);
                    } else {
                        break;
                    }
                }

                std::sort(clustered_indices.begin(), clustered_indices.end());

                // hashing is added in order to prevent duplication of solutions
                if (unique_clusters.insert(clustered_indices).second) {
                  SeedCandidate candidate;
                  candidate.indices = std::move(clustered_indices);
                  output.push_back(std::move(candidate));
                }
              }
            }
        }
    }

    ACTS_DEBUG("===========================================================");
    ACTS_DEBUG(" [Z - CotTheta] process stats ");
    ACTS_DEBUG("===========================================================");
    ACTS_DEBUG("Input sections from Phase 1: " << input.size());
    ACTS_DEBUG("Valid physical clusters found: " << output.size());
    ACTS_DEBUG("===========================================================");
}

void AdaptiveHoughTransformSeeder::processStackZCotThetaSplit(
    std::vector<Acts::Experimental::HoughAccumulatorSection> &input,
    std::vector<Acts::Experimental::HoughAccumulatorSection> &output,
    const std::vector<PreprocessedMeasurement> &measurements) const {
  ExplorationOptions opt;
  opt.xMinBinSize = 101.0f * Acts::UnitConstants::mm;
  opt.yMinBinSize = 10.1f;
  opt.lineFunctor = m_zCotThetaLineParams;
  opt.decisionFunctor =
      [&cfg = m_cfg, &opt](
          const Acts::Experimental::HoughAccumulatorSection &section,
          const std::vector<PreprocessedMeasurement> &) {
        using enum HoughAccumulatorSection::Decision;
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

  exploreHoughParametersSpace(input, measurements, opt, output);
  const unsigned nl = std::accumulate(
      output.begin(), output.end(), 0,
      [](unsigned sum, const Acts::Experimental::HoughAccumulatorSection &s) {
        return sum + s.count();
      });
  ACTS_DEBUG("Z - cotTheta split scan finds "
             << output.size() << " included measurements " << nl);
}

void AdaptiveHoughTransformSeeder::makeSeeds(
    SeedContainer &seeds,
    const std::vector<SeedCandidate> &solutions,
    const std::vector<PreprocessedMeasurement> &measurements) const {
  const SpacePointContainer &spacePoints = seeds.spacePointContainer();

  unsigned skippedCandidatesCounter = 0;

  std::size_t seedIndex = 0;
  for (const SeedCandidate &s : solutions) {
    std::vector<unsigned> sortedIndices = s.indices;
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

    if (spIndex<3){
      skippedCandidatesCounter++;
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

    auto seed = seeds.createSeed();
    seed.assignSpacePointIndices(std::array<SpacePointIndex, 3>{
        sp[0]->index(), sp[1]->index(), sp[2]->index()});

    // for the time the quality is fixed
    // in the future we can use section count for instance
    seed.quality() = 1.0;

    const auto z =
        static_cast<float>(sp[1]->z() - sp[1]->r() * cotThetaEstimate);
    seed.vertexZ() = z;

    seedIndex++;
  }
  if (skippedCandidatesCounter > 0) {
    ACTS_VERBOSE("Diagnostic: Skipped " << skippedCandidatesCounter 
              << " invalid seed candidates (insufficient points) out of " 
              << solutions.size() << " total sections.");
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
  std::ranges::transform(
      input, std::back_inserter(op),
      [](const Acts::Experimental::HoughAccumulatorSection &s) { return &s; });

  auto binaryPredSort =
      [](const Acts::Experimental::HoughAccumulatorSection *a,
         const Acts::Experimental::HoughAccumulatorSection *b) {
        return a->indices() < b->indices();
      };
  auto binaryPredUnique =
      [](const Acts::Experimental::HoughAccumulatorSection *a,
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
