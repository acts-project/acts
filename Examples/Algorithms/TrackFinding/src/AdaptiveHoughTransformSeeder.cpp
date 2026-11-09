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
#include <stdexcept>

namespace ActsExamples {

AdaptiveHoughTransformSeeder::AdaptiveHoughTransformSeeder(
    const Config &cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("AdaptiveHoughTransformSeeder", std::move(logger)),
      m_cfg(cfg),
      m_bFieldA(static_cast<float>(1.0 / cfg.inverseA)) {
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

  std::vector<Phase1Candidate> phase1Candidates;
  {
    Acts::ScopedTimer st("SlidingWindowQPtPhi", logger());
    processQOverPtPhiSlidingWindow(stack1, phase1Candidates, measurements);
    ACTS_DEBUG("Sliding Window Q/pT-Phi generated " << phase1Candidates.size()
                                                    << " candidates.");
  }

  // additional stack for smaller structure as HoughAccumulatorSection is no
  // longer needed
  std::vector<SeedCandidate> stack3;

  // do scan in z_vertex - cot theta
  if (m_cfg.doSecondPhase) {
    Acts::ScopedTimer st("secondPhase", logger());

    processStackZCotTheta(phase1Candidates, stack3, measurements);
    ACTS_DEBUG("past StackZCotTheta phase1Candidates size "
               << phase1Candidates.size() << " stack3 size " << stack3.size());
  }

  // vector for selected candidates with or without second phase
  std::vector<SeedCandidate> selectedCandidates;

  if (m_cfg.doSecondPhase) {
    selectedCandidates = std::move(stack3);
  } else {
    selectedCandidates.reserve(phase1Candidates.size());
    for (auto &candidate : phase1Candidates) {
      selectedCandidates.push_back(
          {{candidate.indices[0], candidate.indices[1], candidate.indices[2]}});
    }
  }

  Acts::ScopedTimer st("seedsMaking", logger());

  // post solutions
  SeedContainer seeds;
  seeds.assignSpacePointContainer(spacePoints);
  seeds.reserve(selectedCandidates.size());

  ACTS_VERBOSE("Number of solutions " << selectedCandidates.size());
  makeSeeds(seeds, selectedCandidates, measurements);

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
    const float rVal = sp.r();
    const float phi = std::atan2(sp.y(), sp.x());
    const float invr = 1.0f / rVal;
    measurements.emplace_back(rVal, invr, phi, sp.z(), sp.index());
    // wrap phi by duplicating some seeds
    if (phi < -std::numbers::pi + m_cfg.phiWrap * std::numbers::pi) {
      measurements.emplace_back(rVal, invr, phi + 2 * std::numbers::pi, sp.z(),
                                sp.index());
    }
    if (phi > std::numbers::pi - m_cfg.phiWrap * std::numbers::pi) {
      measurements.emplace_back(rVal, invr, phi - 2 * std::numbers::pi, sp.z(),
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
  const int nSplits = 215;
  const auto wedgeWidth = static_cast<float>(2.0 * std::numbers::pi / nSplits);

  for (int phiIndex = 0; phiIndex < nSplits; phiIndex++) {
    const auto startPhi = static_cast<float>(
        wedgeWidth * static_cast<float>(phiIndex) - std::numbers::pi);
    stack.emplace_back(1.05f * wedgeWidth, 2.0f * m_cfg.qOverPtMin, startPhi,
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
    int nDrilledSectionsStandard{};
    int nDrilledSectionsExtra{};
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
      [&sStat, &cfg = m_cfg](
          const Acts::Experimental::HoughAccumulatorSection &section,
          const std::vector<PreprocessedMeasurement> &) {
        using enum HoughAccumulatorSection::Decision;

        if (section.count() < cfg.threshold) {
          sStat[section.divisionLevel()].discardedByThresholdCut += 1;
          return Discard;
        }

        if (section.divisionLevel() <= cfg.maxDivision - 1) {
          sStat[section.divisionLevel()].nDrilledSectionsStandard++;
          sStat[section.divisionLevel()].drilledArea +=
              section.xSize() * section.ySize();
          return Drill;
        }

        if (section.count() >= cfg.threshold &&
            section.xSize() <= cfg.phiMinBinSize &&
            section.ySize() <= cfg.qOverPtMinBinSize) {
          sStat[section.divisionLevel()].nAcceptedSections++;
          sStat[section.divisionLevel()].nAcceptedLines += section.count();
          sStat[section.divisionLevel()].acceptedArea +=
              section.xSize() * section.ySize();
          return Accept;
        }
        sStat[section.divisionLevel()].drilledArea +=
            section.xSize() * section.ySize();
        sStat[section.divisionLevel()].nDrilledSectionsExtra++;
        return Drill;
      };

  exploreHoughParametersSpace(input, measurements, opt, output);
  const unsigned nl = std::accumulate(
      output.begin(), output.end(), 0,
      [](unsigned sum, const Acts::Experimental::HoughAccumulatorSection &s) {
        return sum + s.count();
      });
  ACTS_VERBOSE("Phi - qOverPt scan finds "
               << output.size() << " solutions, included measurements " << nl);
  ACTS_VERBOSE("size " << sStat.size());
  ACTS_VERBOSE("===========================================================");
  ACTS_VERBOSE(" [Q/pT - Phi] process stats ");
  ACTS_VERBOSE("===========================================================");
  for (const auto &[div, stats] : sStat) {
    ACTS_VERBOSE(
        "Lvl: " << div << " | Cut(Thresh): " << stats.discardedByThresholdCut
                << " | Cut(Cross): " << stats.discardedByCrossingCut
                << " | Area(Drill): " << stats.drilledArea
                << " | nSections(Drill_standard): "
                << stats.nDrilledSectionsStandard
                << " | nSections(Drill_expand): " << stats.nDrilledSectionsExtra
                << " | Area(Accept): " << stats.acceptedArea
                << " | nSections(Accept): " << stats.nAcceptedSections
                << " | nLines(Accept): " << stats.nAcceptedLines);
  }
  ACTS_VERBOSE("===========================================================");
}

void AdaptiveHoughTransformSeeder::findClustersQOverPtPhi(
    const Acts::Experimental::HoughAccumulatorSection &section,
    std::span<const float> checkPoints,
    const std::vector<PreprocessedMeasurement> &measurements,
    const LineFunctor &lineFunctor, float maxSpread,
    std::unordered_set<std::vector<unsigned>, VectorHash> &uniqueClusters,
    std::vector<Phase1Candidate> &output) const {
  const unsigned cnt = section.count();
  const unsigned reqLines = m_cfg.threshold;

  if (cnt < reqLines) {
    return;
  }

  struct LinePoint {
    float x;
    unsigned idx;
    bool operator<(const LinePoint &other) const { return x < other.x; }
  };

  std::vector<LinePoint> xVals;
  xVals.reserve(cnt);

  // Evaluate line intersections along the primary axis for each check point
  for (float yCheck : checkPoints) {
    xVals.clear();

    for (unsigned i = 0; i < cnt; ++i) {
      unsigned measIdx = section.indices()[i];
      const auto &m = measurements[measIdx];

      xVals.push_back({lineFunctor(m, yCheck), measIdx});
    }

    const unsigned validCnt = xVals.size();
    if (validCnt < reqLines) {
      continue;
    }

    // Sort evaluated line positions along the phi axis to enable 1D sliding
    // window cluster identification in O(N log N) time
    std::sort(xVals.begin(), xVals.end());

    unsigned start = 0;
    unsigned end = reqLines - 1;
    bool candidateFound = false;

    // Two-pointer greedy 1D sliding window scan
    while (end < validCnt) {
      if (xVals[end].x - xVals[start].x <= maxSpread) {
        end++;
      } else {
        if (end - start >= reqLines) {
          std::vector<unsigned> clusteredIndices;
          clusteredIndices.reserve(end - start);
          for (unsigned k = start; k < end; ++k) {
            clusteredIndices.push_back(xVals[k].idx);
          }

          // Sort measurement indices to form a canonical key for hash
          // deduplication
          std::sort(clusteredIndices.begin(), clusteredIndices.end());

          if (uniqueClusters.insert(clusteredIndices).second) {
            Phase1Candidate candidate;
            candidate.indices = std::move(clusteredIndices);
            candidate.zMin = section.history(zSplitMinIndex);
            candidate.zWidth = section.history(zSplitWidthIndex);
            candidate.cotThetaMin = section.history(cotThetaSplitMinIndex);
            candidate.cotThetaWidth = section.history(cotThetaSplitWidthIndex);
            output.push_back(std::move(candidate));
            candidateFound = true;
            break;
          }
        }
        start++;
      }
    }

    if (!candidateFound && end - start >= reqLines) {
      std::vector<unsigned> clusteredIndices;
      clusteredIndices.reserve(end - start);
      for (unsigned k = start; k < end; ++k) {
        clusteredIndices.push_back(xVals[k].idx);
      }

      std::sort(clusteredIndices.begin(), clusteredIndices.end());

      if (uniqueClusters.insert(clusteredIndices).second) {
        Phase1Candidate candidate;
        candidate.indices = std::move(clusteredIndices);
        candidate.zMin = section.history(zSplitMinIndex);
        candidate.zWidth = section.history(zSplitWidthIndex);
        candidate.cotThetaMin = section.history(cotThetaSplitMinIndex);
        candidate.cotThetaWidth = section.history(cotThetaSplitWidthIndex);
        output.push_back(std::move(candidate));
      }
    }
  }
}

void AdaptiveHoughTransformSeeder::processQOverPtPhiSlidingWindow(
    const std::vector<Acts::Experimental::HoughAccumulatorSection> &input,
    std::vector<Phase1Candidate> &output,
    const std::vector<PreprocessedMeasurement> &measurements) const {
  std::unordered_set<std::vector<unsigned>, VectorHash> uniqueClusters;

  const float maxXSpread = m_cfg.phiMinBinSize * 0.35f;
  const unsigned reqLines = m_cfg.threshold;
  constexpr std::size_t kNumCheckPoints = 1;

  for (const auto &section : input) {
    if (section.count() < reqLines) {
      continue;
    }

    const float yBottom = section.yBegin();
    const float yTop = yBottom + section.ySize();

    // Evaluate at the midpoint of the section using std::lerp
    const std::array<float, kNumCheckPoints> checkPoints = {
        std::lerp(yBottom, yTop, 0.5f)};

    findClustersQOverPtPhi(section, checkPoints, measurements,
                           m_phiQOverPtLineParams, maxXSpread, uniqueClusters,
                           output);
  }
}

void AdaptiveHoughTransformSeeder::findClustersZCotTheta(
    const std::vector<unsigned> &indices, std::span<const float> checkPoints,
    const std::vector<PreprocessedMeasurement> &measurements,
    const LineFunctor &lineFunctor, float maxSpread,
    std::unordered_set<std::array<unsigned, 3>, Array3Hash> &uniqueClusters,
    std::vector<SeedCandidate> &output) const {
  const unsigned cnt = indices.size();
  const unsigned reqLines = m_cfg.threshold;

  if (cnt < reqLines) {
    return;
  }

  struct LinePoint {
    float y;
    unsigned idx;
    bool operator<(const LinePoint &other) const { return y < other.y; }
  };

  std::vector<LinePoint> yVals(cnt);

  // Evaluate all lines in the candidate at each discrete z check point
  for (float xCheck : checkPoints) {
    for (unsigned i = 0; i < cnt; ++i) {
      unsigned measIdx = indices[i];
      const auto &m = measurements[measIdx];

      // Evaluate line coordinate using the provided space-specific functor
      yVals[i].y = lineFunctor(m, xCheck);
      yVals[i].idx = measIdx;
    }

    // Sort evaluated line positions along the secondary axis to identify
    // dense line crossings in O(N log N) time
    std::sort(yVals.begin(), yVals.end());

    // 1D sliding window scan with radial triplet verification
    for (unsigned i = 0; i <= cnt - reqLines; ++i) {
      if (yVals[i + reqLines - 1].y - yVals[i].y <= maxSpread) {
        const float baseY = yVals[i].y;

        double minR = measurements[yVals[i].idx].r;
        double maxR = minR;
        unsigned clusterEndIdx = cnt;
        unsigned idxMax = i;
        unsigned idxMin = i;
        unsigned idxMid = i;

        // Find radial bounds (minR and maxR) within the current dense cluster
        for (unsigned j = i; j < cnt; ++j) {
          if (yVals[j].y - baseY <= maxSpread) {
            const double r = measurements[yVals[j].idx].r;
            if (r < minR) {
              minR = r;
              idxMin = j;
            } else if (r > maxR) {
              maxR = r;
              idxMax = j;
            }
          } else {
            clusterEndIdx = j;
            break;
          }
        }

        // Verify that an intermediate space point exists with sufficient
        // radial separation (> 5 mm) from both inner and outer points
        bool validRSpread = false;
        for (unsigned j = i; j < clusterEndIdx; ++j) {
          const double r = measurements[yVals[j].idx].r;
          if (r - minR > 5.0 * Acts::UnitConstants::mm &&
              maxR - r > 5.0 * Acts::UnitConstants::mm) {
            validRSpread = true;
            idxMid = j;
            break;
          }
        }

        if (!validRSpread || maxR - minR < 10.0 * Acts::UnitConstants::mm) {
          continue;
        }

        std::array<unsigned, 3> clusteredIndices = {
            yVals[idxMin].idx, yVals[idxMax].idx, yVals[idxMid].idx};

        // Sort measurement indices to form a canonical key for hash
        // deduplication
        std::sort(clusteredIndices.begin(), clusteredIndices.end());

        if (uniqueClusters.insert(clusteredIndices).second) {
          output.push_back({clusteredIndices});
        }
      }
    }
  }
}

void AdaptiveHoughTransformSeeder::processStackZCotTheta(
    std::vector<Phase1Candidate> &input, std::vector<SeedCandidate> &output,
    const std::vector<PreprocessedMeasurement> &measurements) const {
  Acts::ScopedTimer st("processStackZCotTheta", logger());

  // Hash set to root out duplicate triplets across candidates and check points
  std::unordered_set<std::array<unsigned, 3>, Array3Hash> uniqueClusters;

  const float maxYSpread = m_cfg.cotThetaMinBinSize * 1.5f;
  const unsigned reqLines = m_cfg.threshold;
  constexpr std::size_t kNumCheckPoints = 9;

  for (const auto &candidate : input) {
    if (candidate.indices.size() < reqLines) {
      continue;
    }

    const float xLeft = candidate.zMin;
    const float xRight = xLeft + candidate.zWidth;

    // Generate evenly spaced check points across the z width using std::lerp
    std::array<float, kNumCheckPoints> checkPoints{};
    for (std::size_t i = 0; i < kNumCheckPoints; ++i) {
      checkPoints[i] = std::lerp(xLeft, xRight,
                                 static_cast<float>(i) / (kNumCheckPoints - 1));
    }

    // Perform 1D cluster search in z - cot(theta) space
    findClustersZCotTheta(candidate.indices, checkPoints, measurements,
                          m_zCotThetaLineParams, maxYSpread, uniqueClusters,
                          output);
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
    SeedContainer &seeds, const std::vector<SeedCandidate> &solutions,
    const std::vector<PreprocessedMeasurement> &measurements) const {
  const SpacePointContainer &spacePoints = seeds.spacePointContainer();

  unsigned int skippedCandidatesCounter = 0;

  for (const SeedCandidate &s : solutions) {
    std::array<unsigned, 3> sortedIndices = s.indices;
    std::ranges::sort(sortedIndices, [&measurements](unsigned i1, unsigned i2) {
      return measurements[i1].r < measurements[i2].r;
    });

    unsigned spIndex = 0;
    std::array<std::optional<ConstSpacePointProxy>, 3> sp = {
        std::nullopt, std::nullopt, std::nullopt};
    double lastR = -1000.0;
    for (unsigned sidx : sortedIndices) {
      const PreprocessedMeasurement &m = measurements[sidx];
      if (spIndex == 0 ||
          std::abs(m.r - lastR) > 5.0 * Acts::UnitConstants::mm) {
        sp[spIndex] = spacePoints.at(m.sp);
        lastR = m.r;
        spIndex++;
      }
      if (spIndex >= sp.size()) {
        break;
      }
    }

    if (spIndex < 3) {
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
  }
  if (skippedCandidatesCounter > 0) {
    ACTS_VERBOSE("Diagnostic: Skipped "
                 << skippedCandidatesCounter
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
