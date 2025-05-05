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
    : m_xSize(xw), m_ySize(yw), m_xBegin(xBegin), m_yBegin(yBegin),
      m_divisionLevel(div), m_indices(indices) {}

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

ProcessCode
AdaptiveHoughTransformSeeder::execute(const AlgorithmContext &ctx) const {
  // get inputs
  std::vector<PreprocessedMeasurement> measurements;
  for (const auto &isp : m_inputSpacePoints) {
    const auto &spContainer = (*isp)(ctx);
    ACTS_DEBUG("Inserting " << spContainer.size() << " space points from "
                            << isp->key());
    for (const SimSpacePoint &sp : spContainer) {
      measurements.push_back(PreprocessedMeasurement(
          1.0 / sp.r(), std::atan2(sp.y(), sp.x()),
          Acts::SourceLink(static_cast<const SimSpacePoint *>(&sp))));
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
  while (not sectionsStack.empty()) {
    ACTS_VERBOSE("Processing AccumulatorSection "
                 << sectionsStack.top().xBegin() << " "
                 << sectionsStack.top().xSize() << " "
                 << sectionsStack.top().yBegin() << " "
                 << sectionsStack.top().ySize()
                 << " nlines: " << sectionsStack.top().count());
    processStackHead(sectionsStack, solutions, measurements);
  }

  // post solutions
  SimSeedContainer seeds;
  seeds.reserve(solutions.size());

  ProtoTrackContainer protoTracks;
  protoTracks.reserve(solutions.size());
  ACTS_DEBUG("Number of solutions " << solutions.size());
  for (const AccumulatorSection &s : solutions) {
    ACTS_DEBUG("Solution x: " << s.xBegin() << " " << s.xSize()
                              << " y: " << s.yBegin() << " " << s.ySize()
                              << " nlines: " << s.count());
    const SimSpacePoint *sp0 =
        measurements[s.indices()[0]].link.get<const SimSpacePoint *>();
    const SimSpacePoint *sp1 =
        measurements[s.indices()[1]].link.get<const SimSpacePoint *>();
    const SimSpacePoint *sp2 =
        measurements[s.indices()[2]].link.get<const SimSpacePoint *>();
    seeds.emplace_back(*sp0, *sp1, *sp2);
    float cotThetaEstimate = (sp2->z() - sp0->z()) / (sp2->r() - sp0->r());
    float z = sp1->z() - sp1->r() * cotThetaEstimate;
    seeds.back().setVertexZ(z);
    // for the time the quality is fixed
    // in the future we can use section count for instance
    seeds.back().setQuality(1.0);
  }

  m_outputSeeds(ctx, std::move(seeds));

  return ActsExamples::ProcessCode::SUCCESS;
}

void AdaptiveHoughTransformSeeder::processStackHead(
    std::stack<AccumulatorSection> &sections,
    std::vector<AccumulatorSection> &solutions,
    const std::vector<PreprocessedMeasurement> &measurements) const {
  AccumulatorSection &section = sections.top();

  // check if it is an ~empty section
  if (section.count() < config().threshold) {
    ACTS_VERBOSE("Failed threshold check");
    sections.pop();
    return;
  }
  // check if intersections requirement fails
  if (config().requireIntersections &&
      !passIntersectionsCheck(sections.top(), measurements)) {
    ACTS_VERBOSE("Failed intersection check");
    sections.pop();
    return;
  }
  // maybe we have found solution
  if (section.xSize() <= config().phiMinBinSize &&
      section.ySize() <= config().qOverPtMinBinSize) {
    addSolution(std::move(section), solutions);
    sections.pop();
    return;
  }

  // time to split
  std::vector<AccumulatorSection> divisions;
  if (section.xSize() > config().phiMinBinSize &&
      section.ySize() > config().qOverPtMinBinSize) {
    // need 4 subdivisions
    divisions.push_back(section.topLeft());
    divisions.push_back(section.topRight());
    divisions.push_back(section.bottomLeft());
    divisions.push_back(section.bottomRight());
  } else if (section.xSize() <= config().phiMinBinSize &&
             section.ySize() > config().qOverPtMinBinSize) {
    // only split in q over pT
    divisions.push_back(section.top());
    divisions.push_back(section.bottom());
  } else {
    divisions.push_back(section.left());
    divisions.push_back(section.right());
  }
  ACTS_VERBOSE("This section is split into " << divisions.size());
  sections.pop();
  for (AccumulatorSection &s : divisions) {
    updateSection(s, measurements);
    sections.push(std::move(s));
  }
}

void AdaptiveHoughTransformSeeder::updateSection(
    AccumulatorSection &section,
    const std::vector<PreprocessedMeasurement> &input) const {
  std::vector<unsigned> selectedIndices;
  for (unsigned index : section.indices()) {
    const PreprocessedMeasurement &m = input[index];
    if (section.isLineInside(m.invr * config().inverseA,
                             -m.invr * m.phi * config().inverseA)) {
      selectedIndices.push_back(index);
    }
  }
  section.indices() = std::move(selectedIndices);
}

bool AdaptiveHoughTransformSeeder::passIntersectionsCheck(
    const AccumulatorSection &section,
    const std::vector<PreprocessedMeasurement> &measurements) const {
  unsigned inside = 0;
  for (std::size_t idx1 = 0; idx1 < section.count(); ++idx1) {
    const auto first = section.indices()[idx1];
    const double a1 = measurements[first].invr * config().inverseA;
    const double b1 =
        -measurements[first].invr * measurements[first].phi * config().inverseA;
    for (std::size_t idx2 = idx1 + 1; idx2 < section.count(); ++idx2) {
      const auto second = section.indices()[idx2];
      const double a2 = measurements[second].invr * config().inverseA;
      const double b2 = -measurements[second].invr * measurements[second].phi *
                        config().inverseA;

      const double adif = a1 - a2;
      if (std::abs(adif) < 1e-5) { // Parallel lines, need to skip
        continue;
      }
      const double bdif = b2 - b1;
      const double solX = bdif / adif;
      if (section.xBegin() <= solX &&
          solX <= section.xBegin() + section.xSize()) {
        const double y = std::fma(a1, bdif / adif, b1);
        if (section.yBegin() <= y && y <= section.yBegin() + section.ySize()) {
          inside++;
        }
      }
      if (inside >= config().intersectionsThreshold) {
        return true;
      }
    }
  }
  ACTS_VERBOSE("Number of crossings inside of section " << inside);
  return inside >= config().intersectionsThreshold;
}

void AdaptiveHoughTransformSeeder::addSolution(
    AccumulatorSection &&s, std::vector<AccumulatorSection> &output) const {
  if (config().deduplicate) {
    std::size_t i = output.size() > 4 ? output.size() - 4 : 0;
    for (; i < std::size(output); ++i) {
      if (output[i].indices() == s.indices()) {
        ACTS_VERBOSE(
            "There is a nearby duplicate solution. Skipping this one.");
        return;
      }
    }
  }
  output.push_back(std::move(s));
}
} // namespace ActsExamples
