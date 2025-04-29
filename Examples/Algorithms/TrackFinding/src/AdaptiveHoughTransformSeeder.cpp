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
                                       double yBegin, int div, const std::vector<unsigned>& indices)
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
  return AccumulatorSection(
      m_xSize * xFraction, m_ySize * yFraction, m_xBegin + m_xSize - m_xSize * xFraction,
      m_yBegin + m_ySize - m_ySize * yFraction, m_divisionLevel + 1, m_indices);
}
AccumulatorSection AccumulatorSection::bottomRight(
    float xFraction, float yFraction) const {
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
  for (const auto& spName : config().inputSpacePoints) {
    auto& handle = m_inputSpacePoints.emplace_back(
        std::make_unique<ReadDataHandle<SimSpacePointContainer>>(
            this,
            "InputSpacePoints#" + std::to_string(m_inputSpacePoints.size())));
    handle->initialize(spName);
  }
}

ProcessCode AdaptiveHoughTransformSeeder::execute(
    const AlgorithmContext& ctx) const {
  // get inputs
  std::vector<PreprocessedMeasurement> measurements;

  for (const auto& isp : m_inputSpacePoints) {
    const auto& spContainer = (*isp)(ctx);
    ACTS_DEBUG("Inserting " << spContainer.size() << " space points from "
                            << isp->key());
    for (auto& sp : spContainer) {
      double r = Acts::fastHypot(sp.x(), sp.y());
      double z = sp.z();
    }
  }

  // prepare initial stack
  // will become more complicated with slicing
  std::stack<AccumulatorSection> sectionsStack({AccumulatorSection(
    2.*M_PI, 2.0*config().qOverPtMin, -M_PI, -config().qOverPtMin)});


  // all the measurements should be considered at start
  // therfore indices of all of them are stored
  sectionsStack.top().indices().resize(measurements.size());
  std::iota(std::begin(sectionsStack.top().indices()),
            std::end(sectionsStack.top().indices()), 0);



  std::vector<AccumulatorSection> solutions;
  while (not sectionsStack.empty()) {
    updateSection(sectionsStack.top(), measurements);
    processStackHead(sectionsStack, solutions);
  }
  // post solutions

  return ActsExamples::ProcessCode::SUCCESS;
}

void AdaptiveHoughTransformSeeder::processStackHead(
    std::stack<AccumulatorSection>& sections,
    std::vector<AccumulatorSection>& solutions) const {
  const AccumulatorSection& section = sections.top();

  // check if it is a bad one (most likely)
  if (section.count() < config().threshold) {
    sections.pop();
  } else {
    // we have found solution
    if (section.xSize() <= config().phiMinBinSize &&
        section.ySize() <= config().qOverPtMinBinSize) {
      solutions.push_back(section);
      sections.pop();
    } else {
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
      sections.pop();
      for (const AccumulatorSection& s : divisions) {
        sections.push(std::move(s));
      }
    }
  }
}

void AdaptiveHoughTransformSeeder::updateSection(
    AccumulatorSection& section,
    const std::vector<PreprocessedMeasurement>& input) const {
  std::vector<unsigned> selectedIndices;
  for (unsigned index : section.indices()) {
    const PreprocessedMeasurement& m = input[index];
    if (section.isLineInside(m.r, m.phi)) {
      selectedIndices.push_back(index);
    }
  }
  section.indices() = std::move(selectedIndices);
}

}  // namespace ActsExamples
