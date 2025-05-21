// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// @file AdaptiveHoughTransformSeeder.hpp
// @author Tomasz Bold
// @brief Implements track-seeding using Adaptive Hough Transform.
//

#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <cstddef>
#include <memory>
#include <numbers>
#include <stack>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace ActsExamples {
struct AlgorithmContext;
}

namespace Acts {
class TrackingGeometry;
}

namespace ActsExamples {
// Helper class describing one section of the accumulator space
class AccumulatorSection {
 public:
  AccumulatorSection() = default;

  AccumulatorSection(double xw, double yw, double xBegin, double yBegin,
                     int div = 0, const std::vector<unsigned> &indices = {});

  inline unsigned count() const { return m_indices.size(); }
  inline const std::vector<unsigned> &indices() const { return m_indices; }
  inline std::vector<unsigned> &indices() { return m_indices; }

  // create section that is bottom part this this one
  // +------+
  // |      |
  // +------+
  // |     <|-- this part
  // +------+
  AccumulatorSection bottom(float yFraction = 0.5) const;
  // see @bottom
  AccumulatorSection top(float yFraction = 0.5) const;
  // @see @bottom
  AccumulatorSection left(float xFraction = 0.5) const;
  // @see @bottom
  AccumulatorSection right(float xFraction = 0.5) const;

  // create section that is bottom left corner of this this one
  // by default the section is divided into 4 quadrants,
  // if parameters are provided the quadrants size can be adjusted
  // +---+---+
  // |   |   |
  // +---+---+
  // |   |  <|-- this part
  // +---+---+
  AccumulatorSection bottomRight(float xFraction = 0.5,
                                 float yFraction = 0.5) const;
  AccumulatorSection bottomLeft(float xFraction = 0.5,
                                float yFraction = 0.5) const;

  AccumulatorSection topLeft(float xFraction = 0.5,
                             float yFraction = 0.5) const;
  AccumulatorSection topRight(float xFraction = 0.5,
                              float yFraction = 0.5) const;

  ///  @brief true if the line defined by given parameters passes the section
  /// a and b are line parameters y = ax + b
  inline bool isLineInside(float a, float b) const {
    const float yB = std::fma(a, m_xBegin, b);
    const float yE = std::fma(a, (m_xBegin + m_xSize), b);
    return (a > 0) ? yB < m_yBegin + m_ySize && yE > m_yBegin
                   : yB > m_yBegin && yE < m_yBegin + m_ySize;
  }

  /// @brief check if the lines cross inside the section
  /// @param a1 line 1 parameter a
  /// @param b1 line 1 parameter b
  /// @param a2 line 2 parameter a
  /// @param b2 line 2 parameter b
  /// @return true if the lines cross in the section
  inline bool isCrossingInside(float a1, float b1, float a2, float b2) const {
    const double adif = a1 - a2;
    if (std::abs(adif) < 1e-3) {  // nearly Parallel lines, never cross
      return false;
    }
    const double bdif = b2 - b1;
    const double solX = bdif / adif;
    if (xBegin() <= solX && solX <= xBegin() + xSize()) {
      const double y = std::fma(a1, bdif / adif, b1);
      if (yBegin() <= y && y <= yBegin() + ySize()) {
        return true;
      }
    }
    return false;
  }

  // counter clock wise distance from upper left corner
  // a and b are line parameters y = ax + b
  float distCC(float a, float b) const;
  // anti-counter clock wise distance from upper left corner
  float distACC(float a, float b) const;

  // sizes
  double xSize() const { return m_xSize; }
  double ySize() const { return m_ySize; }
  double xBegin() const { return m_xBegin; }
  double yBegin() const { return m_yBegin; }
  unsigned divisionLevel() const { return m_divisionLevel; }

 private:
  double m_xSize = 0;
  double m_ySize = 0;
  double m_xBegin = 0;
  double m_yBegin = 0;
  unsigned m_divisionLevel =
      0;  // number of times the starting section was already divided
  std::vector<unsigned>
      m_indices;  // indices of measurements contributing to this section
};

// Construct track seeds from space points.
class AdaptiveHoughTransformSeeder final : public IAlgorithm {
 public:
  struct Config {
    /// Input space point collections.
    ///
    /// We allow multiple space point collections to allow different parts of
    /// the detector to use different algorithms for space point construction,
    /// e.g. single-hit space points for pixel-like detectors or double-hit
    /// space points for strip-like detectors.
    /// Note that we don't *need* spacepoints (measurements can be used instead)
    std::vector<std::string> inputSpacePoints;
    /// Output track seed collection.
    std::string outputSeeds;
    /// Output hough track collection.
    std::string outputProtoTracks;
    /// Tracking geometry required to access global-to-local transforms.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    float phiWrap =
        0.1;  // wrap size around angle domain limits (min pT dependent)
    float qOverPtMin = 1.0;  // min q/pt, -1/1 GeV

    float qOverPtMinBinSize = 0.01;  // minimal size of pT bin that the
                                     // algorithm should not go beyond (in GeV)
    float phiMinBinSize = 0.01;  // minimal size of phi bin that the algorithm
                                 // should not go beyond
    float zRange = 200;          // range in z
    float cotThetaRange = 10;    // range in cotTheta

    float zMinBinSize =
        5;  // minimal size of z bin that the algorithm should
            // not go beyond when exploring zvertex-cot(theta) space space
    float cotThetaMinBinSize = 0.2;  // minimal size of cot(theta) bin that
                                     // the algorithm should not go beyond
    unsigned threshold =
        8;  // number of lines passing section for it to be still considered
    bool doSecondPhase = true;  // do the second pass in z-cot(theta) space to
                                // find less solutions
    bool deduplicate = true;    // when adding solutions try avoiding duplicates

    double inverseA =
        1.0 / 3.0e-4;  // Assume B = 2T constant. Can apply corrections to
                       // this with fieldCorrection function
                       // This 3e-4 comes from the 2T field when converted to
                       // units of GeV / (c*mm*e)
  };

  // information that is needed for each measurement
  struct PreprocessedMeasurement {
    /// Construct the measurement
    /// @param inverseR inverse of radius of the SP
    /// @param phiAngle azimuthal angle of the SP
    /// @param zpos z position of the SP
    /// @param l link to space point
    PreprocessedMeasurement(double inverseR, double phiAngle, double zpos,
                            Acts::SourceLink l)
        : invr(inverseR), phi(phiAngle), z(zpos), link(std::move(l)) {}
    double invr;
    double phi;
    double z;
    Acts::SourceLink link;
  };

  template <typename M = PreprocessedMeasurement>
  struct AHTExplorationOptions {
    float xMinBinSize = 1;  // minimum bin size in x direction, beyond that
                            // value the sections are not split
    float yMinBinSize = 1;  // minimum bin size in y direction, beyond that
                            // value the sections are not split
    using LineParamFunctors = std::pair<std::function<float(const M &)>,
                                        std::function<float(const M &)>>;
    LineParamFunctors lineParamFunctors;  // pair of functions needed to obtain
                                          // linear function ax+b parameters,
                                          // first for a, second for b

    enum Decision {
      Discard,  // the section is not to be explored further
      Accept,   // the section should be accepted as solution without further
                // exploration
      Drill,  // the section should be expred further by splitting according to
              // binning definition (split into 4 or 2 left-right or top-bottom)
      Explode,  // the section should be source of 5 subsections as in drill &
                // one in the center (unimplemented)
      Custom    // the custom functor should be used to create subsections
                // (unimplemented)
    };

    using DecisionFunctor = std::function<Decision(
        const AccumulatorSection &section,
        const std::vector<PreprocessedMeasurement> &measurements)>;
    DecisionFunctor decisionFunctor;  // function deciding if the Accumulator
                                      // section should be, discarded, split
                                      // further (and how), or is a solution
  };

  template <typename M>
  void exploreParametersSpace(std::stack<AccumulatorSection> &sectionsStack,
                              const std::vector<M> &measurements,
                              const AHTExplorationOptions<M> &opt,
                              std::vector<AccumulatorSection> &results) const {
    using Decision = AHTExplorationOptions<M>::Decision;
    while (!sectionsStack.empty()) {
      ACTS_VERBOSE("Stack size " << sectionsStack.size());
      AccumulatorSection &thisSection = sectionsStack.top();
      Decision whatNext = opt.decisionFunctor(thisSection, measurements);
      ACTS_VERBOSE("top section "
                   << thisSection.count() << " section " << thisSection.xBegin()
                   << " - " << thisSection.xBegin() + thisSection.xSize() << " "
                   << thisSection.yBegin() << " - "
                   << thisSection.yBegin() + thisSection.ySize()
                   << " nlines: " << thisSection.count() << " div: "
                   << thisSection.divisionLevel() << " decision " << whatNext);

      if (whatNext == Decision::Discard) {
        sectionsStack.pop();
      } else if (whatNext == Decision::Accept) {
        addSolution(std::move(thisSection), results);
        sectionsStack.pop();
      } else {
        // further exploration starts here
        std::vector<AccumulatorSection> divisions;
        if (thisSection.xSize() > opt.xMinBinSize &&
            thisSection.ySize() > opt.yMinBinSize) {
          // need 4 subdivisions
          divisions.push_back(thisSection.topLeft());
          divisions.push_back(thisSection.topRight());
          divisions.push_back(thisSection.bottomLeft());
          divisions.push_back(thisSection.bottomRight());
        } else if (thisSection.xSize() <= opt.xMinBinSize &&
                   thisSection.ySize() > opt.yMinBinSize) {
          // only split in y
          divisions.push_back(thisSection.top());
          divisions.push_back(thisSection.bottom());
        } else {
          // only split in x
          divisions.push_back(thisSection.left());
          divisions.push_back(thisSection.right());
        }
        sectionsStack.pop();  // discard the section that was just split
        for (AccumulatorSection &d : divisions) {
          std::vector<unsigned> selectedIndices;
          for (unsigned index : d.indices()) {
            const PreprocessedMeasurement &m = measurements[index];
            if (d.isLineInside(opt.lineParamFunctors.first(m),
                               opt.lineParamFunctors.second(m))) {
              selectedIndices.push_back(index);
            }
          }
          d.indices() = std::move(selectedIndices);
          sectionsStack.push(std::move(d));
        }
      }
    }
  }

  /// @brief  add solution to the solutions vector
  /// depending on options it may eliminate trivial duplicates
  /// @param s - the solution to be potentially added
  /// @param solutions - the output solutions set
  void addSolution(AccumulatorSection &&s,
                   std::vector<AccumulatorSection> &output) const;

  /// Construct the seeding algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  AdaptiveHoughTransformSeeder(Config cfg, Acts::Logging::Level lvl);

  /// Run the seeding algorithm.
  ///
  /// @param txt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext &ctx) const final;

  /// Const access to the config
  const Config &config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;
  const Acts::Logger &logger() const { return *m_logger; }

  WriteDataHandle<ProtoTrackContainer> m_outputProtoTracks{this,
                                                           "OutputProtoTracks"};

  WriteDataHandle<SimSeedContainer> m_outputSeeds{this, "OutputSeeds"};

  std::vector<std::unique_ptr<ReadDataHandle<SimSpacePointContainer>>>
      m_inputSpacePoints{};

  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};
  /// @brief process sections on the stack
  /// and qualifying them for further division, discarding them or moving to
  /// solutions vector
  /// @param sections is the stack of sectoins to consider
  /// @param solutions is the output set of sections
  /// @param measurements are input measurements
  void processStackHeadQOverPtPhi(
      std::stack<AccumulatorSection> &sections,
      std::vector<AccumulatorSection> &solutions,
      const std::vector<PreprocessedMeasurement> &measurements) const;

  void processStackHeadZCotTheta(
      std::stack<AccumulatorSection> &sections,
      std::vector<AccumulatorSection> &solutions,
      const std::vector<PreprocessedMeasurement> &measurements) const;

  using LineParamFunctors =
      std::pair<std::function<float(const PreprocessedMeasurement &)>,
                std::function<float(const PreprocessedMeasurement &)>>;

  LineParamFunctors m_qOverPtPhiLineParams = {
      [this](const PreprocessedMeasurement &m) {
        return m.invr * config().inverseA;
      },
      [this](const PreprocessedMeasurement &m) {
        return -m.invr * m.phi * config().inverseA;
      }};

  LineParamFunctors m_zCotThetaLineParams = {
      [this](const PreprocessedMeasurement &m) { return -m.invr; },
      [this](const PreprocessedMeasurement &m) { return m.z * m.invr; }};

  /// @brief check if lines intersect in the section
  /// modifies the section leaving only indices of measurements that do so
  /// @param section - the section to check
  /// @param measurements - the measurements that are pointed to by indices in
  /// @param lineParamsAccessor - functions to be used to access line parameters
  /// @param threshold - the number of lines in the section should be at minimum
  bool passIntersectionsCheck(
      const AccumulatorSection &section,
      const std::vector<PreprocessedMeasurement> &measurements,
      const LineParamFunctors &lineParamsAccessor,
      const unsigned threshold) const;
};

}  // namespace ActsExamples
