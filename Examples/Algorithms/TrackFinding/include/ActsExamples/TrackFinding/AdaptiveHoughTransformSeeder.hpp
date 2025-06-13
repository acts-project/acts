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
#include <list>
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

  AccumulatorSection(float xw, float yw, float xBegin, float yBegin,
                     int div = 0, const std::vector<unsigned> &indices = {},
                     const std::vector<float> &history = {});

  /// @brief keep indices and update parameters of the box
  /// This method is useful when changing direction of the search                     
  void updateDimensions(float xw, float yw, float xBegin, float yBegin);

  /// @brief keep indices and update parameters of the box by scalling
  /// @param xs - scale in x direction, if bigger than 1 the size increases
  /// @param ys - scale in y direction
  /// The box is recentred
  void expand(float xs, float ys);


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
  inline bool isLineInside(std::function<float(float)> f) const {
    const float yB = f(m_xBegin);
    const float yE = f(m_xBegin + m_xSize);
    return (yE > yB) ? yB < m_yBegin + m_ySize && yE > m_yBegin
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
  float xSize() const { return m_xSize; }
  float ySize() const { return m_ySize; }
  float xBegin() const { return m_xBegin; }
  float yBegin() const { return m_yBegin; }
  unsigned divisionLevel() const { return m_divisionLevel; }

  /// store additional (arbitrary) info in indexed array
  /// @param index - identifier
  /// @param value - value to store
  void setHistory(unsigned index, float value) { 
    m_history.resize(index+1); m_history[index] = value; 
  }
  /// @brief retrieve history info
  /// @param index - item index
  /// @return value stored by @see setHistory
  float history(unsigned index) const {
    return m_history[index];
  }

  private:
  float m_xSize = 0;
  float m_ySize = 0;
  float m_xBegin = 0;
  float m_yBegin = 0;
  unsigned m_divisionLevel =
      0;  // number of times the starting section was already divided
  std::vector<unsigned>
      m_indices;  // indices of measurements contributing to this section
  std::vector<float> m_history; // additional record where an arbitrary information can be stored
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
    float qOverPtMin = 0.9;  // min q/pt, -1/1 GeV

    float qOverPtMinBinSize = 0.01;  // minimal size of pT bin that the
                                     // algorithm should not go beyond (in GeV)
    float phiMinBinSize = 0.01;  // minimal size of phi bin that the algorithm
                                 // should not go beyond
    float zRange = 200;          // range in z
    float cotThetaRange = 10;    // range in cotTheta

    float zMinBinSize =
        1;  // minimal size of z bin that the algorithm should
            // not go beyond when exploring zvertex-cot(theta) space space
    float cotThetaMinBinSize = 0.1;  // minimal size of cot(theta) bin that
                                     // the algorithm should not go beyond
    unsigned threshold =
        4;  // number of lines passing section for it to be still considered
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
    float expandX = 1.1;    // expand in x direaction (default by 10%) if Expand Decision is made
    float expandY = 1.1;    // expand in y direaction (default by 10%) if Expand Decision is made
    using LineParamFunctor = std::function<float(const M &, float arg)>;
    LineParamFunctor lineParamFunctor;  // pair of functions needed to obtain
                                          // linear function ax+b parameters,
                                          // first for a, second for b

    enum Decision {
      Discard,  // the section is not to be explored further
      Accept,   // the section should be accepted as solution without further
                // exploration
      Drill,  // the section should be expred further by splitting according to
              // binning definition (split into 4 or 2 left-right or top-bottom)
      DrillAndExpand,  // the section should be source of subsections as in the case of drill &
                       // but the sections will be made larger 
                       // size increase is configured in opt by relative factors @see expandX, @see expandY
    };

    using DecisionFunctor = std::function<Decision(
        const AccumulatorSection &section,
        const std::vector<PreprocessedMeasurement> &measurements)>;
    DecisionFunctor decisionFunctor;  // function deciding if the Accumulator
                                      // section should be, discarded, split
                                      // further (and how), or is a solution
  };

  template <typename M>
  void updateSection(AccumulatorSection &section, 
    const std::vector<M> &measurements, 
    const AHTExplorationOptions<M>::LineParamFunctor& lineFunctor) const {
    std::vector<unsigned> selectedIndices;
    for (unsigned index : section.indices()) {
      const PreprocessedMeasurement &m = measurements[index];
      using namespace std::placeholders;
      if (section.isLineInside( std::bind(lineFunctor, m, _1))) {
        selectedIndices.push_back(index);
      }
    }
    section.indices() = std::move(selectedIndices);
  }


  template <typename M>
  void exploreParametersSpace(std::deque<AccumulatorSection> &sectionsStack,
                              const std::vector<M> &measurements,
                              const AHTExplorationOptions<M> &opt,
                              std::deque<AccumulatorSection> &results) const {
    using Decision = AHTExplorationOptions<M>::Decision;
    while (!sectionsStack.empty()) {
      ACTS_VERBOSE("Stack size " << sectionsStack.size());
      AccumulatorSection &thisSection = sectionsStack.back();
      Decision whatNext = opt.decisionFunctor(thisSection, measurements);
      ACTS_VERBOSE("top section "
                   << thisSection.count() << " section " << thisSection.xBegin()
                   << " - " << thisSection.xBegin() + thisSection.xSize() << " "
                   << thisSection.yBegin() << " - "
                   << thisSection.yBegin() + thisSection.ySize()
                   << " nlines: " << thisSection.count() << " div: "
                   << thisSection.divisionLevel() << " decision " << whatNext);

      if (whatNext == Decision::Discard) {
        sectionsStack.pop_back();
      } else if (whatNext == Decision::Accept) {
        results.push_back(std::move(thisSection));
        sectionsStack.pop_back();
      } else {
        // further exploration starts here
        std::vector<AccumulatorSection> divisions;
        if (thisSection.xSize() > opt.xMinBinSize &&
            thisSection.ySize() > opt.yMinBinSize) {
          // need 4 subdivisions
          divisions.push_back(std::move(thisSection.topLeft()));
          divisions.push_back(std::move(thisSection.topRight()));
          divisions.push_back(std::move(thisSection.bottomLeft()));
          divisions.push_back(std::move(thisSection.bottomRight()));
        } else if (thisSection.xSize() <= opt.xMinBinSize &&
                   thisSection.ySize() > opt.yMinBinSize) {
          // only split in y
          divisions.push_back(std::move(thisSection.top()));
          divisions.push_back(std::move(thisSection.bottom()));
        } else {
          // only split in x
          divisions.push_back(std::move(thisSection.left()));
          divisions.push_back(std::move(thisSection.right()));
        }

        if ( whatNext == Decision::DrillAndExpand ) {
          for ( AccumulatorSection& d: divisions ) {
            d.expand(opt.expandX, opt.expandY);
          }
        }
        sectionsStack.pop_back();  // discard the section that was just split
        for (AccumulatorSection &d : divisions) {
          updateSection(d, measurements, opt.lineParamFunctor);
          sectionsStack.push_back(std::move(d));
        }
      }
    }
  }


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
  void processStackQOverPtPhi(
      std::deque<AccumulatorSection> &input,
      std::deque<AccumulatorSection> &output,
      const std::vector<PreprocessedMeasurement> &measurements) const;

  void processStackZCotTheta(
      std::deque<AccumulatorSection> &input,
      std::deque<AccumulatorSection> &output,
      const std::vector<PreprocessedMeasurement> &measurements) const;

  void processStackZCotThetaSplit(
      std::deque<AccumulatorSection> &input,
      std::deque<AccumulatorSection> &output,
      const std::vector<PreprocessedMeasurement> &measurements) const;


  using LineParamFunctor = AHTExplorationOptions<PreprocessedMeasurement>::LineParamFunctor;

  LineParamFunctor m_qOverPtPhiLineParams = 
      [this](const PreprocessedMeasurement &m, float arg) {
            return m.invr * config().inverseA*arg - m.invr * m.phi * config().inverseA;
      };

  LineParamFunctor m_zCotThetaLineParams = 
      [this](const PreprocessedMeasurement &m, float arg) { 
            return -m.invr * arg + m.z * m.invr; };

  /// @brief check if lines intersect in the section
  /// modifies the section leaving only indices of measurements that do so
  /// @param section - the section to check
  /// @param measurements - the measurements that are pointed to by indices in
  /// @param lineParamsAccessor - functions to be used to access line parameters
  /// @param threshold - the number of lines in the section should be at minimum
  bool passIntersectionsCheck(
      const AccumulatorSection &section,
      const std::vector<PreprocessedMeasurement> &measurements,
      const LineParamFunctor &lineParamsAccessor,
      const unsigned threshold) const;

  void deduplicate(std::deque<AccumulatorSection>& input) const;
};

}  // namespace ActsExamples
