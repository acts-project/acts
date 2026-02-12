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
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

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

  inline unsigned int count() const {
    return static_cast<unsigned>(m_indices.size());
  }
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

  /// @brief true if the line defined by given parameters passes the section
  /// @param function is callable used to check crossing at the edges
  template <typename F>
  inline bool isLineInside(F &&function) const &
    requires std::invocable<F, float>
  {
    const float yB = function(m_xBegin);
    const float yE = function(m_xBegin + m_xSize);
    return (yE > yB) ? yB < m_yBegin + m_ySize && yE > m_yBegin
                     : yB > m_yBegin && yE < m_yBegin + m_ySize;
  }

  /// @brief check if the lines cross inside the section
  /// @brief line1 - functional form of line 1
  /// @brief line1 - functional form of line 2
  /// @warning note that this function is assuming that these are lines
  /// It may be incorrect assumption for rapidly changing function or large
  /// sections
  /// @return true if the lines cross in the section
  template <typename F>
  inline bool isCrossingInside(F &&line1, F &&line2) const &
    requires std::invocable<F, float>
  {
    // this microalgorithm idea is illustrated below
    // section left section right
    // example with crossing
    //                                       |            +2
    // line 1 crossing left section edge     +1          _|
    // left edge mid point                   |_           |
    //                                       |            +1
    // line 2crossing left section           +2           |
    //
    // example with no crossing
    //                                       |            +1
    // line 1 crossing left section edge     +1          _|
    // left edge mid point                   |_           |
    //                                       |            +2
    // line 2crossing left section           +2           |
    //
    // if for any of the two lines the condition
    // (line1_left_y-middle_on_the_left_y)*(line1_right_y-middle_on_the_right_y)
    // < 0 means that there is crossing

    float line1_left_y = line1(xBegin());
    float line1_right_y = line1(xBegin() + xSize());
    float line2_left_y = line2(xBegin());
    float line2_right_y = line2(xBegin() + xSize());
    float left_mid = 0.5f * (line1_left_y + line2_left_y);
    float right_mid = 0.5f * (line1_right_y + line2_right_y);
    return (line1_left_y - left_mid) * (line1_right_y - right_mid) < 0;
  }

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
    m_history.resize(index + 1);
    m_history.at(index) = value;
  }
  /// @brief retrieve history info
  /// @param index - item index
  /// @return value stored by @see setHistory
  float history(unsigned index) const { return m_history[index]; }

 private:
  float m_xSize = 0;
  float m_ySize = 0;
  float m_xBegin = 0;
  float m_yBegin = 0;
  unsigned m_divisionLevel =
      0;  // number of times the starting section was already divided
  std::vector<unsigned>
      m_indices;  // indices of measurements contributing to this section
  std::vector<float> m_history;  // additional record where an arbitrary
                                 // information can be stored
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
        0.1f;  // wrap size around angle domain limits (min pT dependent)
    float qOverPtMin = 0.9f;  // min q/pt, -1/1 GeV

    float qOverPtMinBinSize = 0.01f;  // minimal size of pT bin that the
                                      // algorithm should not go beyond (in GeV)
    float phiMinBinSize = 0.01f;  // minimal size of phi bin that the algorithm
                                  // should not go beyond
    float zRange = 200;           // range in z
    float cotThetaRange = 10;     // range in cotTheta

    float zMinBinSize =
        1;  // minimal size of z bin that the algorithm should
            // not go beyond when exploring zvertex-cot(theta) space space
    float cotThetaMinBinSize = 0.1f;  // minimal size of cot(theta) bin that
                                      // the algorithm should not go beyond
    unsigned threshold =
        4;  // number of lines passing section for it to be still considered
    unsigned noiseThreshold = 12;  // number of lines passing section at the
                                   // final split to consider it noise

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

  template <typename measurement_t = PreprocessedMeasurement>
  struct ExplorationOptions {
    float xMinBinSize = 1.0f;  // minimum bin size in x direction, beyond that
                               // value the sections are not split
    float yMinBinSize = 1.0f;  // minimum bin size in y direction, beyond that
                               // value the sections are not split
    float expandX = 1.1f;  // expand in x direaction (default by 10%) if Expand
                           // Decision is made
    float expandY = 1.1f;  // expand in y direaction (default by 10%) if Expand
                           // Decision is made
    using LineParamFunctor =
        std::function<float(const measurement_t &, float arg)>;
    LineParamFunctor lineParamFunctor;  // pair of functions needed to obtain
                                        // linear function ax+b parameters,
                                        // first for a, second for b

    enum class Decision {
      Discard,  // the section is not to be explored further
      Accept,   // the section should be accepted as solution without further
                // exploration
      Drill,  // the section should be expred further by splitting according to
              // binning definition (split into 4 or 2 left-right or top-bottom)
      DrillAndExpand,  // the section should be source of subsections as in the
                       // case of drill & but the sections will be made larger
      // size increase is configured in opt by relative factors @see expandX, @see expandY
    };

    using DecisionFunctor = std::function<Decision(
        const AccumulatorSection &section,
        const std::vector<PreprocessedMeasurement> &measurements)>;
    DecisionFunctor decisionFunctor;  // function deciding if the Accumulator
                                      // section should be, discarded, split
                                      // further (and how), or is a solution
  };

  /// @brief  remove indices pointing to measurements that do not cross this section
  /// @tparam measurement_t - measurements type
  /// @param section - section to update
  /// @param measurements - measurements vector
  /// @param lineFunctor - line definition
  template <typename measurement_t>
  void updateSection(AccumulatorSection &section,
                     const std::vector<measurement_t> &measurements,
                     const ExplorationOptions<measurement_t>::LineParamFunctor
                         &lineFunctor) const {
    std::erase_if(section.indices(), [lineFunctor, &measurements,
                                      &section](unsigned index) {
      const PreprocessedMeasurement &m = measurements[index];
      return !section.isLineInside(std::bind_front(lineFunctor, std::cref(m)));
    });
  }

  template <typename M>
  void exploreParametersSpace(std::vector<AccumulatorSection> &sectionsStack,
                              const std::vector<M> &measurements,
                              const ExplorationOptions<M> &opt,
                              std::vector<AccumulatorSection> &results) const {
    using Decision = ExplorationOptions<M>::Decision;
    while (!sectionsStack.empty()) {
      ACTS_VERBOSE("Stack size " << sectionsStack.size());
      AccumulatorSection &thisSection = sectionsStack.back();
      Decision whatNext = opt.decisionFunctor(thisSection, measurements);
      ACTS_VERBOSE("top section "
                   << thisSection.count() << " section " << thisSection.xBegin()
                   << " - " << thisSection.xBegin() + thisSection.xSize() << " "
                   << thisSection.yBegin() << " - "
                   << thisSection.yBegin() + thisSection.ySize()
                   << " nlines: " << thisSection.count()
                   << " div: " << thisSection.divisionLevel() << " decision "
                   << (whatNext == Decision::Discard ? "Discard"
                                                     : "Drill, Accept"));

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
          divisions.reserve(4);
          divisions.push_back(thisSection.topLeft());
          divisions.push_back(thisSection.topRight());
          divisions.push_back(thisSection.bottomLeft());
          divisions.push_back(thisSection.bottomRight());
        } else if (thisSection.xSize() <= opt.xMinBinSize &&
                   thisSection.ySize() > opt.yMinBinSize) {
          // only split in y
          divisions.reserve(2);
          divisions.push_back(thisSection.top());
          divisions.push_back(thisSection.bottom());
        } else {
          // only split in x
          divisions.reserve(2);
          divisions.push_back(thisSection.left());
          divisions.push_back(thisSection.right());
        }

        if (whatNext == Decision::DrillAndExpand) {
          for (AccumulatorSection &d : divisions) {
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
  AdaptiveHoughTransformSeeder(const Config &cfg, Acts::Logging::Level lvl);

  /// Run the seeding algorithm.
  ///
  /// @param txt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext &ctx) const override;

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

  /// @brief fill vector pf measurements from input space points
  /// @param measurements - vector to fill
  void preparePreprocessedMeasurements(
      const AlgorithmContext &ctx,
      std::vector<PreprocessedMeasurement> &measurements) const;

  /// @brief split the measurements into sections in phi
  /// The pT is not take into account. An overlap in phi is assured.
  /// @param stack - sections stack to fill
  /// @param measurements - measurements to fill the stack
  void fillStackPhiSplit(
      std::vector<AccumulatorSection> &stack,
      const std::vector<PreprocessedMeasurement> &measurements) const;

  /// @brief process sections on the stack
  /// and qualifying them for further division, discarding them or moving to
  /// solutions vector
  /// the search happens in q/pt - phi space
  /// @param sections is the stack of sectoins to consider
  /// @param solutions is the output set of sections
  /// @param measurements are input measurements
  void processStackQOverPtPhi(
      std::vector<AccumulatorSection> &input,
      std::vector<AccumulatorSection> &output,
      const std::vector<PreprocessedMeasurement> &measurements) const;

  /// @brief process sections on the stack
  /// and qualifying them for further division, discarding them or moving to
  /// solutions vector
  /// the search happens in q/pt - phi space
  /// @param sections is the stack of sectoins to consider
  /// @param solutions is the output set of sections
  /// @param measurements are input measurements
  void processStackZCotTheta(
      std::vector<AccumulatorSection> &input,
      std::vector<AccumulatorSection> &output,
      const std::vector<PreprocessedMeasurement> &measurements) const;

  void processStackZCotThetaSplit(
      std::vector<AccumulatorSection> &input,
      std::vector<AccumulatorSection> &output,
      const std::vector<PreprocessedMeasurement> &measurements) const;

  /// @brief produce 3 Sp seeds out of solutions
  /// @param seeds output to fill
  /// @param solutions is the input to be translated
  /// @param measurements are input measurements
  void makeSeeds(
      SimSeedContainer &seeds, const std::vector<AccumulatorSection> &solutions,
      const std::vector<PreprocessedMeasurement> &measurements) const;

  using LineParamFunctor =
      ExplorationOptions<PreprocessedMeasurement>::LineParamFunctor;

  LineParamFunctor m_qOverPtPhiLineParams =
      [this](const PreprocessedMeasurement &m, float arg) {
        return m.invr * config().inverseA * arg -
               m.invr * m.phi * config().inverseA;
      };

  LineParamFunctor m_zCotThetaLineParams = [](const PreprocessedMeasurement &m,
                                              float arg) {
    return -m.invr * arg + m.z * m.invr;
  };

  /// @brief check if lines intersect in the section
  /// modifies the section leaving only indices of measurements that do so
  /// @param section - the section to check
  /// @param measurements - the measurements that are pointed to by indices in
  /// @param lineParamsAccessor - functions to be used to access line parameters
  /// @param threshold - the number of lines in the section should be at minimum
  bool passIntersectionsCheck(
      const AccumulatorSection &section,
      const std::vector<PreprocessedMeasurement> &measurements,
      const LineParamFunctor &lineFunctor, const unsigned threshold) const;

  void deduplicate(std::vector<AccumulatorSection> &input) const;
};

}  // namespace ActsExamples
