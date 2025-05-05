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
} // namespace ActsExamples

using ResultDouble = Acts::Result<double>;
using ResultBool = Acts::Result<bool>;
using ResultUnsigned = Acts::Result<unsigned>;

using FieldCorrector = Acts::Delegate<ResultDouble(
    unsigned, double, double)>; // (unsigned region, double y, double r)
using LayerIDFinder = Acts::Delegate<ResultUnsigned(
    double)>; // (double r) this function will map the r of a measurement to a
              // layer.
using SliceTester = Acts::Delegate<ResultBool(
    double, unsigned, int)>; // (double z,unsigned layer, int slice) returns
                             // true if measurement in slice

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

  // returns true if the line defined by given parameters passes the section
  // a and b are line parameters y = ax + b
  inline bool isLineInside(float a, float b) const {
    const float yB = std::fma(a, m_xBegin, b);
    const float yE = std::fma(a, (m_xBegin + m_xSize), b);
    return yB < m_yBegin + m_ySize && yE > m_yBegin;
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

private:
  double m_xSize;
  double m_ySize;
  double m_xBegin;
  double m_yBegin;
  unsigned m_divisionLevel =
      0; // number of times the starting section was already divided
  std::vector<unsigned>
      m_indices; // indices of measurements contributing to this section
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

    float qOverPtMin = 1.0; // min q/pt, -1/1 GeV

    float qOverPtMinBinSize = 0.01; // minimal size of pT bin that the
                                    // algorithm should not go beyond (in GeV)
    float phiMinBinSize = 0.01; // minimal size of phi bin that the algorithm
                                // should not go beyond

    unsigned threshold =
        8; // number of lines passing section for it to be still considered

    bool deduplicate = true; // when adding solutions try avoiding duplicates

    bool requireIntersections =
        true; // require that lines passing section need to cross inside
              // the count is required to be at lease threshold*(threshold-1):
    unsigned intersectionsThreshold =
        threshold * (threshold - 1) /
        2; // the number of lines in section should be at most this to enable
           // intersection test

    double inverseA =
        1.0 / 3.0e-4; // Assume B = 2T constant. Can apply corrections to
                      // this with fieldCorrection function
                      // This 3e-4 comes from the 2T field when converted to
                      // units of GeV / (c*mm*e)

    // it's up to the user to connect these to the functions they want to use
    FieldCorrector fieldCorrector;
    LayerIDFinder layerIDFinder;
    SliceTester sliceTester;
  };

  // information that is needed for each measurement
  struct PreprocessedMeasurement {
    /// Construct the measurement used internally
    ///
    /// @param ir inverse of radius
    /// @param p azimuthal angle
    /// @param l link to space point
    PreprocessedMeasurement(double ir, double p, Acts::SourceLink l)
        : invr(ir), phi(p), link(l) {}
    double invr;
    double phi;
    Acts::SourceLink link;
  };

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
  void processStackHead(
      std::stack<AccumulatorSection> &sections,
      std::vector<AccumulatorSection> &solutions,
      const std::vector<PreprocessedMeasurement> &measurements) const;

  /// @brief assign measurements to the section
  /// @warning the sections needs to have already indices of measurements to
  /// consider
  /// @warning from previous iteration
  /// @param sections section to be processed
  /// @param vector of measurements - indices in the section need to point to
  /// this vector
  void updateSection(AccumulatorSection &section,
                     const std::vector<PreprocessedMeasurement> &input) const;

  /// @brief check if lines intersect in the section
  /// modifies the section leaving only indices of measurements that do so
  /// @param section - the section to check
  /// @param measurements - the measurements that are pointed to by indices in
  /// section
  bool passIntersectionsCheck(
      const AccumulatorSection &section,
      const std::vector<PreprocessedMeasurement> &measurements) const;

  /// @brief  add solution to the solutions vector
  /// depending on options it may eliminate trivial duplicates
  /// @param s - the solution to be potentially added
  /// @param solutions - the output solutions set
  void addSolution(AccumulatorSection &&s,
                   std::vector<AccumulatorSection> &output) const;
};

} // namespace ActsExamples
