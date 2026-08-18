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

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Seeding/HoughAccumulatorSection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Seed.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace ActsExamples {

// Construct track seeds from space points.
class AdaptiveHoughTransformSeeder final : public IAlgorithm {
 public:
  struct Config {
    /// Input space point collections.
    std::string inputSpacePoints;
    /// Output track seed collection.
    std::string outputSeeds;
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
    /// @param spacePointIndex index of the original space point, needed to construct seeds in the end
    PreprocessedMeasurement(double inverseR, double phiAngle, double zpos,
                            SpacePointIndex spacePointIndex)
        : invr(inverseR), phi(phiAngle), z(zpos), sp(spacePointIndex) {}
    double invr{};
    double phi{};
    double z{};
    SpacePointIndex sp{};
  };
  using HoughAccumulatorSection = Acts::Experimental::HoughAccumulatorSection;
  using ExplorationOptions =
      Acts::Experimental::HoughExplorationOptions<PreprocessedMeasurement>;

  using LineFunctor = ExplorationOptions::LineFunctor;

  /// @brief  remove indices pointing to measurements that do not cross this section
  /// @tparam measurement_t - measurements type
  /// @param section - section to update
  /// @param measurements - measurements vector
  /// @param lineFunctor - line definition
  template <typename measurement_t>
  void updateSection(Acts::Experimental::HoughAccumulatorSection &section,
                     const std::vector<measurement_t> &measurements,
                     const LineFunctor &lineFunctor) const {
    std::erase_if(section.indices(), [lineFunctor, &measurements,
                                      &section](unsigned index) {
      const PreprocessedMeasurement &m = measurements[index];
      return !section.isLineInside(std::bind_front(lineFunctor, std::cref(m)));
    });
  }

  /// Construct the seeding algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  explicit AdaptiveHoughTransformSeeder(
      const Config &cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// Run the seeding algorithm.
  ///
  /// @param txt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext &ctx) const override;

  /// Const access to the config
  const Config &config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<SpacePointContainer> m_inputSpacePoints{this,
                                                         "InputSpacePoints"};

  WriteDataHandle<SeedContainer> m_outputSeeds{this, "OutputSeeds"};

  /// @brief fill vector pf measurements from input space points
  /// @param measurements - vector to fill
  void preparePreprocessedMeasurements(
      const SpacePointContainer &spacePoints,
      std::vector<PreprocessedMeasurement> &measurements) const;

  /// @brief split the measurements into sections in phi
  /// The pT is not take into account. An overlap in phi is assured.
  /// @param stack - sections stack to fill
  /// @param measurements - measurements to fill the stack
  void fillStackPhiSplit(
      std::vector<Acts::Experimental::HoughAccumulatorSection> &stack,
      const std::vector<PreprocessedMeasurement> &measurements) const;

  /// @brief process sections on the stack
  /// and qualifying them for further division, discarding them or moving to
  /// solutions vector
  /// the search happens in q/pt - phi space
  /// @param sections is the stack of sectoins to consider
  /// @param solutions is the output set of sections
  /// @param measurements are input measurements
  void processStackQOverPtPhi(
      std::vector<Acts::Experimental::HoughAccumulatorSection> &input,
      std::vector<Acts::Experimental::HoughAccumulatorSection> &output,
      const std::vector<PreprocessedMeasurement> &measurements) const;

  /// @brief process sections on the stack
  /// and qualifying them for further division, discarding them or moving to
  /// solutions vector
  /// the search happens in q/pt - phi space
  /// @param sections is the stack of sectoins to consider
  /// @param solutions is the output set of sections
  /// @param measurements are input measurements
  void processStackZCotTheta(
      std::vector<Acts::Experimental::HoughAccumulatorSection> &input,
      std::vector<Acts::Experimental::HoughAccumulatorSection> &output,
      const std::vector<PreprocessedMeasurement> &measurements) const;

  void processStackZCotThetaSplit(
      std::vector<Acts::Experimental::HoughAccumulatorSection> &input,
      std::vector<Acts::Experimental::HoughAccumulatorSection> &output,
      const std::vector<PreprocessedMeasurement> &measurements) const;

  /// @brief produce 3 Sp seeds out of solutions
  /// @param seeds output to fill
  /// @param solutions is the input to be translated
  /// @param measurements are input measurements
  void makeSeeds(
      SeedContainer &seeds,
      const std::vector<Acts::Experimental::HoughAccumulatorSection> &solutions,
      const std::vector<PreprocessedMeasurement> &measurements) const;

  LineFunctor m_qOverPtPhiLineParams = [this](const PreprocessedMeasurement &m,
                                              float arg) {
    return m.invr * config().inverseA * arg -
           m.invr * m.phi * config().inverseA;
  };

  LineFunctor m_zCotThetaLineParams = [](const PreprocessedMeasurement &m,
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
      const Acts::Experimental::HoughAccumulatorSection &section,
      const std::vector<PreprocessedMeasurement> &measurements,
      const LineFunctor &lineFunctor, const unsigned threshold) const;

  void deduplicate(
      std::vector<Acts::Experimental::HoughAccumulatorSection> &input) const;
};

}  // namespace ActsExamples
