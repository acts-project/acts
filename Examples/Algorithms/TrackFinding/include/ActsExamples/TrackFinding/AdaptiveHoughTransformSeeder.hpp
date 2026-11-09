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
#include "Acts/Utilities/HashCombine.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Seed.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <unordered_set>
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
    float qOverPtMin = 1.1f;  // min q/pt, -1/1 GeV

    float qOverPtMinBinSize =
        0.055f;                    // minimal size of pT bin that the
                                   // algorithm should not go beyond (in GeV)
    float phiMinBinSize = 0.025f;  // minimal size of phi bin that the algorithm
                                   // should not go beyond
    float zRange = 200;            // range in z
    float cotThetaRange = 10;      // range in cotTheta

    float zMinBinSize =
        1;  // minimal size of z bin that the algorithm should
            // not go beyond when exploring zvertex-cot(theta) space space
    float cotThetaMinBinSize = 0.1f;  // minimal size of cot(theta) bin that
                                      // the algorithm should not go beyond
    unsigned threshold =
        3;  // number of lines passing section for it to be still considered
    unsigned noiseThreshold = 12;  // number of lines passing section at the
                                   // final split to consider it noise

    bool doSecondPhase = true;  // do the second pass in z-cot(theta) space to
                                // find less solutions
    bool deduplicate = true;    // when adding solutions try avoiding duplicates

    unsigned maxDivision = 6;  // Level of maximal standard section division

    double inverseA =
        1.0 / 3.0e-4;  // Assume B = 2T constant. Can apply corrections to
                       // this with fieldCorrection function
                       // This 3e-4 comes from the 2T field when converted to
                       // units of GeV / (c*mm*e)
  };

  /// @brief Lightweight representation of a 3-spacepoint seed candidate
  struct SeedCandidate {
    std::array<unsigned, 3> indices{};
  };

  /// @brief Intermediate candidate extracted from the first phase (q/pT - phi)
  /// carrying measurement indices and historical z - cot(theta) bounds
  struct Phase1Candidate {
    std::vector<unsigned> indices;
    float zMin{};
    float zWidth{};
    float cotThetaMin{};
    float cotThetaWidth{};
  };

  /// @brief Hash functor using hashMixAndCombine utility to enable O(1) deduplication
  struct VectorHash {
    std::size_t operator()(const std::vector<unsigned> &v) const {
      std::size_t seed = v.size();
      for (unsigned idx : v) {
        seed = Acts::hashMixAndCombine(seed, idx);
      }
      return static_cast<std::size_t>(seed);
    }
  };

  /// @brief Hash functor for std::array<unsigned, 3> using hashMixAndCombine
  struct Array3Hash {
    std::size_t operator()(const std::array<unsigned, 3> &a) const {
      return static_cast<std::size_t>(
          Acts::hashMixAndCombine(a[0], a[1], a[2]));
    }
  };

  // information that is needed for each measurement
  struct PreprocessedMeasurement {
    /// Construct the measurement
    /// @param radius radius of the SP
    /// @param inverseR inverse of radius of the SP
    /// @param phiAngle azimuthal angle of the SP
    /// @param zpos z position of the SP
    /// @param spacePointIndex index of the original space point, needed to construct seeds in the end
    PreprocessedMeasurement(double radius, double inverseR, double phiAngle,
                            double zpos, SpacePointIndex spacePointIndex)
        : r(radius),
          invr(inverseR),
          phi(phiAngle),
          z(zpos),
          sp(spacePointIndex) {}
    double r{};
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
  /// @param logger is the logging instance
  explicit AdaptiveHoughTransformSeeder(
      const Config &cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// Run the seeding algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext &ctx) const override;

  /// Const access to the config
  const Config &config() const { return m_cfg; }

 private:
  Config m_cfg;
  float m_bFieldA{};

  ReadDataHandle<SpacePointContainer> m_inputSpacePoints{this,
                                                         "InputSpacePoints"};

  WriteDataHandle<SeedContainer> m_outputSeeds{this, "OutputSeeds"};

  /// @brief fill vector of measurements from input space points
  /// @param spacePoints - input space points container
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
  /// @param input is the stack of sections to consider
  /// @param output is the output set of sections
  /// @param measurements are input measurements
  void processStackQOverPtPhi(
      std::vector<Acts::Experimental::HoughAccumulatorSection> &input,
      std::vector<Acts::Experimental::HoughAccumulatorSection> &output,
      const std::vector<PreprocessedMeasurement> &measurements) const;

  /// @brief 1D sliding window search for clusters in q/pt - phi space
  /// @param input is the stack of sections to consider in 1D cluster search
  /// @param output is the output set of candidates after the first phase
  /// @param measurements are input measurements
  void processQOverPtPhiSlidingWindow(
      const std::vector<Acts::Experimental::HoughAccumulatorSection> &input,
      std::vector<Phase1Candidate> &output,
      const std::vector<PreprocessedMeasurement> &measurements) const;

  /// @brief Evaluates lines at discrete check points in q/pT - phi space and
  /// extracts 1D clusters into Phase1Candidates.
  /// @param section current accumulator section being evaluated
  /// @param checkPoints discrete q/pT coordinates at which phi intersections are evaluated
  /// @param measurements global collection of preprocessed measurements
  /// @param lineFunctor functor evaluating the phi coordinate for a given q/pT coordinate
  /// @param maxSpread maximum allowed spread along the phi axis to form a cluster
  /// @param uniqueClusters hash set tracking already discovered index combinations
  /// @param output output vector of unique Phase1Candidates
  void findClustersQOverPtPhi(
      const Acts::Experimental::HoughAccumulatorSection &section,
      std::span<const float> checkPoints,
      const std::vector<PreprocessedMeasurement> &measurements,
      const LineFunctor &lineFunctor, float maxSpread,
      std::unordered_set<std::vector<unsigned>, VectorHash> &uniqueClusters,
      std::vector<Phase1Candidate> &output) const;

  /// @brief process candidates from Phase 1 in z - cot(theta) space
  /// @param input is the vector of Phase 1 candidates to consider
  /// @param output is the output set of 3-point seed candidates
  /// @param measurements are input measurements
  void processStackZCotTheta(
      std::vector<Phase1Candidate> &input, std::vector<SeedCandidate> &output,
      const std::vector<PreprocessedMeasurement> &measurements) const;

  /// @brief Evaluates lines at discrete check points in z - cot(theta) space
  /// and extracts radially separated 3-point seed candidates.
  /// @param indices measurement indices present in the current candidate
  /// @param checkPoints discrete z coordinates to evaluate
  /// @param measurements global collection of preprocessed measurements
  /// @param lineFunctor functor evaluating cot(theta) for a given z coordinate
  /// @param maxSpread maximum allowed spread along the cot(theta) axis to form a cluster
  /// @param uniqueClusters hash set tracking already discovered triplets
  /// @param output output vector of unique 3-point SeedCandidates
  void findClustersZCotTheta(
      const std::vector<unsigned> &indices, std::span<const float> checkPoints,
      const std::vector<PreprocessedMeasurement> &measurements,
      const LineFunctor &lineFunctor, float maxSpread,
      std::unordered_set<std::array<unsigned, 3>, Array3Hash> &uniqueClusters,
      std::vector<SeedCandidate> &output) const;

  void processStackZCotThetaSplit(
      std::vector<Acts::Experimental::HoughAccumulatorSection> &input,
      std::vector<Acts::Experimental::HoughAccumulatorSection> &output,
      const std::vector<PreprocessedMeasurement> &measurements) const;

  /// @brief produce 3 Sp seeds out of solutions
  /// @param seeds output to fill
  /// @param solutions is the input to be translated
  /// @param measurements are input measurements
  void makeSeeds(
      SeedContainer &seeds, const std::vector<SeedCandidate> &solutions,
      const std::vector<PreprocessedMeasurement> &measurements) const;

  /// @brief line functor to access q/pT - phi line parameters
  LineFunctor m_qOverPtPhiLineParams = [this](const PreprocessedMeasurement &m,
                                              float arg) {
    return m.invr * config().inverseA * arg -
           m.invr * m.phi * config().inverseA;
  };

  /// @brief line functor to access phi - q/pT line parameters
  LineFunctor m_phiQOverPtLineParams = [this](const PreprocessedMeasurement &m,
                                              float arg) {
    return static_cast<float>(m.phi + arg * m.r * m_bFieldA);
  };

  /// @brief line functor to access cot(theta) - z line parameters
  LineFunctor m_zCotThetaLineParams = [](const PreprocessedMeasurement &m,
                                         float arg) {
    return -m.invr * arg + m.z * m.invr;
  };

  /// @brief check if lines intersect in the section
  /// modifies the section leaving only indices of measurements that do so
  /// @param section - the section to check
  /// @param measurements - the measurements that are pointed to by indices in the section
  /// @param lineFunctor - functions to be used to access line parameters
  /// @param threshold - the number of lines in the section should be at minimum
  bool passIntersectionsCheck(
      const Acts::Experimental::HoughAccumulatorSection &section,
      const std::vector<PreprocessedMeasurement> &measurements,
      const LineFunctor &lineFunctor, const unsigned threshold) const;

  void deduplicate(
      std::vector<Acts::Experimental::HoughAccumulatorSection> &input) const;
};

}  // namespace ActsExamples
