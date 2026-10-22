// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <vector>

namespace Acts {

/// One doublet of a middle space point, found by the doublet seed finder.
struct DoubletForMiddleSp {
  /// Index of the doublet's other space point
  SpacePointIndex spacePointIndex{};
  /// Cotangent of the doublet's polar angle
  float cotTheta{};
  /// Inverse of the two space points' distance in the transverse plane
  float iDeltaR{};
  /// Squared uncertainty of cotTheta
  float er{};
  /// Conformal coordinate u of the other space point
  float u{};
  /// Conformal coordinate v of the other space point
  float v{};
  /// Position x of the other space point relative to the middle one, in the
  /// frame rotated to the middle one's direction
  float x{};
  /// Position y of the other space point, as x
  float y{};
};

/// The doublets of one middle space point.
using DoubletsForMiddleSp = std::vector<DoubletForMiddleSp>;

/// Pair of doublet index and cotTheta value.
struct IndexAndCotTheta {
  /// Doublet index
  std::uint32_t index{};
  /// Cotangent of theta
  float cotTheta{};
};

/// Storage for sortByCotTheta, kept by the caller so that sorting does not
/// allocate
struct SortStorage {
  /// Doublets being sorted
  std::vector<IndexAndCotTheta> scratch;
  /// Bucket counts and indices
  std::vector<std::uint32_t> buckets;
};

/// Sort doublets by cotTheta. Doublets with equal cotTheta keep their index
/// order.
/// @param doublets The doublets to sort
/// @param indexAndCotTheta Output vector containing sorted index and cotTheta pairs
void sortByCotTheta(const DoubletsForMiddleSp& doublets,
                    std::vector<IndexAndCotTheta>& indexAndCotTheta);

/// Sort doublets by cotTheta, with storage the caller keeps between calls.
/// Doublets with equal cotTheta keep their index order.
/// @param doublets The doublets to sort
/// @param indexAndCotTheta Output vector containing sorted index and cotTheta pairs
/// @param storage Storage for the sort, of any content
void sortByCotTheta(const DoubletsForMiddleSp& doublets,
                    std::vector<IndexAndCotTheta>& indexAndCotTheta,
                    SortStorage& storage);

/// The doublets of a middle space point in the order of a list of their
/// indices, as sortByCotTheta gives it.
class Subspan {
 public:
  /// @param doublets The doublets
  /// @param indexAndCotTheta Their index and cotTheta pairs, in the order to
  ///   visit them
  Subspan(const DoubletsForMiddleSp& doublets,
          std::span<const IndexAndCotTheta> indexAndCotTheta)
      : m_doublets(&doublets), m_indexAndCotTheta(indexAndCotTheta) {}

  /// @return The number of doublets visited
  [[nodiscard]] std::size_t size() const { return m_indexAndCotTheta.size(); }
  /// @return Whether no doublet is visited
  [[nodiscard]] bool empty() const { return m_indexAndCotTheta.empty(); }
  /// @param i Position in the subspan
  /// @return The doublet at that position
  const DoubletForMiddleSp& operator[](std::size_t i) const {
    return (*m_doublets)[m_indexAndCotTheta[i].index];
  }
  /// @param offset Number of doublets to skip
  /// @return The doublets from position offset on
  Subspan subspan(std::size_t offset) const {
    return {*m_doublets, m_indexAndCotTheta.subspan(offset)};
  }

 private:
  const DoubletsForMiddleSp* m_doublets{};
  std::span<const IndexAndCotTheta> m_indexAndCotTheta;
};

/// Derived quantities for the middle space point in a doublet.
struct MiddleSpInfo {
  /// minus one over radius of middle SP
  float uIP{};
  /// square of uIP
  float uIP2{};
  /// ratio between middle SP x position and radius
  float cosPhiM{};
  /// ratio between middle SP y position and radius
  float sinPhiM{};
};

/// Interface and a collection of standard implementations for a doublet seed
/// finder. Given a starting space point and a collection of candidates, it
/// finds all doublets that satisfy the selection criteria. For the standard
/// implementations the criteria are given by interaction point cuts.
///
/// @note The standard implementations rely on virtual function dispatch which
/// did not turn out to affect the performance after measurement.
class DoubletSeedFinder {
 public:
  /// Collection of configuration parameters for the doublet seed finder. This
  /// includes doublet cuts, steering switches, and assumptions about the space
  /// points.
  struct Config {
    /// Whether the input space points are sorted by radius
    bool spacePointsSortedByRadius = false;

    /// Direction of the doublet candidate space points. Either forward, also
    /// called top doublet, or backward, also called bottom doublet.
    Direction candidateDirection = Direction::Forward();

    /// Minimum radial distance between two doublet components
    float deltaRMin = 5 * UnitConstants::mm;
    /// Maximum radial distance between two doublet components
    float deltaRMax = 270 * UnitConstants::mm;

    /// Minimal z distance between two doublet components
    float deltaZMin = -std::numeric_limits<float>::infinity();
    /// Maximum z distance between two doublet components
    float deltaZMax = std::numeric_limits<float>::infinity();

    /// Maximum value of impact parameter estimation of the seed candidates
    float impactMax = 20 * UnitConstants::mm;

    /// Enable cut on the compatibility between interaction point and doublet,
    /// this is an useful approximation to speed up the seeding
    bool interactionPointCut = false;
    /// Enable the cut on the time compatibility
    /// The space point container has to provide the
    /// `SpacePointColumns::Time` and `SpacePointColumns::VarianceT` columns
    bool useTime = false;

    /// Limiting location of collision region in z-axis used to check if doublet
    /// origin is within reasonable bounds
    float collisionRegionMin = -150 * UnitConstants::mm;
    /// Maximum collision region boundary in z-axis for doublet origin checks
    float collisionRegionMax = +150 * UnitConstants::mm;

    /// Maximum allowed cotTheta between two space-points in doublet, used to
    /// check if forward angle is within bounds
    float cotThetaMax = 10.01788;  // equivalent to eta = 3 (pseudorapidity)

    /// Minimum transverse momentum (pT) used to check the r-z slope
    /// compatibility of triplets with maximum multiple scattering effect
    /// (produced by the minimum allowed pT particle) + a certain uncertainty
    /// term. Check the documentation for more information
    /// https://acts.readthedocs.io/en/latest/core/reconstruction/pattern_recognition/seeding.html
    float minPt = 400 * UnitConstants::MeV;
    /// Parameter which can loosen the tolerance of the track seed to form a
    /// helix. This is useful for e.g. misaligned seeding.
    float helixCutTolerance = 1;
    /// Maximum allowed difference of the times of the two space points,
    /// corrected for the time of flight between them, in units of their
    /// combined time resolution. Only used when `useTime` is enabled.
    float timeCutNSigma = 5;

    /// If set, calibrate both strip space points against their chord after the
    /// z origin cut, with this tolerance in strip half lengths. The following
    /// cuts and the doublet use the calibrated positions.
    /// Requires `SpacePointColumns::StripCalibrationDetails`.
    std::optional<float> stripCalibrationTolerance;

    /// Type alias for delegate to apply experiment specific cuts during doublet
    /// finding
    using ExperimentCuts =
        Delegate<bool(const ConstSpacePointProxy& /*middle*/,
                      const ConstSpacePointProxy& /*other*/, float /*cotTheta*/,
                      bool /*isBottomCandidate*/)>;

    /// Delegate to apply experiment specific cuts during doublet finding
    ExperimentCuts experimentCuts;
  };

  /// Derived configuration for the doublet seed finder using a magnetic field.
  struct DerivedConfig : public Config {
    /// Constructor for derived configuration with magnetic field
    /// @param config Base configuration to derive from
    /// @param bFieldInZ Magnetic field strength in z-direction
    DerivedConfig(const Config& config, float bFieldInZ);

    /// Magnetic field strength in z-direction for helix calculation
    float bFieldInZ = std::numeric_limits<float>::quiet_NaN();
    /// Squared minimum helix diameter derived from magnetic field and minimum
    /// pT
    float minHelixDiameter2 = std::numeric_limits<float>::quiet_NaN();
    /// Squared time cut in units of sigma, allowing the cut to be applied on
    /// the combined time variance
    float timeCutNSigma2 = std::numeric_limits<float>::quiet_NaN();
  };

  /// Computes additional quantities from the middle space point which can be
  /// reused during doublet finding.
  /// @param spM Middle space point for doublet computation
  /// @return MiddleSpInfo structure with computed quantities
  static MiddleSpInfo computeMiddleSpInfo(const ConstSpacePointProxy& spM);

  /// Creates a new doublet seed finder instance given the configuration.
  /// @param config Configuration for the doublet seed finder
  /// @return Unique pointer to new DoubletSeedFinder instance
  static std::unique_ptr<DoubletSeedFinder> create(const DerivedConfig& config);

  virtual ~DoubletSeedFinder() = default;

  /// Returns the configuration of the doublet seed finder.
  /// @return Reference to the configuration object
  virtual const DerivedConfig& config() const = 0;

  /// Creates compatible dublets by applying a series of cuts that can be
  /// tested with only two SPs.
  ///
  /// @param middleSp Space point candidate to be used as middle SP in a seed
  /// @param middleSpInfo Information about the middle space point
  /// @param candidateSps Subset of space points to be used as candidates for
  ///   middle SP in a seed
  /// @param compatibleDoublets Output container for compatible doublets
  virtual void createDoublets(
      const ConstSpacePointProxy& middleSp, const MiddleSpInfo& middleSpInfo,
      SpacePointContainer::ConstSubset candidateSps,
      DoubletsForMiddleSp& compatibleDoublets) const = 0;

  /// Creates compatible dublets by applying a series of cuts that can be
  /// tested with only two SPs.
  ///
  /// @param middleSp Space point candidate to be used as middle SP in a seed
  /// @param middleSpInfo Information about the middle space point
  /// @param candidateSps Range of space points to be used as candidates for
  ///   middle SP in a seed
  /// @param compatibleDoublets Output container for compatible doublets
  virtual void createDoublets(
      const ConstSpacePointProxy& middleSp, const MiddleSpInfo& middleSpInfo,
      SpacePointContainer::ConstRange candidateSps,
      DoubletsForMiddleSp& compatibleDoublets) const = 0;
};

}  // namespace Acts
