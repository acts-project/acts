// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/EventData/StripSpacePointCalibrationDetails.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Seeding/CartesianSpacePointGrid.hpp"
#include "Acts/Seeding/CosmicSeedCuts.hpp"
#include "Acts/Seeding/DoubletSeedFinder.hpp"
#include "Acts/Seeding/TripletSeedFinder.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

namespace Acts {

// finds triplet seeds from cosmic ray space points, applying the configured
// cuts
class CosmicTripletSeedFinder {
 public:
  explicit CosmicTripletSeedFinder(const CosmicSeedCuts& cuts = {});

  void createSeeds(
      const SpacePointContainer& spacePoints, SeedContainer& outputSeeds,
      const std::vector<CosmicSpacePointClass>* spClasses = nullptr,
      SpacePointContainer* calibratedSpacePoints = nullptr,
      const std::vector<GeometryIdentifier>* spSurfaces = nullptr) const;

 private:
  // position of a SP in a triplet: the one from the container,
  // or the recalculated one with the direction of the triplet if it is a strip
  // SP
  struct Point {
    double x = 0., y = 0., z = 0.;
    std::size_t index = 0;
  };

  // information about the event, derived from the space points and their
  // classes
  struct EventData {
    const SpacePointContainer* spacePoints = nullptr;
    // what the doublet and triplet finders read, with the indices of
    // spacePoints: x, y, z and the variances (scaled with covarianceScale if
    // useCovarianceScale). The "r" of its zr column is y (see prepareEvent)
    SpacePointContainer seeding;
    // the same with -y as "r", only without requireTopMiddleBottom: a
    // candidate below the middle SP is then also a top and one above it also
    // a bottom
    SpacePointContainer reflected;
    // positions in double precision, and the z bin (kNoZBin: not used)
    std::vector<Point> points;
    std::vector<std::uint32_t> zBin;
    std::vector<std::uint8_t> isStrip;
    std::vector<std::uint8_t> isCalibratable;
    std::vector<OuterStripSpacePointCalibrationDetailsDerived> stripDetails;
    const GeometryIdentifier* surfaces = nullptr;
    // the extent to which the y-coordinate of a strip SP changes upon
    // recalculation [mm]
    double stripShiftY = 0.0;
    bool anyCalibratable = false;

    bool strip(std::size_t index) const {
      return !isStrip.empty() && isStrip[index] != 0;
    }
    bool calibratable(std::size_t index) const {
      return anyCalibratable && isCalibratable[index] != 0;
    }
  };

  // doublet cuts of the cosmic seeding, given to Acts::DoubletSeedFinder as
  // its experimentCuts
  struct DoubletCuts {
    const CosmicTripletSeedFinder* finder = nullptr;
    const EventData* ev = nullptr;
    // the y window of the candidates is wider than [minDeltaY, maxDeltaY]
    // because of the recalculated strip SPs: cut it again on the others
    bool exactDeltaY = false;

    // doublets of the container positions
    bool operator()(const ConstSpacePointProxy& middle,
                    const ConstSpacePointProxy& other, float cotTheta,
                    bool isBottomCandidate) const;
    // doublets of a triplet with recalculated strip SPs, whose other cuts
    // were done on the triplet
    bool recalculated(const ConstSpacePointProxy& middle,
                      const ConstSpacePointProxy& other, float cotTheta,
                      bool isBottomCandidate) const;
  };

  // fill ev with the necessary information for each space point, depending on
  // the cuts that are enabled
  void prepareEvent(const SpacePointContainer& spacePoints,
                    const std::vector<CosmicSpacePointClass>* spClasses,
                    const std::vector<GeometryIdentifier>* spSurfaces,
                    EventData& ev) const;

  // index of the z bin of z (see CosmicSeedCuts::zBinEdges)
  std::size_t zBin(double z) const;

  // a doublet finder of the cosmic seeding: candidates ordered in y,
  // |dy| in [deltaYMin, deltaYMax], no collision region or cotTheta cut
  // (they are cut on dz/dr by the experimentCuts)
  static std::unique_ptr<DoubletSeedFinder> makeDoubletFinder(
      bool bottom, double deltaYMin, double deltaYMax,
      const DoubletSeedFinder::Config::ExperimentCuts& experimentCuts);

  // the Athena quality (dT weighted by the segment lengths) of a triplet that
  // passed the TripletSeedFinder, from its two doublets and its curvature
  float seedQuality(const DoubletsForMiddleSp::Proxy& bottom,
                    const DoubletsForMiddleSp::Proxy& top,
                    const ConstSpacePointProxy& middle, float curvature) const;

  bool calibrateStripSpacePoints(Point& top, Point& middle, Point& bottom,
                                 const EventData& ev) const;

  bool tripletDirections(const Point& top, const Point& middle,
                         const Point& bottom,
                         std::array<std::array<float, 3>, 3>& dirs) const;

  // a segment of a triplet: d1 = middle - top or d2 = bottom - middle, with
  // what the triplet cuts compute from it alone. A doublet is a segment of
  // every triplet it is part of, so this is computed once per doublet (the
  // transverse length and the polar angle only if a cut asks for them)
  struct Segment {
    double dx = 0., dy = 0., dz = 0.;
    double length = 0.;  // |d|

    // |d| in xy (pT-dependent and pixel cuts)
    double transverse() const {
      if (std::isnan(m_transverse)) {
        m_transverse = std::hypot(dx, dy);
      }
      return m_transverse;
    }
    // atan2(|d| in xy, dz) (pT-dependent cuts)
    double polar() const {
      if (std::isnan(m_polar)) {
        m_polar = std::atan2(transverse(), dz);
      }
      return m_polar;
    }

   private:
    mutable double m_transverse = std::numeric_limits<double>::quiet_NaN();
    mutable double m_polar = std::numeric_limits<double>::quiet_NaN();
  };
  // the segment d = to - from
  static Segment segment(const Point& from, const Point& to);

  // cuts of the triplet that the TripletSeedFinder does not have: order in y,
  // z range, kink (or cos(dphi) for pixels), |dy| per segment, cosThetaY of
  // the top-bottom chord and residual of the middle SP to that chord
  bool passGeometricCuts(const Point& top, const Point& middle,
                         const Point& bottom, const EventData& ev) const;
  // the same with the two segments already computed:
  // top = segment(top, middle), bottom = segment(middle, bottom)
  bool passGeometricCuts(const Point& top, const Point& middle,
                         const Point& bottom, const Segment& topSegment,
                         const Segment& bottomSegment,
                         const EventData& ev) const;

  struct PtDependentTerms {
    double polarKink = 0.;         // |theta1 - theta2| [rad]
    double kinkSlack = 0.;         // multiple scattering in angle [rad]
    double parallelResidual = 0.;  // residual out of the transverse plane [mm]
    double residualSlack = 0.;     // multiple scattering in residual [mm]
  };
  // residual = distance of the middle SP to the top-bottom chord
  bool ptDependentTerms(const Segment& topSegment, const Segment& bottomSegment,
                        double residual, PtDependentTerms& out) const;

  CosmicSeedCuts m_cuts;
  double m_dzdrMax = 0.0;
  // the grid of the space points of an event (built in createSeeds)
  // with CosmicSeedCuts::useCartesianGridBinning (the default), equidistant
  // bins in x, y and z [mm], the space point at (x, y, z)
  CartesianSpacePointGrid::Config m_gridConfig;

  // the triplet compatibility and curvature cuts (see makeTripletFinderConfig
  // in the cpp)
  std::shared_ptr<const TripletSeedFinder> m_tripletFinder;
  // the TripletSeedFinder has a cut to make (dT with
  // useLongitudinalCompatibility, or the pT cut) without them it would let
  // every triplet through and is not called
  bool m_useTripletFinder = false;
  // multiple scattering terms of the quality, the same as in the dT of Athena:
  // without pT cut, (1 + cot^2) * scatteringMinPt2; with pT cut,
  // (1 + cot^2) * curvature^2 * scatteringPerCurvature2
  double m_scatteringMinPt2 = 0.0;
  double m_scatteringPerCurvature2 = 0.0;
  bool m_ptCut = false;

  bool m_ptDependent = false;
  // constants of the pT-dependent kink and residual:
  // sqrt(multipleScatteringCof), pT per mm of radius, 1 / (minPtFactor minPt)
  // and 1 / maxPtForScattering
  double m_sqrtScatteringCof = 0.0;
  double m_ptPerRadius = 0.0;
  double m_maxInvPt = 0.0;
  double m_minInvPtScattering = 0.0;
  // the kink cut acos(cos(kink)) > maxKinkAngle without the acos: certain below
  // cos(maxKinkAngle) - kKinkCosBand (rejected) and above cos(maxKinkAngle) +
  // kKinkCosBand (kept)
  double m_cosKinkReject = -2.0;
  double m_cosKinkAccept = -2.0;
  // lower limit of the cosine of the angle between each segment and the y-axis
  // for the prefilter
  double m_segCosThetaYMin = 0.0;
};

}  // namespace Acts
