// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// cosmictripletseedfinder.cpp

#include "Acts/Seeding/CosmicTripletSeedFinder.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/SpacePointFormation/detail/StripSpacePointCalibrationImpl.hpp"
#include "Acts/Utilities/GridBinFinder.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Acts {

namespace {

// varianceR() and varianceZ():
double varianceR(const ConstSpacePointProxy& sp) {
  return static_cast<double>(sp.varianceR());
}
double varianceZ(const ConstSpacePointProxy& sp) {
  return static_cast<double>(sp.varianceZ());
}

// kPrefilterTol avoids numerical issues in the prefiltering cuts (it is never
// stricter than the exact cuts)
constexpr double kPrefilterTol = 1e-9;

// half width of the band around cos(maxKinkAngle) where the kink cut computes
// the acos (see m_cosKinkReject in the .hpp)
constexpr double kKinkCosBand = 1e-12;

// a term of the doublet or triplet finder that has to let everything through
constexpr float kNoLimit = 1e30f;

// a global variance (r or z) smaller than this fraction of the other one is
// rounding noise of the projection of the local covariance, not a resolution:
// for a pixel module perpendicular to z (endcap disk) the variance in z comes
// out as cos(pi/2)^2 * local variance, ~1e-37 mm^2
// kept as is, a doublet with dz = 0 (two space points on the same disk)
constexpr double kVarianceRoundingFraction = 1e-12;

// z bin of a space point that does not take part (coordinates not finite or z
// outside [zMin, zMax])
constexpr std::uint32_t kNoZBin = std::numeric_limits<std::uint32_t>::max();

// frame of the conformal transformation of every doublet:
// xNewFrame = dx cosPhiM + dy sinPhiM = dy, yNewFrame = dy cosPhiM - dx sinPhiM
// = -dx and it saves u = xNewFrame/|dxy|^2, v = yNewFrame/|dxy|^2.
// MiddleSpInfo's order: {uIP, uIP2, cosPhiM, sinPhiM}
constexpr MiddleSpInfo kFrameAlongY{0.f, 0.f, 0.f, 1.f};

TripletSeedFinder::DerivedConfig makeTripletFinderConfig(
    const CosmicSeedCuts& cuts, bool ptCut) {
  const float minPt =
      static_cast<float>(cuts.minPtFactor * cuts.minPt * UnitConstants::MeV);
  const float cof = static_cast<float>(cuts.multipleScatteringCof *
                                       UnitConstants::MeV * UnitConstants::MeV);

  TripletSeedFinder::Config config;
  // the strip SPs are recalculated by CosmicTripletSeedFinder itself, with
  // the direction of the triplet instead of the one of a helix from the beam
  // line
  config.useStripInfo = false;
  // every top doublet is tested, without the early break of the sorted mode
  config.sortedByCotTheta = false;
  config.impactMax = std::numeric_limits<float>::infinity();
  config.helixCutTolerance = 1.0f;
  if (cuts.minPt > 0.0) {
    config.minPt = minPt;
  }
  // bFieldInZ is pT per helix radius in ACTS units, with ptPerFieldRadius
  // the DerivedConfig needs a field != 0 to stay finite
  // without pT cut every term that depends on it is replaced below
  const double ptPerRadiusPerTesla =
      cuts.ptPerFieldRadius * UnitConstants::MeV / UnitConstants::mm;
  const float bFieldInZ = static_cast<float>(
      (cuts.magneticFieldBz != 0.0 ? cuts.magneticFieldBz : 1.0) *
      ptPerRadiusPerTesla);

  TripletSeedFinder::DerivedConfig derived(config, bFieldInZ);
  if (!cuts.useLongitudinalCompatibility) {
    derived.multipleScattering2 = kNoLimit;
    derived.sigmapT2perRadius = kNoLimit;
  } else {
    derived.multipleScattering2 =
        cuts.minPt > 0.0 ? cof / (minPt * minPt) : 0.0f;
    derived.sigmapT2perRadius =
        ptCut ? 4.0f * cof / (bFieldInZ * bFieldInZ) : kNoLimit;
  }
  if (!ptCut) {
    derived.minHelixDiameter2 = 0.0f;
  }
  return derived;
}

}  // namespace

// CosmicTripletSeedFinder constructor: checks the validity of the cuts and
// precomputes some derived quantities
CosmicTripletSeedFinder::CosmicTripletSeedFinder(const CosmicSeedCuts& cuts)
    : m_cuts(cuts) {
  // N
  if (!(m_cuts.zMin < m_cuts.zMax)) {
    throw std::invalid_argument("CosmicSeedCuts: zMin must be < zMax");
  }
  // kink should be >= 0
  if (m_cuts.maxKinkAngle < 0.0) {
    throw std::invalid_argument("CosmicSeedCuts: maxKinkAngle must be >= 0");
  }
  // |dy| window: min >= 0 and max > 0
  if (!(m_cuts.minDeltaY >= 0.0) || !(m_cuts.maxDeltaY > m_cuts.minDeltaY)) {
    throw std::invalid_argument(
        "CosmicSeedCuts: need 0 <= minDeltaY < maxDeltaY");
  }
  // STEP 1: validation of constants
  if (m_cuts.useStripCalibration) {
    // strip crossing tolerance must be > 0
    if (!(m_cuts.stripCalibrationTolerance > 0.0)) {
      throw std::invalid_argument(
          "CosmicSeedCuts: stripCalibrationTolerance must be > 0");
    }
    // at least one iteration of the strip calculation
    if (m_cuts.stripCalibrationMaxIterations < 1) {
      throw std::invalid_argument(
          "CosmicSeedCuts: stripCalibrationMaxIterations must be >= 1");
    }
    // convergence threshold in mm must be > 0
    if (!(m_cuts.stripCalibrationConvergence > 0.0)) {
      throw std::invalid_argument(
          "CosmicSeedCuts: stripCalibrationConvergence must be > 0");
    }
    // minimum cos with the normal in [0,1) (1 rejects all directions which are
    // not normal to the module)
    if (!(m_cuts.stripCalibrationMinCosToNormal >= 0.0 &&
          m_cuts.stripCalibrationMinCosToNormal < 1.0)) {
      throw std::invalid_argument(
          "CosmicSeedCuts: stripCalibrationMinCosToNormal must be in [0, 1)");
    }
  }
  // cosThetaY is |dy|/|d| always in [0,1]
  if (m_cuts.minCosThetaY < 0.0 || m_cuts.minCosThetaY > 1.0) {
    throw std::invalid_argument(
        "CosmicSeedCuts: minCosThetaY must be in [0, 1]");
  }
  // |eta| max > 0
  if (!(m_cuts.maxAbsEta > 0.0)) {
    throw std::invalid_argument("CosmicSeedCuts: maxAbsEta must be > 0");
  }
  // pT min >= 0 (if pT = 0: without cut in pT)
  if (m_cuts.minPt < 0.0) {
    throw std::invalid_argument("CosmicSeedCuts: minPt must be >= 0");
  }
  // variance COF >= 0 (MeV^2)
  if (m_cuts.multipleScatteringCof < 0.0) {
    throw std::invalid_argument(
        "CosmicSeedCuts: multipleScatteringCof must be >= 0");
  }
  // minPtFactor and ptPerFieldRadius divide/multiply pT; they should be finite
  // and positive
  if (!(m_cuts.minPtFactor > 0.0) || !std::isfinite(m_cuts.minPtFactor)) {
    throw std::invalid_argument(
        "CosmicSeedCuts: minPtFactor must be finite and > 0");
  }
  if (!(m_cuts.ptPerFieldRadius > 0.0) ||
      !std::isfinite(m_cuts.ptPerFieldRadius)) {
    throw std::invalid_argument(
        "CosmicSeedCuts: ptPerFieldRadius must be finite and > 0");
  }

  // acceptance window in dz/dr = dz / sqrt(dx^2 + dy^2)
  const double etaMax = std::abs(m_cuts.maxAbsEta);
  // |dz/dxy| max: theta = 2 atan(exp(-eta)) and dz/dxy = cot(theta) = sinh(eta)
  m_dzdrMax = 1.0 / std::tan(2.0 * std::atan(std::exp(-etaMax)));

  // STEP 2: z bins and the neighbours of each one
  // where zBinNeighbors gives (first, second) bins
  // range where their tops and bottoms are looked for
  const std::vector<double>& edges = m_cuts.zBinEdges;
  // edges should be finite and strictly increasing
  for (std::size_t i = 0; i < edges.size(); ++i) {
    if (!std::isfinite(edges[i]) || (i > 0 && !(edges[i] > edges[i - 1]))) {
      throw std::invalid_argument(
          "CosmicSeedCuts: zBinEdges must be finite and strictly increasing");
    }
  }
  // N interior edges -> N + 1 bins
  const std::size_t nZBins = edges.size() + 1;
  if (!m_cuts.zBinNeighbors.empty() && m_cuts.zBinNeighbors.size() != nZBins) {
    throw std::invalid_argument(
        "CosmicSeedCuts: zBinNeighbors needs one pair per z bin "
        "(zBinEdges.size() + 1), or none");
  }
  // the grid of the space points (see createSeeds): its axes are equidistant
  // and the z bins are not, so its z axis is the index of the z bin, one bin
  // of width 1 per z bin; x and y are a single bin each
  m_gridConfig.nXbins = 1;
  m_gridConfig.xMin = -kNoLimit;
  m_gridConfig.xMax = kNoLimit;
  m_gridConfig.nYbins = 1;
  m_gridConfig.yMin = -kNoLimit;
  m_gridConfig.yMax = kNoLimit;
  m_gridConfig.nZbins = static_cast<std::uint32_t>(nZBins);
  m_gridConfig.zMin = 0.0f;
  m_gridConfig.zMax = static_cast<float>(nZBins);
  // bins list of neighbors per bin
  std::vector<std::pair<int, int>> zNeighbors;
  zNeighbors.reserve(nZBins);
  for (std::size_t bin = 0; bin < nZBins; ++bin) {
    // (first, second) displacement related to the bin
    // (-1,+1) = the previous bin, the bin itself, and the next one
    const auto [first, second] = m_cuts.zBinNeighbors.empty()
                                     ? std::pair<int, int>{-1, 1}
                                     : m_cuts.zBinNeighbors[bin];
    // the range cannot be empty
    if (first > second) {
      throw std::invalid_argument(
          "CosmicSeedCuts: each zBinNeighbors pair needs first <= second");
    }
    zNeighbors.emplace_back(first, second);
  }
  // the range is restricted to the existing bins by the grid: its underflow
  // and overflow bins stay empty
  m_gridConfig.bottomBinFinder = GridBinFinder<3>(0, 0, zNeighbors);
  m_gridConfig.topBinFinder = m_gridConfig.bottomBinFinder;

  // STEP 2b: equidistant bins in x, y and z [mm]
  if (m_cuts.useCartesianGridBinning) {
    for (std::size_t axis = 0; axis < 3; ++axis) {
      if (m_cuts.gridBins[axis] == 0 || m_cuts.gridNeighbors[axis] < 0) {
        throw std::invalid_argument(
            "CosmicSeedCuts: gridBins must be >= 1 and gridNeighbors >= 0");
      }
    }
    auto validRange = [](double lo, double hi) {
      return std::isfinite(lo) && std::isfinite(hi) && lo < hi;
    };
    if (!validRange(m_cuts.gridXMin, m_cuts.gridXMax) ||
        !validRange(m_cuts.gridYMin, m_cuts.gridYMax) ||
        !validRange(m_cuts.zMin, m_cuts.zMax)) {
      throw std::invalid_argument(
          "CosmicSeedCuts: gridXMin/gridXMax, gridYMin/gridYMax and zMin/zMax "
          "must be finite with min < max");
    }
    // the axes are [min, max): the next float above max keeps max inside
    auto upper = [](double max) {
      return std::nextafter(static_cast<float>(max),
                            std::numeric_limits<float>::infinity());
    };
    m_gridConfig.nXbins = m_cuts.gridBins[0];
    m_gridConfig.xMin = static_cast<float>(m_cuts.gridXMin);
    m_gridConfig.xMax = upper(m_cuts.gridXMax);
    m_gridConfig.nYbins = m_cuts.gridBins[1];
    m_gridConfig.yMin = static_cast<float>(m_cuts.gridYMin);
    m_gridConfig.yMax = upper(m_cuts.gridYMax);
    m_gridConfig.nZbins = m_cuts.gridBins[2];
    m_gridConfig.zMin = static_cast<float>(m_cuts.zMin);
    m_gridConfig.zMax = upper(m_cuts.zMax);
    m_gridConfig.bottomBinFinder =
        GridBinFinder<3>(m_cuts.gridNeighbors[0], m_cuts.gridNeighbors[1],
                         m_cuts.gridNeighbors[2]);
    m_gridConfig.topBinFinder = m_gridConfig.bottomBinFinder;
  }

  // check covarianceScale entries if used
  if (m_cuts.useCovarianceScale) {
    for (double s : m_cuts.covarianceScale) {
      if (!(s > 0.0) || !std::isfinite(s)) {
        throw std::invalid_argument(
            "CosmicSeedCuts: covarianceScale entries must be finite and > 0");
      }
    }
  }

  // check pixelMinCosDeltaPhi if used, it should be in the [-1, 1] range
  if (m_cuts.useSubdetectorAngularCuts &&
      !(m_cuts.pixelMinCosDeltaPhi >= -1.0 &&
        m_cuts.pixelMinCosDeltaPhi <= 1.0)) {
    throw std::invalid_argument(
        "CosmicSeedCuts: pixelMinCosDeltaPhi must be in [-1, 1]");
  }

  // check pt-dependent kink residual cuts if used
  // field should be != 0 and minPt > 0
  if (m_cuts.usePtDependentKinkResidual) {
    if (!m_cuts.usePtCut || !(m_cuts.minPt > 0.0) ||
        m_cuts.magneticFieldBz == 0.0) {
      throw std::invalid_argument(
          "CosmicSeedCuts: usePtDependentKinkResidual requires usePtCut, "
          "minPt > 0 and magneticFieldBz != 0");
    }
    // maxPtForScattering > 0
    if (!(m_cuts.maxPtForScattering > 0.0)) {
      throw std::invalid_argument(
          "CosmicSeedCuts: maxPtForScattering must be > 0");
    }
  }

  // STEP 3: triplet compatibility and curvature: Acts::TripletSeedFinder
  m_ptCut =
      m_cuts.usePtCut && m_cuts.minPt > 0.0 && m_cuts.magneticFieldBz != 0.0;
  // without dT and pT cut every term of the finder is kNoLimit or 0 and it
  // would accept every triplet: it is not called, and the quality does not
  // need its curvature (without pT cut the scattering uses minPt)
  m_useTripletFinder = m_cuts.useLongitudinalCompatibility || m_ptCut;
  // finder is created
  const TripletSeedFinder::DerivedConfig finderConfig =
      makeTripletFinderConfig(m_cuts, m_ptCut);
  // create() returns a shared_ptr to the implementation
  m_tripletFinder = TripletSeedFinder::create(finderConfig);
  // COF in ACTS units
  const double cof =
      m_cuts.multipleScatteringCof * UnitConstants::MeV * UnitConstants::MeV;
  // scattering with min pT: COF / (0.9 minPt)^2
  if (m_cuts.minPt > 0.0) {
    const double minPt = m_cuts.minPtFactor * m_cuts.minPt * UnitConstants::MeV;
    m_scatteringMinPt2 = cof / (minPt * minPt);
  }
  // scattering per curvature^2 : 4 COF / bFieldInZ^2
  if (m_ptCut) {
    const double bFieldInZ = finderConfig.bFieldInZ;
    m_scatteringPerCurvature2 = 4.0 * cof / (bFieldInZ * bFieldInZ);
  }

  // STEP 4: derived quantities for the prefilters
  m_ptDependent = m_cuts.usePtDependentKinkResidual;
  m_sqrtScatteringCof = std::sqrt(m_cuts.multipleScatteringCof);
  if (m_ptDependent) {
    // pT [MeV] = 0.3 * |Bz|
    m_ptPerRadius = m_cuts.ptPerFieldRadius * std::abs(m_cuts.magneticFieldBz);
    // 1/pT max permitted = 1/(0.9 minPt)
    m_maxInvPt = 1.0 / (m_cuts.minPtFactor * m_cuts.minPt);
    // 1/pT min of the scattering: pT_eff = min(pT, maxPtForScattering)
    m_minInvPtScattering = 1.0 / m_cuts.maxPtForScattering;
  }
  if (!m_ptDependent) {
    // k = maxKinkAngle [rad]
    const double k = m_cuts.maxKinkAngle;
    // if cos(kink) < cos(k) - 1e-12 it is rejected
    m_cosKinkReject = std::cos(k) - kKinkCosBand;
    // if cos(kink) > cos(k) + 1e-12 it is accepted
    m_cosKinkAccept = std::cos(k) + kKinkCosBand;
    // sum of the two segments, and each segment forms at most "k" with it
    // the angle of the chord with the y-axis is at most acos(minCosThetaY);
    // by the triangle inequality on the sphere,
    // each segment forms at most acos(minCosThetaY) + k with the y-axis
    const double maxAngleY = std::acos(m_cuts.minCosThetaY) + k;
    // if the angle comes to 90°, the bound is zero and the prefilter is
    // deactivated
    m_segCosThetaYMin =
        maxAngleY < 0.5 * std::numbers::pi ? std::cos(maxAngleY) : 0.0;
  }
}

// STEP 5: prepareEvent(): what is calculated once per SP in each event
// - ev.seeding: container read by the two ACTS finders, with the same indices
// as SPs: xy, (z,"r") and the variances
// - variances: with covarianceScale it's applied to y
// - ev.points: positions in double; ev.zBin bin in z of each SP
// - ev.reflected: copy with -y as "r", only if requireTopMiddleBottom is off,
// so that an SP from the bottom can also be a "top", and one from the top can
// also be a "bottom", only for doublets (using DoubletSeedFinder)
void CosmicTripletSeedFinder::prepareEvent(
    const SpacePointContainer& spacePoints,
    const std::vector<CosmicSpacePointClass>* spClasses,
    const std::vector<GeometryIdentifier>* spSurfaces, EventData& ev) const {
  // to prepare the event, we need to fill ev with the necessary information for
  // each space point, depending on the cuts that are enabled, this includes the
  // container for the doublet and triplet finders (with the covariance scaling)
  // and type information, as well as any surface identifiers if
  // rejectSameSurface is enabled pointer to the original container
  ev.spacePoints = &spacePoints;
  // number of SPs in the event
  const std::size_t n = spacePoints.size();

  // the veto compares ids: a GeometryIdentifier is needed per SP
  if (m_cuts.rejectSameSurface) {
    if (spSurfaces == nullptr || spSurfaces->size() != n) {
      throw std::invalid_argument(
          "CosmicTripletSeedFinder: rejectSameSurface requires one surface "
          "GeometryIdentifier per space point");
    }
    // each id should be from a sensitive module: sensitive() != 0
    for (const GeometryIdentifier& surface : *spSurfaces) {
      if (surface.sensitive() == 0) {
        throw std::invalid_argument(
            "CosmicTripletSeedFinder: rejectSameSurface requires sensitive "
            "surface ids, got one with sensitive() == 0");
      }
    }
    // pointer to the vector's data
    ev.surfaces = spSurfaces->data();
  }

  // covariance scaling and type information (pixel/strip) for each space point
  const bool needType =
      m_cuts.requireSameSubdetector || m_cuts.useSubdetectorAngularCuts;
  // class (barrel/endcap, pixel/strip) is needed to scale the covariance
  const bool needClasses = m_cuts.useCovarianceScale || needType;
  // at least one class per SP
  if (needClasses && (spClasses == nullptr || spClasses->size() != n)) {
    throw std::invalid_argument(
        "CosmicTripletSeedFinder: useCovarianceScale / "
        "requireSameSubdetector / useSubdetectorAngularCuts require one "
        "CosmicSpacePointClass per space point");
  }
  // what the doublet and triplet finders read, with the indices of SPs:  x, y,
  // z and the variances the space points of a cosmic are ordered in y, not in r
  // (r decreases and then grows again along the track),
  // the "r" of the zr column is y: the deltaR of Acts::DoubletSeedFinder
  // (r_middle - r_bottom, r_top - r_middle) is then the y separation
  ev.seeding = SpacePointContainer(
      SpacePointColumns::PackedXY | SpacePointColumns::PackedZR |
      SpacePointColumns::VarianceZ | SpacePointColumns::VarianceR);
  // memory for n SP
  ev.seeding.reserve(static_cast<std::uint32_t>(n));
  // positions in double, one per SP, the same index
  ev.points.resize(n);
  // all the bins start as 'not participating'; they are assigned to a bin below
  // if they pass
  ev.zBin.assign(n, kNoZBin);
  // pixel/strip vector, only if it is used
  if (needType) {
    ev.isStrip.resize(n);
  }
  // fill the buffers with the necessary information for each space point
  for (ConstSpacePointProxy sp : spacePoints) {
    const std::size_t i = sp.index();
    if (i >= n) {
      throw std::invalid_argument(
          "CosmicTripletSeedFinder: space point index out of range");
    }
    // class by default 0 (BarrelPixel) when it is not assigned
    // if this is the case, covariance is not scaled and it's not used
    std::size_t cls = 0;
    if (needClasses) {
      //
      cls = static_cast<std::size_t>((*spClasses)[i]);
      if (cls >= m_cuts.covarianceScale.size()) {
        throw std::invalid_argument(
            "CosmicTripletSeedFinder: invalid CosmicSpacePointClass");
      }
    }
    // BarrelPixel or BarrelStrip class
    const bool barrel =
        cls == static_cast<std::size_t>(CosmicSpacePointClass::BarrelPixel) ||
        cls == static_cast<std::size_t>(CosmicSpacePointClass::BarrelStrip);
    // BarrelStrip or EndcapStrip class
    const bool strip =
        cls == static_cast<std::size_t>(CosmicSpacePointClass::BarrelStrip) ||
        cls == static_cast<std::size_t>(CosmicSpacePointClass::EndcapStrip);
    // errorsc: only the dominant coordinate is scaled (z in barrel, r in
    // endcap); the other remains with the variance of the container.
    double covR = varianceR(sp);
    double covZ = varianceZ(sp);
    // error: covariance is multiplied by the worst measured coordinate
    if (m_cuts.useCovarianceScale) {
      const double scale = m_cuts.covarianceScale[cls];
      (barrel ? covZ : covR) *= scale;
    }
    // rounding noise of the projection -> 0 (see kVarianceRoundingFraction)
    if (covZ < kVarianceRoundingFraction * covR) {
      covZ = 0.0;
    }
    if (covR < kVarianceRoundingFraction * covZ) {
      covR = 0.0;
    }
    // SP coordinates in float
    const float x = sp.x();
    const float y = sp.y();
    const float z = sp.z();
    // new SP in seeding; their index is the same as i
    MutableSpacePointProxy out = ev.seeding.createSpacePoint();
    out.xy() = {x, y};
    // zr = (z,y), the "r" that is seen by DoubletSeedFinder is the height "y"
    out.zr() = {z, y};
    // scaled variances in float
    out.varianceZ() = static_cast<float>(covZ);
    out.varianceR() = static_cast<float>(covR);
    // SP copy with the index, for geometric cuts
    ev.points[i] = Point{x, y, z, i};
    // the space points with finite coordinates in [zMin, zMax]
    if (std::isfinite(x) && std::isfinite(y) && std::isfinite(z) &&
        z >= m_cuts.zMin && z <= m_cuts.zMax) {
      ev.zBin[i] = static_cast<std::uint32_t>(zBin(z));
    }
    // 1 = strip, 0 = pixel
    if (needType) {
      ev.isStrip[i] = strip ? 1 : 0;
    }
  }
  // copy with r = -y
  if (!m_cuts.requireTopMiddleBottom) {
    ev.reflected = ev.seeding;
    for (MutableSpacePointProxy sp : ev.reflected) {
      sp.zr()[1] = -sp.zr()[1];
    }
  }

  // two strip SPs that are recalculated based on the direction of each seed
  // (useStripCalibration) stripShiftY limits how much the y coordinate of an SP
  // can change when it is recalculated (shift from the inner strip to the outer
  // strip + movement along the outer strip within the tolerance)
  // -> this expands the window in the y direction.

  if (m_cuts.useStripCalibration &&
      spacePoints.hasColumns(SpacePointColumns::StripCalibrationDetails)) {
    // none that can be recalculated by default
    ev.isCalibratable.assign(n, 0);
    // details derived from strips, one per space point
    ev.stripDetails.resize(n);
    for (ConstSpacePointProxy sp : spacePoints) {
      // saves the info about the two centers of the strips which make the space
      // point
      const OuterStripSpacePointCalibrationDetails& d =
          sp.outerStripCalibrationDetails();
      // squared norm of a vector of 3 floats
      auto norm2 = [](const std::array<float, 3>& v) {
        return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
      };
      // skip SPs that are not recalibratable (e.g. pixel SPs or strip SPs with
      // zero vectors)
      if (!(norm2(d.outerHalfVector) > 0.f) ||
          !(norm2(d.innerHalfVector) > 0.f)) {
        continue;
      }
      // mark the SP as recalibratable and store its calibration details
      ev.isCalibratable[sp.index()] = 1;
      ev.anyCalibratable = true;
      // cross products are precomputed (ihv x ohv, iosv x ohv, iosv x ihv)
      ev.stripDetails[sp.index()] =
          detail::deriveOuterStripSpacePointCalibrationDetails(d);
      // Y shift for the recalibration: the maximum shift in y when
      // recalibrating a strip SP this is the sum of the absolute values of the
      // inner to outer separation in y, the inner half vector in y, and the
      // outer half vector in y scaled by the strip calibration tolerance.
      const double shift =
          std::abs(static_cast<double>(d.innerToOuterSeparation[1])) +
          std::abs(static_cast<double>(d.innerHalfVector[1])) +
          m_cuts.stripCalibrationTolerance *
              std::abs(static_cast<double>(d.outerHalfVector[1]));
      // the worst case scenario of the event is saved, thus widening the window
      // in y
      ev.stripShiftY = std::max(ev.stripShiftY, shift);
    }
  }
}

// STEP 2: zBin(): index of the z bin, binary search in zBinEdges
std::size_t CosmicTripletSeedFinder::zBin(double z) const {
  // lower_bound: first edge >= z ; upper_bound: first edge > z
  const std::vector<double>& edges = m_cuts.zBinEdges;
  const auto it = z > 0.0 ? std::upper_bound(edges.begin(), edges.end(), z)
                          : std::lower_bound(edges.begin(), edges.end(), z);
  // returns the index: how many edges there are before the first edge above z
  // it is useful to group per bin and decide if a SP enters the neighbor pool
  // of a bin
  return static_cast<std::size_t>(it - edges.begin());
}

// STEP 6: the doublet finder of the cosmic seeding
std::unique_ptr<DoubletSeedFinder> CosmicTripletSeedFinder::makeDoubletFinder(
    bool bottom, double deltaYMin, double deltaYMax,
    const DoubletSeedFinder::Config::ExperimentCuts& experimentCuts) {
  // config by default
  DoubletSeedFinder::Config config;
  // in the reflected container (r = -y)
  config.spacePointsSortedByRadius = false;
  // Backward() looks for candidates with lower 'r'
  // Forward() looks for candidates with higher 'r'
  config.candidateDirection =
      bottom ? Direction::Backward() : Direction::Forward();
  // deltaR window = window in |dy|
  config.deltaRMin = static_cast<float>(deltaYMin);
  config.deltaRMax = static_cast<float>(deltaYMax);
  config.interactionPointCut = false;
  config.collisionRegionMin = -kNoLimit;
  config.collisionRegionMax = kNoLimit;
  //
  config.cotThetaMax = kNoLimit;
  config.experimentCuts = experimentCuts;
  // the field only enters the interaction point cut, which is off
  return DoubletSeedFinder::create(
      DoubletSeedFinder::DerivedConfig(config, 1.0f));
}

// STEP 6: doublet cuts on the positions of the container, same subdetector, not
// on the same surface, the pseudorapidity window and the segment cosThetaY
// prefilter
bool CosmicTripletSeedFinder::DoubletCuts::operator()(
    const ConstSpacePointProxy& middle, const ConstSpacePointProxy& other,
    float cotTheta, bool /*isBottomCandidate*/) const {
  const CosmicSeedCuts& cuts = finder->m_cuts;
  const std::size_t indexM = middle.index();
  const std::size_t index = other.index();
  if (cuts.requireSameSubdetector &&
      ev->isStrip[index] != ev->isStrip[indexM]) {
    return false;
  }
  // surface veto: if the two SPs are on the same surface, they cannot be part
  // of the same triplet
  if (ev->surfaces != nullptr && ev->surfaces[index] == ev->surfaces[indexM]) {
    return false;
  }
  if (ev->calibratable(index) || ev->calibratable(indexM)) {
    return true;
  }
  const Point& pM = ev->points[indexM];
  const Point& p = ev->points[index];
  const double absDy = std::abs(p.y - pM.y);
  if (exactDeltaY && (absDy < cuts.minDeltaY || absDy > cuts.maxDeltaY)) {
    return false;
  }
  // pseudorapidity window for the segments, symmetric in Tz = dz / rT
  if (std::abs(cotTheta) > finder->m_dzdrMax) {
    return false;
  }
  // bound per segment derived from minCosThetaY and maxKinkAngle (see the
  // .hpp): the top-bottom vector is the sum of the two segments and each
  // segment forms at most `maxKinkAngle` with it with middle pixel and
  // useSubdetectorAngularCuts there is no 3D kink cut, so there is no bound
  const bool pixelMid = cuts.useSubdetectorAngularCuts && !ev->strip(indexM);
  if (!pixelMid && finder->m_segCosThetaYMin > 0.0) {
    const double dx = p.x - pM.x;
    const double dz = p.z - pM.z;
    const double length = std::sqrt(dx * dx + absDy * absDy + dz * dz);
    if (absDy < finder->m_segCosThetaYMin * length * (1.0 - kPrefilterTol)) {
      return false;
    }
  }
  return true;
}

// STEP 11: the doublets of a triplet with recalculated strip SPs: the triplet
// already passed the geometric cuts, only the pseudorapidity window is left
bool CosmicTripletSeedFinder::DoubletCuts::recalculated(
    const ConstSpacePointProxy& /*middle*/,
    const ConstSpacePointProxy& /*other*/, float cotTheta,
    bool /*isBottomCandidate*/) const {
  return !(std::abs(cotTheta) > finder->m_dzdrMax);
}

// STEP 12: quality = dT * invDistTop * invDistBottom / (invDistTop +
// invDistBottom), as in Athena, with dT = (cotB - cotT)^2 - error2 -
// scattering: the margin of the cut of the TripletSeedFinder, negative for a
// compatible triplet the errors are those of the finder; the scattering is
// Athena's, also when useLongitudinalCompatibility is off and the finder does
// not cut on dT if useSeedQuality is false, the quality is left as the fixed
// value m_cuts.seedQuality
float CosmicTripletSeedFinder::seedQuality(
    const DoubletsForMiddleSp::Proxy& bottom,
    const DoubletsForMiddleSp::Proxy& top, const ConstSpacePointProxy& middle,
    float curvature) const {
  if (!m_cuts.useSeedQuality) {
    return m_cuts.seedQuality;
  }
  const double cotB = bottom.cotTheta();
  const double cotT = top.cotTheta();
  const double iDeltaRB = bottom.iDeltaR();
  const double iDeltaRT = top.iDeltaR();
  const double error2 =
      top.er() + bottom.er() +
      2.0 * (cotB * cotT * varianceR(middle) + varianceZ(middle)) * iDeltaRB *
          iDeltaRT;
  const double iSinTheta2 = 1.0 + cotB * cotB;
  const double scattering = m_ptCut
                                ? iSinTheta2 * static_cast<double>(curvature) *
                                      curvature * m_scatteringPerCurvature2
                                : iSinTheta2 * m_scatteringMinPt2;
  const double deltaCotTheta = cotB - cotT;
  const double dT = deltaCotTheta * deltaCotTheta - error2 - scattering;
  const double denom = iDeltaRB + iDeltaRT;
  return denom > 0.0 ? static_cast<float>(dT * iDeltaRB * iDeltaRT / denom)
                     : static_cast<float>(dT);
}

// createSeeds() generates triplet seeds from the given space points, applying
// the configured cuts and storing the results in outputSeeds. It also
// optionally recalibrates strip space points and stores them in
// calibratedSpacePoints.
void CosmicTripletSeedFinder::createSeeds(
    const SpacePointContainer& spacePoints, SeedContainer& outputSeeds,
    const std::vector<CosmicSpacePointClass>* spClasses,
    SpacePointContainer* calibratedSpacePoints,
    const std::vector<GeometryIdentifier>* spSurfaces) const {
  if (!outputSeeds.hasSpacePointContainer()) {
    outputSeeds.assignSpacePointContainer(spacePoints);
  }

  // a cosmic gives many more seeds than space points (every triplet of its
  // space points in the window): ~n^2 for n space points, capped by 2^16
  const std::size_t expectedSeeds =
      std::min<std::size_t>(spacePoints.size() * spacePoints.size(), 1u << 16);
  outputSeeds.reserve(
      static_cast<SeedContainer::Index>(outputSeeds.size() + expectedSeeds),
      3.0f);

  // STEP 5: event data
  EventData ev;
  prepareEvent(spacePoints, spClasses, spSurfaces, ev);
  const SpacePointContainer& seeding = ev.seeding;
  const SpacePointContainer& reflected = ev.reflected;
  // true = top > middle > bottom order
  const bool ordered = m_cuts.requireTopMiddleBottom;

  // Best quality is used to determine if a recalibrated strip space point
  // should be updated in the calibratedSpacePoints container lower quality =
  // better quality
  std::vector<float> bestQuality;
  if (calibratedSpacePoints != nullptr && !ev.isCalibratable.empty()) {
    if (calibratedSpacePoints->size() != spacePoints.size()) {
      throw std::invalid_argument(
          "CosmicTripletSeedFinder: calibratedSpacePoints must be a copy of "
          "spacePoints");
    }
    // give the best quality vector the same size as the number of space points,
    // initialized to infinity
    bestQuality.assign(spacePoints.size(),
                       std::numeric_limits<float>::infinity());
  }

  // STEP 6: window in y for the segments, |dy| in [minDeltaY, maxDeltaY]
  // if the segment has a recalibratable strip SP, each end can move up to
  // stripShiftY in y: the window is enlarged by that and the exact cut is done
  // by the doublet cuts (the other doublets) and by passGeometricCuts (the
  // recalculated triplets)
  const double slack = ev.anyCalibratable ? 2.0 * ev.stripShiftY + 1e-6 : 0.0;
  const double deltaYMin = m_cuts.minDeltaY - slack;
  const double deltaYMax = m_cuts.maxDeltaY + slack;

  // the doublet finders of this event, with its doublet cuts
  const DoubletCuts doubletCuts{this, &ev, ev.anyCalibratable};
  DoubletSeedFinder::Config::ExperimentCuts eventCuts;
  // // the operator() method is bound to the doubletCuts object
  eventCuts.connect<&DoubletCuts::operator()>(&doubletCuts);
  // finder for event bottoms (backward) and tops (forward)
  const std::unique_ptr<DoubletSeedFinder> bottomFinder =
      makeDoubletFinder(true, deltaYMin, deltaYMax, eventCuts);
  const std::unique_ptr<DoubletSeedFinder> topFinder =
      makeDoubletFinder(false, deltaYMin, deltaYMax, eventCuts);
  // and the ones of the triplets with recalculated strip SPs, which read the
  // recalculated positions of a triplet from `scratch`:
  // bottom, middle, middle, top (see tryRecalculated)
  std::unique_ptr<DoubletSeedFinder> recalculatedBottomFinder;
  std::unique_ptr<DoubletSeedFinder> recalculatedTopFinder;
  // 4 SP container with the same columns as 'seeding'
  SpacePointContainer scratch(
      SpacePointColumns::PackedXY | SpacePointColumns::PackedZR |
      SpacePointColumns::VarianceZ | SpacePointColumns::VarianceR);
  if (ev.anyCalibratable) {
    // saves with the eta cut |cotTheta| ≤ sinh(maxAbsEta)
    DoubletSeedFinder::Config::ExperimentCuts recalculatedCuts;
    recalculatedCuts.connect<&DoubletCuts::recalculated>(&doubletCuts);
    // it uses [minDeltaY, maxDeltaY]
    recalculatedBottomFinder = makeDoubletFinder(
        true, m_cuts.minDeltaY, m_cuts.maxDeltaY, recalculatedCuts);
    recalculatedTopFinder = makeDoubletFinder(
        false, m_cuts.minDeltaY, m_cuts.maxDeltaY, recalculatedCuts);
    // the four empty slots are created once, then tryRecalculated overwrites
    // them
    scratch.reserve(4);
    for (int i = 0; i < 4; ++i) {
      // to ignore the proxy
      static_cast<void>(scratch.createSpacePoint());
    }
  }
  const SpacePointContainer& scratchView = scratch;

  // emit(): saves the seed (order: bottom, middle, top; vertexZ = 0) and, if it
  // has recalculated strips, copies their new positions to
  // calibratedSpacePoints—provided that this seed is of better quality than
  // previous ones that used that same SP
  auto emit = [&](const Point& pB, const Point& pM, const Point& pT,
                  float quality, bool recalculated) {
    // bottom, middle, top indices
    const std::array<SpacePointIndex, 3> seedSpacePoints{
        static_cast<SpacePointIndex>(pB.index),
        static_cast<SpacePointIndex>(pM.index),
        static_cast<SpacePointIndex>(pT.index)};
    // new seed in output
    MutableSeedProxy seed = outputSeeds.createSeed();
    seed.assignSpacePointIndices(seedSpacePoints);
    // quality: lower is better
    seed.quality() = quality;
    // for cosmic rays there is no primary vertex: as in Athena, zVertex = 0
    seed.vertexZ() = 0.0f;

    // if any of the space points are recalibrated strip space points, update
    // their positions in the calibratedSpacePoints container if the new quality
    // is better than the previous best quality
    if (!recalculated || bestQuality.empty()) {
      return;
    }
    // for each of the 3 SPs
    for (const Point* p : {&pB, &pM, &pT}) {
      // only recalculable items pass, and only if this seed is better than the
      // best previous one used by that SP.
      // !(a < b) also skips NaNs.
      if (!ev.calibratable(p->index) || !(quality < bestQuality[p->index])) {
        continue;
      }
      // the best quality and the SP's position in the output are updated
      bestQuality[p->index] = quality;
      MutableSpacePointProxy out =
          calibratedSpacePoints->at(static_cast<SpacePointIndex>(p->index));
      out.x() = static_cast<float>(p->x);
      out.y() = static_cast<float>(p->y);
      out.z() = static_cast<float>(p->z);
      // if the container has columns r and phi, they are recalculated so that
      // they remain consistent with the new x and y
      if (calibratedSpacePoints->hasColumns(SpacePointColumns::R)) {
        out.r() = static_cast<float>(std::hypot(p->x, p->y));
      }
      if (calibratedSpacePoints->hasColumns(SpacePointColumns::Phi)) {
        out.phi() = static_cast<float>(std::atan2(p->y, p->x));
      }
    }
  };

  // caches, reused for every middle space point
  std::vector<SpacePointIndex> pool;
  DoubletsForMiddleSp bottomDoublets;
  DoubletsForMiddleSp topDoublets;
  TripletTopCandidates candidates;
  // the segment of each doublet (see Segment in the .hpp)
  std::vector<Segment> bottomSegments;
  std::vector<Segment> topSegments;
  // top doublet of each space point, for the candidates of the
  // TripletSeedFinder
  std::vector<DoubletsForMiddleSp::Index> topOfSpacePoint(
      m_useTripletFinder ? spacePoints.size() : 0);
  // the same for one triplet with recalculated strip SPs
  DoubletsForMiddleSp recalculatedBottom;
  DoubletsForMiddleSp recalculatedTop;
  TripletTopCandidates recalculatedCandidates;

  // STEP 11: a triplet with a recalculated strip SP: the strip SPs are
  // recalculated with the direction of the triplet, and the cuts (geometric,
  // then the doublet and triplet finders) are done with the new positions
  auto tryRecalculated = [&](Point pB, Point pM, Point pT) {
    // first the iterative recalculous and then the geometric cuts
    if (!calibrateStripSpacePoints(pT, pM, pB, ev) ||
        !passGeometricCuts(pT, pM, pB, ev)) {
      return;
    }
    // "r" (= y, see prepareEvent) of each pair with the sign that makes it grow
    // from the bottom to the middle and from the middle to the top: always +
    // with requireTopMiddleBottom, which passGeometricCuts checked
    const double signB = pM.y >= pB.y ? 1.0 : -1.0;
    const double signT = pT.y >= pM.y ? 1.0 : -1.0;
    // write a recalculated point in slot, with r = sign * y
    auto place = [&](SpacePointIndex slot, const Point& p, double sign) {
      MutableSpacePointProxy sp = scratch.at(slot);
      const ConstSpacePointProxy original =
          seeding[static_cast<SpacePointIndex>(p.index)];
      sp.xy() = {static_cast<float>(p.x), static_cast<float>(p.y)};
      sp.zr() = {static_cast<float>(p.z), static_cast<float>(sign * p.y)};
      sp.varianceZ() = original.varianceZ();
      sp.varianceR() = original.varianceR();
    };
    // indexes: 0 = bottom, 1 = middle, 2 = middle, 3 = top
    place(0, pB, signB);
    place(1, pM, signB);
    place(2, pM, signT);
    place(3, pT, signT);
    // cleaning the doublets
    recalculatedBottom.clear();
    recalculatedTop.clear();
    // bottom doublet: middle = index 1, candidate range [0,1)
    recalculatedBottomFinder->createDoublets(scratchView[1], kFrameAlongY,
                                             scratchView.range({0, 1}),
                                             recalculatedBottom);
    // top doublet: middle = slot 2, candidate range [3, 4)
    recalculatedTopFinder->createDoublets(scratchView[2], kFrameAlongY,
                                          scratchView.range({3, 4}),
                                          recalculatedTop);
    //
    if (recalculatedBottom.empty() || recalculatedTop.empty()) {
      return;
    }
    const ConstSpacePointProxy spM =
        seeding[static_cast<SpacePointIndex>(pM.index)];
    float curvature = 0.0f;
    if (m_useTripletFinder) {
      recalculatedCandidates.clear();
      static_cast<void>(m_tripletFinder->createTripletTopCandidates(
          seeding, spM, recalculatedBottom[0], recalculatedTop.range(),
          recalculatedCandidates));
      if (recalculatedCandidates.size() == 0) {
        return;
      }
      curvature = recalculatedCandidates[0].curvature();
    }
    emit(pB, pM, pT,
         seedQuality(recalculatedBottom[0], recalculatedTop[0], spM, curvature),
         true);
  };

  // STEP 7: the space points that take part, in the grid: by z bin (inserted in
  // container order) and, within each bin, by y (the candidate pools)
  CartesianSpacePointGrid grid(m_gridConfig);
  for (std::size_t i = 0; i < ev.zBin.size(); ++i) {
    if (ev.zBin[i] != kNoZBin) {
      const Point& p = ev.points[i];
      const Vector3 position = m_cuts.useCartesianGridBinning
                                   ? Vector3(p.x, p.y, p.z)
                                   : Vector3(p.x, p.y, ev.zBin[i] + 0.5);
      grid.insert(static_cast<SpacePointIndex>(i), position);
    }
  }
  // ascending y, ties by index (so that the order does not depend on the sort)
  const auto byY = [&ev](SpacePointIndex l, SpacePointIndex r) {
    const double yl = ev.points[l].y;
    const double yr = ev.points[r].y;
    return yl < yr || (yl == yr && l < r);
  };
  for (std::size_t bin = 0; bin < grid.numberOfBins(); ++bin) {
    std::ranges::sort(grid.at(bin), byY);
  }
  // the middle SPs of a bin, in container order (the order of the seeds)
  std::vector<SpacePointIndex> middles;

  // STEP 8: for each z bin with activity, the central point ("middle"/r0) is
  // taken only from that bin, while the candidates for "top" and "bottom" are
  // taken from the union of their neighboring bins (buildFrameWork() /
  // production3Sp: rb/rt are filled from m_rfz_ib/m_rfz_it) this reduces the
  // combinatorial complexity from N^3 to something much smaller when there are
  // many space points per event.
  for (const auto [bottomBins, middleBin, topBins] : grid.binnedGroup()) {
    // pool: union of the neighboring bins of the current bin, ordered by y
    // (the bottom and top bins are the same, see m_gridConfig)
    pool.clear();
    for (const std::size_t bin : topBins) {
      const std::vector<SpacePointIndex>& binSpacePoints = grid.at(bin);
      const auto merged =
          pool.insert(pool.end(), binSpacePoints.begin(), binSpacePoints.end());
      std::inplace_merge(pool.begin(), merged, pool.end(), byY);
    }
    if (pool.empty()) {
      continue;
    }
    middles.assign(grid.at(middleBin).begin(), grid.at(middleBin).end());
    std::ranges::sort(middles);

    for (const SpacePointIndex indexM : middles) {
      const Point& pM = ev.points[indexM];
      const double yM = pM.y;
      const bool middleRecalculated = ev.calibratable(indexM);
      // the cheap kink rejection of passGeometricCuts, also done here: most
      // combinations of space points of different muons stop there, and
      // here they do not pay the call
      const bool kinkPrefilter =
          !m_ptDependent &&
          !(m_cuts.useSubdetectorAngularCuts && !ev.strip(indexM));
      auto kinkRejects = [&](const Segment& s1, const Segment& s2) {
        return kinkPrefilter && s1.dx * s2.dx + s1.dy * s2.dy + s1.dz * s2.dz <
                                    m_cosKinkReject * (s1.length * s2.length);
      };

      // candidates in the y window: below the middle SP [b0, b1) and above it
      // [t0, t1); a bottom must be below and a top above the middle SP, unless
      // requireTopMiddleBottom is off, where both are candidates of both
      auto yOf = [&ev](SpacePointIndex i) { return ev.points[i].y; };
      // b0: first SP in the pool with yM - y <= deltaYMax
      const auto b0 = std::partition_point(
          pool.begin(), pool.end(),
          [&](SpacePointIndex i) { return yM - yOf(i) > deltaYMax; });
      // b1: starting from b0, the first SP that is no longer underneath, with a
      // separation >= deltaYMin
      //  [b0, b1) = candidates for bottom
      const auto b1 =
          std::partition_point(b0, pool.end(), [&](SpacePointIndex i) {
            const double y = yOf(i);
            return y < yM && yM - y >= deltaYMin;
          });
      // t0: starting from b1, skip those that are too close to the middle
      const auto t0 =
          std::partition_point(b1, pool.end(), [&](SpacePointIndex i) {
            const double y = yOf(i);
            return !(y > yM && y - yM >= deltaYMin);
          });
      // t1: from t0, those above are at <= deltaYMax
      // [t0, t1) = top candidates
      const auto t1 = std::partition_point(
          t0, pool.end(),
          [&](SpacePointIndex i) { return yOf(i) - yM <= deltaYMax; });
      // subset of indices which are passed on to DoubletSeedFinder
      const SpacePointContainer::IndexSubset below(b0, b1);
      const SpacePointContainer::IndexSubset above(t0, t1);

      // doublets, both lists in ascending y
      const ConstSpacePointProxy spM = seeding[indexM];
      bottomDoublets.clear();
      bottomFinder->createDoublets(spM, kFrameAlongY, seeding.subset(below),
                                   bottomDoublets);
      if (!ordered) {
        bottomFinder->createDoublets(reflected[indexM], kFrameAlongY,
                                     reflected.subset(above), bottomDoublets);
      }
      if (bottomDoublets.empty()) {
        continue;
      }
      topDoublets.clear();
      if (!ordered) {
        topFinder->createDoublets(reflected[indexM], kFrameAlongY,
                                  reflected.subset(below), topDoublets);
      }
      topFinder->createDoublets(spM, kFrameAlongY, seeding.subset(above),
                                topDoublets);
      if (topDoublets.empty()) {
        continue;
      }

      const DoubletsForMiddleSp::Index nBottoms = bottomDoublets.size();
      const DoubletsForMiddleSp::Index nTops = topDoublets.size();
      bottomSegments.resize(nBottoms);
      for (DoubletsForMiddleSp::Index b = 0; b < nBottoms; ++b) {
        bottomSegments[b] =
            segment(pM, ev.points[bottomDoublets[b].spacePointIndex()]);
      }
      // d1 = middle - top segment of each doublet
      topSegments.resize(nTops);
      for (DoubletsForMiddleSp::Index t = 0; t < nTops; ++t) {
        const SpacePointIndex indexT = topDoublets[t].spacePointIndex();
        topSegments[t] = segment(ev.points[indexT], pM);
        // index map of SP's top
        if (m_useTripletFinder) {
          topOfSpacePoint[indexT] = t;
        }
      }
      for (DoubletsForMiddleSp::Index n = 0; n < nBottoms; ++n) {
        // with requireTopMiddleBottom in a vertical order, from the lowest to
        // the one furthest away
        const DoubletsForMiddleSp::Index b = ordered ? nBottoms - 1 - n : n;
        const SpacePointIndex indexB = bottomDoublets[b].spacePointIndex();
        const Point& pB = ev.points[indexB];
        // If the middle or the bottom are recalculable, all the triplets in
        // this bottom are processed via tryRecalculated
        const bool bottomRecalculated =
            middleRecalculated || ev.calibratable(indexB);

        // STEP 9: clear candidates from the previous bottom
        candidates.clear();
        if (m_useTripletFinder && !bottomRecalculated) {
          static_cast<void>(m_tripletFinder->createTripletTopCandidates(
              seeding, spM, bottomDoublets[b], topDoublets.range(),
              candidates));
          // case without recalculable strips in the event
          if (!ev.anyCalibratable) {
            // only the tops that passed, in ascending y
            for (TripletTopCandidates::Index c = 0; c < candidates.size();
                 ++c) {
              // top doublet of candidate c
              const DoubletsForMiddleSp::Index t =
                  topOfSpacePoint[candidates[c].spacePoint()];
              const SpacePointIndex indexT = topDoublets[t].spacePointIndex();
              // other cuts are symmetric under top<->bottom, without the
              // ordering each unordered triplet would appear twice
              if (!ordered && indexT <= indexB) {
                continue;
              }
              const Point& pT = ev.points[indexT];
              // kink prefilter and geometric cuts then
              if (!kinkRejects(topSegments[t], bottomSegments[b]) &&
                  passGeometricCuts(pT, pM, pB, topSegments[t],
                                    bottomSegments[b], ev)) {
                // returns the seed
                emit(pB, pM, pT,
                     seedQuality(bottomDoublets[b], topDoublets[t], spM,
                                 candidates[c].curvature()),
                     false);
              }
            }
            continue;
          }
        }

        // STEP 10: every top in ascending y (without TripletSeedFinder, or with
        // recalculated strip SPs, whose triplets are interleaved with the
        // others so that the seeds keep the order of the previous version);
        // the candidates of the TripletSeedFinder come in the order of the tops
        TripletTopCandidates::Index c = 0;
        for (DoubletsForMiddleSp::Index t = 0; t < nTops; ++t) {
          const SpacePointIndex indexT = topDoublets[t].spacePointIndex();
          bool passed = !m_useTripletFinder;
          float curvature = 0.0f;
          if (c < candidates.size() && candidates[c].spacePoint() == indexT) {
            passed = true;
            curvature = candidates[c].curvature();
            ++c;
          }
          // other cuts are symmetric under top<->bottom, without the ordering
          // each unordered triplet would appear twice
          if (!ordered && indexT <= indexB) {
            continue;
          }
          const Point& pT = ev.points[indexT];
          // recalculate strips
          if (bottomRecalculated || ev.calibratable(indexT)) {
            tryRecalculated(pB, pM, pT);
            continue;
          }
          // if it does not pass the finder, kink or geometric cuts, there is no
          // seed
          if (!passed || kinkRejects(topSegments[t], bottomSegments[b]) ||
              !passGeometricCuts(pT, pM, pB, topSegments[t], bottomSegments[b],
                                 ev)) {
            continue;
          }
          emit(pB, pM, pT,
               seedQuality(bottomDoublets[b], topDoublets[t], spM, curvature),
               false);
        }
      }
    }
  }
}

// the segment from `from` to `to`
// segment(): the vector between two points and its length
CosmicTripletSeedFinder::Segment CosmicTripletSeedFinder::segment(
    const Point& from, const Point& to) {
  Segment s;
  // vector components: to - from
  s.dx = to.x - from.x;
  s.dy = to.y - from.y;
  s.dz = to.z - from.z;
  // 3D length
  s.length = std::sqrt(s.dx * s.dx + s.dy * s.dy + s.dz * s.dz);
  return s;
}

// passGeometricCuts applies vertical order, kink angle, vertical separation,
// residual relative to the top-bottom line and z envelope
// the curvature and the longitudinal compatibility (dT) are cut by the
// TripletSeedFinder
bool CosmicTripletSeedFinder::passGeometricCuts(const Point& top,
                                                const Point& middle,
                                                const Point& bottom,
                                                const EventData& ev) const {
  // recalculated position should be finite
  if (!std::isfinite(top.x) || !std::isfinite(top.y) || !std::isfinite(top.z) ||
      !std::isfinite(middle.x) || !std::isfinite(middle.y) ||
      !std::isfinite(middle.z) || !std::isfinite(bottom.x) ||
      !std::isfinite(bottom.y) || !std::isfinite(bottom.z)) {
    return false;
  }
  // calculates d1 = middle - top and d2 = bottom - middle
  return passGeometricCuts(top, middle, bottom, segment(top, middle),
                           segment(middle, bottom), ev);
}

// notation: d1 = middle - top (topSegment), d2 = bottom - middle
// (bottomSegment), chord = bottom - top
// the cuts are ordered from the computationally cheapest and most selective to
// the most expensive ones
bool CosmicTripletSeedFinder::passGeometricCuts(const Point& top,
                                                const Point& middle,
                                                const Point& bottom,
                                                const Segment& topSegment,
                                                const Segment& bottomSegment,
                                                const EventData& ev) const {
  // s1 = d1 (top), s2 = d2 (bottom)
  const Segment& s1 = topSegment;
  const Segment& s2 = bottomSegment;
  // pixel middle with cuts per subdetector
  const bool pixelAngularCut =
      m_cuts.useSubdetectorAngularCuts && ev.isStrip[middle.index] == 0;
  // the 3D kink is cut unless it is pT dependent (polar kink) or the middle
  // SP is a pixel with useSubdetectorAngularCuts (cos(dphi))
  const bool kinkCut = !m_ptDependent && !pixelAngularCut;

  // the kink certainly above maxKinkAngle: cos(kink) = d1.d2 / (|d1| |d2|)
  const double dot = s1.dx * s2.dx + s1.dy * s2.dy + s1.dz * s2.dz;
  // rejects the kink without division and acos
  if (kinkCut && dot < m_cosKinkReject * (s1.length * s2.length)) {
    return false;
  }
  // segment with length 0 -> two SPs at the same point, so it's rejected
  if (s1.length == 0.0 || s2.length == 0.0) {
    return false;
  }
  if (pixelAngularCut) {
    // transverse lengths, dphi = phi2 - phi1, with phi_i = atan2(dy_i, dx_i)
    const double c1 = s1.transverse();
    const double c2 = s2.transverse();
    if (c1 == 0.0 || c2 == 0.0) {
      return false;
    }
    // cos(phi2 - phi1) = (d1_xy . d2_xy)/(|d1_xy||d2_xy|)
    const double cosDphi = (s1.dx * s2.dx + s1.dy * s2.dy) / (c1 * c2);
    // cos(phi2 - phi1) should be higher than pixelMinCosDeltaPhi
    if (!(cosDphi > m_cuts.pixelMinCosDeltaPhi)) {
      return false;
    }
  }

  // all the points must be of the same type (pixel or strip) if
  // requireSameSubdetector is enabled
  if (m_cuts.requireSameSubdetector) {
    const std::uint8_t kind = ev.isStrip[middle.index];
    if (ev.isStrip[top.index] != kind || ev.isStrip[bottom.index] != kind) {
      return false;
    }
  }

  // vertical order: top > middle > bottom in y
  if (m_cuts.requireTopMiddleBottom &&
      !(top.y > middle.y && middle.y > bottom.y)) {
    return false;
  }

  // 3 SPs in [zMin, zMax]
  if (top.z < m_cuts.zMin || top.z > m_cuts.zMax || middle.z < m_cuts.zMin ||
      middle.z > m_cuts.zMax || bottom.z < m_cuts.zMin ||
      bottom.z > m_cuts.zMax) {
    return false;
  }

  // Δy per segment in [minDeltaY, maxDeltaY]
  // |dy| of each element in the direct window
  const double dyTop = std::abs(top.y - middle.y);
  const double dyBottom = std::abs(middle.y - bottom.y);
  if (dyTop < m_cuts.minDeltaY || dyTop > m_cuts.maxDeltaY ||
      dyBottom < m_cuts.minDeltaY || dyBottom > m_cuts.maxDeltaY) {
    return false;
  }

  // chord = bottom - top
  const double dx = bottom.x - top.x;
  const double dy = bottom.y - top.y;
  const double dz = bottom.z - top.z;
  // chord length
  const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
  // rejects overlapping points
  if (distance == 0.0) {
    return false;
  }
  // vertical extension of triplets
  const double verticalSpan = std::abs(dy);
  // cosThetaY = cosine of the angle between the chord and the y-axis
  const double cosThetaY = verticalSpan / distance;
  if (cosThetaY < m_cuts.minCosThetaY) {
    return false;
  }
  // distance from the middle point to the chord, |d1 x chord| / |chord|
  auto residualToChord = [&]() {
    const double cx = s1.dy * dz - s1.dz * dy;
    const double cy = s1.dz * dx - s1.dx * dz;
    const double cz = s1.dx * dy - s1.dy * dx;
    return std::sqrt(cx * cx + cy * cy + cz * cz) / distance;
  };

  // kink and residual
  double kinkLimit = m_cuts.maxKinkAngle;
  double residualLimit = m_cuts.maxResidual;
  PtDependentTerms terms;
  if (m_ptDependent) {
    // if pT cannot be estimated (lengths = 0) or if pT < 0.9 minPt, it is
    // rejected
    if (!ptDependentTerms(s1, s2, residualToChord(), terms)) {
      return false;
    }
    kinkLimit += terms.kinkSlack;
    residualLimit += terms.residualSlack;
  }
  // with pT, the kink that is cut off is the polar |theta1 - theta2|
  if (m_ptDependent && !pixelAngularCut && terms.polarKink > kinkLimit) {
    return false;
  }
  if (kinkCut) {
    // the acos only near maxKinkAngle (see kKinkCosBand)
    const double cosAngle =
        std::clamp(dot / (s1.length * s2.length), -1.0, 1.0);
    // (a separate unlikely branch: otherwise the compiler computes the acos
    // of every triplet in advance)
    if (cosAngle < m_cosKinkAccept) [[unlikely]] {
      if (std::acos(cosAngle) > kinkLimit) {
        return false;
      }
    }
  }

  // maxResidual < 0 = without residual cutoff
  if (m_cuts.maxResidual >= 0.0) {
    const double residual =
        m_ptDependent ? terms.parallelResidual : residualToChord();
    if (residual > residualLimit) {
      return false;
    }
  }
  return true;
}

// pT-dependent terms: curvature, kink, residual, and scattering tolerance
// Notation: d1 = middle - top, d2 = bottom - middle
// c_i = transverse length (xy plane) of d_i, L_i = 3D length, chord = bottom -
// top
bool CosmicTripletSeedFinder::ptDependentTerms(const Segment& topSegment,
                                               const Segment& bottomSegment,
                                               double residual,
                                               PtDependentTerms& out) const {
  // d1 = middle - top and d2 = bottom - middle components
  const double x1 = topSegment.dx;
  const double y1 = topSegment.dy;
  const double z1 = topSegment.dz;
  const double x2 = bottomSegment.dx;
  const double y2 = bottomSegment.dy;
  const double z2 = bottomSegment.dz;

  // c1, c2: lengths in xy of each segment.
  const double c1 = topSegment.transverse();
  const double c2 = bottomSegment.transverse();
  // cc: length in xy of the chord d1 + d2
  const double cc = std::hypot(x1 + x2, y1 + y2);
  // L1, L2: 3D lengths of the segments
  const double L1 = topSegment.length;
  const double L2 = bottomSegment.length;
  // Lc: 3D length of the chord
  const double Lc = std::sqrt((x1 + x2) * (x1 + x2) + (y1 + y2) * (y1 + y2) +
                              (z1 + z2) * (z1 + z2));
  // any length of 0 would result in division by 0
  if (c1 == 0.0 || c2 == 0.0 || cc == 0.0 || Lc == 0.0) {
    return false;
  }

  // pT estimated: the circle through the three points in the xy plane has the
  // radius of the Menger formula, and pT cut
  const double cross = x1 * y2 - y1 * x2;
  // Menger's formula: 1/R = 4 Area/(a b c) = 2|cross|/(c1 c2 cc)
  const double invRho = 2.0 * std::abs(cross) / (c1 * c2 * cc);  // 1/mm
  // 1/pT = (1/R)/(0.3 |Bz|)
  const double invPt = invRho / m_ptPerRadius;  // 1/MeV
  // pT >= 0.9 minPt cut
  if (invPt > m_maxInvPt) {
    return false;
  }
  // effective scattering pT = min(pT, maxPtForScattering)
  const double invPtEff = std::max(invPt, m_minInvPtScattering);
  // sin(theta) of the chord with respect to z
  const double sinThetaChord = cc / Lc;
  // theta_MS = sqrt(COF) sin(theta)/pT = sqrt(COF)/p
  const double thetaMS = m_sqrtScatteringCof * sinThetaChord * invPtEff;
  // kink = theta_MS tolerance
  out.kinkSlack = thetaMS;
  // L2/(L1 + L2) height of the triangle
  out.residualSlack = thetaMS * L1 * L2 / (L1 + L2);

  // the transverse turn is already bounded by the pT cut
  // the kink compares the polar deviation |theta1 - theta2|, with theta_i =
  // atan2(c_i, dz_i)
  out.polarKink = std::abs(topSegment.polar() - bottomSegment.polar());
  // the residual is the distance from the middle point to the line connecting
  // top and bottom which can be decomposed into a transverse component (h_T)
  // and a parallel component (h_par) residual out of the transverse plane:
  // h3D^2 = h_T^2 + h_par^2
  const double hT = std::abs(cross) / cc;
  out.parallelResidual =
      std::sqrt(std::max(0.0, residual * residual - hT * hT));
  return true;
}

// track direction at bottom, middle and top (dirs[0..2])
// without magnetic field, the direction is the bottom-top chord
// with magnetic field, each SP uses the tangent at its own point:
// in xy, perpendicular to the radius of the circle through the three points
// in z, the slope of the chord
bool CosmicTripletSeedFinder::tripletDirections(
    const Point& top, const Point& middle, const Point& bottom,
    std::array<std::array<float, 3>, 3>& dirs) const {
  // top - bottom chord
  const double cx = top.x - bottom.x;
  const double cy = top.y - bottom.y;
  const double cz = top.z - bottom.z;
  const double chord = std::sqrt(cx * cx + cy * cy + cz * cz);
  if (!(chord > 0.0)) {
    return false;
  }
  // normalized straight-line direction, as a float (as required by detail::).
  const std::array<float, 3> straight = {static_cast<float>(cx / chord),
                                         static_cast<float>(cy / chord),
                                         static_cast<float>(cz / chord)};
  // the three SPs use the straight direction by default
  dirs = {straight, straight, straight};
  // chord length in xy
  const double chordXY = std::hypot(cx, cy);
  // without field -> straight direction
  if (m_cuts.magneticFieldBz == 0.0 || chordXY == 0.0) {
    return true;
  }

  // positions relative to the middle a = bottom - middle, q = top - middle
  const double ax = bottom.x - middle.x;
  const double ay = bottom.y - middle.y;
  const double qx = top.x - middle.x;
  const double qy = top.y - middle.y;
  // 2x2 system determinant
  const double det = ax * qy - ay * qx;
  const double a2 = ax * ax + ay * ay;
  const double q2 = qx * qx + qy * qy;
  if (det == 0.0) {
    return true;
  }
  // center of the circle by Cramer's rule: 2 u.a = |a|^2, 2 u.q = |q|^2
  const double ux = 0.5 * (a2 * qy - q2 * ay) / det;
  const double uy = 0.5 * (ax * q2 - qx * a2) / det;
  // radius > 100 km: practically straight; the direction is left straight
  constexpr double kStraightRadius = 1e8;  // mm
  if (!(std::hypot(ux, uy) < kStraightRadius)) {
    return true;
  }
  // dz/ds slope in xy
  const double slope = cz / chordXY;  // dz / ds_xy
  // relative positions of bottom, middle and top
  const std::array<std::array<double, 2>, 3> rel = {
      {{ax, ay}, {0.0, 0.0}, {qx, qy}}};
  for (std::size_t i = 0; i < 3; ++i) {
    // tangent in xy: perpendicular to the radius, oriented like the chord
    double tx = -(rel[i][1] - uy);
    double ty = rel[i][0] - ux;
    if (tx * cx + ty * cy < 0.0) {
      tx = -tx;
      ty = -ty;
    }
    // xy normalization
    const double tn = std::hypot(tx, ty);
    if (!(tn > 0.0)) {
      return false;
    }
    const double dx = tx / tn;
    const double dy = ty / tn;
    // unit 3D vectors (dx, dy, slope)/sqrt(1 + slope^2).
    const double n = std::sqrt(1.0 + slope * slope);
    dirs[i] = {static_cast<float>(dx / n), static_cast<float>(dy / n),
               static_cast<float>(slope / n)};
  }
  return true;
}

// recalculates the strip SPs of the seed with the direction of the track
// with the direction of the seed the crossing is calculated, with the new
// crossings a new direction, and so on (at most stripCalibrationMaxIterations
// times) until no point moves more than stripCalibrationConvergence then, with
// the final direction, the crossing must fall within the two strips with
// stripCalibrationTolerance (if not, the seed is not a track)
bool CosmicTripletSeedFinder::calibrateStripSpacePoints(
    Point& top, Point& middle, Point& bottom, const EventData& ev) const {
  // pointers to the three points in the order: bottom, middle, top
  const std::array<Point*, 3> points = {&bottom, &middle, &top};
  // identify which of the three are recalculable strips
  std::array<bool, 3> strip{};
  for (std::size_t i = 0; i < 3; ++i) {
    strip[i] = ev.calibratable(points[i]->index);
  }
  // direction of the triplet in bottom, middle and top
  auto directions = [&](std::array<std::array<float, 3>, 3>& dirs) {
    return tripletDirections(top, middle, bottom, dirs);
  };
  // note: a direction (almost) parallel to the module plane crosses the strips
  // at infinity, so its crossing is not valid
  const double minCosToNormal = m_cuts.stripCalibrationMinCosToNormal;
  // crossesModule: n = ihv x ohv is normal to the plane of the strips
  auto crossesModule = [&ev, minCosToNormal](const std::array<float, 3>& d,
                                             std::size_t index) {
    const auto& n = ev.stripDetails[index].innerCrossOuterHalfVector;
    const double scale = static_cast<double>(d[0]) * n[0] +
                         static_cast<double>(d[1]) * n[1] +
                         static_cast<double>(d[2]) * n[2];
    const double norm = std::sqrt(static_cast<double>(n[0]) * n[0] +
                                  static_cast<double>(n[1]) * n[1] +
                                  static_cast<double>(n[2]) * n[2]);
    return std::abs(scale) > minCosToNormal * norm;
  };
  // infinite tolerance during iterations: no rejection occurs while the
  // direction is still being adjusted
  constexpr float kUnboundedCrossing = std::numeric_limits<float>::infinity();
  std::array<std::array<float, 3>, 3> dirs{};
  bool converged = false;
  for (int iteration = 0; iteration < m_cuts.stripCalibrationMaxIterations;
       ++iteration) {
    // direction with current positions
    if (!directions(dirs)) {
      return false;
    }
    // largest displacement in this iteration
    double maxShift = 0.0;
    for (std::size_t i = 0; i < 3; ++i) {
      if (!strip[i]) {
        continue;
      }
      // rejects if direction is almost parallel to the module
      if (!crossesModule(dirs[i], points[i]->index)) {
        return false;
      }
      std::array<float, 3> r{};
      // intersection of the direction with the outer strip: new r position
      detail::calibrateOuterStripSpacePoint(
          dirs[i], ev.stripDetails[points[i]->index], r, kUnboundedCrossing);
      // a non-finite crossing invalidates the triplet
      if (!std::isfinite(r[0]) || !std::isfinite(r[1]) ||
          !std::isfinite(r[2])) {
        return false;
      }
      // displacement relative to the previous position
      const double sx = r[0] - points[i]->x;
      const double sy = r[1] - points[i]->y;
      const double sz = r[2] - points[i]->z;
      maxShift = std::max(maxShift, std::sqrt(sx * sx + sy * sy + sz * sz));
      // new position is saved
      points[i]->x = r[0];
      points[i]->y = r[1];
      points[i]->z = r[2];
    }
    // if the maximum shift is less than the convergence threshold, the
    // calibration is considered converged
    if (maxShift < m_cuts.stripCalibrationConvergence) {
      converged = true;
      break;
    }
  }
  // without convergence or without a valid final direction -> no seed
  if (!converged || !directions(dirs)) {
    return false;
  }
  // tolerance for the crossing of the strips with the final direction
  // if the crossing is outside the two strips by more than this tolerance, the
  // seed is rejected
  const auto tolerance = static_cast<float>(m_cuts.stripCalibrationTolerance);
  // With the final orientation, each recalculable strip must cross the module
  // and fall within the two strips, accounting for the tolerance
  for (std::size_t i = 0; i < 3; ++i) {
    std::array<float, 3> r{};
    if (strip[i] &&
        (!crossesModule(dirs[i], points[i]->index) ||
         !detail::calibrateOuterStripSpacePoint(
             dirs[i], ev.stripDetails[points[i]->index], r, tolerance))) {
      return false;
    }
  }
  return true;
}

}  // namespace Acts
