//cosmictripletseedfinder.cpp

#include "Acts/Seeding/CosmicTripletSeedFinder.hpp"
#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/SpacePointFormation/detail/StripSpacePointCalibrationImpl.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include <utility>

namespace Acts {

namespace {

double varianceR(const ConstSpacePointProxy& sp) {
  return static_cast<double>(sp.varianceR());
}
double varianceZ(const ConstSpacePointProxy& sp) {
  return static_cast<double>(sp.varianceZ());
}

constexpr int kNumZBins = 11;

// zBin() returns the index of the z bin for a given z coordinate, based on predefined ranges
int zBin(double z) {
  if (z > 0.0) {
    if (z < 250.0) return 5;
    if (z < 450.0) return 6;
    if (z < 925.0) return 7;
    if (z < 1400.0) return 8;
    if (z < 2500.0) return 9;
    return 10;
  }
  if (z > -250.0) return 5;
  if (z > -450.0) return 4;
  if (z > -925.0) return 3;
  if (z > -1400.0) return 2;
  if (z > -2500.0) return 1;
  return 0;
}

// neighbors are used to compare near space points from possible trayectories
const std::array<std::vector<int>, kNumZBins> kZNeighbors = {{
    {0, 1},
    {0, 1, 2},
    {1, 2, 3},
    {2, 3, 4, 5},
    {3, 4, 5},
    {4, 5, 6},
    {5, 6, 7},
    {5, 6, 7, 8},
    {7, 8, 9},
    {8, 9, 10},
    {9, 10},
}};


// kPrefilterTol avoid numerical issues in the prefiltering cuts (its never stricter than passCosmicCuts)
constexpr double kPrefilterTol = 1e-9;

// Local copy (double) of a space point + the original proxy
// which is still used for exact cuts and for the indices of the seed
struct Hit {
  ConstSpacePointProxy sp;
  double x, y, z;
};

using ZBins = std::array<std::vector<Hit>, kNumZBins>;

// fillLists() reorganizes in 11 bins the space points of z
ZBins fillZBins(const SpacePointContainer& spacePoints, double zMin,
                double zMax, const std::vector<std::uint8_t>* skip) {
  ZBins bins;
  for (ConstSpacePointProxy sp : spacePoints) {
    if (skip != nullptr && (*skip)[sp.index()] != 0) {
      continue;
    }
    const double x = sp.x();
    const double y = sp.y();
    const double z = sp.z();
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
      continue;
    }
    if (z < zMin || z > zMax) {
      continue;
    }
    bins[zBin(z)].push_back(Hit{sp, x, y, z});
  }
  return bins;
}

// Candidate (top or bottom) relative to a middle space point
struct Candidate {
  const Hit* hit;
  double dx, dy;      // distances middle -> candidate
  double cT2;         // dx^2 + dy^2  (transverse distance squared)
  double y;           // absolute y of the candidate (for sorting)
  double ux, uy, uz;  // 3D unit vector middle -> candidate
  bool calibratable;  // true if the candidate or the middle space point is a strip that can be recalibrated
};

// calibrateStripSpacePoints recalibrates the y of the space points of a triplet
// if any of them is a strip that can be recalibrated
constexpr int kMaxCalibrationIterations = 50;
constexpr double kCalibrationConvergence = 1e-3;
// |cos| minimum between the direction of the triplet and the normal of the module,
// below which the recalibration is not attempted
constexpr double kMinCosToModuleNormal = 1e-6;
}



// CosmicTripletSeedFinder constructor: checks the validity of the cuts and precomputes some derived quantities
CosmicTripletSeedFinder::CosmicTripletSeedFinder(const CosmicSeedCuts& cuts)
    : m_cuts(cuts) {
  if (!(m_cuts.zMin < m_cuts.zMax)) {
    throw std::invalid_argument("CosmicSeedCuts: zMin must be < zMax");
  }
  if (m_cuts.maxKinkAngle < 0.0) {
    throw std::invalid_argument("CosmicSeedCuts: maxKinkAngle must be >= 0");
  }
  if (!(m_cuts.minDeltaY >= 0.0) || !(m_cuts.maxDeltaY > m_cuts.minDeltaY)) {
    throw std::invalid_argument(
        "CosmicSeedCuts: need 0 <= minDeltaY < maxDeltaY");
  }
  if (m_cuts.useStripCalibration && !(m_cuts.stripCalibrationTolerance > 0.0)) {
    throw std::invalid_argument(
        "CosmicSeedCuts: stripCalibrationTolerance must be > 0");
  }
  if (m_cuts.minCosThetaY < 0.0 || m_cuts.minCosThetaY > 1.0) {
    throw std::invalid_argument(
        "CosmicSeedCuts: minCosThetaY must be in [0, 1]");
  }
  if (!(m_cuts.maxAbsEta > 0.0)) {
    throw std::invalid_argument("CosmicSeedCuts: maxAbsEta must be > 0");
  }
  if (m_cuts.minPt < 0.0) {
    throw std::invalid_argument("CosmicSeedCuts: minPt must be >= 0");
  }
  if (m_cuts.multipleScatteringCof < 0.0) {
    throw std::invalid_argument(
        "CosmicSeedCuts: multipleScatteringCof must be >= 0");
  }

  // acceptance window in dz/dr = dz / sqrt(dx^2 + dy^2)
  const double etaMax = std::abs(m_cuts.maxAbsEta);
  m_dzdrMax = 1.0 / std::tan(2.0 * std::atan(std::exp(-etaMax)));

  // check covarianceScale entries if used
  if (m_cuts.useCovarianceScale) {
    for (double s : m_cuts.covarianceScale) {
      if (!(s > 0.0) || !std::isfinite(s)) {
        throw std::invalid_argument(
            "CosmicSeedCuts: covarianceScale entries must be finite and > 0");
      }
    }
  }

  // check pixelMinCosDeltaPhi if used
  if (m_cuts.useSubdetectorAngularCuts &&
      !(m_cuts.pixelMinCosDeltaPhi >= -1.0 &&
        m_cuts.pixelMinCosDeltaPhi <= 1.0)) {
    throw std::invalid_argument(
        "CosmicSeedCuts: pixelMinCosDeltaPhi must be in [-1, 1]");
  }

  // check pt-dependent kink residual cuts if used
  if (m_cuts.usePtDependentKinkResidual) {
    if (!m_cuts.usePtCut || !(m_cuts.minPt > 0.0) ||
        m_cuts.magneticFieldBz == 0.0) {
      throw std::invalid_argument(
          "CosmicSeedCuts: usePtDependentKinkResidual requires usePtCut, "
          "minPt > 0 and magneticFieldBz != 0");
    }
    if (!(m_cuts.maxPtForScattering > 0.0)) {
      throw std::invalid_argument(
          "CosmicSeedCuts: maxPtForScattering must be > 0");
    }
  }

  // Derived quantities for prefiltering and passCosmicCuts
  m_ptDependent = m_cuts.usePtDependentKinkResidual;
  m_ptPrefilter = m_cuts.usePtCut && m_cuts.minPt > 0.0 &&
                  m_cuts.magneticFieldBz != 0.0;
  if (m_ptPrefilter) {
    // pT_est >= 0.9 minPt  <=>  rho >= 0.9 minPt / (0.3 |Bz|)  [mm]
    m_rhoMin = 0.9 * m_cuts.minPt / (0.3 * std::abs(m_cuts.magneticFieldBz));
  }
  if (!m_ptDependent) {
    // Kink: 1 - k^2/2 <= cos(k) for every k in [0, maxKinkAngle]
    const double k = m_cuts.maxKinkAngle;
    m_kinkCosThreshold = 1.0 - 0.5 * k * k;
    // add of the both segments, and each segment forms at most "k" with it
    // the angle of the chord with the y-axis is at most acos(minCosThetaY);
    // by the triangular inequality on the sphere,
    // each segment forms at most acos(minCosThetaY) + k with the y-axis.
    const double maxAngleY = std::acos(m_cuts.minCosThetaY) + k;
    m_segCosThetaYMin = maxAngleY < 0.5 * M_PI ? std::cos(maxAngleY) : 0.0;
  }
}

CosmicTripletSeedFinder::EventData CosmicTripletSeedFinder::prepareEvent(
    const SpacePointContainer& spacePoints,
    const std::vector<CosmicSpacePointClass>* spClasses,
    const std::vector<GeometryIdentifier>* spSurfaces, bool needCovariance,
    bool needType, EventBuffers& buffers) const {
  // to prepare the event, we need to fill the buffers with the necessary information for each space point,
  // depending on the cuts that are enabled. This includes covariance scaling and type information,
  // as well as any surface identifiers if rejectSameSurface is enabled
  EventData ev;
  ev.spacePoints = &spacePoints;

  // the veto compares ids: one without the sensitive part (sensitive() == 0) is a layer or a volume,
  // and would veto all pairs within it
  if (m_cuts.rejectSameSurface) {
    if (spSurfaces == nullptr || spSurfaces->size() != spacePoints.size()) {
      throw std::invalid_argument(
          "CosmicTripletSeedFinder: rejectSameSurface requires one surface "
          "GeometryIdentifier per space point");
    }
    for (const GeometryIdentifier& surface : *spSurfaces) {
      if (surface.sensitive() == 0) {
        throw std::invalid_argument(
            "CosmicTripletSeedFinder: rejectSameSurface requires sensitive "
            "surface ids, got one with sensitive() == 0");
      }
    }
    ev.surfaces = spSurfaces->data();
  }

  // covariance scaling and type information (pixel/strip) for each space point
  if (needCovariance || needType) {
    const std::size_t n = spacePoints.size();
    if (spClasses == nullptr || spClasses->size() != n) {
      throw std::invalid_argument(
          "CosmicTripletSeedFinder: useCovarianceScale / "
          "requireSameSubdetector / useSubdetectorAngularCuts require one "
          "CosmicSpacePointClass per space point");
    }
    // resize buffers only if needed, to avoid unnecessary allocations
    if (needCovariance) {
      buffers.covR.resize(n);
      buffers.covZ.resize(n);
    }
    // resize buffers only if needed, to avoid unnecessary allocations
    if (needType) {
      buffers.isStrip.resize(n);
    }
    // fill the buffers with the necessary information for each space point
    for (ConstSpacePointProxy sp : spacePoints) {
      const std::size_t i = sp.index();
      if (i >= n) {
        throw std::invalid_argument(
            "CosmicTripletSeedFinder: space point index out of range");
      }
      const auto cls = static_cast<std::size_t>((*spClasses)[i]);
      if (cls >= m_cuts.covarianceScale.size()) {
        throw std::invalid_argument(
            "CosmicTripletSeedFinder: invalid CosmicSpacePointClass");
      }
      const bool barrel = cls == static_cast<std::size_t>(
                                     CosmicSpacePointClass::BarrelPixel) ||
                          cls == static_cast<std::size_t>(
                                     CosmicSpacePointClass::BarrelStrip);
      const bool strip = cls == static_cast<std::size_t>(
                                    CosmicSpacePointClass::BarrelStrip) ||
                         cls == static_cast<std::size_t>(
                                    CosmicSpacePointClass::EndcapStrip);
      // fill the covariance buffers with the scaled variances if needed
      // errorsc: only the dominant coordinate is scaled (z in barrel, r in endcap);
      // the other remains with the variance of the container.
      if (needCovariance) {
        const double scale = m_cuts.covarianceScale[cls];
        buffers.covR[i] = barrel ? varianceR(sp) : varianceR(sp) * scale;
        buffers.covZ[i] = barrel ? varianceZ(sp) * scale : varianceZ(sp);
      }
      if (needType) {
        buffers.isStrip[i] = strip ? 1 : 0;
      }
    }
    // assign the covariance buffers to the event data
    if (needCovariance) {
      ev.covR = buffers.covR.data();
      ev.covZ = buffers.covZ.data();
    }
    // assign the type buffer to the event data
    if (needType) {
      ev.isStrip = buffers.isStrip.data();
    }
  }

  // Two strip SPs that are recalculated based on the direction of each seed
  // (useStripCalibration): those with non-zero values in StripCalibrationDetails
  // (pixel SPs retain their default value, all zeros).
  // stripShiftY limits how much the y coordinate of an SP can change when it is recalculated
  // (shift from the inner strip to the outer strip + movement along the outer strip
  // within the tolerance); this expands the window in the y direction.
  if (m_cuts.useStripCalibration &&
      spacePoints.hasColumns(SpacePointColumns::StripCalibrationDetails)) {
    const std::size_t n = spacePoints.size();
    buffers.isCalibratable.assign(n, 0);
    buffers.stripDetails.resize(n);
    for (ConstSpacePointProxy sp : spacePoints) {
      const OuterStripSpacePointCalibrationDetails& d =
          sp.outerStripCalibrationDetails();
      auto norm2 = [](const std::array<float, 3>& v) {
        return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
      };
      // skip SPs that are not recalibratable (e.g., pixel SPs or strip SPs with zero vectors)
      if (!(norm2(d.outerHalfVector) > 0.f) ||
          !(norm2(d.innerHalfVector) > 0.f)) {
        continue;
      }
      // mark the SP as recalibratable and store its calibration details
      buffers.isCalibratable[sp.index()] = 1;
      buffers.stripDetails[sp.index()] =
          detail::deriveOuterStripSpacePointCalibrationDetails(d);
      // Y shift for the recalibration: the maximum shift in y when recalibrating a strip SP
      // this is the sum of the absolute values of the inner to outer separation in y,
      // the inner half vector in y, and the outer half vector in y scaled by the strip calibration tolerance.
      const double shift =
          std::abs(static_cast<double>(d.innerToOuterSeparation[1])) +
          std::abs(static_cast<double>(d.innerHalfVector[1])) +
          m_cuts.stripCalibrationTolerance *
              std::abs(static_cast<double>(d.outerHalfVector[1]));
      buffers.stripShiftY = std::max(buffers.stripShiftY, shift);
    }
    ev.isCalibratable = buffers.isCalibratable.data();
    ev.stripDetails = buffers.stripDetails.data();
  }
  return ev;
}

// createSeeds() generates triplet seeds from the given space points, applying the configured cuts and storing the results in outputSeeds.
// It also optionally recalibrates strip space points and stores them in calibratedSpacePoints.
void CosmicTripletSeedFinder::createSeeds(
    const SpacePointContainer& spacePoints, SeedContainer& outputSeeds,
    const std::vector<CosmicSpacePointClass>* spClasses,
    SpacePointContainer* calibratedSpacePoints,
    const std::vector<GeometryIdentifier>* spSurfaces) const {
  if (!outputSeeds.hasSpacePointContainer()) {
    outputSeeds.assignSpacePointContainer(spacePoints);
  }

  outputSeeds.reserve(outputSeeds.size() + spacePoints.size(), 3.0f);

  const bool needType =
      m_cuts.requireSameSubdetector || m_cuts.useSubdetectorAngularCuts;
  EventBuffers buffers;
  const EventData ev = prepareEvent(spacePoints, spClasses, spSurfaces,
                                    m_cuts.useCovarianceScale, needType,
                                    buffers);
  const double stripShiftY = buffers.stripShiftY;
  auto calibratable = [&ev](std::size_t index) {
    return ev.isCalibratable != nullptr && ev.isCalibratable[index] != 0;
  };

  // Best quality is used to determine if a recalibrated strip space point should be updated in the calibratedSpacePoints container.
  // minor quality = better quality
  std::vector<float> bestQuality;
  if (calibratedSpacePoints != nullptr && ev.isCalibratable != nullptr) {
    if (calibratedSpacePoints->size() != spacePoints.size()) {
      throw std::invalid_argument(
          "CosmicTripletSeedFinder: calibratedSpacePoints must be a copy of "
          "spacePoints");
    }
    // assign the best quality vector to have the same size as the number of space points, initialized to infinity
    bestQuality.assign(spacePoints.size(), std::numeric_limits<float>::infinity());
  }

  const ZBins bins = fillZBins(spacePoints, m_cuts.zMin, m_cuts.zMax, nullptr);

  // tryEmit() attempts to create a seed from a triplet of space points (top, middle, bottom)
  auto tryEmit = [&](const Hit& top, const Hit& middle, const Hit& bottom) {
    Point pT{top.x, top.y, top.z, top.sp.index()};
    Point pM{middle.x, middle.y, middle.z, middle.sp.index()};
    Point pB{bottom.x, bottom.y, bottom.z, bottom.sp.index()};
    const bool anyStrip =
        calibratable(pT.index) || calibratable(pM.index) ||
        calibratable(pB.index);
    if (anyStrip && !calibrateStripSpacePoints(pT, &pM, pB, ev)) {
      return;
    }
    float quality = m_cuts.seedQuality;
    if (!passCosmicCuts(pT, pM, pB, quality, ev)) {
      return;
    }
    // acts convention:bottom-middle-top 
    std::array<SpacePointIndex, 3> seedSpacePoints{bottom.sp.index(), middle.sp.index(), top.sp.index()};

    // create a new seed in the outputSeeds container and assign the space point indices, quality, and estimated vertex Z position
    MutableSeedProxy seed = outputSeeds.createSeed();
    seed.assignSpacePointIndices(seedSpacePoints);
    seed.quality() = quality;
    seed.vertexZ() = estimateVertexZ(pT, pM, pB);

    // if any of the space points are recalibratable strip space points, update their positions
    // in the calibratedSpacePoints container if the new quality is better than the previous best quality
    if (anyStrip && !bestQuality.empty()) {
      for (const Point* p : {&pB, &pM, &pT}) {
        if (!calibratable(p->index) || !(quality < bestQuality[p->index])) {
          continue;
        }
        bestQuality[p->index] = quality;
        MutableSpacePointProxy out = calibratedSpacePoints->at(static_cast<SpacePointIndex>(p->index));
        out.x() = static_cast<float>(p->x);
        out.y() = static_cast<float>(p->y);
        out.z() = static_cast<float>(p->z);
        if (calibratedSpacePoints->hasColumns(SpacePointColumns::R)) {
          out.r() = static_cast<float>(std::hypot(p->x, p->y));
        }
        if (calibratedSpacePoints->hasColumns(SpacePointColumns::Phi)) {
          out.phi() = static_cast<float>(std::atan2(p->y, p->x));
        }
      }
    }
  };

  // preselection cuts: same subdetector, not in the same surface, dz/dr, kink/cosThetaY, pT (rho)
  // if any of the two is a recalibratable strip SP, the position is not yet the final
  // one and the angular cuts are left for passCosmicCuts.
  auto makeCandidate = [&](const Hit& h, const Hit& mid, bool pixelMid, Candidate& c) {
    if (m_cuts.requireSameSubdetector &&
        ev.isStrip[h.sp.index()] != ev.isStrip[mid.sp.index()]) {
      return false;
    }
    // surface veto: if the two SPs are in the same surface, they cannot be part of the same triplet
    if (ev.surfaces != nullptr &&
        ev.surfaces[h.sp.index()] == ev.surfaces[mid.sp.index()]) {
      return false;
    }
    const double dx = h.x - mid.x;
    const double dy = h.y - mid.y;
    const double dz = h.z - mid.z;
    const double cT2 = dx * dx + dy * dy;
    const bool strip = calibratable(h.sp.index()) || calibratable(mid.sp.index());
    // computeUV: r2 = dx^2 + dy^2 = 0 -> cannot compute the unit vector in the transverse plane
    if (cT2 == 0.0 && !strip) {
      return false;
    }
    // compute the 3D distance between the middle and candidate space points
    const double L = std::sqrt(cT2 + dz * dz);
    if (L == 0.0) {
      return false;
    }
    if (!strip) {
      const double cT = std::sqrt(cT2);
      // Tz = dz / rT : symmetric window
      if (std::abs(dz) > m_dzdrMax * cT * (1.0 + kPrefilterTol)) {
        return false;
      }
      // cote by segment derived from minCosThetaY and maxKinkAngle (see the .hpp):
      // the top-bottom vector is the sum of the two segments and each segment forms at most `maxKinkAngle` with it
      // with middle pixel and useSubdetectorAngularCuts there is no 3D kink cut, so there is no bound
      if (!pixelMid &&
          std::abs(dy) < m_segCosThetaYMin * L * (1.0 - kPrefilterTol)) {
        return false;
      }
    }
    // candidate structure: pointer to the hit, distances, transverse distance squared,
    // absolute y, unit vector components, and whether it is calibratable
    c = Candidate{&h, dx, dy, cT2, h.y, dx / L, dy / L, dz / L, strip};
    return true;
  };

  // prefilter cuts: same subdetector, not in the same surface, dz/dr, kink/cosThetaY, pT (rho)
  auto pairPrefilter = [&](const Candidate& t, const Candidate& b, bool pixelMid) {
    if (t.calibratable || b.calibratable) {
      return true;  // container position is not yet final, leave the cuts for passCosmicCuts
    }
    if (pixelMid) {
      // cos(dphi) between the xy proyections (middle - top) and (bottom - middle) = -(T . B) / (|T| |B|)
      const double cosDphi = -(t.dx * b.dx + t.dy * b.dy);
      if (cosDphi < (m_cuts.pixelMinCosDeltaPhi - kPrefilterTol) * std::sqrt(t.cT2 * b.cT2)) {
        return false;
      }
    } else if (!m_ptDependent) {
      // cos(kink) between (middle - top) and (bottom - middle) = -(u_t . u_b)
      const double cosKink = -(t.ux * b.ux + t.uy * b.uy + t.uz * b.uz);
      if (cosKink < m_kinkCosThreshold - kPrefilterTol) {
        return false;
      }
    }
    if (m_ptPrefilter) {
      // pT cut: rho = 1 / |k| = pT / (0.3 |Bz|) >= rhoMin
      const double cross = t.dx * b.dy - t.dy * b.dx;
      const double ex = t.dx - b.dx;
      const double ey = t.dy - b.dy;
      const double chord2 = ex * ex + ey * ey;
      if (4.0 * cross * cross * m_rhoMin * m_rhoMin >
          t.cT2 * b.cT2 * chord2 * (1.0 + 1e-6)) {
        return false;
      }
    }
    return true;
  };

  // pool: union of the neighboring bins of the current bin, ordered by y
  std::vector<const Hit*> pool;
  std::vector<Candidate> bottoms;
  std::vector<Candidate> tops;

  // for each z bin with activity, the central point ("middle"/r0) is taken only from that bin,
  // while the candidates for "top" and "bottom" are taken from the union of their neighboring
  // bins (buildFrameWork() / production3Sp: rb/rt are filled from m_rfz_ib/m_rfz_it)
  // this reduces the combinatorial complexity from N^3 to something much smaller when there are many space points per event.
  for (int a = 0; a < kNumZBins; ++a) {
    if (bins[a].empty()) {
      continue;
    }

    pool.clear();
    for (int neighbor : kZNeighbors[a]) {
      for (const Hit& h : bins[neighbor]) {
        pool.push_back(&h);
      }
    }
    if (!m_cuts.requireTopMiddleBottom) {
      // if the top-middle-bottom ordering is not required, the pool is not sorted by y
      // and the candidates are generated in a forced manner
      for (const Hit& middle : bins[a]) {
        for (const Hit* top : pool) {
          if (top->sp.index() == middle.sp.index()) {
            continue;
          }
          if (m_cuts.requireSameSubdetector &&
              ev.isStrip[top->sp.index()] != ev.isStrip[middle.sp.index()]) {
            continue;
          }
          for (const Hit* bottom : pool) {
            if (bottom->sp.index() == middle.sp.index() ||
                bottom->sp.index() == top->sp.index()) {
              continue;
            }
            if (m_cuts.requireSameSubdetector &&
                ev.isStrip[bottom->sp.index()] !=
                    ev.isStrip[middle.sp.index()]) {
              continue;
            }

            // other cuts are symmetric under top<->bottom, each unordered triplet would appear twice.
            if (top->sp.index() < bottom->sp.index()) {
              continue;
            }
            // try to emit a seed from the triplet (top, middle, bottom)
            tryEmit(*top, middle, *bottom);
          }
        }
      }
      continue;
    }

    // neighbors pool sorted by y (once per bin, not once per middle point)
    std::sort(pool.begin(), pool.end(), [](const Hit* l, const Hit* r) {
      return l->y < r->y || (l->y == r->y && l->sp.index() < r->sp.index());
    });
    // angular cuts are applied in makeCandidate() and pairPrefilter(), which are called for each middle point and candidate pair
    for (const Hit& middle : bins[a]) {
      const bool isPixelMid = ev.isStrip != nullptr && ev.isStrip[middle.sp.index()] == 0;
      const bool pixelMid = m_cuts.useSubdetectorAngularCuts && isPixelMid;
      // window in y for the segments, |Δy| in [minDeltaY, maxDeltaY]
      // iff the segment has a recalibratable strip SP, each end can move up to stripShiftY in y:
      // the window is enlarged by that and the exact cut is done by passCosmicCuts
      const double margin = 2.0 * stripShiftY;
      // bottom: y < y_mid
      // top: y > y_mid
      auto lo = std::lower_bound(
          pool.begin(), pool.end(), middle.y,
          [](const Hit* h, double v) { return h->y < v; });
      auto hi = std::upper_bound(
          lo, pool.end(), middle.y,
          [](double v, const Hit* h) { return v < h->y; });

      tops.clear();
      bottoms.clear();
      Candidate c{};
      // generate candidates for the top and bottom space points relative to the middle point
      for (auto it = hi; it != pool.end(); ++it) {
        if (makeCandidate(**it, middle, pixelMid, c)) {
          tops.push_back(c);
        }
      }
      for (auto it = lo; it != pool.begin();) {
        --it;
        if (makeCandidate(**it, middle, pixelMid, c)) {
          bottoms.push_back(c);
        }
      }
      if (tops.empty() || bottoms.empty()) {
        continue;
      }
      
      // for each bottom candidate, find the top candidates that are within the y window defined by minDeltaY and maxDeltaY
      // taking into account the margin for recalibratable strip SPs
      for (const Candidate& b : bottoms) {
        const double dyB = middle.y - b.y;
        const double slackB = b.calibratable ? margin : 0.0;
        if (dyB < m_cuts.minDeltaY - slackB - 1e-6 || dyB > m_cuts.maxDeltaY + slackB + 1e-6) {
          continue;
        }
        // binary search for the first top candidate whose y is greater than middle.y + minDeltaY - margin
        auto first = std::lower_bound(
            tops.begin(), tops.end(), middle.y + m_cuts.minDeltaY - margin - 1e-6,
            [](const Candidate& t, double v) { return t.y < v; });
        auto last = std::upper_bound(
            first, tops.end(), middle.y + m_cuts.maxDeltaY + margin + 1e-6,
            [](double v, const Candidate& t) { return v < t.y; });
        for (auto it = first; it != last; ++it) {
          const double dyT = it->y - middle.y;
          if (!it->calibratable && (dyT < m_cuts.minDeltaY - 1e-6 ||
                                    dyT > m_cuts.maxDeltaY + 1e-6)) {
            continue;
          }
          if (!pairPrefilter(*it, b, pixelMid)) {
            continue;
          }
          // try to emit a seed from the triplet (top, middle, bottom)
          tryEmit(*it->hit, middle, *b.hit);
        }
      }
    }
  }
}

// createTwoSpacePointSeeds() generates doublet seeds from the given space points, applying the configured cuts and storing the results in outputSeeds.
// its optional recalibration of strip space points is not performed, as doublets do not have a third point to determine the direction
void CosmicTripletSeedFinder::createTwoSpacePointSeeds(
    const SpacePointContainer& spacePoints, SeedContainer& outputSeeds,
    const std::vector<CosmicSpacePointClass>* spClasses,
    const std::vector<std::uint8_t>* skipSpacePoints,
    const std::vector<GeometryIdentifier>* spSurfaces) const {
  if (!outputSeeds.hasSpacePointContainer()) {
    outputSeeds.assignSpacePointContainer(spacePoints);
  }
  if (skipSpacePoints != nullptr &&
      skipSpacePoints->size() != spacePoints.size()) {
    throw std::invalid_argument(
        "CosmicTripletSeedFinder: skipSpacePoints needs one entry per space "
        "point");
  }
  outputSeeds.reserve(outputSeeds.size() + spacePoints.size(), 2.0f);

  // variances only enter in dT, which a doublet does not have:
  // only the pixel/strip type (type restriction) and strip details are needed.
  EventBuffers buffers;
  const EventData ev = prepareEvent(spacePoints, spClasses, spSurfaces, false,
                                    m_cuts.requireSameSubdetector, buffers);
  auto calibratable = [&ev](std::size_t index) {
    return ev.isCalibratable != nullptr && ev.isCalibratable[index] != 0;
  };

  // the margin is twice the maximum shift in y for recalibratable strip space points
  const double margin = 2.0 * buffers.stripShiftY;

  const ZBins bins = fillZBins(spacePoints, m_cuts.zMin, m_cuts.zMax, skipSpacePoints);

  std::vector<const Hit*> pool;
  for (int a = 0; a < kNumZBins; ++a) {
    if (bins[a].empty()) {
      continue;
    }

    // the same neighbors as production3Sp, ordered by y
    pool.clear();
    for (int neighbor : kZNeighbors[a]) {
      for (const Hit& h : bins[neighbor]) {
        pool.push_back(&h);
      }
    }
    std::sort(pool.begin(), pool.end(), [](const Hit* l, const Hit* r) {
      return l->y < r->y || (l->y == r->y && l->sp.index() < r->sp.index());
    });

    // top is used as the central point: each pair is considered only once, from the top bin
    for (const Hit& top : bins[a]) {
      auto first = std::lower_bound(
          pool.begin(), pool.end(), top.y - m_cuts.maxDeltaY - margin - 1e-6,
          [](const Hit* h, double v) { return h->y < v; });
      auto last = std::upper_bound(
          first, pool.end(), top.y - m_cuts.minDeltaY + margin + 1e-6,
          [](double v, const Hit* h) { return v < h->y; });
      for (auto it = first; it != last; ++it) {
        const Hit& bottom = **it;
        if (!(bottom.y < top.y)) {
          continue;
        }
        if (m_cuts.requireSameSubdetector &&
            ev.isStrip[bottom.sp.index()] != ev.isStrip[top.sp.index()]) {
          continue;
        }
        if (ev.surfaces != nullptr &&
            ev.surfaces[bottom.sp.index()] == ev.surfaces[top.sp.index()]) {
          continue;
        }
        Point pT{top.x, top.y, top.z, top.sp.index()};
        Point pB{bottom.x, bottom.y, bottom.z, bottom.sp.index()};
        if ((calibratable(pT.index) || calibratable(pB.index)) &&
            !calibrateStripSpacePoints(pT, nullptr, pB, ev)) {
          continue;
        }
        if (!passDoubletCuts(pT, pB)) {
          continue;
        }
        std::array<SpacePointIndex, 2> seedSpacePoints{bottom.sp.index(), top.sp.index()};
        MutableSeedProxy seed = outputSeeds.createSeed();
        seed.assignSpacePointIndices(seedSpacePoints);
        // without a third point there is no dT, so all doublets have the same quality 
        seed.quality() = m_cuts.seedQuality;
        seed.vertexZ() = 0.0f;
      }
    }
  }
}
// bottom-middle segment cuts: the same as production3Sp, but without the kink cut (no third point to determine the direction)
bool CosmicTripletSeedFinder::passDoubletCuts(const Point& top,
                                              const Point& bottom) const {
  if (!std::isfinite(top.x) || !std::isfinite(top.y) ||
      !std::isfinite(top.z) || !std::isfinite(bottom.x) ||
      !std::isfinite(bottom.y) || !std::isfinite(bottom.z)) {
    return false;
  }
  if (top.z < m_cuts.zMin || top.z > m_cuts.zMax || bottom.z < m_cuts.zMin ||
      bottom.z > m_cuts.zMax) {
    return false;
  }
  const double dx = top.x - bottom.x;
  const double dy = top.y - bottom.y;
  const double dz = top.z - bottom.z;
  // Δy in [minDeltaY, maxDeltaY] (mindRadius/maxdRadius)
  if (dy < m_cuts.minDeltaY || dy > m_cuts.maxDeltaY) {
    return false;
  }
  // window in Tz = dz / rT
  const double rT = std::hypot(dx, dy);
  if (!(rT > 0.0) || std::abs(dz) > m_dzdrMax * rT) {
    return false;
  }
  // y aligment
  const double distance = std::sqrt(rT * rT + dz * dz);
  return dy / distance >= m_cuts.minCosThetaY;
}

// passCosmicCuts applies vertical order, kink angle, vertical separation/alignment,
// residual relative to the top-bottom line, Z envelope, and finally the
// curvature/longitudinal cut (which also calculates the quality).
bool CosmicTripletSeedFinder::passCosmicCuts(const Point& top,
                                             const Point& middle,
                                             const Point& bottom,
                                             float& quality,
                                             const EventData& ev) const {

  if (!std::isfinite(top.x) || !std::isfinite(top.y) ||
      !std::isfinite(top.z) || !std::isfinite(middle.x) ||
      !std::isfinite(middle.y) || !std::isfinite(middle.z) ||
      !std::isfinite(bottom.x) || !std::isfinite(bottom.y) ||
      !std::isfinite(bottom.z)) {
    return false;
  }

  // all the points must be of the same type (pixel or strip) if requireSameSubdetector is enabled
  if (m_cuts.requireSameSubdetector) {
    const std::uint8_t kind = ev.isStrip[middle.index];
    if (ev.isStrip[top.index] != kind || ev.isStrip[bottom.index] != kind) {
      return false;
    }
  }

 // vertical order: top > middle > bottom in y
  if (m_cuts.requireTopMiddleBottom && !(top.y > middle.y && middle.y > bottom.y)) {
    return false;
  }

  if (top.z < m_cuts.zMin || top.z > m_cuts.zMax || middle.z < m_cuts.zMin ||
      middle.z > m_cuts.zMax || bottom.z < m_cuts.zMin ||
      bottom.z > m_cuts.zMax) {
    return false;
  }

  const double dx1 = middle.x - top.x;
  const double dy1 = middle.y - top.y;
  const double dz1 = middle.z - top.z;

  const double dx2 = bottom.x - middle.x;
  const double dy2 = bottom.y - middle.y;
  const double dz2 = bottom.z - middle.z;

  const double n1 = std::sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
  const double n2 = std::sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
  if (n1 == 0.0 || n2 == 0.0) {
    return false;
  }

  const double dot = dx1 * dx2 + dy1 * dy2 + dz1 * dz2;
  const double cosAngle = std::clamp(dot / (n1 * n2), -1.0, 1.0);
  const double angle = std::acos(cosAngle);

  // Kink and residual
  double kinkValue = angle;
  double kinkLimit = m_cuts.maxKinkAngle;
  double residualLimit = m_cuts.maxResidual;
  double parallelResidual = 0.0;
  if (m_ptDependent) {
    PtDependentTerms terms;
    if (!ptDependentTerms(top, middle, bottom, terms)) {
      return false;
    }
    kinkValue = terms.polarKink;
    kinkLimit += terms.kinkSlack;
    residualLimit += terms.residualSlack;
    parallelResidual = terms.parallelResidual;
  }
  if (m_cuts.useSubdetectorAngularCuts && ev.isStrip[middle.index] == 0) {
    // dphi = phi2 - phi1, with phi_i = atan2(dy_i, dx_i)
    const double c1 = std::hypot(dx1, dy1);
    const double c2 = std::hypot(dx2, dy2);
    if (c1 == 0.0 || c2 == 0.0) {
      return false;
    }
    const double cosDphi = (dx1 * dx2 + dy1 * dy2) / (c1 * c2);
    if (!(cosDphi > m_cuts.pixelMinCosDeltaPhi)) {
      return false;
    }
  } else if (kinkValue > kinkLimit) {
    return false;
  }

  const double dx = bottom.x - top.x;
  const double dy = bottom.y - top.y;
  const double dz = bottom.z - top.z;
  const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
  if (distance == 0.0) {
    return false;
  }

  // Δy per segment in [minDeltaY, maxDeltaY]
  const double dyTop = std::abs(top.y - middle.y);
  const double dyBottom = std::abs(middle.y - bottom.y);
  if (dyTop < m_cuts.minDeltaY || dyTop > m_cuts.maxDeltaY || dyBottom < m_cuts.minDeltaY || dyBottom > m_cuts.maxDeltaY) {
    return false;
  }

  const double verticalSpan = std::abs(dy);
  const double cosThetaY = verticalSpan / distance;
  if (cosThetaY < m_cuts.minCosThetaY) {
    return false;
  }

  if (m_cuts.maxResidual >= 0.0) {
    const double residual = m_ptDependent
                                ? parallelResidual
                                : residualToTopBottomLine(top, middle, bottom);
    if (residual > residualLimit) {
      return false;
    }
  }
  // return true if all cuts are passed, and update the quality parameter with the calculated quality
  return passCurvatureAndQuality(top, middle, bottom, quality, ev);
}

// pT-dependent terms: curvature, kink, residual, and scattering tolerance
// Notation: d1 = middle - top, d2 = bottom - middle
// c_i = transverse length (xy plane) of d_i, L_i = 3D length, chord = bottom - top.



bool CosmicTripletSeedFinder::ptDependentTerms(const Point& top,
                                               const Point& middle,
                                               const Point& bottom,
                                               PtDependentTerms& out) const {
  const double x1 = middle.x - top.x;
  const double y1 = middle.y - top.y;
  const double z1 = middle.z - top.z;
  const double x2 = bottom.x - middle.x;
  const double y2 = bottom.y - middle.y;
  const double z2 = bottom.z - middle.z;

  const double c1 = std::hypot(x1, y1);
  const double c2 = std::hypot(x2, y2);
  const double cc = std::hypot(x1 + x2, y1 + y2);
  const double L1 = std::sqrt(x1 * x1 + y1 * y1 + z1 * z1);
  const double L2 = std::sqrt(x2 * x2 + y2 * y2 + z2 * z2);
  const double Lc = std::sqrt((x1 + x2) * (x1 + x2) + (y1 + y2) * (y1 + y2) +
                              (z1 + z2) * (z1 + z2));
  if (c1 == 0.0 || c2 == 0.0 || cc == 0.0 || Lc == 0.0) {
    return false;
  }

  // pT estimated: the circunference for the three points in the xy plane has radius rho = c1 c2 c_chord / (2 |d1 x d2|_z) (Menger formula),
  // and pT[MeV] = 0.3 Bz[T] rho[mm]
  // the cut pT >= 0.9 minPt is applied
  const double cross = x1 * y2 - y1 * x2;
  const double invRho = 2.0 * std::abs(cross) / (c1 * c2 * cc);       // 1/mm
  const double invPt = invRho / (0.3 * std::abs(m_cuts.magneticFieldBz));  // 1/MeV
  if (invPt > 1.0 / (0.9 * m_cuts.minPt)) {
    return false;
  }
  const double invPtEff = std::max(invPt, 1.0 / m_cuts.maxPtForScattering);

  // 4) tolerancia de scattering múltiple
  const double sinThetaChord = cc / Lc;
  const double thetaMS =
      std::sqrt(m_cuts.multipleScatteringCof) * sinThetaChord * invPtEff;
  out.kinkSlack = thetaMS;
  out.residualSlack = thetaMS * L1 * L2 / (L1 + L2);

  // the transverse turn is already bounded by the pT cut
  // the kink compare the polar deviation |theta1 - theta2|, with theta_i = atan2(c_i, dz_i)
  out.polarKink = std::abs(std::atan2(c1, z1) - std::atan2(c2, z2));
  // the residual is the distance from the middle point to the line connecting top and bottom
  // which can be decomposed into a transverse component (h_T) and a parallel component (h_par)
  // residual out of transversal plane: h3D^2 = h_T^2 + h_par^2
  const double hT = std::abs(cross) / cc;
  const double h3 = residualToTopBottomLine(top, middle, bottom);
  out.parallelResidual = std::sqrt(std::max(0.0, h3 * h3 - hT * hT));
  return true;
}

// built the local frame (U, V) centered in "middle", calculates the slopes Tz of each segment,
// applies the pseudorapidity window, the longitudinal discriminant dT
// and (depending on B) the momentum cut by curvature and its scattering correction.
bool CosmicTripletSeedFinder::passCurvatureAndQuality(const Point& top, const Point& middle,
                                                      const Point& bottom, float& quality,
                                                      const EventData& ev) const {
  const double X = middle.x;
  const double Y = middle.y;
  const double dirX = top.x - X;
  const double dirY = top.y - Y;
  const double R = std::sqrt(dirX * dirX + dirY * dirY);
  if (R == 0.0) {
    return false;
  }
  const double ax = dirX / R;
  const double ay = dirY / R;

  struct UV {
    double U = 0.0, V = 0.0;
    double invDist = 0.0;  // inverse of the distance from the middle point to the other point
    double dz = 0.0;
  };

  auto computeUV = [&](const Point& sp, UV& out) -> bool {
    const double dx = sp.x - X;
    const double dy = sp.y - Y;
    const double x = dx * ax + dy * ay;
    const double y = -dx * ay + dy * ax;
    const double r2 = x * x + y * y;
    if (r2 == 0.0) {
      return false;
    }
    const double invR2 = 1.0 / r2;
    out.U = x * invR2;
    out.V = y * invR2;
    out.invDist = std::sqrt(invR2);
    out.dz = sp.z - middle.z;
    return true;
  };

  UV b, t;
  if (!computeUV(bottom, b) || !computeUV(top, t)) {
    return false;
  }

  // pseudorapidity window for the segments, symmetric in Tz = dz / rT
  const double TzB = -b.dz * b.invDist;
  const double TzT = t.dz * t.invDist;
  if (TzB < -m_dzdrMax || TzB > m_dzdrMax || TzT < -m_dzdrMax ||
      TzT > m_dzdrMax) {
    return false;
  }

  // Varianzas: las del container, o las efectivas (con escala de clase) si
  // useCovarianceScale está activo.
  // container variances are used for doublets and effective variances for triplets
  auto covRof = [&ev](const Point& sp) {
    return ev.covR != nullptr
               ? ev.covR[sp.index]
               : varianceR(ev.spacePoints->at(
                     static_cast<SpacePointIndex>(sp.index)));
  };
  auto covZof = [&ev](const Point& sp) {
    return ev.covZ != nullptr
               ? ev.covZ[sp.index]
               : varianceZ(ev.spacePoints->at(
                     static_cast<SpacePointIndex>(sp.index)));
  };

  const double covr0 = covRof(middle);
  const double covz0 = covZof(middle);

  const double ErB = (covz0 + covZof(bottom) +
                       TzB * TzB * (covr0 + covRof(bottom))) *
                      (b.invDist * b.invDist);
  const double ErT = (covz0 + covZof(top) +
                       TzT * TzT * (covr0 + covRof(top))) *
                      (t.invDist * t.invDist);

  const double Ts = 0.5 * (TzB + TzT);
  const double dTz = TzB - TzT;

  // longitudinal discriminant dT, which is negative for a compatible triplet
  double dT = dTz * dTz - ErB - ErT - 2.0 * b.invDist * t.invDist * (Ts * Ts * covr0 + covz0);

  // momentum cut by curvature and its scattering correction
  if (m_cuts.minPt > 0.0) {
    const bool hasField = m_cuts.usePtCut && m_cuts.magneticFieldBz != 0.0;
    if (hasField) {
      const double dU = t.U - b.U;
      if (dU == 0.0) {
        return false;
      }
      const double A = (t.V - b.V) / dU;
      const double Bc = t.V - A * t.U;
      const double S2 = 1.0 + A * A;
      const double S = std::sqrt(S2);
      const double K = 2.0 / (0.3 * m_cuts.magneticFieldBz);
      const double ipt = 1.0 / std::abs(0.9 * m_cuts.minPt);
      const double BK = std::abs(Bc * K);

      // momentum cut by curvature
      if (BK > ipt * S) {
        return false;
      }

      // scattering correction: dT -= (BK^2) * multipleScatteringCof * (1 + TzB^2) / S^2
      const double SA = 1.0 + TzB * TzB;
      dT -= (BK * BK) * m_cuts.multipleScatteringCof * SA / S2;
    } else {
      // if there is no magnetic field, the scattering correction is applied with the effective pT cut
      const double ipt = 1.0 / std::abs(0.9 * m_cuts.minPt);
      const double SA = 1.0 + TzB * TzB;
      dT -= (ipt * ipt) * m_cuts.multipleScatteringCof * SA;
    }
  }

  if (!std::isfinite(dT)) {
    return false;
  }

  // if dT > 0, the triplet is incompatible with a cosmic track
  // if useLongitudinalCompatibility is true, the triplet is rejected
  if (m_cuts.useLongitudinalCompatibility && dT > 0.0) {
    return false;
  }

  // if useSeedQuality is true, the quality is calculated as dT * invDistTop * invDistBottom / (invDistTop + invDistBottom)
  // else, the quality is left as the fixed value m_cuts.seedQuality
  if (m_cuts.useSeedQuality) {
    const double denom = b.invDist + t.invDist;
    quality = denom > 0.0
                  ? static_cast<float>(dT * b.invDist * t.invDist / denom)
                  : static_cast<float>(dT);
  }
  return true;
}

// trace direction in bottom, middle and top (dirs[0..2])
// without magnetic field, the direction is the bottom-top chord
// with magnetic field, each SP uses the tangent at its own point:
// in xy, perpendicular to the radius of the circle through the three points
// in z, the slope of the chord
bool CosmicTripletSeedFinder::tripletDirections(const Point& top,
                                                const Point& middle,
                                                const Point& bottom,
                                                std::array<std::array<float, 3>, 3>& dirs) const {
  const double cx = top.x - bottom.x;
  const double cy = top.y - bottom.y;
  const double cz = top.z - bottom.z;
  const double chord = std::sqrt(cx * cx + cy * cy + cz * cz);
  if (!(chord > 0.0)) {
    return false;
  }
  const std::array<float, 3> straight = {static_cast<float>(cx / chord),
                                         static_cast<float>(cy / chord),
                                         static_cast<float>(cz / chord)};
  dirs = {straight, straight, straight};
  const double chordXY = std::hypot(cx, cy);
  if (m_cuts.magneticFieldBz == 0.0 || chordXY == 0.0) {
    return true;
  }

  // circunderence through the three points in xy
  // 2 u.a = |a|^2 and 2 u.q = |q|^2 for a = bottom - middle, q = top - middle
  const double ax = bottom.x - middle.x;
  const double ay = bottom.y - middle.y;
  const double qx = top.x - middle.x;
  const double qy = top.y - middle.y;
  const double det = ax * qy - ay * qx;
  const double a2 = ax * ax + ay * ay;
  const double q2 = qx * qx + qy * qy;
  if (det == 0.0) {
    return true;
  }
  const double ux = 0.5 * (a2 * qy - q2 * ay) / det;
  const double uy = 0.5 * (ax * q2 - qx * a2) / det;
  constexpr double kStraightRadius = 1e8;  // mm
  if (!(std::hypot(ux, uy) < kStraightRadius)) {
    return true;
  }
  const double slope = cz / chordXY;  // dz / ds_xy
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
    const double tn = std::hypot(tx, ty);
    if (!(tn > 0.0)) {
      return false;
    }
    const double dx = tx / tn;
    const double dy = ty / tn;
    const double n = std::sqrt(1.0 + slope * slope);
    dirs[i] = {static_cast<float>(dx / n),
                static_cast<float>(dy / n),
                static_cast<float>(slope / n)};
  }
  return true;
}

// recalculates the SPs of strips of the seed with the direction of the track
// with the direction of the seed the crossing is calculated, with the new crossings a new direction, and
// so on until no point moves more than kCalibrationConvergence
// then, with the final direction, the crossing must fall within the two strips with stripCalibrationTolerance
// (if not, the seed is not a track)

bool CosmicTripletSeedFinder::calibrateStripSpacePoints(Point& top,
                                                        Point* middle,
                                                        Point& bottom,
                                                        const EventData& ev) const {
  const std::array<Point*, 3> points = {&bottom, middle, &top};
  std::array<bool, 3> strip{};
  for (std::size_t i = 0; i < 3; ++i) {
    strip[i] = points[i] != nullptr && ev.isCalibratable != nullptr &&
               ev.isCalibratable[points[i]->index] != 0;
  }
  // direction in bottom, middle and top: that of the triplet or, in a doublet, the chord
  auto directions = [&](std::array<std::array<float, 3>, 3>& dirs) {
    if (middle != nullptr) {
      return tripletDirections(top, *middle, bottom, dirs);
    }
    const double cx = top.x - bottom.x;
    const double cy = top.y - bottom.y;
    const double cz = top.z - bottom.z;
    const double chord = std::sqrt(cx * cx + cy * cy + cz * cz);
    if (!(chord > 0.0)) {
      return false;
    }
    const std::array<float, 3> d = {static_cast<float>(cx / chord),
                                    static_cast<float>(cy / chord),
                                    static_cast<float>(cz / chord)};
    dirs = {d, d, d};
    return true;
  };
  // note: a direction (almost) parallel to the module plane crosses the strips at infinity, so its crossing is not valid
  auto crossesModule = [&ev](const std::array<float, 3>& d,
                             std::size_t index) {
    const auto& n = ev.stripDetails[index].innerCrossOuterHalfVector;
    const double scale = static_cast<double>(d[0]) * n[0] +
                         static_cast<double>(d[1]) * n[1] +
                         static_cast<double>(d[2]) * n[2];
    const double norm = std::sqrt(static_cast<double>(n[0]) * n[0] +
                                  static_cast<double>(n[1]) * n[1] +
                                  static_cast<double>(n[2]) * n[2]);
    return std::abs(scale) > kMinCosToModuleNormal * norm;
  };
  constexpr float kNoLimit = std::numeric_limits<float>::infinity();
  std::array<std::array<float, 3>, 3> dirs{};
  bool converged = false;
  for (int iteration = 0; iteration < kMaxCalibrationIterations;
       ++iteration) {
    if (!directions(dirs)) {
      return false;
    }
    double maxShift = 0.0;
    for (std::size_t i = 0; i < 3; ++i) {
      if (!strip[i]) {
        continue;
      }
      if (!crossesModule(dirs[i], points[i]->index)) {
        return false;
      }
      std::array<float, 3> r{};
      detail::calibrateOuterStripSpacePoint(
          dirs[i], ev.stripDetails[points[i]->index], r, kNoLimit);
      if (!std::isfinite(r[0]) || !std::isfinite(r[1]) ||
          !std::isfinite(r[2])) {
        return false;
      }
      const double sx = r[0] - points[i]->x;
      const double sy = r[1] - points[i]->y;
      const double sz = r[2] - points[i]->z;
      maxShift = std::max(maxShift, std::sqrt(sx * sx + sy * sy + sz * sz));
      points[i]->x = r[0];
      points[i]->y = r[1];
      points[i]->z = r[2];
    }
    // if the maximum shift is less than the convergence threshold, its considered that the calibration converged
    if (maxShift < kCalibrationConvergence) {
      converged = true;
      break;
    }
  }
  if (!converged || !directions(dirs)) {
    return false;
  }
  // tolerance for the crossing of the strips with the final direction
  // if the crossing is outside the two strips by more than this tolerance, the seed is rejected
  const auto tolerance = static_cast<float>(m_cuts.stripCalibrationTolerance);
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

// calculates the distance from the middle point to the line connecting top and bottom
double CosmicTripletSeedFinder::residualToTopBottomLine(const Point& top,
                                                       const Point& middle,
                                                       const Point& bottom) {
  const double ax = middle.x - top.x;
  const double ay = middle.y - top.y;
  const double az = middle.z - top.z;

  const double bx = bottom.x - top.x;
  const double by = bottom.y - top.y;
  const double bz = bottom.z - top.z;

  const double bNorm = std::sqrt(bx * bx + by * by + bz * bz);
  if (bNorm == 0.0) {
    return std::numeric_limits<double>::infinity();
  }

  const double cx = ay * bz - az * by;
  const double cy = az * bx - ax * bz;
  const double cz = ax * by - ay * bx;

  return std::sqrt(cx * cx + cy * cy + cz * cz) / bNorm;
}

// estimates the vertex z position from the three space points, but for cosmic rays,
// there is no primary vertex, therefore the reference implementation always records zVertex = 0 instead of estimating it from the SPs.
float CosmicTripletSeedFinder::estimateVertexZ(const Point& top,
                                               const Point& middle,
                                               const Point& bottom) {
  (void)top;
  (void)middle;
  (void)bottom;
  return 0.0f;
}

}
