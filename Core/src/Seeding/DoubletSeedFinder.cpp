// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/DoubletSeedFinder.hpp"

#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Seeding/detail/BlockPass.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <type_traits>

#include <boost/mp11.hpp>
#include <boost/mp11/algorithm.hpp>

namespace Acts {

namespace {

/// How the cuts on deltaZ and on the collision region are made, chosen in
/// create(). Where a range is symmetric about zero, it is tested with one
/// comparison on the absolute value; an unbounded deltaZ range rejects no
/// candidate, so its cut is not made.
enum class ZCuts {
  /// two comparisons each
  eGeneral,
  /// symmetric deltaZ range and collision region
  eSymmetric,
  /// symmetric collision region, unbounded deltaZ range
  eSymmetricNoDeltaZ,
};

template <bool isBottomCandidate, bool interactionPointCut, bool sortedByR,
          bool experimentCuts, bool useTime, ZCuts zCuts>
class Impl final : public DoubletSeedFinder {
 public:
  explicit Impl(const DerivedConfig& config) : m_cfg(config) {}

  const DerivedConfig& config() const override { return m_cfg; }

  /// Iterates over dublets and tests the compatibility by applying a series of
  /// cuts that can be tested with only two SPs.
  ///
  /// @param config Doublet cuts that define the compatibility of space points
  /// @param middleSp Space point candidate to be used as middle SP in a seed
  /// @param middleSpInfo Information about the middle space point
  /// @param candidateSps Range or subet of space points to be used as candidates
  ///   for middle SP in a seed
  /// @param compatibleDoublets Output container for compatible doublets
  template <typename CandidateSps>
  void createDoubletsImpl(const ConstSpacePointProxy& middleSp,
                          const MiddleSpInfo& middleSpInfo,
                          CandidateSps candidateSps,
                          DoubletsForMiddleSp& compatibleDoublets) const {
    const float impactMax =
        isBottomCandidate ? -m_cfg.impactMax : m_cfg.impactMax;

    const float xM = middleSp.xy()[0];
    const float yM = middleSp.xy()[1];
    const float zM = middleSp.zr()[0];
    const float rM = middleSp.zr()[1];
    const float varianceZM = middleSp.varianceZ();
    const float varianceRM = middleSp.varianceR();

    // time of the middle space point and its variance. only filled when the
    // time cut is enabled
    [[maybe_unused]] float tM = 0;
    [[maybe_unused]] float varianceTM = 0;
    if constexpr (useTime) {
      tM = middleSp.time();
      varianceTM = middleSp.varianceT();
    }

    // equivalent to impactMax / (rM * rM);
    const float vIPAbs = impactMax * middleSpInfo.uIP2;

    const auto outsideRangeCheck = [](const float value, const float min,
                                      const float max) {
      // intentionally using `|` after profiling. faster due to better branch
      // prediction
      return static_cast<bool>(static_cast<int>(value < min) |
                               static_cast<int>(value > max));
    };
    // Same result as outsideRangeCheck(value, -max, max) for every input, but
    // one comparison instead of two.
    const auto outsideSymmetricRange = [](const float value, const float max) {
      return std::abs(value) > max;
    };

    const auto calculateError = [&](float varianceZO, float varianceRO,
                                    float iDeltaR2, float cotTheta) {
      return iDeltaR2 * ((varianceZM + varianceZO) +
                         (cotTheta * cotTheta) * (varianceRM + varianceRO));
    };

    // The cut values the cuts on z and r read.
    const float deltaRMin = m_cfg.deltaRMin;
    const float deltaRMax = m_cfg.deltaRMax;
    [[maybe_unused]] const float deltaZMin = m_cfg.deltaZMin;
    [[maybe_unused]] const float deltaZMax = m_cfg.deltaZMax;
    [[maybe_unused]] const float collisionRegionMin = m_cfg.collisionRegionMin;
    const float collisionRegionMax = m_cfg.collisionRegionMax;

    // The cuts that need only the candidate's z and r: deltaR, deltaZ and the
    // collision region.
    const auto passesZAndRCuts = [=](float zO, float rO) {
      const float deltaR = isBottomCandidate ? rM - rO : rO - rM;
      const float deltaZ = isBottomCandidate ? zM - zO : zO - zM;
      bool outside = false;
      if constexpr (!sortedByR) {
        outside = outsideRangeCheck(deltaR, deltaRMin, deltaRMax);
      }
      if constexpr (zCuts == ZCuts::eGeneral) {
        outside |= outsideRangeCheck(deltaZ, deltaZMin, deltaZMax);
      } else if constexpr (zCuts == ZCuts::eSymmetric) {
        // deltaZMin is -deltaZMax
        outside |= outsideSymmetricRange(deltaZ, deltaZMax);
      }
      // the longitudinal impact parameter zOrigin is defined as (zM - rM *
      // cotTheta) where cotTheta is the ratio Z/R (forward angle) of space
      // point duplet but instead we calculate (zOrigin * deltaR) and multiply
      // collisionRegion by deltaR to avoid divisions
      const float zOriginTimesDeltaR = zM * deltaR - rM * deltaZ;
      // check if duplet origin on z axis within collision region
      if constexpr (zCuts == ZCuts::eGeneral) {
        outside |=
            outsideRangeCheck(zOriginTimesDeltaR, collisionRegionMin * deltaR,
                              collisionRegionMax * deltaR);
      } else {
        // collisionRegionMin * deltaR is exactly -(collisionRegionMax *
        // deltaR), since negation is exact
        outside |= outsideSymmetricRange(zOriginTimesDeltaR,
                                         collisionRegionMax * deltaR);
      }
      return static_cast<std::uint32_t>(!outside);
    };

    const SpacePointContainer& container = candidateSps.container();
    const auto xy = container.xyColumn().data();
    const auto zr = container.zrColumn().data();
    const auto varianceZ = container.varianceZColumn().data();
    const auto varianceR = container.varianceRColumn().data();

    // The candidates as indices into the container: a range of them is
    // contiguous, a subset lists them.
    constexpr bool kContiguous =
        std::is_same_v<CandidateSps, SpacePointContainer::ConstRange>;
    const auto candidateIndex = [&](std::uint32_t k) -> SpacePointIndex {
      if constexpr (kContiguous) {
        return candidateSps.range().first + k;
      } else {
        return candidateSps.subset()[k];
      }
    };

    // Sorted by radius, the candidates inside deltaR end at the first one
    // outside it.
    std::uint32_t end = static_cast<std::uint32_t>(candidateSps.size());
    if constexpr (sortedByR) {
      assert(std::ranges::is_sorted(
          std::views::iota(std::uint32_t{0}, end), {},
          [&](std::uint32_t k) { return zr[candidateIndex(k)][1]; }));
      end = *std::ranges::partition_point(
          std::views::iota(std::uint32_t{0}, end), [&](std::uint32_t k) {
            const float rO = zr[candidateIndex(k)][1];
            return isBottomCandidate ? !(rM - rO < deltaRMin)
                                     : !(rO - rM > deltaRMax);
          });
    }

    // The candidate in the u-v reference frame of the middle space point.
    struct UvFrame {
      float xNewFrame;
      float yNewFrame;
      float uT;
      float vT;
      float iDeltaR2;
    };
    const auto transform = [&](SpacePointIndex indexO) {
      const float deltaX = xy[indexO][0] - xM;
      const float deltaY = xy[indexO][1] - yM;

      const float xNewFrame =
          deltaX * middleSpInfo.cosPhiM + deltaY * middleSpInfo.sinPhiM;
      const float yNewFrame =
          deltaY * middleSpInfo.cosPhiM - deltaX * middleSpInfo.sinPhiM;

      const float deltaR2 = deltaX * deltaX + deltaY * deltaY;
      const float iDeltaR2 = 1 / deltaR2;

      return UvFrame{xNewFrame, yNewFrame, xNewFrame * iDeltaR2,
                     yNewFrame * iDeltaR2, iDeltaR2};
    };

    // The experiment cuts, and the doublet if it passes them, for a candidate
    // that passed all other cuts.
    const auto addDoublet = [&](SpacePointIndex indexO, float deltaZ,
                                const UvFrame& uv) {
      const float iDeltaR = std::sqrt(uv.iDeltaR2);
      const float cotTheta = deltaZ * iDeltaR;

      // discard doublets based on experiment specific cuts
      if constexpr (experimentCuts) {
        if (!m_cfg.experimentCuts(middleSp, container[indexO], cotTheta,
                                  isBottomCandidate)) {
          return;
        }
      }

      const float er = calculateError(varianceZ[indexO], varianceR[indexO],
                                      uv.iDeltaR2, cotTheta);

      // fill output vectors
      compatibleDoublets.emplace_back(indexO, cotTheta, iDeltaR, er, uv.uT,
                                      uv.vT, uv.xNewFrame, uv.yNewFrame);
    };

    // The curvature cut, for a candidate in the u-v reference frame.
    [[maybe_unused]] const auto failsCurvatureCut = [&](float uT, float vT,
                                                        float yNewFrame) {
      // in the rotated frame the interaction point is positioned at x = -rM
      // and y ~= impactParam
      const float vIP = (yNewFrame > 0) ? -vIPAbs : vIPAbs;

      // we can obtain aCoef as the slope dv/du of the linear function,
      // estimated using du and dv between the two SP bCoef is obtained by
      // inserting aCoef into the linear equation
      const float aCoef = (vT - vIP) / (uT - middleSpInfo.uIP);
      const float bCoef = vIP - aCoef * middleSpInfo.uIP;
      // the distance of the straight line from the origin (radius of the
      // circle) is related to aCoef and bCoef by d^2 = bCoef^2 / (1 +
      // aCoef^2) = 1 / (radius^2) and we can apply the cut on the curvature
      return (bCoef * bCoef) * m_cfg.minHelixDiameter2 > 1 + aCoef * aCoef;
    };

    // The candidates are processed in blocks. For each block, the cuts on z and
    // r are evaluated for every candidate and the candidates that pass are
    // collected. With the interaction point cut, the remaining cuts are applied
    // to the survivors of the block in the same way; without it, the survivors
    // are completed one at a time.
    constexpr std::uint32_t kBlockSize = 256;
    // Not cleared: each entry is written before it is read.
    // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init)
    std::array<std::uint32_t, kBlockSize> passedZAndR;
    std::array<SpacePointIndex, kBlockSize> survivors;
    // NOLINTEND(cppcoreguidelines-pro-type-member-init)
    for (std::uint32_t block = 0; block < end; block += kBlockSize) {
      const std::uint32_t n = std::min(kBlockSize, end - block);

      // the cuts on z and r
      std::uint32_t nPassed = 0;
      if constexpr (kContiguous) {
        // the z and r of the block's candidates, which are contiguous
        const std::array<float, 2>* const zrBlock =
            zr.data() + candidateIndex(block);
        nPassed = detail::evaluateCut<8>(
            n,
            [&](std::size_t i) {
              return passesZAndRCuts(zrBlock[i][0], zrBlock[i][1]);
            },
            passedZAndR.data());
      } else {
        // the candidates of a subset, each looked up by its index
        nPassed = detail::evaluateCut<1>(
            n,
            [&](std::size_t i) {
              const std::array<float, 2>& zrO =
                  zr[candidateIndex(block + static_cast<std::uint32_t>(i))];
              return passesZAndRCuts(zrO[0], zrO[1]);
            },
            passedZAndR.data());
      }
      if (nPassed == 0) {
        continue;
      }
      std::uint32_t nSurvivors = detail::collectFlagged(
          n, passedZAndR.data(),
          [&](std::uint32_t i) { return candidateIndex(block + i); },
          survivors.data());

      // check the time compatibility of the two space points. the time
      // difference is corrected for the time of flight along the straight line
      // between them, assuming an outgoing particle at the speed of light
      // (which is 1 in ACTS units). placed after the cheap range checks above
      // to avoid the square root for candidates that are rejected anyway
      if constexpr (useTime) {
        std::uint32_t nInTime = 0;
        for (std::uint32_t k = 0; k < nSurvivors; ++k) {
          const SpacePointIndex indexO = survivors[k];
          const ConstSpacePointProxy otherSp = container[indexO];
          const float distance = fastHypot(
              xy[indexO][0] - xM, xy[indexO][1] - yM, zr[indexO][0] - zM);
          float dt = 0;
          if constexpr (isBottomCandidate) {
            dt = tM - otherSp.time() - distance;
          } else {
            dt = otherSp.time() - tM - distance;
          }
          const float varianceDt = varianceTM + otherSp.varianceT();
          survivors[nInTime] = indexO;
          nInTime += static_cast<std::uint32_t>(
              !(dt * dt > m_cfg.timeCutNSigma2 * varianceDt));
        }
        nSurvivors = nInTime;
      }

      if constexpr (interactionPointCut) {
        // The values the remaining cuts compute for each survivor. Not cleared,
        // as above.
        struct PerSurvivor {
          std::array<float, kBlockSize> xNewFrame;
          std::array<float, kBlockSize> yNewFrame;
          std::array<float, kBlockSize> uT;
          std::array<float, kBlockSize> vT;
          std::array<float, kBlockSize> iDeltaR2;
          std::array<float, kBlockSize> deltaZ;
          std::array<std::uint32_t, kBlockSize> kept;
          std::array<std::uint32_t, kBlockSize> keepMask;
          std::array<std::uint32_t, kBlockSize> curvatureIndex;
        };
        // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init)
        PerSurvivor perSurvivor;
        // NOLINTEND(cppcoreguidelines-pro-type-member-init)

        // The transformation and the cuts on the impact parameter and cotTheta
        // for each survivor; those that pass the cotTheta cut and fail the
        // impact parameter test are listed for the curvature cut.
        std::uint32_t nCurvature = 0;
        for (std::uint32_t k = 0; k < nSurvivors; ++k) {
          const SpacePointIndex indexO = survivors[k];
          const float xO = xy[indexO][0];
          const float yO = xy[indexO][1];
          const float zO = zr[indexO][0];
          const float rO = zr[indexO][1];
          const float deltaR = isBottomCandidate ? rM - rO : rO - rM;
          const float deltaZ = isBottomCandidate ? zM - zO : zO - zM;

          // transform SP coordinates to the u-v reference frame
          const float deltaX = xO - xM;
          const float deltaY = yO - yM;

          const float xNewFrame =
              deltaX * middleSpInfo.cosPhiM + deltaY * middleSpInfo.sinPhiM;
          const float yNewFrame =
              deltaY * middleSpInfo.cosPhiM - deltaX * middleSpInfo.sinPhiM;

          const float deltaR2 = deltaX * deltaX + deltaY * deltaY;
          const float iDeltaR2 = 1 / deltaR2;

          const float uT = xNewFrame * iDeltaR2;
          const float vT = yNewFrame * iDeltaR2;

          // We check the interaction point by evaluating the minimal distance
          // between the origin and the straight line connecting the two points
          // in the doublets. Using a geometric similarity, the Im is given by
          // yNewFrame * rM / deltaR > config.impactMax
          // However, we make here an approximation of the impact parameter
          // which is valid under the assumption yNewFrame / xNewFrame is small
          // The correct computation would be:
          // yNewFrame * yNewFrame * rM * rM > config.impactMax *
          // config.impactMax * deltaR2
          const bool ipTest = std::abs(rM * yNewFrame) > impactMax * xNewFrame;

          // check if duplet cotTheta is within the region of interest
          // cotTheta is defined as (deltaZ / deltaR) but instead we multiply
          // cotThetaMax by deltaR to avoid division
          const bool cutCotTheta =
              outsideSymmetricRange(deltaZ, m_cfg.cotThetaMax * deltaR);

          perSurvivor.xNewFrame[k] = xNewFrame;
          perSurvivor.yNewFrame[k] = yNewFrame;
          perSurvivor.uT[k] = uT;
          perSurvivor.vT[k] = vT;
          perSurvivor.iDeltaR2[k] = iDeltaR2;
          perSurvivor.deltaZ[k] = deltaZ;

          // Kept unless cut on cotTheta; where ipTest holds, the curvature cut
          // below decides.
          perSurvivor.keepMask[k] =
              static_cast<std::uint32_t>(!cutCotTheta && !ipTest);
          perSurvivor.curvatureIndex[nCurvature] = k;
          nCurvature += static_cast<std::uint32_t>(ipTest && !cutCotTheta);
        }

        // The curvature cut, for the survivors listed above.
        for (std::uint32_t j = 0; j < nCurvature; ++j) {
          const std::uint32_t k = perSurvivor.curvatureIndex[j];
          perSurvivor.keepMask[k] = static_cast<std::uint32_t>(
              !failsCurvatureCut(perSurvivor.uT[k], perSurvivor.vT[k],
                                 perSurvivor.yNewFrame[k]));
        }

        // the survivors kept
        const std::uint32_t nKept = detail::collectFlagged(
            nSurvivors, perSurvivor.keepMask.data(),
            [](std::uint32_t k) { return k; }, perSurvivor.kept.data());

        for (std::uint32_t j = 0; j < nKept; ++j) {
          const std::uint32_t k = perSurvivor.kept[j];
          const SpacePointIndex indexO = survivors[k];

          const float iDeltaR = std::sqrt(perSurvivor.iDeltaR2[k]);
          const float cotTheta = perSurvivor.deltaZ[k] * iDeltaR;

          // discard doublets based on experiment specific cuts
          if constexpr (experimentCuts) {
            if (!m_cfg.experimentCuts(middleSp, container[indexO], cotTheta,
                                      isBottomCandidate)) {
              continue;
            }
          }

          const float er = calculateError(varianceZ[indexO], varianceR[indexO],
                                          perSurvivor.iDeltaR2[k], cotTheta);

          // fill output vectors
          compatibleDoublets.emplace_back(indexO, cotTheta, iDeltaR, er,
                                          perSurvivor.uT[k], perSurvivor.vT[k],
                                          perSurvivor.xNewFrame[k],
                                          perSurvivor.yNewFrame[k]);
        }
      } else {
        // without the interaction point cut: the cotTheta cut, then the
        // transformation and the doublet
        for (std::uint32_t k = 0; k < nSurvivors; ++k) {
          const SpacePointIndex indexO = survivors[k];
          const float zO = zr[indexO][0];
          const float rO = zr[indexO][1];
          const float deltaR = isBottomCandidate ? rM - rO : rO - rM;
          const float deltaZ = isBottomCandidate ? zM - zO : zO - zM;
          if (outsideSymmetricRange(deltaZ, m_cfg.cotThetaMax * deltaR)) {
            continue;
          }
          addDoublet(indexO, deltaZ, transform(indexO));
        }
      }
    }
  }

  void createDoublets(const ConstSpacePointProxy& middleSp,
                      const MiddleSpInfo& middleSpInfo,
                      SpacePointContainer::ConstSubset candidateSps,
                      DoubletsForMiddleSp& compatibleDoublets) const override {
    createDoubletsImpl(middleSp, middleSpInfo, candidateSps,
                       compatibleDoublets);
  }

  void createDoublets(const ConstSpacePointProxy& middleSp,
                      const MiddleSpInfo& middleSpInfo,
                      SpacePointContainer::ConstRange candidateSps,
                      DoubletsForMiddleSp& compatibleDoublets) const override {
    createDoubletsImpl(middleSp, middleSpInfo, candidateSps,
                       compatibleDoublets);
  }

 private:
  DerivedConfig m_cfg;
};

}  // namespace

std::unique_ptr<DoubletSeedFinder> DoubletSeedFinder::create(
    const DerivedConfig& config) {
  using BooleanOptions =
      boost::mp11::mp_list<std::bool_constant<false>, std::bool_constant<true>>;

  using IsBottomCandidateOptions = BooleanOptions;
  using InteractionPointCutOptions = BooleanOptions;
  using SortedByROptions = BooleanOptions;
  using ExperimentCutsOptions = BooleanOptions;
  using UseTimeOptions = BooleanOptions;
  using ZCutsOptions = boost::mp11::mp_list<
      std::integral_constant<ZCuts, ZCuts::eGeneral>,
      std::integral_constant<ZCuts, ZCuts::eSymmetric>,
      std::integral_constant<ZCuts, ZCuts::eSymmetricNoDeltaZ>>;

  using DoubletOptions =
      boost::mp11::mp_product<boost::mp11::mp_list, IsBottomCandidateOptions,
                              InteractionPointCutOptions, SortedByROptions,
                              ExperimentCutsOptions, UseTimeOptions,
                              ZCutsOptions>;

  // Compared exactly: only an exactly symmetric range gives the same cut.
  const bool symmetricCollisionRegion =
      config.collisionRegionMin == -config.collisionRegionMax;
  const bool unboundedDeltaZ =
      config.deltaZMin == -std::numeric_limits<float>::infinity() &&
      config.deltaZMax == std::numeric_limits<float>::infinity();
  ZCuts configZCuts = ZCuts::eGeneral;
  if (symmetricCollisionRegion && unboundedDeltaZ) {
    configZCuts = ZCuts::eSymmetricNoDeltaZ;
  } else if (symmetricCollisionRegion &&
             config.deltaZMin == -config.deltaZMax) {
    configZCuts = ZCuts::eSymmetric;
  }

  std::unique_ptr<DoubletSeedFinder> result;
  boost::mp11::mp_for_each<DoubletOptions>([&](auto option) {
    using OptionType = decltype(option);

    using IsBottomCandidate = boost::mp11::mp_at_c<OptionType, 0>;
    using InteractionPointCut = boost::mp11::mp_at_c<OptionType, 1>;
    using SortedByR = boost::mp11::mp_at_c<OptionType, 2>;
    using ExperimentCuts = boost::mp11::mp_at_c<OptionType, 3>;
    using UseTime = boost::mp11::mp_at_c<OptionType, 4>;
    using ZCutsOption = boost::mp11::mp_at_c<OptionType, 5>;

    const bool configIsBottomCandidate =
        config.candidateDirection == Direction::Backward();

    if (configIsBottomCandidate != IsBottomCandidate::value ||
        config.interactionPointCut != InteractionPointCut::value ||
        config.spacePointsSortedByRadius != SortedByR::value ||
        config.experimentCuts.connected() != ExperimentCuts::value ||
        config.useTime != UseTime::value || configZCuts != ZCutsOption::value) {
      return;  // skip if the configuration does not match
    }

    // check if we already have an implementation for this configuration
    if (result != nullptr) {
      throw std::runtime_error(
          "DoubletSeedFinder: Multiple implementations found for one "
          "configuration");
    }

    // create the implementation for the given configuration
    result = std::make_unique<Impl<
        IsBottomCandidate::value, InteractionPointCut::value, SortedByR::value,
        ExperimentCuts::value, UseTime::value, ZCutsOption::value>>(config);
  });
  if (result == nullptr) {
    throw std::runtime_error(
        "DoubletSeedFinder: No implementation found for the given "
        "configuration");
  }
  return result;
}

DoubletSeedFinder::DerivedConfig::DerivedConfig(const Config& config,
                                                float bFieldInZ_)
    : Config(config), bFieldInZ(bFieldInZ_) {
  // bFieldInZ is in (pT/radius) natively, no need for conversion
  const float pTPerHelixRadius = bFieldInZ;
  minHelixDiameter2 = square(minPt * 2 / pTPerHelixRadius) * helixCutTolerance;
  timeCutNSigma2 = square(timeCutNSigma);
}

MiddleSpInfo DoubletSeedFinder::computeMiddleSpInfo(
    const ConstSpacePointProxy& spM) {
  const float rM = spM.zr()[1];
  const float uIP = -1 / rM;
  const float cosPhiM = -spM.xy()[0] * uIP;
  const float sinPhiM = -spM.xy()[1] * uIP;
  const float uIP2 = uIP * uIP;

  return {uIP, uIP2, cosPhiM, sinPhiM};
}

}  // namespace Acts
