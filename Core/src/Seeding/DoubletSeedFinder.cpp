// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/DoubletSeedFinder.hpp"

#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/SpacePointFormation/detail/StripSpacePointCalibrationImpl.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <type_traits>

#include <boost/mp11.hpp>
#include <boost/mp11/algorithm.hpp>

#include "detail/SortByKey.hpp"

namespace Acts {

void DoubletsForMiddleSp::sortByCotTheta(
    const IndexRange& range,
    std::vector<IndexAndCotTheta>& indexAndCotTheta) const {
  SortStorage storage;
  sortByCotTheta(range, indexAndCotTheta, storage);
}

void DoubletsForMiddleSp::sortByCotTheta(
    const IndexRange& range, std::vector<IndexAndCotTheta>& indexAndCotTheta,
    SortStorage& storage) const {
  indexAndCotTheta.clear();
  indexAndCotTheta.reserve(range.second - range.first);
  for (Index i = range.first; i < range.second; ++i) {
    indexAndCotTheta.emplace_back(i, m_doublets[i].cotTheta);
  }

  const auto cotThetaOf = [](const IndexAndCotTheta& item) {
    return item.cotTheta;
  };

  detail::sortByKey(std::span(indexAndCotTheta), storage, cotThetaOf);
}

namespace {

template <bool isBottomCandidate, bool interactionPointCut, bool sortedByR,
          bool useTime, bool stripCalibration, bool symmetricCollisionRegion>
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
    [[maybe_unused]] const float vIPAbs = impactMax * middleSpInfo.uIP2;

    // only filled with the strip calibration
    [[maybe_unused]] float calibrationTolerance = 0;
    [[maybe_unused]] OuterStripSpacePointCalibrationDetailsDerived
        calibrationM{};
    if constexpr (stripCalibration) {
      calibrationTolerance = *m_cfg.stripCalibrationTolerance;
      calibrationM = detail::deriveOuterStripSpacePointCalibrationDetails(
          middleSp.outerStripCalibrationDetails());
    }

    const auto outsideRangeCheck = [](const float value, const float min,
                                      const float max) {
      // intentionally using `|` after profiling. faster due to better branch
      // prediction
      return static_cast<bool>(static_cast<int>(value < min) |
                               static_cast<int>(value > max));
    };

    const auto calculateError = [&](float varianceZO, float varianceRO,
                                    float iDeltaR2, float cotTheta) {
      return iDeltaR2 * ((varianceZM + varianceZO) +
                         (cotTheta * cotTheta) * (varianceRM + varianceRO));
    };

    const SpacePointContainer& container = candidateSps.container();
    const auto xyColumn = container.xyColumn().data();
    const auto zrColumn = container.zrColumn().data();
    const auto varianceZColumn = container.varianceZColumn().data();
    const auto varianceRColumn = container.varianceRColumn().data();
    // only filled when the time cut is enabled
    [[maybe_unused]] std::span<const float> timeColumn;
    [[maybe_unused]] std::span<const float> varianceTColumn;
    if constexpr (useTime) {
      timeColumn = container.timeColumn().data();
      varianceTColumn = container.varianceTColumn().data();
    }

    // the index in the container of the k-th candidate
    const auto candidateIndex = [&](std::size_t k) -> std::size_t {
      if constexpr (std::is_same_v<CandidateSps,
                                   SpacePointContainer::ConstRange>) {
        return candidateSps.range().first + k;
      } else {
        return candidateSps.subset()[k];
      }
    };

    // the difference of a coordinate between the outer and the inner space
    // point of the doublet
    const auto outerMinusInner = [](float middle, float other) {
      return isBottomCandidate ? middle - other : other - middle;
    };

    // the number of candidates to look at; if they are sorted by radius, they
    // end at the first candidate outside the radial range
    std::size_t nCandidates = candidateSps.size();
    if constexpr (sortedByR) {
      const auto insideRadialRange = [&](const ConstSpacePointProxy& sp) {
        const float deltaR = outerMinusInner(rM, sp.zr()[1]);
        return isBottomCandidate ? deltaR >= m_cfg.deltaRMin
                                 : deltaR <= m_cfg.deltaRMax;
      };
      nCandidates =
          std::ranges::partition_point(candidateSps, insideRadialRange) -
          candidateSps.begin();
    }

    // the candidates go through in blocks: in each block the cuts that need
    // only z and r are made for every candidate first, then the cuts in the
    // u-v frame for those that passed, then the other cuts for those that
    // passed them too
    constexpr std::size_t blockSize = 256;
    // measured difference with initialization
    // please don't remove this without profiling
    // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init)
    std::array<std::uint32_t, blockSize> passesZRTimeCuts;
    std::array<SpacePointIndex, blockSize> passedIndices;
    // the x-y difference of each candidate that passed, and its u-v frame with
    // the cuts made in it; only filled without the strip calibration
    [[maybe_unused]] std::array<float, blockSize> deltaXs;
    [[maybe_unused]] std::array<float, blockSize> deltaYs;
    [[maybe_unused]] std::array<float, blockSize> xNewFrames;
    [[maybe_unused]] std::array<float, blockSize> yNewFrames;
    [[maybe_unused]] std::array<float, blockSize> iDeltaR2s;
    [[maybe_unused]] std::array<std::uint32_t, blockSize> passesFrameCuts;
    [[maybe_unused]] std::array<std::uint32_t, blockSize> framePassedPositions;
    // NOLINTEND(cppcoreguidelines-pro-type-member-init)

    // the experiment specific cuts, which need the doublet's cotTheta
    const auto passesExperimentCuts = [&](const SpacePointIndex indexO,
                                          const float cotTheta) {
      return !m_cfg.experimentCuts.connected() ||
             m_cfg.experimentCuts(middleSp, container[indexO], cotTheta,
                                  isBottomCandidate);
    };

    // the doublet's error and its position in the u-v frame, added to the
    // output
    const auto addDoublet = [&](const SpacePointIndex indexO,
                                const float cotTheta, const float iDeltaR,
                                const float iDeltaR2, const float xNewFrame,
                                const float yNewFrame) {
      const float uT = xNewFrame * iDeltaR2;
      const float vT = yNewFrame * iDeltaR2;
      const float er = calculateError(
          varianceZColumn[indexO], varianceRColumn[indexO], iDeltaR2, cotTheta);
      compatibleDoublets.emplace_back(indexO, cotTheta, iDeltaR, er, uT, vT,
                                      xNewFrame, yNewFrame);
    };

    for (std::size_t blockBegin = 0; blockBegin < nCandidates;
         blockBegin += blockSize) {
      const std::size_t blockCount =
          std::min(blockSize, nCandidates - blockBegin);

      // the z-r and time cuts for every candidate of the block, one flag each
      std::uint32_t nPassed = 0;
      for (std::size_t i = 0; i < blockCount; ++i) {
        const std::size_t indexO = candidateIndex(blockBegin + i);
        const auto& zrO = zrColumn[indexO];
        const float deltaR = outerMinusInner(rM, zrO[1]);
        const float deltaZ = outerMinusInner(zM, zrO[0]);

        bool outside = false;
        if constexpr (!sortedByR) {
          outside = outsideRangeCheck(deltaR, m_cfg.deltaRMin, m_cfg.deltaRMax);
        }

        outside |= std::abs(deltaZ) > m_cfg.absDeltaZMax;

        // the longitudinal impact parameter zOrigin is defined as (zM - rM *
        // cotTheta) where cotTheta is the ratio Z/R (forward angle) of space
        // point duplet but instead we calculate (zOrigin * deltaR) and
        // multiply collisionRegion by deltaR to avoid divisions
        const float zOriginTimesDeltaR = zM * deltaR - rM * deltaZ;
        // check if duplet origin on z axis within collision region
        if constexpr (symmetricCollisionRegion) {
          outside |=
              std::abs(zOriginTimesDeltaR) > m_cfg.collisionRegionMax * deltaR;
        } else {
          outside |= outsideRangeCheck(zOriginTimesDeltaR,
                                       m_cfg.collisionRegionMin * deltaR,
                                       m_cfg.collisionRegionMax * deltaR);
        }

        // cotTheta = deltaZ / deltaR, compared without the division. The strip
        // calibration moves deltaZ and deltaR, so there the cut comes after it
        if constexpr (!stripCalibration) {
          outside |= std::abs(deltaZ) > m_cfg.cotThetaMax * deltaR;
        }

        // check the time compatibility of the two space points. the time
        // difference is corrected for the time of flight along the straight
        // line between them, assuming an outgoing particle at the speed of
        // light (which is 1 in ACTS units)
        if constexpr (useTime) {
          const auto& xyO = xyColumn[indexO];
          const float distance =
              fastHypot(xyO[0] - xM, xyO[1] - yM, zrO[0] - zM);
          const float dt = outerMinusInner(tM, timeColumn[indexO]) - distance;
          const float varianceDt = varianceTM + varianceTColumn[indexO];
          outside |= dt * dt > m_cfg.timeCutNSigma2 * varianceDt;
        }

        passesZRTimeCuts[i] = outside ? 0 : 1;
        nPassed += passesZRTimeCuts[i];
      }

      if (nPassed == 0) {
        continue;
      }

      // collect the indices of the candidates that passed, in order
      std::size_t next = 0;
      for (std::size_t i = 0; i < blockCount; ++i) {
        passedIndices[next] =
            static_cast<SpacePointIndex>(candidateIndex(blockBegin + i));
        next += passesZRTimeCuts[i];
      }

      if constexpr (!stripCalibration) {
        // the x-y difference of every candidate that passed, read in a loop of
        // its own: the u-v pass is not vectorised with these loads in it
        for (std::uint32_t k = 0; k < nPassed; ++k) {
          const auto& xyO = xyColumn[passedIndices[k]];
          deltaXs[k] = xyO[0] - xM;
          deltaYs[k] = xyO[1] - yM;
        }

        // the u-v frame of every candidate that passed, and the interaction
        // point and curvature cuts in it, one flag each, with no exit between
        // candidates so that the arithmetic is vectorised
        for (std::uint32_t k = 0; k < nPassed; ++k) {
          // transform SP coordinates to the u-v reference frame
          const float deltaX = deltaXs[k];
          const float deltaY = deltaYs[k];

          const float xNewFrame =
              deltaX * middleSpInfo.cosPhiM + deltaY * middleSpInfo.sinPhiM;
          const float yNewFrame =
              deltaY * middleSpInfo.cosPhiM - deltaX * middleSpInfo.sinPhiM;

          const float deltaR2 = deltaX * deltaX + deltaY * deltaY;
          const float iDeltaR2 = 1 / deltaR2;

          xNewFrames[k] = xNewFrame;
          yNewFrames[k] = yNewFrame;
          iDeltaR2s[k] = iDeltaR2;

          if constexpr (interactionPointCut) {
            const float uT = xNewFrame * iDeltaR2;
            const float vT = yNewFrame * iDeltaR2;

            // We check the interaction point by evaluating the minimal
            // distance between the origin and the straight line connecting the
            // two points in the doublets. Using a geometric similarity, the Im
            // is given by yNewFrame * rM / deltaR > config.impactMax
            // However, we make here an approximation of the impact parameter
            // which is valid under the assumption yNewFrame / xNewFrame is
            // small. The correct computation would be:
            // yNewFrame * yNewFrame * rM * rM > config.impactMax *
            // config.impactMax * deltaR2
            const bool farFromIP =
                std::abs(rM * yNewFrame) > impactMax * xNewFrame;

            // in the rotated frame the interaction point is positioned at
            // x = -rM and y ~= impactParam
            const float vIP = (yNewFrame > 0) ? -vIPAbs : vIPAbs;

            // we can obtain aCoef as the slope dv/du of the linear function,
            // estimated using du and dv between the two SP bCoef is obtained
            // by inserting aCoef into the linear equation.
            // the slope is made for every candidate so that the loop
            // vectorises, where before it was made only for those far from the
            // interaction point. the others pass within impactMax of the
            // origin, and du vanishes only for a line through the origin at the
            // other space point's radius, so for them du is of the order of
            // 1 / r, as for the others: no division by zero, and no overflow
            // in aCoef or bCoef, is added
            const float aCoef = (vT - vIP) / (uT - middleSpInfo.uIP);
            const float bCoef = vIP - aCoef * middleSpInfo.uIP;
            // the distance of the straight line from the origin (radius of the
            // circle) is related to aCoef and bCoef by d^2 = bCoef^2 / (1 +
            // aCoef^2) = 1 / (radius^2) and we can apply the cut on the
            // curvature
            const bool tooCurved =
                (bCoef * bCoef) * m_cfg.minHelixDiameter2 > 1 + aCoef * aCoef;

            // `&`, not `&&`: no branch
            passesFrameCuts[k] = (farFromIP & tooCurved) ? 0 : 1;
          } else {
            // with the cut off every candidate is kept
            passesFrameCuts[k] = 1;
          }
        }

        // collect the positions, in the arrays of the candidates that passed
        // the z-r cuts, of those that also passed the frame cuts, in order
        std::uint32_t nFramePassed = 0;
        for (std::uint32_t k = 0; k < nPassed; ++k) {
          framePassedPositions[nFramePassed] = k;
          nFramePassed += passesFrameCuts[k];
        }

        for (const std::uint32_t k :
             std::span(framePassedPositions).first(nFramePassed)) {
          const SpacePointIndex indexO = passedIndices[k];
          const float deltaZ = outerMinusInner(zM, zrColumn[indexO][0]);
          const float iDeltaR = std::sqrt(iDeltaR2s[k]);
          const float cotTheta = deltaZ * iDeltaR;
          if (!passesExperimentCuts(indexO, cotTheta)) {
            continue;
          }
          addDoublet(indexO, cotTheta, iDeltaR, iDeltaR2s[k], xNewFrames[k],
                     yNewFrames[k]);
        }
      } else {
        for (const SpacePointIndex indexO :
             std::span(passedIndices).first(nPassed)) {
          float xO = xyColumn[indexO][0];
          float yO = xyColumn[indexO][1];
          float zO = zrColumn[indexO][0];
          const float rO = zrColumn[indexO][1];

          float deltaR = outerMinusInner(rM, rO);
          float deltaZ = outerMinusInner(zM, zO);

          // calibrate both space points against their chord
          const std::array<float, 3> chord{xO - xM, yO - yM, zO - zM};
          std::array<float, 3> calibratedM{};
          std::array<float, 3> calibratedO{};
          if (!detail::calibrateOuterStripSpacePoint(
                  chord, calibrationM, calibratedM, calibrationTolerance) ||
              !detail::calibrateOuterStripSpacePoint(
                  chord,
                  detail::deriveOuterStripSpacePointCalibrationDetails(
                      container[indexO].outerStripCalibrationDetails()),
                  calibratedO, calibrationTolerance)) {
            continue;
          }

          // middle space point of this doublet, calibrated
          const float xMd = calibratedM[0];
          const float yMd = calibratedM[1];
          const float rMd = fastHypot(xMd, yMd);
          const float uIP = -1 / rMd;
          const float vIPAbsD = impactMax * uIP * uIP;

          xO = calibratedO[0];
          yO = calibratedO[1];
          zO = calibratedO[2];
          const float rOd = fastHypot(xO, yO);

          deltaR = outerMinusInner(rMd, rOd);
          deltaZ = outerMinusInner(calibratedM[2], zO);

          // if interactionPointCut is false we apply z cuts before coordinate
          // transformation to avoid unnecessary calculations. If
          // interactionPointCut is true we apply the curvature cut first
          // because it is more frequent but requires the coordinate
          // transformation
          if constexpr (!interactionPointCut) {
            // the cotTheta cut, after the calibration moved deltaZ and deltaR
            if (std::abs(deltaZ) > m_cfg.cotThetaMax * deltaR) {
              continue;
            }
          }

          // transform SP coordinates to the u-v reference frame
          const float deltaX = xO - xMd;
          const float deltaY = yO - yMd;

          const float xNewFrame =
              deltaX * middleSpInfo.cosPhiM + deltaY * middleSpInfo.sinPhiM;
          const float yNewFrame =
              deltaY * middleSpInfo.cosPhiM - deltaX * middleSpInfo.sinPhiM;

          const float deltaR2 = deltaX * deltaX + deltaY * deltaY;
          const float iDeltaR2 = 1 / deltaR2;

          if constexpr (interactionPointCut) {
            const float uT = xNewFrame * iDeltaR2;
            const float vT = yNewFrame * iDeltaR2;

            // the interaction point and curvature cuts as above, with the
            // calibrated middle space point
            if (std::abs(rMd * yNewFrame) > impactMax * xNewFrame) {
              const float vIP = (yNewFrame > 0) ? -vIPAbsD : vIPAbsD;
              const float aCoef = (vT - vIP) / (uT - uIP);
              const float bCoef = vIP - aCoef * uIP;
              if ((bCoef * bCoef) * m_cfg.minHelixDiameter2 >
                  1 + aCoef * aCoef) {
                continue;
              }
            }

            // the cotTheta cut, after the calibration moved deltaZ and deltaR
            if (std::abs(deltaZ) > m_cfg.cotThetaMax * deltaR) {
              continue;
            }
          }

          const float iDeltaR = std::sqrt(iDeltaR2);
          const float cotTheta = deltaZ * iDeltaR;
          if (!passesExperimentCuts(indexO, cotTheta)) {
            continue;
          }
          addDoublet(indexO, cotTheta, iDeltaR, iDeltaR2, xNewFrame, yNewFrame);
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
  using UseTimeOptions = BooleanOptions;
  using StripCalibrationOptions = BooleanOptions;
  using SymmetricCollisionRegionOptions = BooleanOptions;

  using DoubletOptions =
      boost::mp11::mp_product<boost::mp11::mp_list, IsBottomCandidateOptions,
                              InteractionPointCutOptions, SortedByROptions,
                              UseTimeOptions, StripCalibrationOptions,
                              SymmetricCollisionRegionOptions>;

  std::unique_ptr<DoubletSeedFinder> result;
  boost::mp11::mp_for_each<DoubletOptions>([&](auto option) {
    using OptionType = decltype(option);

    using IsBottomCandidate = boost::mp11::mp_at_c<OptionType, 0>;
    using InteractionPointCut = boost::mp11::mp_at_c<OptionType, 1>;
    using SortedByR = boost::mp11::mp_at_c<OptionType, 2>;
    using UseTime = boost::mp11::mp_at_c<OptionType, 3>;
    using StripCalibration = boost::mp11::mp_at_c<OptionType, 4>;
    using SymmetricCollisionRegion = boost::mp11::mp_at_c<OptionType, 5>;

    const bool configIsBottomCandidate =
        config.candidateDirection == Direction::Backward();

    if (configIsBottomCandidate != IsBottomCandidate::value ||
        config.interactionPointCut != InteractionPointCut::value ||
        config.spacePointsSortedByRadius != SortedByR::value ||
        config.useTime != UseTime::value ||
        config.stripCalibrationTolerance.has_value() !=
            StripCalibration::value ||
        (config.collisionRegionMin == -config.collisionRegionMax) !=
            SymmetricCollisionRegion::value) {
      return;  // skip if the configuration does not match
    }

    // check if we already have an implementation for this configuration
    if (result != nullptr) {
      throw std::runtime_error(
          "DoubletSeedFinder: Multiple implementations found for one "
          "configuration");
    }

    // create the implementation for the given configuration
    result = std::make_unique<
        Impl<IsBottomCandidate::value, InteractionPointCut::value,
             SortedByR::value, UseTime::value, StripCalibration::value,
             SymmetricCollisionRegion::value>>(config);
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
  // the deltaZ cut is symmetric, deltaZMin is only checked
  if (deltaZMin != -std::numeric_limits<float>::infinity() &&
      deltaZMin != -deltaZMax) {
    throw std::invalid_argument(
        "DoubletSeedFinder: the deltaZ cut is symmetric: deltaZMin must be "
        "-deltaZMax");
  }

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
