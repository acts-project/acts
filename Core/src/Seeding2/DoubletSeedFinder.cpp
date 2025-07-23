// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/DoubletSeedFinder.hpp"

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <memory>
#include <stdexcept>
#include <type_traits>

#include <boost/mp11.hpp>

namespace Acts::Experimental {

template <bool isBottomCandidate, bool interactionPointCut, bool sortedByR>
class DoubletSeedFinder::Impl final : public DoubletSeedFinder::ImplBase {
 public:
  explicit Impl(const DerivedConfig& config) : ImplBase(config) {}

  /// Iterates over dublets and tests the compatibility by applying a series of
  /// cuts that can be tested with only two SPs.
  ///
  /// @param config Doublet cuts that define the compatibility of space points
  /// @param middleSp Space point candidate to be used as middle SP in a seed
  /// @param middleSpInfo Information about the middle space point
  /// @param candidateSps Range or subet of space points to be used as candidates
  ///   for middle SP in a seed. In case of `sortedByR` - an offset will be
  ///   applied based on the middle SP radius.
  /// @param compatibleDoublets Output container for compatible doublets
  template <typename CandidateSps>
  void createDoubletsImpl(
      const ConstSpacePointProxy2& middleSp,
      const DoubletSeedFinder::MiddleSpInfo& middleSpInfo,
      const CandidateSps& candidateSps,
      DoubletSeedFinder::DoubletsForMiddleSp& compatibleDoublets) const {
    const float impactMax =
        isBottomCandidate ? -m_cfg.impactMax : m_cfg.impactMax;

    const float xM = middleSp.x();
    const float yM = middleSp.y();
    const float zM = middleSp.z();
    const float rM = middleSp.r();
    const float varianceZM = middleSp.varianceZ();
    const float varianceRM = middleSp.varianceR();

    // equivalent to impactMax / (rM * rM);
    const float vIPAbs = impactMax * middleSpInfo.uIP2;

    float deltaR = 0.;
    float deltaZ = 0.;

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

    const SpacePointContainer2& container = candidateSps.container();
    for (auto [indexO, xO, yO, zO, rO, varianceZO, varianceRO] :
         candidateSps.zip(container.xColumn(), container.yColumn(),
                          container.zColumn(), container.rColumn(),
                          container.varianceZColumn(),
                          container.varianceRColumn())) {
      if constexpr (isBottomCandidate) {
        deltaR = rM - rO;

        if constexpr (sortedByR) {
          // if r-distance is too small we are done
          if (deltaR < m_cfg.deltaRMin) {
            break;
          }
        }
      } else {
        deltaR = rO - rM;

        if constexpr (sortedByR) {
          // if r-distance is too big we are done
          if (deltaR > m_cfg.deltaRMax) {
            break;
          }
        }
      }

      if constexpr (!sortedByR) {
        if (outsideRangeCheck(deltaR, m_cfg.deltaRMin, m_cfg.deltaRMax)) {
          continue;
        }
      }

      if constexpr (isBottomCandidate) {
        deltaZ = zM - zO;
      } else {
        deltaZ = zO - zM;
      }

      if (outsideRangeCheck(deltaZ, m_cfg.deltaZMin, m_cfg.deltaZMax)) {
        continue;
      }

      // the longitudinal impact parameter zOrigin is defined as (zM - rM *
      // cotTheta) where cotTheta is the ratio Z/R (forward angle) of space
      // point duplet but instead we calculate (zOrigin * deltaR) and multiply
      // collisionRegion by deltaR to avoid divisions
      const float zOriginTimesDeltaR = zM * deltaR - rM * deltaZ;
      // check if duplet origin on z axis within collision region
      if (outsideRangeCheck(zOriginTimesDeltaR,
                            m_cfg.collisionRegionMin * deltaR,
                            m_cfg.collisionRegionMax * deltaR)) {
        continue;
      }

      // if interactionPointCut is false we apply z cuts before coordinate
      // transformation to avoid unnecessary calculations. If
      // interactionPointCut is true we apply the curvature cut first because it
      // is more frequent but requires the coordinate transformation
      if constexpr (!interactionPointCut) {
        // check if duplet cotTheta is within the region of interest
        // cotTheta is defined as (deltaZ / deltaR) but instead we multiply
        // cotThetaMax by deltaR to avoid division
        if (outsideRangeCheck(deltaZ, -m_cfg.cotThetaMax * deltaR,
                              m_cfg.cotThetaMax * deltaR)) {
          continue;
        }

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

        const float iDeltaR = std::sqrt(iDeltaR2);
        const float cotTheta = deltaZ * iDeltaR;

        const float er =
            calculateError(varianceZO, varianceRO, iDeltaR2, cotTheta);

        // fill output vectors
        compatibleDoublets.emplace_back(
            indexO, {cotTheta, iDeltaR, er, uT, vT, xNewFrame, yNewFrame});
        continue;
      }

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
      // between the origin and the straight line connecting the two points in
      // the doublets. Using a geometric similarity, the Im is given by
      // yNewFrame * rM / deltaR > config.impactMax
      // However, we make here an approximation of the impact parameter
      // which is valid under the assumption yNewFrame / xNewFrame is small
      // The correct computation would be:
      // yNewFrame * yNewFrame * rM * rM > config.impactMax *
      // config.impactMax * deltaR2
      if (std::abs(rM * yNewFrame) > impactMax * xNewFrame) {
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
        if ((bCoef * bCoef) * m_cfg.minHelixDiameter2 > 1 + aCoef * aCoef) {
          continue;
        }
      }

      // check if duplet cotTheta is within the region of interest
      // cotTheta is defined as (deltaZ / deltaR) but instead we multiply
      // cotThetaMax by deltaR to avoid division
      if (outsideRangeCheck(deltaZ, -m_cfg.cotThetaMax * deltaR,
                            m_cfg.cotThetaMax * deltaR)) {
        continue;
      }

      const float iDeltaR = std::sqrt(iDeltaR2);
      const float cotTheta = deltaZ * iDeltaR;

      // discard bottom-middle doublets in a certain (r, eta) region according
      // to detector specific cuts
      if constexpr (isBottomCandidate) {
        if (m_cfg.experimentCuts.connected() &&
            !m_cfg.experimentCuts(rO, cotTheta)) {
          continue;
        }
      }

      const float er =
          calculateError(varianceZO, varianceRO, iDeltaR2, cotTheta);

      // fill output vectors
      compatibleDoublets.emplace_back(
          indexO, {cotTheta, iDeltaR, er, uT, vT, xNewFrame, yNewFrame});
    }
  }

  void createDoublets(const ConstSpacePointProxy2& middleSp,
                      const MiddleSpInfo& middleSpInfo,
                      SpacePointContainer2::ConstSubset& candidateSps,
                      DoubletsForMiddleSp& compatibleDoublets) const override {
    if constexpr (sortedByR) {
      const float rM = middleSp.r();

      // find the first SP inside the radius region of interest and update
      // the iterator so we don't need to look at the other SPs again
      std::uint32_t offset = 0;
      for (ConstSpacePointProxy2 otherSp : candidateSps) {
        if constexpr (isBottomCandidate) {
          // if r-distance is too big, try next SP in bin
          if (rM - otherSp.r() <= m_cfg.deltaRMax) {
            break;
          }
        } else {
          // if r-distance is too small, try next SP in bin
          if (otherSp.r() - rM >= m_cfg.deltaRMin) {
            break;
          }
        }

        ++offset;
      }
      candidateSps = candidateSps.container().subset(
          candidateSps.subset().subspan(offset));
    }

    createDoubletsImpl(middleSp, middleSpInfo, candidateSps,
                       compatibleDoublets);
  }

  void createDoublets(const ConstSpacePointProxy2& middleSp,
                      const MiddleSpInfo& middleSpInfo,
                      SpacePointContainer2::ConstRange& candidateSps,
                      DoubletsForMiddleSp& compatibleDoublets) const override {
    if constexpr (sortedByR) {
      const float rM = middleSp.r();

      // find the first SP inside the radius region of interest and update
      // the iterator so we don't need to look at the other SPs again
      std::uint32_t offset = 0;
      for (ConstSpacePointProxy2 otherSp : candidateSps) {
        if constexpr (isBottomCandidate) {
          // if r-distance is too big, try next SP in bin
          if (rM - otherSp.r() <= m_cfg.deltaRMax) {
            break;
          }
        } else {
          // if r-distance is too small, try next SP in bin
          if (otherSp.r() - rM >= m_cfg.deltaRMin) {
            break;
          }
        }

        ++offset;
      }
      candidateSps = candidateSps.subrange(offset);
    }

    createDoubletsImpl(middleSp, middleSpInfo, candidateSps,
                       compatibleDoublets);
  }
};

std::shared_ptr<DoubletSeedFinder::ImplBase> DoubletSeedFinder::makeImpl(
    const DerivedConfig& config) {
  using BooleanOptions =
      boost::mp11::mp_list<std::bool_constant<false>, std::bool_constant<true>>;
  using IsBottomCandidateOptions = BooleanOptions;
  using InteractionPointCutOptions = BooleanOptions;
  using SortedByROptions = BooleanOptions;

  using DoubletOptions =
      boost::mp11::mp_product<boost::mp11::mp_list, IsBottomCandidateOptions,
                              InteractionPointCutOptions, SortedByROptions>;

  std::shared_ptr<DoubletSeedFinder::ImplBase> result;
  boost::mp11::mp_for_each<DoubletOptions>([&](auto option) {
    using OptionType = decltype(option);
    using IsBottomCandidate = boost::mp11::mp_first<OptionType>;
    using InteractionPointCut = boost::mp11::mp_second<OptionType>;
    using SortedByR = boost::mp11::mp_third<OptionType>;

    const bool configIsBottomCandidate =
        config.candidateDirection == Direction::Backward();

    if (configIsBottomCandidate != IsBottomCandidate::value ||
        config.interactionPointCut != InteractionPointCut::value ||
        config.spacePointsSortedByRadius != SortedByR::value) {
      return;  // skip if the configuration does not match
    }

    // check if we already have an implementation for this configuration
    if (result != nullptr) {
      throw std::runtime_error(
          "DoubletSeedFinder: Multiple implementations found for one "
          "configuration");
    }

    // create the implementation for the given configuration
    result = std::make_shared<
        DoubletSeedFinder::Impl<IsBottomCandidate::value,
                                InteractionPointCut::value, SortedByR::value>>(
        config);
  });
  return result;
}

DoubletSeedFinder::DerivedConfig::DerivedConfig(const Config& config,
                                                float bFieldInZ_)
    : Config(config), bFieldInZ(bFieldInZ_) {
  // bFieldInZ is in (pT/radius) natively, no need for conversion
  const float pTPerHelixRadius = bFieldInZ;
  minHelixDiameter2 = square(minPt * 2 / pTPerHelixRadius) * helixCutTolerance;
}

DoubletSeedFinder::MiddleSpInfo DoubletSeedFinder::computeMiddleSpInfo(
    const ConstSpacePointProxy2& spM) {
  const float rM = spM.r();
  const float uIP = -1 / rM;
  const float cosPhiM = -spM.x() * uIP;
  const float sinPhiM = -spM.y() * uIP;
  const float uIP2 = uIP * uIP;

  return {uIP, uIP2, cosPhiM, sinPhiM};
}

DoubletSeedFinder::DoubletSeedFinder(const DerivedConfig& config)
    : m_impl(makeImpl(config)) {
  if (m_impl == nullptr) {
    throw std::runtime_error(
        "DoubletSeedFinder: No implementation found for the given "
        "configuration");
  }
}

}  // namespace Acts::Experimental
