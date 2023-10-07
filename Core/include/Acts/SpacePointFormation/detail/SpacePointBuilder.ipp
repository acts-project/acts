// This file is part of the Acts project.
//
// Copyright (C) 2018-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
namespace Acts {

template <typename spacepoint_t>
SpacePointBuilder<spacepoint_t>::SpacePointBuilder(
    const SpacePointBuilderConfig& cfg,
    std::function<spacepoint_t(Acts::Vector3, Acts::Vector2,
                               boost::container::static_vector<SourceLink, 2>)>
        func,
    std::unique_ptr<const Logger> logger)
    : m_config(cfg), m_spConstructor(func), m_logger(std::move(logger)) {
  m_spUtility = std::make_shared<SpacePointUtility>(cfg);
}

template <typename spacepoint_t>
template <template <typename...> typename container_t>
void SpacePointBuilder<spacepoint_t>::buildSpacePoint(
    const GeometryContext& gctx, const std::vector<SourceLink>& sourceLinks,
    const SpacePointBuilderOptions& opt,
    std::back_insert_iterator<container_t<spacepoint_t>> spacePointIt) const {
  const unsigned int num_slinks = sourceLinks.size();

  Acts::Vector3 gPos = Acts::Vector3::Zero();
  Acts::Vector2 gCov = Acts::Vector2::Zero();

  if (num_slinks == 1) {  // pixel SP formation
    auto slink = sourceLinks.at(0);
    auto [param, cov] = opt.paramCovAccessor(sourceLinks.at(0));
    auto gPosCov = m_spUtility->globalCoords(gctx, slink, param, cov);
    gPos = gPosCov.first;
    gCov = gPosCov.second;
  } else if (num_slinks == 2) {  // strip SP formation

    const auto& ends1 = opt.stripEndsPair.first;
    const auto& ends2 = opt.stripEndsPair.second;

    Acts::SpacePointParameters spParams;

    if (!m_config.usePerpProj) {  // default strip SP building

      auto spFound = m_spUtility->calculateStripSPPosition(
          ends1, ends2, opt.vertex, spParams, opt.stripLengthTolerance);

      if (!spFound.ok()) {
        spFound = m_spUtility->recoverSpacePoint(spParams,
                                                 opt.stripLengthGapTolerance);
      }

      if (!spFound.ok()) {
        return;
      }

      gPos = 0.5 *
             (ends1.first + ends1.second + spParams.m * spParams.firstBtmToTop);

    } else {  // for cosmic without vertex constraint

      auto resultPerpProj =
          m_spUtility->calcPerpendicularProjection(ends1, ends2, spParams);

      if (!resultPerpProj.ok()) {
        return;
      }
      gPos = ends1.first + resultPerpProj.value() * spParams.firstBtmToTop;
    }

    double theta =
        acos(spParams.firstBtmToTop.dot(spParams.secondBtmToTop) /
             (spParams.firstBtmToTop.norm() * spParams.secondBtmToTop.norm()));

    gCov = m_spUtility->calcRhoZVars(gctx, sourceLinks.at(0), sourceLinks.at(1),
                                     opt.paramCovAccessor, gPos, theta);

  } else {
    ACTS_ERROR("More than 2 sourceLinks are given for a space point.");
  }
  boost::container::static_vector<SourceLink, 2> slinks(sourceLinks.begin(),
                                                        sourceLinks.end());

  spacePointIt = m_spConstructor(gPos, gCov, std::move(slinks));
}

template <typename spacepoint_t>
void SpacePointBuilder<spacepoint_t>::makeSourceLinkPairs(
    const GeometryContext& gctx, const std::vector<SourceLink>& slinksFront,
    const std::vector<SourceLink>& slinksBack,
    std::vector<std::pair<SourceLink, SourceLink>>& slinkPairs,
    const StripPairOptions& pairOpt) const {
  if (slinksFront.empty() || slinksBack.empty()) {
    return;
  }
  double minDistance = 0;
  unsigned int closestIndex = 0;

  for (unsigned int i = 0; i < slinksFront.size(); i++) {
    const auto& slinkFront = slinksFront[i];
    minDistance = std::numeric_limits<double>::max();
    closestIndex = slinksBack.size();
    for (unsigned int j = 0; j < slinksBack.size(); j++) {
      const auto& slinkBack = slinksBack[j];

      const auto [paramFront, covFront] = pairOpt.paramCovAccessor(slinkFront);
      const auto [gposFront, gcovFront] =
          m_spUtility->globalCoords(gctx, slinkFront, paramFront, covFront);

      const auto [paramBack, covBack] = pairOpt.paramCovAccessor(slinkBack);
      const auto [gposBack, gcovBack] =
          m_spUtility->globalCoords(gctx, slinkBack, paramBack, covBack);

      auto res = m_spUtility->differenceOfMeasurementsChecked(
          gposFront, gposBack, pairOpt.vertex, pairOpt.diffDist,
          pairOpt.diffPhi2, pairOpt.diffTheta2);
      if (!res.ok()) {
        continue;
      }
      const auto distance = res.value();
      if (distance >= 0. && distance < minDistance) {
        minDistance = distance;
        closestIndex = j;
      }
    }
    if (closestIndex < slinksBack.size()) {
      slinkPairs.emplace_back(slinksFront[i], slinksBack[closestIndex]);
    }
  }
}

}  // namespace Acts
