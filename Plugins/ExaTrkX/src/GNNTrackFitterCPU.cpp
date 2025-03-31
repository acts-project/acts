// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/GNNTrackFitterCPU.hpp"

#include "Acts/EventData/TrackParameterHelpers.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"

namespace Acts {

std::optional<BoundTrackParameters> GNNParametersBuilderCPU::buildParameters(
    const std::vector<float> &spacepointFeatures,
    const std::vector<Acts::GeometryIdentifier> &geoIds,
    const std::vector<int> &candidate,
    Acts::MagneticFieldProvider::Cache &bCache,
    const Acts::GeometryContext &gctx) const {
  ACTS_VERBOSE("Try to get seed from prototrack with " << candidate.size());

  // in this case we cannot seed properly
  if (candidate.size() < 3) {
    ACTS_VERBOSE(
        "Cannot seed because less then three hits with unique (layer, "
        "volume)");
    return {};
  }

  auto getR = [&](int sp) {
    return spacepointFeatures.at(sp * m_cfg.nFeatures + m_cfg.rIdx) *
           m_cfg.rScale;
  };
  auto getZ = [&](int sp) {
    return spacepointFeatures.at(sp * m_cfg.nFeatures + m_cfg.zIdx) *
           m_cfg.zScale;
  };
  auto getPhi = [&](int sp) {
    return spacepointFeatures.at(sp * m_cfg.nFeatures + m_cfg.phiIdx) *
           m_cfg.phiScale;
  };
  auto getXYZ = [&](int sp) {
    return Acts::Vector3{getR(sp) * std::cos(getPhi(sp)),
                         getR(sp) * std::sin(getPhi(sp)), getZ(sp)};
  };

  ACTS_VERBOSE("Cand idx: " << [&]() {
    std::stringstream ss;
    for (auto t : candidate) {
      ss << t << " ";
    }
    return ss.str();
  }());
  ACTS_VERBOSE("Cand r:   " << [&]() {
    std::stringstream ss;
    for (auto t : candidate) {
      ss << getR(t) << " ";
    }
    return ss.str();
  }());
  ACTS_VERBOSE("Cand phi: " << [&]() {
    std::stringstream ss;
    for (auto t : candidate) {
      ss << getPhi(t) << " ";
    }
    return ss.str();
  }());
  ACTS_VERBOSE("Cand z:   " << [&]() {
    std::stringstream ss;
    for (auto t : candidate) {
      ss << getZ(t) << " ";
    }
    return ss.str();
  }());

  auto tmpCand = candidate;
  std::ranges::sort(
      tmpCand, {}, [&](const auto &t) { return std::hypot(getR(t), getZ(t)); });

  tmpCand.erase(
      std::unique(tmpCand.begin(), tmpCand.end(),
                  [&](auto &a, auto &b) { return getR(a) == getR(b); }),
      tmpCand.end());

  if (tmpCand.size() < 3) {
    ACTS_WARNING("Not more then 3 spacepoints unique in R, skip!");
    return {};
  }

  Acts::Vector2 prevZR{getZ(tmpCand.front()), getR(tmpCand.front())};
  tmpCand.erase(std::remove_if(std::next(tmpCand.begin()), tmpCand.end(),
                               [&](auto &a) {
                                 Acts::Vector2 currentZR{getZ(a), getR(a)};
                                 if ((currentZR - prevZR).norm() <
                                     m_cfg.minSpacepointDist) {
                                   return true;
                                 } else {
                                   prevZR = currentZR;
                                   return false;
                                 }
                               }),
                tmpCand.end());

  if (tmpCand.size() < 3) {
    ACTS_WARNING(
        "Not more then 3 spacepoints remaining after minimum distance "
        "check!");
    return {};
  }

  const auto bottomGeoId = geoIds.at(tmpCand.at(0));
  if (m_cfg.stripVolumes.contains(bottomGeoId.volume())) {
    ACTS_VERBOSE("Bottom spacepoint is in strips, skip it!");
    return {};
  }

  const auto s = tmpCand.size();
  auto seed = m_cfg.buildTightSeeds
                  ? std::array{getXYZ(tmpCand.at(0)), getXYZ(tmpCand.at(1)),
                               getXYZ(tmpCand.at(2))}
                  : std::array{getXYZ(tmpCand.at(0)), getXYZ(tmpCand.at(s / 2)),
                               getXYZ(tmpCand.at(s - 1))};

  Acts::Vector3 field = Acts::Vector3::Zero();
  try {
    auto fieldRes = m_cfg.bField->getField(seed.at(0), bCache);
    if (!fieldRes.ok()) {
      ACTS_ERROR("Field lookup error: " << fieldRes.error());
      return {};
    }
    field = *fieldRes;
  } catch (std::exception &e) {
    ACTS_ERROR("Field lookup exception: " << e.what());
    return {};
  }

  if (field.norm() < m_cfg.bFieldMin) {
    ACTS_WARNING("Magnetic field at seed is too small " << field.norm());
    return {};
  }

  auto freePars =
      Acts::estimateTrackParamsFromSeed(seed[0], seed[1], seed[2], field);

  auto surface = m_cfg.tGeometry->findSurface(bottomGeoId);
  auto boundRes =
      Acts::transformFreeToBoundParameters(freePars, *surface, gctx);

  if (!boundRes.ok()) {
    ACTS_DEBUG("Skip track because of bad parameters");
    return {};
  }

  if (!Acts::isBoundVectorValid(*boundRes, true)) {
    ACTS_WARNING("Skipped seed because bound params not valid");
    return {};
  }

  auto params = Acts::BoundTrackParameters(
      surface->getSharedPtr(), *boundRes,
      Acts::estimateTrackParamCovariance(m_cfg.covCfg, *boundRes, false),
      m_cfg.partHypot);

  if (params.absoluteMomentum() > 1.e5) {
    ACTS_WARNING("Momentum estimate is " << params.absoluteMomentum());
    return {};
  }

  return params;
}

}  // namespace Acts
