// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/HashedString.hpp"

#include <concepts>

namespace Acts {

namespace detail_tpc {
using ParametersMap = Eigen::Map<BoundVector>;
using ConstParametersMap = Eigen::Map<const BoundVector>;
using CovarianceMap = Eigen::Map<BoundMatrix>;
using ConstCovarianceMap = Eigen::Map<const BoundMatrix>;

constexpr HashedString kConceptKey = hashString("TrackProxyConceptKey");
}  // namespace detail_tpc

template <typename T>
concept TrackProxyConcept = requires(const T& cv, HashedString key) {
  { T::ReadOnly } -> std::same_as<const bool&>;

  { cv.index() } -> std::same_as<TrackIndexType>;
  { cv.tipIndex() } -> std::convertible_to<TrackIndexType>;
  { cv.stemIndex() } -> std::convertible_to<TrackIndexType>;

  { cv.referenceSurface() } -> std::same_as<const Surface&>;
  { cv.hasReferenceSurface() } -> std::same_as<bool>;

  { cv.parameters() } -> std::same_as<detail_tpc::ConstParametersMap>;
  { cv.covariance() } -> std::same_as<detail_tpc::ConstCovarianceMap>;

  { cv.particleHypothesis() } -> std::same_as<ParticleHypothesis>;

  { cv.loc0() } -> std::same_as<double>;
  { cv.loc1() } -> std::same_as<double>;
  { cv.time() } -> std::same_as<double>;
  { cv.theta() } -> std::same_as<double>;
  { cv.phi() } -> std::same_as<double>;
  { cv.qOverP() } -> std::same_as<double>;
  { cv.absoluteMomentum() } -> std::same_as<double>;
  { cv.transverseMomentum() } -> std::same_as<double>;
  { cv.charge() } -> std::same_as<double>;

  { cv.nMeasurements() } -> std::convertible_to<unsigned int>;
  { cv.nHoles() } -> std::convertible_to<unsigned int>;
  { cv.nOutliers() } -> std::convertible_to<unsigned int>;
  { cv.nSharedHits() } -> std::convertible_to<unsigned int>;
  { cv.chi2() } -> std::convertible_to<float>;
  { cv.nDoF() } -> std::convertible_to<unsigned int>;

  { cv.nTrackStates() } -> std::same_as<unsigned int>;

  { cv.hasColumn(key) } -> std::same_as<bool>;
  { cv.hasColumn(detail_tpc::kConceptKey) } -> std::same_as<bool>;

  { cv.template component<int>(key) } -> std::same_as<const int&>;
  {
    cv.template component<int, detail_tpc::kConceptKey>()
  } -> std::same_as<const int&>;
};

template <typename T>
concept ConstTrackProxyConcept = TrackProxyConcept<T>;

template <typename T>
concept MutableTrackProxyConcept =
    TrackProxyConcept<T> && requires(T v, HashedString key) {
      { v.tipIndex() } -> std::same_as<TrackIndexType&>;
      { v.stemIndex() } -> std::same_as<TrackIndexType&>;

      { v.parameters() } -> std::same_as<detail_tpc::ParametersMap>;
      { v.covariance() } -> std::same_as<detail_tpc::CovarianceMap>;

      { v.nMeasurements() } -> std::same_as<unsigned int&>;
      { v.nHoles() } -> std::same_as<unsigned int&>;
      { v.nOutliers() } -> std::same_as<unsigned int&>;
      { v.nSharedHits() } -> std::same_as<unsigned int&>;
      { v.chi2() } -> std::same_as<float&>;
      { v.nDoF() } -> std::same_as<unsigned int&>;

      { v.template component<int>(key) } -> std::same_as<int&>;
      {
        v.template component<int, detail_tpc::kConceptKey>()
      } -> std::same_as<int&>;
    };

}  // namespace Acts
