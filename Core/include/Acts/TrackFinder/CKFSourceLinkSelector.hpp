// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <limits>

#include "Acts/EventData/SourceLinkConcept.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/TrackFinder/CombinatorialKalmanFilterError.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts {

/// @brief This struct selects those source links compatible with the given
/// track parameter against provided criteria
///
/// The selection criteria could be allowed maximum chi2
/// and allowed maximum number of source links on one surface
///
/// If there is no compatible source link, it will return the source link with
/// min chi2 and tag it as an outlier.

struct CKFSourceLinkSelector {
 public:
  /// @brief nested config struct
  ///
  /// Criteria at different detector level.
  /// The layer-level criteria has the highest priority, then volume-level
  /// criteria, and the global criteria has the lowest priority
  ///
  struct Config {
    // Criteria type down to detector volume
    using VolumeChisq = std::map<GeometryID::Value, double>;
    using VolumeNumMeas = std::map<GeometryID::Value, size_t>;
    // Criteria type down to detector layer
    using LayerChisq = std::map<GeometryID::Value, VolumeChisq>;
    using LayerNumMeas = std::map<GeometryID::Value, VolumeNumMeas>;

    // Global maximum chi2
    double maxChi2 = 10;

    // Volume-level maximum chi2
    VolumeChisq volumeMaxChi2;

    // Layer-level maximum chi2
    LayerChisq layerMaxChi2;

    // Global maximum number of source links on surface
    double maxNumSourcelinksOnSurface = std::numeric_limits<size_t>::max();

    // Volume-level maximum number of source links on surface
    VolumeNumMeas volumeMaxNumSourcelinksOnSurface;

    // Layer-level maximum number of source links on surface
    LayerNumMeas layerMaxNumSourcelinksOnSurface;
  };

  /// @brief Default constructor
  CKFSourceLinkSelector() = default;

  /// @brief Constructor with config and (non-owning) logger
  ///
  /// @param config a config instance
  /// @param logger a logger instance
  CKFSourceLinkSelector(
      Config cfg,
      std::shared_ptr<const Logger> logger = std::shared_ptr<const Logger>(
          getDefaultLogger("CKFSourceLinkSelector", Logging::INFO).release()))
      : m_config(std::move(cfg)), m_logger(std::move(logger)) {}

  /// @brief Operater that select the source links compatible with
  /// the given track parameter on a surface
  ///
  /// @tparam calibrator_t The type of calibrator
  /// @tparam source_link_t The type of source link
  ///
  /// @param calibrator The measurement calibrator
  /// @param predictedParams The predicted track parameter on a surface
  /// @param sourcelinks The pool of source links
  ///
  /// @return the compatible or outlier source link indice(s)
  template <typename calibrator_t, typename source_link_t>
  Result<std::pair<std::vector<size_t>, bool>> operator()(
      const calibrator_t& calibrator, const BoundParameters& predictedParams,
      const std::vector<source_link_t>& sourcelinks) const {
    ACTS_VERBOSE("Invoked CKFSourceLinkSelector");

    using CovMatrix_t = typename BoundParameters::CovMatrix_t;

    auto surface = &predictedParams.referenceSurface();
    // Get volume and layer ID
    auto geoID = surface->geoID();
    auto volumeID = geoID.volume();
    auto layerID = geoID.layer();

    // First check if the layer-level criteria is configured.
    // If not, check if the volume-level criteria is configured.
    // Otherwise, use world criteria
    double chi2Cutoff = std::numeric_limits<double>::max();
    size_t numSlsCutoff = std::numeric_limits<size_t>::max();
    // Get the allowed maximum chi2 on this surface
    while (true) {
      // layer-level criteria
      auto lvMaxChi2 = m_config.layerMaxChi2.find(volumeID);
      if (lvMaxChi2 != m_config.layerMaxChi2.end()) {
        auto lvlMaxChi2 = lvMaxChi2->second.find(layerID);
        if (lvlMaxChi2 != lvMaxChi2->second.end()) {
          chi2Cutoff = lvlMaxChi2->second;
          break;
        }
      }
      // volume-level criteria
      auto vvMaxChi2 = m_config.volumeMaxChi2.find(volumeID);
      if (vvMaxChi2 != m_config.volumeMaxChi2.end()) {
        chi2Cutoff = vvMaxChi2->second;
        break;
      }
      // world-level criteria
      chi2Cutoff = m_config.maxChi2;
      break;
    }

    // Get the allowed maximum number of source link candidates on this surface
    while (true) {
      // layer-level criteria
      auto lvMaxNumSls =
          m_config.layerMaxNumSourcelinksOnSurface.find(volumeID);
      if (lvMaxNumSls != m_config.layerMaxNumSourcelinksOnSurface.end()) {
        auto lvlMaxNumSls = lvMaxNumSls->second.find(layerID);
        if (lvlMaxNumSls != lvMaxNumSls->second.end()) {
          numSlsCutoff = lvlMaxNumSls->second;
          break;
        }
      }
      // volume-level criteria
      auto vvMaxNumSls =
          m_config.volumeMaxNumSourcelinksOnSurface.find(volumeID);
      if (vvMaxNumSls != m_config.volumeMaxNumSourcelinksOnSurface.end()) {
        numSlsCutoff = vvMaxNumSls->second;
        break;
      }
      // world-level criteria
      numSlsCutoff = m_config.maxNumSourcelinksOnSurface;
      break;
    }

    std::vector<std::pair<size_t, double>> candidateChi2;
    double minChi2 = std::numeric_limits<double>::max();
    size_t minIndex = 0;
    size_t index = 0;
    // Loop over all source links to select the compatible source links
    for (const auto& sourcelink : sourcelinks) {
      std::visit(
          [&](const auto& calibrated) {
            // The measurement surface should be the same as parameter surface
            assert(&calibrated.referenceSurface() == surface);

            // type of measurement
            using meas_t =
                typename std::remove_const<typename std::remove_reference<
                    decltype(calibrated)>::type>::type;
            // measurement (local) parameter vector
            using meas_par_t = typename meas_t::ParVector_t;
            // type of projection
            using projection_t = typename meas_t::Projection_t;

            // Take the projector (measurement mapping function)
            const projection_t& H = calibrated.projector();
            // Take the parameter covariance
            const CovMatrix_t& predicted_covariance =
                *predictedParams.covariance();
            // Get the residual
            meas_par_t residual = calibrated.residual(predictedParams);
            // Get the chi2
            double chi2 = (residual.transpose() *
                           ((calibrated.covariance() +
                             H * predicted_covariance * H.transpose()))
                               .inverse() *
                           residual)
                              .eval()(0, 0);

            ACTS_VERBOSE("Chi2: " << chi2
                                  << " and Chi2 criteria: " << chi2Cutoff);
            // Push the source link and tag it as measurement if satisfying
            // criteria
            if (chi2 < chi2Cutoff) {
              candidateChi2.push_back({index, chi2});
            }
            // To search for the source link with the min chisq
            if (chi2 < minChi2) {
              minChi2 = chi2;
              minIndex = index;
            }
          },
          calibrator(sourcelink, predictedParams));
      index++;
    }

    // Check the number of source links against provided criteria
    // Sort the source link candidates based on chi2
    sort(candidateChi2.begin(), candidateChi2.end(),
         [=](const std::pair<size_t, double>& achi2,
             const std::pair<size_t, double>& bchi2) {
           return achi2.second < bchi2.second;
         });
    // Only store the allowed number of source link candidates
    std::vector<size_t> candidateIndices;
    size_t nCandidates = 0;
    for (const auto& [id, chi2] : candidateChi2) {
      if (numSlsCutoff <= nCandidates) {
        break;
      }
      candidateIndices.push_back(id);
      nCandidates++;
    }

    ACTS_VERBOSE(
        "Number of measurement candidates: " << candidateIndices.size());

    // If there is no selected source link, return the source link with the best
    // chisq and tag it as an outlier
    bool isOutlier = false;
    if (index > 0 and candidateIndices.empty()) {
      candidateIndices.push_back(minIndex);
      isOutlier = true;
      ACTS_DEBUG("No measurement candidate. Return an outlier source link.");
    }

    // Return error if neither source link candidates nor outlier
    if (candidateIndices.empty()) {
      return CombinatorialKalmanFilterError::SourcelinkSelectionFailed;
    }

    return std::make_pair(std::move(candidateIndices), isOutlier);
  }

  /// The config
  Config m_config;

  /// Pointer to a logger that is owned by the parent, track finder
  std::shared_ptr<const Logger> m_logger{nullptr};

  /// Getter for the logger, to support logging macros
  const Logger& logger() const {
    assert(m_logger);
    return *m_logger;
  }
};

}  // namespace Acts
