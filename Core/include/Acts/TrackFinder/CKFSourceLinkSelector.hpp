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
/// track parameter against provided criteria.
///
/// The selection criteria could be allowed maximum chi2
/// and allowed maximum number of source links on one surface
///
/// If there is no compatible source link, the source link with the mininum chi2
/// will be selected and the status will be tagged as an outlier
///
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
    using VolumeChi2 = std::map<GeometryID::Value, double>;
    using VolumeNumMeasurements = std::map<GeometryID::Value, size_t>;
    // Criteria type down to detector layer
    using LayerChi2 = std::map<GeometryID::Value, VolumeChi2>;
    using LayerNumMeasurements =
        std::map<GeometryID::Value, VolumeNumMeasurements>;

    // Global maximum chi2
    double maxChi2 = std::numeric_limits<double>::max();

    // Volume-level maximum chi2
    VolumeChi2 volumeMaxChi2;

    // Layer-level maximum chi2
    LayerChi2 layerMaxChi2;

    // Global maximum number of source links on surface (value 1 means selecting
    // only the most compatible source link)
    size_t maxNumSourcelinksOnSurface = 1;

    // Volume-level maximum number of source links on surface
    VolumeNumMeasurements volumeMaxNumSourcelinksOnSurface;

    // Layer-level maximum number of source links on surface
    LayerNumMeasurements layerMaxNumSourcelinksOnSurface;
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
  /// @param sourcelinkChi2 The container for index and chi2 of intermediate
  /// source link candidates
  /// @param sourcelinkCandidateIndices The container for index of final source
  /// link candidates
  /// @param isOutlier The indicator for outlier or not
  ///
  template <typename calibrator_t, typename source_link_t>
  Result<void> operator()(
      const calibrator_t& calibrator, const BoundParameters& predictedParams,
      const std::vector<source_link_t>& sourcelinks,
      std::vector<std::pair<size_t, double>>& sourcelinkChi2,
      std::vector<size_t>& sourcelinkCandidateIndices, bool& isOutlier) const {
    ACTS_VERBOSE("Invoked CKFSourceLinkSelector");

    using CovMatrix_t = typename BoundParameters::CovMatrix_t;

    // Return error if no source link
    if (sourcelinks.empty()) {
      return CombinatorialKalmanFilterError::SourcelinkSelectionFailed;
    }

    auto surface = &predictedParams.referenceSurface();
    // Get volume and layer ID
    auto geoID = surface->geoID();
    auto volumeID = geoID.volume();
    auto layerID = geoID.layer();

    // First check if the layer-level criteria is configured.
    // If not, check if the volume-level criteria is configured.
    // Otherwise, use world criteria
    double chi2Cutoff = std::numeric_limits<double>::max();
    size_t numSourcelinksCutoff = std::numeric_limits<size_t>::max();
    // Get the allowed maximum chi2 on this surface
    while (true) {
      // layer-level criteria
      auto layerMaxChi2_volume_it = m_config.layerMaxChi2.find(volumeID);
      if (layerMaxChi2_volume_it != m_config.layerMaxChi2.end()) {
        auto layerMaxChi2_layer_it =
            layerMaxChi2_volume_it->second.find(layerID);
        if (layerMaxChi2_layer_it != layerMaxChi2_volume_it->second.end()) {
          chi2Cutoff = layerMaxChi2_layer_it->second;
          break;
        }
      }
      // volume-level criteria
      auto volumeMaxChi2_volume_it = m_config.volumeMaxChi2.find(volumeID);
      if (volumeMaxChi2_volume_it != m_config.volumeMaxChi2.end()) {
        chi2Cutoff = volumeMaxChi2_volume_it->second;
        break;
      }
      // world-level criteria
      chi2Cutoff = m_config.maxChi2;
      break;
    }
    ACTS_VERBOSE("Allowed maximum chi2: " << chi2Cutoff);

    // Get the allowed maximum number of source link candidates on this surface
    while (true) {
      // layer-level criteria
      auto layerMaxNumSourcelinks_volume_it =
          m_config.layerMaxNumSourcelinksOnSurface.find(volumeID);
      if (layerMaxNumSourcelinks_volume_it !=
          m_config.layerMaxNumSourcelinksOnSurface.end()) {
        auto layerMaxNumSourcelinks_layer_it =
            layerMaxNumSourcelinks_volume_it->second.find(layerID);
        if (layerMaxNumSourcelinks_layer_it !=
            layerMaxNumSourcelinks_volume_it->second.end()) {
          numSourcelinksCutoff = layerMaxNumSourcelinks_layer_it->second;
          break;
        }
      }
      // volume-level criteria
      auto volumeMaxNumSourcelinks_volume_it =
          m_config.volumeMaxNumSourcelinksOnSurface.find(volumeID);
      if (volumeMaxNumSourcelinks_volume_it !=
          m_config.volumeMaxNumSourcelinksOnSurface.end()) {
        numSourcelinksCutoff = volumeMaxNumSourcelinks_volume_it->second;
        break;
      }
      // world-level criteria
      numSourcelinksCutoff = m_config.maxNumSourcelinksOnSurface;
      break;
    }
    ACTS_VERBOSE(
        "Allowed maximum number of source links: " << numSourcelinksCutoff);

    sourcelinkChi2.resize(sourcelinks.size());
    double minChi2 = std::numeric_limits<double>::max();
    size_t minIndex = 0;
    size_t index = 0;
    size_t nInitialCandidates = 0;
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

            ACTS_VERBOSE("Chi2: " << chi2);
            // Push the source link index and chi2 if satisfying the criteria
            if (chi2 < chi2Cutoff) {
              sourcelinkChi2.at(nInitialCandidates) = {index, chi2};
              nInitialCandidates++;
            }
            // Search for the source link with the min chi2
            if (chi2 < minChi2) {
              minChi2 = chi2;
              minIndex = index;
            }
          },
          calibrator(sourcelink, predictedParams));
      index++;
    }

    // Get the number of source link candidates with provided constraint
    // considered
    size_t nFinalCandidates =
        std::min(nInitialCandidates, numSourcelinksCutoff);

    // Size of sourcelinkCandidates
    size_t containerSize = sourcelinkCandidateIndices.size();

    // If there is no selected source link, return the source link with the best
    // chi2 and tag it as an outlier
    if (nFinalCandidates == 0) {
      sourcelinkCandidateIndices.resize(1);
      sourcelinkCandidateIndices.at(0) = minIndex;
      ACTS_DEBUG("No measurement candidate. Return an outlier source link.");
      isOutlier = true;
      return Result<void>::success();
    }

    ACTS_VERBOSE("Number of measurement candidates: " << nFinalCandidates);
    sourcelinkCandidateIndices.resize(nFinalCandidates);
    // Sort the initial source link candidates based on chi2 in ascending order
    std::sort(sourcelinkChi2.begin(),
              sourcelinkChi2.begin() + nInitialCandidates,
              [](const std::pair<size_t, double>& lchi2,
                 const std::pair<size_t, double>& rchi2) {
                return lchi2.second < rchi2.second;
              });
    // Get only allowed number of source link candidates, i.e. nFinalCandidates,
    // from the front and reset the values in the container
    size_t nRecorded = 0;
    for (const auto& [id, chi2] : sourcelinkChi2) {
      if (nRecorded >= nFinalCandidates) {
        break;
      }
      sourcelinkCandidateIndices.at(nRecorded) = id;
      nRecorded++;
    }
    isOutlier = false;
    return Result<void>::success();
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
