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
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts {

/// @brief This struct selects those source links compatible with the given
/// track parameter against provided criteria
///
/// The selection criteria could be allowed maximum chi2
/// and allowed maximum number of source links on one surface
///
struct CKFSourceLinkSelector {
 public:
  /// @brief nested config struct
  ///
  struct Config {
    // Criteria type down to detector volume
    using VolumeChisq = std::map<GeometryID::Value, double>;
    using VolumeNumMeas = std::map<GeometryID::Value, size_t>;
    // Criteria type down to detector layer
    using DetectorChisq = std::map<GeometryID::Value, VolumeChisq>;
    using DetectorNumMeas = std::map<GeometryID::Value, VolumeNumMeas>;

    // Global maximum chi2
    double maxChi2 = 10;

    // Volume-level maximum chi2 for detector
    DetectorChisq detectorMaxChi2;

    // Layer-level maximum chi2 for detector volume
    VolumeChisq volumeMaxChi2;

    // Global maximum number of source links on surface
    double maxNumSourcelinksOnSurface = std::numeric_limits<size_t>::max();

    // Volume-level maximum number of source links on surface for detector
    DetectorNumMeas detectorMaxNumSourcelinks;

    // Layer-level maximum number of source links on surface for detector volume
    VolumeNumMeas volumeMaxNumSourcelinks;
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
          getDefaultLogger("CKFSourceLinkSelector", Logging::WARNING)
              .release()))
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
  /// @return the compatible source links
  template <typename calibrator_t, typename source_link_t>
  std::vector<source_link_t> operator()(
      const calibrator_t& calibrator, const BoundParameters& predictedParams,
      const std::vector<source_link_t>& sourcelinks) const {
    ACTS_VERBOSE("Invoked CKFSourceLinkSelector");

    using CovMatrix_t = typename BoundParameters::CovMatrix_t;

    std::vector<source_link_t> sourcelinkCandidates;
    sourcelinkCandidates.reserve(sourcelinks.size());

    for (const auto& sourcelink : sourcelinks) {
      std::visit(
          [&](const auto& calibrated) {
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

            // Get the surface info
            auto surface = &calibrated.referenceSurface();
            if (surface and surface->associatedDetectorElement()) {
              auto geoID = surface->geoID();
              auto volumeID = geoID.volume();
              auto layerID = geoID.layer();
              // First check if the layer-level criteria is configured.
              // If not, check if the volume-level criteria is configured.
              // Otherwise, use world criteria
              double maxChi2;
              while (true) {
                auto dvMaxChi2 = m_config.detectorMaxChi2.find(volumeID);
                if (dvMaxChi2 != m_config.detectorMaxChi2.end()) {
                  auto dvlMaxChi2 = dvMaxChi2->second.find(layerID);
                  if (dvlMaxChi2 != dvMaxChi2->second.end()) {
                    maxChi2 = dvlMaxChi2->second;
                    break;
                  }
                }

                auto vvMaxChi2 = m_config.volumeMaxChi2.find(volumeID);
                if (vvMaxChi2 != m_config.volumeMaxChi2.end()) {
                  maxChi2 = vvMaxChi2->second;
                  break;
                }

                maxChi2 = m_config.maxChi2;
                break;
              }

              // Push the source link if satisfying criteria
              if (chi2 < maxChi2) {
                ACTS_VERBOSE("Chi2: " << chi2 << " is within chi2 criteria: "
                                      << maxChi2);
                sourcelinkCandidates.push_back(sourcelink);
              }
            }
          },
          calibrator(sourcelink, predictedParams));
    }

    // TODO: check the number of source links against provided criteria
    //-> sort the chi2 for source link candidates
    //-> remove the 'overflow' source links

    ACTS_VERBOSE(
        "Number of source link candidates: " << sourcelinkCandidates.size());

    return std::move(sourcelinkCandidates);
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
