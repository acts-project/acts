// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/variant.hpp>
#include <memory>
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @brief Update step of Kalman Filter using gain matrix formalism
///
/// This is implemented as a boost vistor pattern for use of the
/// boost variant container
class GainMatrixUpdator
{
private:
  using return_type = std::unique_ptr<const BoundParameters>;

public:
  /// @brief Public call operator for the boost visitor pattern
  ///
  /// @tparam measurement_t Type of the measurement to be used
  /// @param m The measurement
  /// @param pars The predicted parameters
  ///
  /// @todo Include the incremental chi-2 here or do some outlier steering ?
  ///
  /// @return The updated parameters
  template <typename measurement_t>
  return_type
  operator()(const measurement_t& m, const BoundParameters& pars) const
  {
    GainMatrixUpdatorImpl impl(pars);
    return boost::apply_visitor(impl, m);
  }

private:
  struct GainMatrixUpdatorImpl : public boost::static_visitor<return_type>
  {
  public:
    /// @brief Explicit constructor of the GainMatrix updator
    ///
    /// @param pars The predicted parameters
    explicit GainMatrixUpdatorImpl(const BoundParameters& pars)
      : m_pParameters(&pars)
    {
    }

    /// @brief Private call operator for the visitor pattern
    ///
    /// @tparam measurement_t Type of the measurement to be used
    /// @param m The measurement
    ///
    /// @todo Include the incremental chi-2 here or do some outlier steering ?
    ///
    /// @return The updated parameters
    template <typename measurement_t>
    return_type
    operator()(const measurement_t& m) const
    {
      static const ActsSymMatrixD<Acts::NGlobalPars> unit
          = ActsSymMatrixD<Acts::NGlobalPars>::Identity();

      const auto* pCov_trk = m_pParameters->covariance();
      if (!pCov_trk) return nullptr;

      // Take the projector (measurement mapping function)
      const auto& H = m.projector();
      // The Kalman gain matrix
      ActsMatrixD<Acts::NGlobalPars, measurement_t::size()> K = (*pCov_trk)
          * H.transpose()
          * (H * (*pCov_trk) * H.transpose() + m.covariance()).inverse();
      // New parameters after update
      BoundParameters::ParVector_t newParValues
          = m_pParameters->parameters() + K * m.residual(*m_pParameters);
      // New covaraincd after update
      BoundParameters::CovMatrix_t newCov = (unit - K * H) * (*pCov_trk);

      // Create a new measurement with the updated parameters and covariance
      return std::make_unique<const BoundParameters>(
          std::make_unique<const BoundParameters::CovMatrix_t>(
              std::move(newCov)),
          newParValues,
          m_pParameters->referenceSurface());
    }

  private:
    const BoundParameters* m_pParameters;
  };
};

}  // namespace Acts