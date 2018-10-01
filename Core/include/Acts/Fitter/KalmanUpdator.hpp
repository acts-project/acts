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
/// @tparam parameters_t Type of the parameters to be used
///
/// This is implemented as a boost vistor pattern for use of the
/// boost variant container
template <typename parameters_t>
class GainMatrixUpdator
{

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
  parameters_t
  operator()(const measurement_t& m, const parameters_t& pars) const
  {
    GainMatrixUpdatorImpl impl(pars);
    return boost::apply_visitor(impl, m);
  }

private:
  struct GainMatrixUpdatorImpl : public boost::static_visitor<parameters_t>
  {
  public:
    /// @brief Explicit constructor of the GainMatrix updator
    ///
    /// @param pars The predicted parameters
    explicit GainMatrixUpdatorImpl(const parameters_t& pars)
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
    parameters_t
    operator()(const measurement_t& m) const
    {
      static const ActsSymMatrixD<Acts::NGlobalPars> unit
          = ActsSymMatrixD<Acts::NGlobalPars>::Identity();

      const auto* pCov_trk = m_pParameters->covariance();
      // Take the projector (measurement mapping function)
      const auto& H = m.projector();

      // The Kalman gain matrix
      ActsMatrixD<Acts::NGlobalPars, measurement_t::size()> K = (*pCov_trk)
          * H.transpose()
          * (H * (*pCov_trk) * H.transpose() + m.covariance()).inverse();
      // New parameters after update
      typename parameters_t::ParVector_t newParValues
          = m_pParameters->parameters() + K * m.residual(*m_pParameters);
      // New covaraincd after update
      typename parameters_t::CovMatrix_t newCov = (unit - K * H) * (*pCov_trk);

      // Create a new measurement with the updated parameters and covariance
      return parameters_t(
          std::make_unique<const typename parameters_t::CovMatrix_t>(
              std::move(newCov)),
          newParValues,
          m_pParameters->referenceSurface());
    }

  private:
    const parameters_t* m_pParameters;
  };
};

}  // namespace Acts