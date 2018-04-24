// This file is part of the ACTS project.
//
// Copyright (C) 2016-2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_KALMANUPDATOR_H
#define ACTS_KALMANUPDATOR_H 1

// STL include(s)
#include <memory>

// BOOST include(s)
#include <boost/variant.hpp>

// ATS include(s)
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
/**
 * @brief update step of Kalman Filter using gain matrix formalism
 */
class GainMatrixUpdator
{
private:
  typedef std::unique_ptr<const BoundParameters> return_type;

public:
  template <typename Meas_t>
  return_type
  operator()(const Meas_t& m, const BoundParameters& pars) const
  {
    GainMatrixUpdatorImpl impl(pars);

    return boost::apply_visitor(impl, m);
  }

private:
  struct GainMatrixUpdatorImpl : public boost::static_visitor<return_type>
  {
  public:
    explicit GainMatrixUpdatorImpl(const BoundParameters& pars)
      : m_pParameters(&pars)
    {
    }

    template <typename Meas_t>
    return_type
    operator()(const Meas_t& m) const
    {
      static const ActsSymMatrixD<Acts::NGlobalPars> unit
          = ActsSymMatrixD<Acts::NGlobalPars>::Identity();

      const auto* pCov_trk = m_pParameters->covariance();
      if (!pCov_trk) return nullptr;

      const auto& H = m.projector();
      ActsMatrixD<Acts::NGlobalPars, Meas_t::size()> K = (*pCov_trk)
          * H.transpose()
          * (H * (*pCov_trk) * H.transpose() + m.covariance()).inverse();

      BoundParameters::ParVector_t newParValues
          = m_pParameters->parameters() + K * m.residual(*m_pParameters);
      BoundParameters::CovMatrix_t newCov = (unit - K * H) * (*pCov_trk);

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

}  // end of namespace Acts

#endif  // ACTS_KALMANUPDATOR_H
