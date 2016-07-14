#ifndef ATS_KALMANUPDATOR_H
#define ATS_KALMANUPDATOR_H 1

// STL include(s)
#include <memory>

// BOOST include(s)
#include <boost/variant.hpp>

// ATS include(s)
#include "ACTS/EventData/Measurement.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {
/**
 * @brief update step of Kalman Filter using gain matrix formalism
 */
class GainMatrixUpdator
{
private:
  typedef std::unique_ptr<BoundParameters> return_type;

public:
  template <typename Meas_t>
  return_type
  operator()(const Meas_t& m, const BoundParameters& pars) const
  {
    static GainMatrixUpdatorImpl impl(pars);
    impl.setParameters(pars);

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

    void
    setParameters(const BoundParameters& pars)
    {
      m_pParameters = &pars;
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

      return std::make_unique<BoundParameters>(
          std::make_unique<BoundParameters::CovMatrix_t>(std::move(newCov)),
          newParValues,
          m_pParameters->associatedSurface());
    }

  private:
    const BoundParameters* m_pParameters;
  };
};

}  // end of namespace Acts

#endif  // ATS_KALMANUPDATOR_H
