#ifndef ACTS_SINGLETRACKPARAMETERS_H
#define ACTS_SINGLETRACKPARAMETERS_H 1

// STL include(s)
#include <type_traits>

#include "ACTS/Utilities/Definitions.h"
// ACTS includes
#include "ACTS/EventData/TrackParametersBase.h"
#include "ACTS/EventData/detail/CoordinateTransformations.h"

namespace Acts
{
  template<class ChargePolicy>
  class SingleTrackParameters : public TrackParametersBase
  {
  public:
    typedef typename TrackParametersBase::ParVector_t ParVector_t;
    typedef typename TrackParametersBase::CovMatrix_t CovMatrix_t;
    typedef std::unique_ptr<CovMatrix_t> Cov_uptr;                   ///< type for unique pointer to covariance matrix

    template<typename T = ChargePolicy,std::enable_if_t<std::is_same<T,ChargedPolicy>::value,int> = 0>
    SingleTrackParameters(Cov_uptr cov,
                          const ParVector_t& parValues,
                          const ActsVectorD<3>& position,
                          const ActsVectorD<3>& momentum):
      TrackParametersBase(),
      m_oChargePolicy(coordinate_transformation::parameters2charge(parValues)),
      m_oParameters(std::move(cov),parValues),
      m_vPosition(position),
      m_vMomentum(momentum)
    {}

    template<typename T = ChargePolicy,std::enable_if_t<std::is_same<T,NeutralPolicy>::value,int> = 0>
    SingleTrackParameters(Cov_uptr cov,
                          const ParVector_t& parValues,
                          const ActsVectorD<3>& position,
                          const ActsVectorD<3>& momentum):
      TrackParametersBase(),
      m_oChargePolicy(),
      m_oParameters(std::move(cov),parValues),
      m_vPosition(position),
      m_vMomentum(momentum)
    {}

    /**
     * @brief copy constructor
     */
    SingleTrackParameters(const SingleTrackParameters<ChargePolicy>& copy):
      TrackParametersBase(copy),
      m_oChargePolicy(copy.m_oChargePolicy),
      m_oParameters(copy.m_oParameters),
      m_vPosition(copy.m_vPosition),
      m_vMomentum(copy.m_vMomentum)
    {}

    /**
     * @brief move constructor
     */
    SingleTrackParameters(SingleTrackParameters<ChargePolicy>&& copy):
      TrackParametersBase(std::move(copy)),
      m_oChargePolicy(std::move(copy.m_oChargePolicy)),
      m_oParameters(std::move(copy.m_oParameters)),
      m_vPosition(std::move(copy.m_vPosition)),
      m_vMomentum(std::move(copy.m_vMomentum))
    {}

    virtual ~SingleTrackParameters() = default;

    virtual SingleTrackParameters<ChargePolicy>* clone() const override = 0;

    /** Access method for the position */
    virtual ActsVectorD<3> position() const final
    {
      return m_vPosition;
    }

    /** Access method for the momentum */
    virtual ActsVectorD<3> momentum() const final
    {
      return m_vMomentum;
    }

    /**
     * @brief copy assignment operator
     */
    SingleTrackParameters<ChargePolicy>& operator=(const SingleTrackParameters<ChargePolicy>& rhs)
    {
      // check for self-assignment
      if(this != &rhs)
      {
        TrackParametersBase::operator=(rhs);

        m_oChargePolicy = rhs.m_oChargePolicy;
        m_oParameters   = rhs.m_oParameters;
        m_vPosition     = rhs.m_vPosition;
        m_vMomentum     = rhs.m_vMomentum;
      }

      return *this;
    }

    /**
     * @brief move assignment operator
     */
    SingleTrackParameters<ChargePolicy>& operator=(SingleTrackParameters<ChargePolicy>&& rhs)
    {
      // check for self-assignment
      if(this != &rhs)
      {
        TrackParametersBase::operator=(std::move(rhs));

        m_oChargePolicy = std::move(rhs.m_oChargePolicy);
        m_oParameters   = std::move(rhs.m_oParameters);
        m_vPosition     = std::move(rhs.m_vPosition);
        m_vMomentum     = std::move(rhs.m_vMomentum);
      }

      return *this;
    }

    virtual bool operator==(const TrackParametersBase& rhs) const override
    {
      auto casted = dynamic_cast<decltype(this)>(&rhs);
      if(!casted)
        return false;

      return (m_oChargePolicy == casted->m_oChargePolicy && m_oParameters == casted->m_oParameters && m_vPosition == casted->m_vPosition && m_vMomentum == casted->m_vMomentum);
    }

    virtual double charge() const final
    {
      return m_oChargePolicy.getCharge();
    }

  protected:
    template<typename T>
    void updateGlobalCoordinates(const T&)
    {
      m_vMomentum = coordinate_transformation::parameters2globalMomentum(getParameterSet().getParameters());
    }

    void updateGlobalCoordinates(const local_parameter&)
    {
      m_vPosition = coordinate_transformation::parameters2globalPosition(getParameterSet().getParameters(),this->associatedSurface());
    }

    virtual FullParameterSet& getParameterSet() final
    {
      return m_oParameters;
    }

    virtual const FullParameterSet& getParameterSet() const final
    {
      return m_oParameters;
    }

    ChargePolicy      m_oChargePolicy;
    FullParameterSet  m_oParameters;
    ActsVectorD<3>     m_vPosition;
    ActsVectorD<3>     m_vMomentum;
  };
}  // end of namespace Acts

#endif // ACTS_SINGLETRACKPARAMETERS_h
