// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ATS_KALMANFITTER_H
#define ATS_KALMANFITTER_H 1

// STL include(s)
#include <list>
#include <memory>
#include <type_traits>

// ATS include(s)
#include "ACTS/EventData/Measurement.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {
class Track;

namespace KF {
  template <typename ID>
  struct Step
  {
  public:
    typedef ActsMatrixD<Acts::NGlobalPars, Acts::NGlobalPars> JacobianMatrix;
    const BoundParameters*
    getPredictedState() const
    {
      return m_pPredicted.get();
    }
    const BoundParameters*
    getFilteredState() const
    {
      return m_pFiltered.get();
    }
    const BoundParameters*
    getSmoothedState() const
    {
      return m_pSmoothed.get();
    }
    const FittableMeasurement<ID>*
    getCalibratedMeasurement() const
    {
      return m_pCalibratedMeasurement.get();
    }
    const JacobianMatrix*
    getJacobian() const
    {
      return m_pJacobian.get();
    }

    void
    setPredictedState(std::unique_ptr<const BoundParameters> newPars)
    {
      m_pPredicted = std::move(newPars);
    }
    void
    setFilteredState(std::unique_ptr<const BoundParameters> newPars)
    {
      m_pFiltered = std::move(newPars);
    }
    void
    setSmoothedState(std::unique_ptr<const BoundParameters> newPars)
    {
      m_pSmoothed = std::move(newPars);
    }
    void
    setCalibratedMeasurement(
        std::unique_ptr<const FittableMeasurement<ID>> newMeasurement)
    {
      m_pCalibratedMeasurement = std::move(newMeasurement);
    }
    void
    setJacobian(std::unique_ptr<const JacobianMatrix> newJacobian)
    {
      m_pJacobian = std::move(newJacobian);
    }

  private:
    std::unique_ptr<const BoundParameters>         m_pPredicted;
    std::unique_ptr<const BoundParameters>         m_pFiltered;
    std::unique_ptr<const BoundParameters>         m_pSmoothed;
    std::unique_ptr<const JacobianMatrix>          m_pJacobian;
    std::unique_ptr<const FittableMeasurement<ID>> m_pCalibratedMeasurement;
  };
}

template <typename Extrapolator,
          typename CacheGenerator,
          typename Calibrator,
          typename Updator>
class KalmanFitter
{
public:
  template <typename MeasurementContainer>
  auto
  fit(const MeasurementContainer&            vMeasurements,
      std::unique_ptr<const BoundParameters> pInitialPars = nullptr) const
  {
    typedef typename MeasurementContainer::value_type Meas_t;
    typedef std::result_of_t<Extrapolator(const Meas_t&,
                                          const TrackParameters&)>
                                                       ExResult;
    typedef std::result_of_t<CacheGenerator(ExResult)> StepCache;
    typedef std::list<StepCache>                       Cache;

    Cache c = forwardFilter(vMeasurements, std::move(pInitialPars));
    applySmoothing(c);

    return convertCacheToTrack(std::move(c));
  }

  template <typename MeasurementContainer>
  auto
  forwardFilter(const MeasurementContainer&            vMeasurements,
                std::unique_ptr<const BoundParameters> pInitialPars) const
  {
    // typedef to actual measurement type
    typedef typename MeasurementContainer::value_type Meas_t;
    typedef std::result_of_t<Extrapolator(const Meas_t&,
                                          const TrackParameters&)>
                                                       ExResult;
    typedef std::result_of_t<CacheGenerator(ExResult)> StepCache;
    typedef std::list<StepCache>                       Cache;

    if (not pInitialPars) {
      ActsSymMatrixD<Acts::NGlobalPars> cov;
      cov << 100, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 10, 0,
          0, 0, 0, 0, 1;
      ActsVectorD<5> parValues;
      parValues << 0, 0, 0, 0, 0.001;
      std::cout << *std::begin(vMeasurements) << std::endl;
      getSurface(*std::begin(vMeasurements));
      pInitialPars = std::make_unique<const BoundParameters>(
          std::make_unique<const ActsSymMatrixD<Acts::NGlobalPars>>(
              std::move(cov)),
          parValues,
          getSurface(*std::begin(vMeasurements)));
    }

    Cache                  c;
    const BoundParameters* pPredicted = nullptr;
    const TrackParameters* pUpdated   = pInitialPars.get();
    for (const Meas_t& m : vMeasurements) {
      StepCache step = m_oCacheGenerator(m_oExtrapolator(m, *pUpdated));

      pPredicted = step->getPredictedState();
      step->setCalibratedMeasurement(m_oCalibrator(m, *pPredicted));
      step->setFilteredState(m_oUpdator(m, *pPredicted));
      pUpdated = step->getFilteredState();
      c.push_back(std::move(step));
    }

    return c;
  }

  template <typename StepCache>
  void
  applySmoothing(std::list<StepCache>& cache) const
  {
    typedef ActsMatrixD<Acts::NGlobalPars, Acts::NGlobalPars> GMatrix;
    // smoothing update matrix
    GMatrix G;
    // smoothed parameter vector and covariance matrix
    BoundParameters::ParVector_t smoothedPars;
    BoundParameters::CovMatrix_t smoothedCov;
    // smoothed track parameters
    std::unique_ptr<const BoundParameters> pSmoothed = nullptr;

    auto it = cache.rbegin();

    // for the last measurement the filtered state and the smoothed state are
    // equal
    (*it)->setSmoothedState(std::unique_ptr<const BoundParameters>(
        (*it)->getFilteredState()->clone()));
    // remember the previous step cache and move on
    decltype(it) pLast = it++;
    // loop over the remaining caches
    for (; it != cache.rend(); ++it, ++pLast) {
      G = (*(*it)->getFilteredState()->covariance())
          * (*(*it)->getJacobian()).transpose()
          * (*(*pLast)->getPredictedState()->covariance()).inverse();
      smoothedPars = (*it)->getFilteredState()->parameters()
          + G * ((*pLast)->getSmoothedState()->parameters()
                 - (*pLast)->getPredictedState()->parameters());
      smoothedCov = *(*it)->getFilteredState()->covariance()
          - G * (*(*pLast)->getPredictedState()->covariance()
                 - *(*pLast)->getSmoothedState()->covariance())
              * G.transpose();

      // create smoothed track parameters
      pSmoothed = std::make_unique<const BoundParameters>(
          std::make_unique<const decltype(smoothedCov)>(std::move(smoothedCov)),
          smoothedPars,
          (*it)->getFilteredState()->referenceSurface());
      (*it)->setSmoothedState(std::move(pSmoothed));
    }
  }

  template <typename Cache>
  auto
  convertCacheToTrack(Cache c) const
  {
    return c;
  }

  Extrapolator   m_oExtrapolator;
  CacheGenerator m_oCacheGenerator;
  Calibrator     m_oCalibrator;
  Updator        m_oUpdator;
};

}  // end of namespace Acts

#endif  // ATS_KALMANFITTER_H
