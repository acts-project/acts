#include <iostream>
#include <math.h>
#include <memory>
#include <random>
#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/EventData/Measurement.hpp"
#include "ACTS/Examples/BuildGenericDetector.hpp"
#include "ACTS/Extrapolation/ExtrapolationCell.hpp"
#include "ACTS/Extrapolation/ExtrapolationEngine.hpp"
#include "ACTS/Extrapolation/IExtrapolationEngine.hpp"
#include "ACTS/Extrapolation/MaterialEffectsEngine.hpp"
#include "ACTS/Extrapolation/RungeKuttaEngine.hpp"
#include "ACTS/Extrapolation/StaticEngine.hpp"
#include "ACTS/Extrapolation/StaticNavigationEngine.hpp"
#include "ACTS/Fitter/KalmanFitter.hpp"
#include "ACTS/Fitter/KalmanUpdator.hpp"
#include "ACTS/MagneticField/IMagneticFieldSvc.hpp"
#include "ACTS/Surfaces/PerigeeSurface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Logger.hpp"

using namespace Acts;

typedef FittableMeasurement<long int> FitMeas_t;
template <ParID_t... pars>
using Meas_t = Measurement<long int, pars...>;

class ConstantField : public IMagneticFieldSvc
{
public:
  ConstantField(double bx, double by, double bz) : m_field()
  {
    m_field[0] = bx;
    m_field[1] = by;
    m_field[2] = bz;
  }

  void
  getField(const double*, double* bxyz, double* = nullptr) const override
  {
    bxyz[0] = m_field[0];
    bxyz[1] = m_field[1];
    bxyz[2] = m_field[2];
  }

  void
  getFieldZR(const double* xyz,
             double*       bxyz,
             double*       deriv = nullptr) const override
  {
    return getField(xyz, bxyz, deriv);
  }

private:
  double m_field[3];
};

struct MyCache
{
  std::unique_ptr<const KF::Step<long int>::JacobianMatrix> jacobian;
  std::unique_ptr<const BoundParameters>                    parameters;

  MyCache()               = default;
  MyCache(const MyCache&) = delete;
  MyCache(MyCache&&)      = default;
};

class MyExtrapolator
{
public:
  MyExtrapolator(std::shared_ptr<const IExtrapolationEngine> exEngine = nullptr)
    : m_exEngine(std::move(exEngine))
  {
  }

  MyCache
  operator()(const FitMeas_t& m, const TrackParameters& tp) const
  {
    auto exCell = std::make_unique<ExtrapolationCell<TrackParameters>>(tp);
    exCell->addConfigurationMode(ExtrapolationMode::Destination);
    exCell->addConfigurationMode(ExtrapolationMode::FATRAS);
    exCell->addConfigurationMode(ExtrapolationMode::CollectJacobians);
    const Surface& sf = getSurface(m);

    ExtrapolationCode eCode = m_exEngine->extrapolate(*exCell, &sf);
    MyCache           c;
    auto j = exCell->extrapolationSteps.back().transportJacobian.release();
    c.jacobian.reset(new KF::Step<long int>::JacobianMatrix(*j));
    auto pars
        = static_cast<const BoundParameters*>(exCell->leadParameters->clone());
    c.parameters.reset(pars);

    return c;
  }

private:
  std::shared_ptr<const IExtrapolationEngine> m_exEngine;
};

class NoCalibration
{
public:
  std::unique_ptr<const FitMeas_t>
  operator()(const FitMeas_t& m, const BoundParameters&) const
  {
    return std::make_unique<const FitMeas_t>(m);
  }
};

class CacheGenerator
{
public:
  std::unique_ptr<KF::Step<long int>>
  operator()(MyCache m) const
  {
    auto step = std::make_unique<KF::Step<long int>>();
    step->setPredictedState(std::move(m.parameters));
    step->setJacobian(std::move(m.jacobian));

    return step;
  }
};

std::shared_ptr<IExtrapolationEngine>
initExtrapolator(const std::shared_ptr<const TrackingGeometry>& geo)
{
  auto propConfig         = RungeKuttaEngine::Config();
  propConfig.fieldService = std::make_shared<ConstantField>(0, 0, 0.002);
  auto propEngine         = std::make_shared<RungeKuttaEngine>(propConfig);

  auto matConfig      = MaterialEffectsEngine::Config();
  auto materialEngine = std::make_shared<MaterialEffectsEngine>(matConfig);

  auto navConfig                  = StaticNavigationEngine::Config();
  navConfig.propagationEngine     = propEngine;
  navConfig.materialEffectsEngine = materialEngine;
  navConfig.trackingGeometry      = geo;
  auto navEngine = std::make_shared<StaticNavigationEngine>(navConfig);

  auto statConfig                  = StaticEngine::Config();
  statConfig.propagationEngine     = propEngine;
  statConfig.navigationEngine      = navEngine;
  statConfig.materialEffectsEngine = materialEngine;
  auto statEngine                  = std::make_shared<StaticEngine>(statConfig);

  auto exEngineConfig                 = ExtrapolationEngine::Config();
  exEngineConfig.trackingGeometry     = geo;
  exEngineConfig.propagationEngine    = propEngine;
  exEngineConfig.navigationEngine     = navEngine;
  exEngineConfig.extrapolationEngines = {statEngine};
  auto exEngine = std::make_shared<ExtrapolationEngine>(exEngineConfig);
  exEngine->setLogger(
      getDefaultLogger("ExtrapolationEngine", Logging::VERBOSE));

  return exEngine;
}

template <typename T>
void
fit(const std::vector<FitMeas_t>& measurements, const T& fitter)
{
}

int
main()
{
  // options are stage = 0, 1, 2
  std::shared_ptr<const Acts::TrackingGeometry> geo
    = Acts::buildGenericDetector(Acts::Logging::VERBOSE,
                                 Acts::Logging::VERBOSE,
                                 Acts::Logging::VERBOSE,
                                 0);
  ActsVector<ParValue_t, NGlobalPars> pars;
  pars << 0, 0, M_PI / 2, M_PI / 2, 0.0001;
  auto startCov = std::make_unique<ActsSymMatrix<ParValue_t, NGlobalPars>>(
      ActsSymMatrix<ParValue_t, NGlobalPars>::Identity());

  const Surface* pSurf   = geo->getBeamline();
  auto           startTP = std::make_unique<BoundParameters>(
      std::move(startCov), std::move(pars), *pSurf);

  ExtrapolationCell<TrackParameters> exCell(*startTP);
  exCell.addConfigurationMode(ExtrapolationMode::CollectSensitive);
  exCell.addConfigurationMode(ExtrapolationMode::StopAtBoundary);

  auto exEngine = initExtrapolator(geo);
  exEngine->extrapolate(exCell);

  std::cout << "got " << exCell.extrapolationSteps.size()
            << " extrapolation steps" << std::endl;

  std::vector<FitMeas_t> vMeasurements;
  vMeasurements.reserve(exCell.extrapolationSteps.size());

  // identifier
  long int id = 0;
  // random numbers for smearing measurements
  std::default_random_engine             e;
  std::uniform_real_distribution<double> std_loc1(1, 5);
  std::uniform_real_distribution<double> std_loc2(0.1, 2);
  std::normal_distribution<double>       g(0, 1);

  double std1, std2, l1, l2;
  for (const auto& step : exCell.extrapolationSteps) {
    const auto& tp = step.parameters;
    if (tp->associatedSurface().type() != Surface::Cylinder) continue;

    std1 = std_loc1(e);
    std2 = std_loc2(e);
    l1   = tp->get<eLOC_1>() + std1 * g(e);
    l2   = tp->get<eLOC_2>() + std2 * g(e);
    ActsSymMatrixD<2> cov;
    cov << std1 * std1, 0, 0, std2 * std2;
    vMeasurements.push_back(Meas_t<eLOC_1, eLOC_2>(
        tp->associatedSurface(), id, std::move(cov), l1, l2));
    ++id;
  }

  std::cout << "created " << vMeasurements.size() << " pseudo-measurements"
            << std::endl;
  for (const auto& m : vMeasurements) std::cout << m << std::endl << std::endl;

  KalmanFitter<MyExtrapolator, CacheGenerator, NoCalibration, GainMatrixUpdator>
      KF;
  KF.m_oCacheGenerator = CacheGenerator();
  KF.m_oCalibrator     = NoCalibration();
  KF.m_oExtrapolator   = MyExtrapolator(exEngine);
  KF.m_oUpdator        = GainMatrixUpdator();

  std::cout << "start fit" << std::endl;
  auto track
      = KF.fit(vMeasurements, std::make_unique<BoundParameters>(*startTP));

  // dump track
  for (const auto& p : track) {
    std::cout << *p->getCalibratedMeasurement() << std::endl;
    std::cout << *p->getSmoothedState() << std::endl;
  }
    
  return 0;
}
