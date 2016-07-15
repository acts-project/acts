#include <iostream>
#include <math.h>
#include <memory>
#include <random>
#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/EventData/Measurement.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Examples/BuildGenericDetector.hpp"
#include "ACTS/Extrapolation/ExtrapolationCell.hpp"
#include "ACTS/Extrapolation/ExtrapolationEngine.hpp"
#include "ACTS/Extrapolation/IExtrapolationEngine.hpp"
#include "ACTS/Extrapolation/MaterialEffectsEngine.hpp"
#include "ACTS/Extrapolation/RungeKuttaEngine.hpp"
#include "ACTS/Extrapolation/StaticEngine.hpp"
#include "ACTS/Extrapolation/StaticNavigationEngine.hpp"
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

std::unique_ptr<IExtrapolationEngine>
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

  return std::make_unique<ExtrapolationEngine>(exEngineConfig);
}

int
main()
{
  std::shared_ptr<const TrackingGeometry> geo(trackingGeometry(Logging::DEBUG));
  ActsVector<ParValue_t, NGlobalPars> pars;
  pars << 0, 0, M_PI / 2, M_PI / 2, 0.0001;
  const Surface* pSurf = geo->getBeamline();
  auto           startTP
      = std::make_unique<BoundParameters>(nullptr, std::move(pars), *pSurf);

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

  return 0;
}
