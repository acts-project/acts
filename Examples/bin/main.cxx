#include <iostream>
#include <math.h>
#include <memory>
#include "ACTS/Detector/TrackingGeometry.hpp"
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
//  exCell.addConfigurationMode(ExtrapolationMode::CollectSensitive);
  exCell.addConfigurationMode(ExtrapolationMode::StopAtBoundary);

  auto exEngine = initExtrapolator(geo);
  exEngine->extrapolate(exCell);

  std::cout << "crossed " << exCell.extrapolationSteps.size()
            << " active sensor surfaces" << std::endl;
  for (const auto& step : exCell.extrapolationSteps)
    std::cout << *step.parameters << std::endl;

  return 0;
}
