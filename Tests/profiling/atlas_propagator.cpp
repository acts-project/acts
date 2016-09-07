#include <iostream>
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/RungeKuttaEngine.hpp"
#include "ACTS/MagneticField/ConstantFieldSvc.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Utilities/Logger.hpp"
#include "ACTS/Utilities/Units.hpp"

using namespace Acts;

int
main()
{
  typedef ConstantFieldSvc BField_type;
  BField_type::Config      b_config;
  b_config.field = {0, 0, 0.002};
  std::unique_ptr<BField_type> magnetic_field
      = std::make_unique<BField_type>(b_config);

  const CylinderSurface cs(nullptr, 1 * units::_m, 3 * units::_m);

  RungeKuttaEngine::Config c;
  c.fieldService  = std::move(magnetic_field);
  c.maxPathLength = 35 * units::_cm;

  RungeKuttaEngine propagator(c);

  Vector3D                           pos(0, 0, 0);
  Vector3D                           mom(1000, 0, 0);
  CurvilinearParameters              pars(nullptr, pos, mom, +1);
  ExtrapolationCell<TrackParameters> exCell(pars);
  exCell.addConfigurationMode(ExtrapolationMode::StopWithPathLimit);

  double totalPathLength = 0;
  for (unsigned int i = 0; i < 10000; ++i) {
    exCell.pathLength         = 0;
    exCell.leadParameters     = &pars;
    exCell.lastLeadParameters = 0;
    propagator.propagate(exCell, cs);
    totalPathLength += exCell.pathLength;
  }

  std::cout << totalPathLength << std::endl;
  return 0;
}
