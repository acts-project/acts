#include <iostream>
#include "../Extrapolation/atlas_propagator_fixture.hpp"
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
  Test::atlas_propagator_fixture f_propagator;

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
    f_propagator.propagator->propagate(exCell, *f_propagator.surface);
    totalPathLength += exCell.pathLength;
  }

  std::cout << totalPathLength << std::endl;
  return 0;
}
