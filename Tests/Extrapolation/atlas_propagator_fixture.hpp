#ifndef ACTS_ATLAS_PROPAGATOR_FIXTURE_HPP
#define ACTS_ATLAS_PROPAGATOR_FIXTURE_HPP 1

#include <memory>
#include "ACTS/Extrapolation/RungeKuttaEngine.hpp"
#include "ACTS/MagneticField/ConstantFieldSvc.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace Acts {

namespace Test {

  struct atlas_propagator_fixture
  {
    atlas_propagator_fixture()
    {
      typedef ConstantFieldSvc BField_type;
      BField_type::Config      b_config;
      b_config.field = {0, 0, 0.002};
      std::unique_ptr<BField_type> magnetic_field
          = std::make_unique<BField_type>(b_config);

      RungeKuttaEngine::Config c;
      c.fieldService  = std::move(magnetic_field);
      c.maxPathLength = 5 * units::_m;

      propagator = std::make_unique<RungeKuttaEngine>(c);

      surface = std::make_unique<CylinderSurface>(
          nullptr, 100 * units::_m, 30 * units::_m);
    }

    ~atlas_propagator_fixture() = default;

    std::unique_ptr<RungeKuttaEngine> propagator = nullptr;
    std::unique_ptr<const Surface>    surface    = nullptr;
  };

}  // namespace Test

}  // namespace Acts

#endif  // ACTS_ATLAS_PROPAGATOR_FIXTURE_HPP
