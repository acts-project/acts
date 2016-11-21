// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_ATLAS_PROPAGATOR_FIXTURE_HPP
#define ACTS_ATLAS_PROPAGATOR_FIXTURE_HPP 1

#include <memory>
#include "ACTS/Extrapolation/RungeKuttaEngine.hpp"
#include "ACTS/MagneticField/ConstantBField.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace Acts {

namespace Test {

  struct atlas_propagator_fixture
  {
    atlas_propagator_fixture()
    {
      typedef ConstantBField       BField_type;
      std::unique_ptr<BField_type> magnetic_field
          = std::make_unique<BField_type>(0, 0, 2 * units::_T);

      RungeKuttaEngine<>::Config c;
      c.fieldService  = std::move(magnetic_field);
      c.maxPathLength = 5 * units::_m;

      propagator = std::make_unique<RungeKuttaEngine<>>(c);

      surface = std::make_unique<CylinderSurface>(
          nullptr, 100 * units::_m, 30 * units::_m);
    }

    ~atlas_propagator_fixture() = default;

    std::unique_ptr<RungeKuttaEngine<>> propagator = nullptr;
    std::unique_ptr<const Surface>      surface    = nullptr;
  };

}  // namespace Test

}  // namespace Acts

#endif  // ACTS_ATLAS_PROPAGATOR_FIXTURE_HPP
