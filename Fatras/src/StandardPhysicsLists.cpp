// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/Physics/StandardPhysicsLists.hpp"

ActsFatras::ChargedElectroMagneticPhysicsList
ActsFatras::makeChargedElectroMagneticPhysicsList(double minimumAbsMomentum) {
  ChargedElectroMagneticPhysicsList pl;
  pl.get<detail::StandardBetheBloch>().selectOutputParticle.valMin =
      minimumAbsMomentum;
  pl.get<detail::StandardBetheHeitler>().selectOutputParticle.valMin =
      minimumAbsMomentum;
  return pl;
}
