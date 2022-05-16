// This file is part of the Acts project.
//
// Copyright (C) 2017-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/EDM4hep/EDM4hepConverter.hpp"

#include "Acts/Definitions/Units.hpp"

#include "edm4hep/SimTrackerHitCollection.h"

#include <fstream> // TODO remove

namespace Acts {

ActsFatras::Hit convertEDM4hepSimHit(const edm4hep::SimTrackerHit& sth) {
  std::ofstream debug("/home/andreas/debug.txt", std::ios::app);
  debug << "EDM4hep convert cell id " << sth.getCellID() << std::endl;

  const auto geometryId = sth.getCellID();
  ActsFatras::Barcode particleId;
  particleId.setParticle(sth.getMCParticle().getPDG()); // TODO or id()?

  const auto mass = sth.getMCParticle().getMass();
  const ActsVector<3> momentum{
      sth.getMomentum().x * Acts::UnitConstants::GeV,
      sth.getMomentum().y * Acts::UnitConstants::GeV,
      sth.getMomentum().z * Acts::UnitConstants::GeV,
  };
  const auto energy = std::sqrt(momentum.squaredNorm() + mass * mass);

  ActsFatras::Hit::Vector4 pos4{
      sth.getPosition().x * Acts::UnitConstants::mm,
      sth.getPosition().y * Acts::UnitConstants::mm,
      sth.getPosition().z * Acts::UnitConstants::mm,
      sth.getTime() * Acts::UnitConstants::ns,
  };
  ActsFatras::Hit::Vector4 mom4{
      momentum.x(),
      momentum.y(),
      momentum.z(),
      energy,
  };
  ActsFatras::Hit::Vector4 delta4{
      0 * Acts::UnitConstants::GeV, 0 * Acts::UnitConstants::GeV,
      0 * Acts::UnitConstants::GeV,
      0 * Acts::UnitConstants::GeV,  // sth.getEDep()
  };
  int32_t index = -1;

  return ActsFatras::Hit(geometryId, particleId, pos4, mom4, mom4 + delta4,
                         index);
}

}  // namespace Acts
