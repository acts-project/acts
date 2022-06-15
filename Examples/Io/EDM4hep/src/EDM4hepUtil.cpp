// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/EventData/SimHit.hpp"

namespace ActsExamples {

ActsFatras::Hit EDM4hepUtil::fromSimHit(const edm4hep::SimTrackerHit& from,
                                        MapParticleIdFrom particleMapper,
                                        MapGeometryIdFrom geometryMapper) {
  ActsFatras::Barcode particleId = particleMapper(from.getMCParticle());
  Acts::GeometryIdentifier geometryId = geometryMapper(from.getCellID());

  const auto mass = from.getMCParticle().getMass();
  const Acts::ActsVector<3> momentum{
      from.getMomentum().x * Acts::UnitConstants::GeV,
      from.getMomentum().y * Acts::UnitConstants::GeV,
      from.getMomentum().z * Acts::UnitConstants::GeV,
  };
  const auto energy = std::hypot(momentum.norm(), mass);

  ActsFatras::Hit::Vector4 pos4{
      from.getPosition().x * Acts::UnitConstants::mm,
      from.getPosition().y * Acts::UnitConstants::mm,
      from.getPosition().z * Acts::UnitConstants::mm,
      from.getTime() * Acts::UnitConstants::ns,
  };
  ActsFatras::Hit::Vector4 mom4{
      momentum.x(),
      momentum.y(),
      momentum.z(),
      energy,
  };

  // TODO no EDM4hep equivalent?
  ActsFatras::Hit::Vector4 delta4{
      0 * Acts::UnitConstants::GeV, 0 * Acts::UnitConstants::GeV,
      0 * Acts::UnitConstants::GeV,
      0 * Acts::UnitConstants::GeV,  // sth.getEDep()
  };
  // TODO no EDM4hep equivalent?
  int32_t index = -1;

  return ActsFatras::Hit(geometryId, particleId, pos4, mom4, mom4 + delta4,
                         index);
}

ActsFatras::Particle EDM4hepUtil::fromParticle(
    edm4hep::MCParticle from, MapParticleIdFrom particleMapper) {
  ActsFatras::Barcode particleId = particleMapper(from);

  ActsFatras::Particle to(particleId, (Acts::PdgParticle)from.getPDG(),
                          from.getCharge() * Acts::UnitConstants::e,
                          from.getMass() * Acts::UnitConstants::GeV);

  // TODO do we have that in EDM4hep?
  // particle.setProcess(static_cast<ActsFatras::ProcessType>(data.process));

  to.setPosition4(from.getVertex()[0] * Acts::UnitConstants::mm,
                  from.getVertex()[1] * Acts::UnitConstants::mm,
                  from.getVertex()[2] * Acts::UnitConstants::mm,
                  from.getTime() * Acts::UnitConstants::ns);

  // Only used for direction; normalization/units do not matter
  to.setDirection(from.getMomentum()[0], from.getMomentum()[1],
                  from.getMomentum()[2]);
  to.setAbsoluteMomentum(std::hypot(from.getMomentum()[0],
                                    from.getMomentum()[1],
                                    from.getMomentum()[2]) *
                         Acts::UnitConstants::GeV);

  return to;
}

void EDM4hepUtil::toParticle(const ActsFatras::Particle& from,
                             edm4hep::MutableMCParticle to) {
  // TODO what about particleId?

  to.setPDG(from.pdg());
  to.setCharge(from.charge() / Acts::UnitConstants::e);
  to.setMass(from.mass() / Acts::UnitConstants::GeV);
  to.setVertex({from.position().x(), from.position().y(), from.position().z()});
  to.setMomentum({(float)from.fourMomentum().x(),
                  (float)from.fourMomentum().y(),
                  (float)from.fourMomentum().z()});
}

}  // namespace ActsExamples
