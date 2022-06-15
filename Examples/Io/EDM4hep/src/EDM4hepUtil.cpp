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

ActsFatras::Hit EDM4hepUtil::fromSimHit(
    const edm4hep::SimTrackerHit& simTrackerHit,
    MapParticleIdFrom particleMapper, MapGeometryIdFrom geometryMapper) {
  ActsFatras::Barcode particleId =
      particleMapper(simTrackerHit.getMCParticle());
  Acts::GeometryIdentifier geometryId =
      geometryMapper(simTrackerHit.getCellID());

  const auto mass = simTrackerHit.getMCParticle().getMass();
  const Acts::ActsVector<3> momentum{
      simTrackerHit.getMomentum().x * Acts::UnitConstants::GeV,
      simTrackerHit.getMomentum().y * Acts::UnitConstants::GeV,
      simTrackerHit.getMomentum().z * Acts::UnitConstants::GeV,
  };
  const auto energy = std::hypot(momentum.norm(), mass);

  ActsFatras::Hit::Vector4 pos4{
      simTrackerHit.getPosition().x * Acts::UnitConstants::mm,
      simTrackerHit.getPosition().y * Acts::UnitConstants::mm,
      simTrackerHit.getPosition().z * Acts::UnitConstants::mm,
      simTrackerHit.getTime() * Acts::UnitConstants::ns,
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
    edm4hep::MCParticle particle, MapParticleIdFrom particleMapper) {
  ActsFatras::Barcode particleId = particleMapper(particle);

  ActsFatras::Particle result(particleId, Acts::PdgParticle(particle.getPDG()),
                              particle.getCharge() * Acts::UnitConstants::e,
                              particle.getMass() * Acts::UnitConstants::GeV);

  // TODO do we have that in EDM4hep?
  // particle.setProcess(static_cast<ActsFatras::ProcessType>(data.process));

  result.setPosition4(particle.getVertex()[0] * Acts::UnitConstants::mm,
                      particle.getVertex()[1] * Acts::UnitConstants::mm,
                      particle.getVertex()[2] * Acts::UnitConstants::mm,
                      particle.getTime() * Acts::UnitConstants::ns);

  // Only used for direction; normalization/units do not matter
  result.setDirection(particle.getMomentum()[0], particle.getMomentum()[1],
                      particle.getMomentum()[2]);
  result.setAbsoluteMomentum(std::hypot(particle.getMomentum()[0],
                                        particle.getMomentum()[1],
                                        particle.getMomentum()[2]) *
                             Acts::UnitConstants::GeV);

  return result;
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
