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

void EDM4hepUtil::toSimHit(const ActsFatras::Hit& from,
                           edm4hep::MutableSimTrackerHit to,
                           MapParticleIdTo particleMapper,
                           MapGeometryIdTo geometryMapper) {
  const Acts::Vector4& globalPos4 = from.fourPosition();
  const Acts::Vector4& momentum4Before = from.momentum4Before();
  const auto delta4 = from.momentum4After() - momentum4Before;

  to.setMCParticle(particleMapper(from.particleId()));

  // TODO what about the digitization?
  to.setCellID(geometryMapper(from.geometryId()));

  to.setTime(globalPos4[Acts::eTime] / Acts::UnitConstants::ns);
  to.setPosition({
      globalPos4[Acts::ePos0] / Acts::UnitConstants::mm,
      globalPos4[Acts::ePos1] / Acts::UnitConstants::mm,
      globalPos4[Acts::ePos2] / Acts::UnitConstants::mm,
  });

  to.setMomentum({
      (float)(momentum4Before[Acts::eMom0] / Acts::UnitConstants::GeV),
      (float)(momentum4Before[Acts::eMom1] / Acts::UnitConstants::GeV),
      (float)(momentum4Before[Acts::eMom2] / Acts::UnitConstants::GeV),
  });

  to.setEDep(delta4[Acts::eEnergy] / Acts::UnitConstants::GeV);
}

void EDM4hepUtil::toMeasurement(const ActsExamples::Measurement& from,
                                edm4hep::MutableTrackerHitPlane to,
                                const Cluster* fromCluster,
                                edm4hep::TrackerHitCollection& toClusters) {
  std::visit(
      [&](const auto& m) {
        Acts::GeometryIdentifier geoId = m.sourceLink().geometryId();

        auto parameters = (m.expander() * m.parameters()).eval();

        // TODO map to DD4hep?
        to.setCellID(geoId.value());

        to.setTime(parameters[Acts::eBoundTime] / Acts::UnitConstants::ns);

        to.setU({(float)parameters[Acts::eBoundLoc0],
                 (float)parameters[Acts::eBoundLoc1]});

        // auto covariance = (m.expander() * m.covariance() *
        // m.expander().transpose()).eval();

        if (fromCluster != nullptr) {
          for (auto& c : fromCluster->channels) {
            auto toChannel = toClusters.create();
            to.addToRawHits(toChannel.getObjectID());

            // TODO map to DD4hep?
            toChannel.setCellID(to.getCellID());

            // TODO get EDM4hep fixed
            // misusing some fields to store ACTS specific information
            // don't ask ...
            toChannel.setType(c.bin[0]);
            toChannel.setQuality(c.bin[1]);
            toChannel.setTime(c.activation);
          }
        }
      },
      from);
}

}  // namespace ActsExamples
