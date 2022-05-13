#include "Acts/Plugins/EDM4hep/EDM4hepConverter.hpp"

#include "Acts/Definitions/Units.hpp"

#include "edm4hep/SimTrackerHitCollection.h"

namespace Acts {

ActsFatras::Hit convertEDM4hepSimHit(const edm4hep::SimTrackerHit& sth) {
  const auto geometryId =
      0;  // Acts::GeometryIdentifier(sth.getCellID());  // TODO
  const auto particleId =
      0;  // ActsFatras::Barcode(sth.getMCParticle().getPDG());  // TODO

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
      0 * Acts::UnitConstants::GeV,  // TODO
      0 * Acts::UnitConstants::GeV,  // TODO
      0 * Acts::UnitConstants::GeV,  // TODO
      0 * Acts::UnitConstants::GeV,  // TODO sth.getEDep()
  };
  int32_t index = -1;  // TODO

  return ActsFatras::Hit(geometryId, particleId, pos4, mom4, mom4 + delta4,
                         index);
}

}  // namespace Acts
