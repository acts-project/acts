// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>

#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"
#include "Acts/Utilities/PdgParticle.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ActsFatras/EventData/Particle.hpp"

namespace Dataset {

namespace data = boost::unit_test::data;
using namespace Acts::UnitLiterals;

// default test material
const auto material = Acts::Test::makeBeryllium();
const Acts::MaterialProperties thinSlab(material, 1_mm);
const Acts::MaterialProperties thickSlab(material, 15_cm);

// particle identity
const auto particlePdg = data::make({
    Acts::PdgParticle::eElectron,
    Acts::PdgParticle::ePositron,
    Acts::PdgParticle::eMuon,
    Acts::PdgParticle::eAntiMuon,
});

// kinematic particle parameters
const auto momentumPhi = data::xrange(0_degree, 360_degree, 60_degree);
const auto momentumLambda = data::xrange(-45_degree, 45_degree, 15_degree);
const auto momentumAbs = data::xrange(500_MeV, 10_GeV, 500_MeV);

// combined particle dataset:
// cartesian grid over kinematic parameters, join over identity parameters
const auto particleParameters =
    particlePdg * momentumPhi * momentumLambda * momentumAbs;

// utility function to build a particle from the dataset parameters
inline ActsFatras::Particle makeParticle(Acts::PdgParticle pdg, double phi,
                                         double lambda, double p) {
  const auto id = ActsFatras::Barcode().setVertexPrimary(1).setParticle(1);
  return ActsFatras::Particle(id, pdg)
      .setPosition4(0, 0, 0, 0)
      .setDirection(std::cos(lambda) * std::cos(phi),
                    std::cos(lambda) * std::sin(phi), std::sin(lambda))
      .setMomentum(p);
}

}  // namespace Dataset
