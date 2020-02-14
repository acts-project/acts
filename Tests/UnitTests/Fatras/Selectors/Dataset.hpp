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
namespace {

using namespace Acts::UnitLiterals;

Acts::Material material = Acts::Test::makeBeryllium();
Acts::MaterialProperties thinSlab(material, 1_mm);
Acts::MaterialProperties thickSlab(material, 15_cm);

constexpr auto massElectron = 0.51099891_MeV;
constexpr auto massMuon = 105.658367_MeV;
constexpr auto massPion = 134.9766_MeV;

ActsFatras::Particle makeParticle(double z, double eta, Acts::PdgParticle pdg,
                                  double m, double q) {
  const auto id = ActsFatras::Barcode().setVertexPrimary(1).setParticle(1);
  return ActsFatras::Particle(id, pdg, m, q)
      .setPosition4(0.0, 0.0, z, 0.0)
      .setDirection(1.0 / std::cosh(eta), 0.0, std::tanh(eta))
      .setMomentum(1.5_GeV);
}

const auto centralElectron =
    makeParticle(0_mm, 0.0, Acts::PdgParticle::eElectron, massElectron, -1_e);
const auto centralPositron =
    makeParticle(0_mm, 0.0, Acts::PdgParticle::ePositron, massElectron, 1_e);
const auto centralMuon =
    makeParticle(0_mm, 0.0, Acts::PdgParticle::eMuon, massMuon, -1_e);
const auto centralAntiMuon =
    makeParticle(0_mm, 0.0, Acts::PdgParticle::eAntiMuon, massMuon, 1_e);
const auto backwardPion =
    makeParticle(-100_mm, -4.5, Acts::PdgParticle::ePionMinus, massPion, -1_e);
const auto centralPion =
    makeParticle(0_mm, 0.0, Acts::PdgParticle::ePionMinus, massPion, -1_e);
const auto forwardPion =
    makeParticle(100_mm, 4.5, Acts::PdgParticle::ePionMinus, massPion, -1_e);
const auto centralNeutron =
    makeParticle(1_mm, 0.1, Acts::PdgParticle::eNeutron, 1_GeV, 0_e);

}  // namespace
}  // namespace Dataset
