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

ActsFatras::Particle centralElectron({0_mm, 0_mm, 0_mm},
                                     {0_GeV, 1.5_GeV, 0_GeV}, massElectron,
                                     -1_e, Acts::PdgParticle::eElectron);
ActsFatras::Particle centralPositron({0_mm, 0_mm, 0_mm},
                                     {0_GeV, 1.5_GeV, 0_GeV}, massElectron, 1_e,
                                     Acts::PdgParticle::ePositron);
ActsFatras::Particle centralMuon({0_mm, 0_mm, 0_mm}, {0_GeV, 1.5_GeV, 0_GeV},
                                 massMuon, -1_e, Acts::PdgParticle::eMuon);
ActsFatras::Particle centralAntiMuon({0_mm, 0_mm, 0_mm},
                                     {0_GeV, 1.5_GeV, 0_GeV}, massMuon, 1_e,
                                     Acts::PdgParticle::eAntiMuon);
ActsFatras::Particle backwardPion({0_mm, 0_mm, -100_mm},
                                  {10_MeV, 10_MeV, -1.5_GeV}, massPion, -1_e,
                                  Acts::PdgParticle::ePionMinus);
ActsFatras::Particle centralPion({0_mm, 0_mm, 0_mm}, {0_GeV, 1.5_GeV, 0_GeV},
                                 massPion, -1_e, Acts::PdgParticle::ePionMinus);
ActsFatras::Particle forwardPion({0_mm, 0_mm, 100_mm},
                                 {10_MeV, 10_MeV, 1.5_GeV}, massPion, -1_e,
                                 Acts::PdgParticle::ePionMinus);

}  // namespace
}  // namespace Dataset
