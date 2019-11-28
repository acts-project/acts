// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @brief compute energy loss tables using the Acts implementation

#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "Acts/Material/Interactions.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Utilities/PdgParticle.hpp"
#include "Acts/Utilities/Units.hpp"

#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"

using namespace Acts::UnitLiterals;

/// Return the particle mass given a PDG particle identifier.
static float findMass(int pdg) {
  switch (pdg) {
    case Acts::PdgParticle::eElectron:
    case Acts::PdgParticle::eAntiElectron:
      return 511_keV;
    case Acts::PdgParticle::eMuon:
    case Acts::PdgParticle::eAntiMuon:
      return 105.7_MeV;
    case Acts::PdgParticle::ePionPlus:
    case Acts::PdgParticle::ePionMinus:
      return 139.6_MeV;
    case Acts::PdgParticle::eProton:
    case Acts::PdgParticle::eAntiProton:
      return 938.3_MeV;
    default:
      std::cerr << "Use zero mass for unknown pdg number " << pdg << '\n';
      return 0.0f;
  }
}

/// Return the particle charge given a PDG particle identifier.
static float findCharge(int pdg) {
  switch (pdg) {
    case Acts::PdgParticle::eElectron:
    case Acts::PdgParticle::eMuon:
    case Acts::PdgParticle::eTau:
    case Acts::PdgParticle::ePionMinus:
    case Acts::PdgParticle::eAntiProton:
      return -1_e;
    case Acts::PdgParticle::eAntiElectron:
    case Acts::PdgParticle::eAntiMuon:
    case Acts::PdgParticle::eAntiTau:
    case Acts::PdgParticle::ePionPlus:
    case Acts::PdgParticle::eProton:
      return 1_e;
    default:
      std::cerr << "Use 1e charge for unknown pdg number " << pdg << '\n';
      return 1_e;
  }
}

static constexpr int width = 11;
static constexpr int precision = 3;
static constexpr char separator = ' ';

static void printHeader(std::ostream& os, const Acts::MaterialProperties& slab,
                        int pdg, float mass, float charge) {
  os << "# material " << slab << '\n';
  os << "# particle pdg number " << pdg << '\n';
  os << "# particle mass " << mass / 1_MeV << "MeV\n";
  os << "# particle charge " << charge / 1_e << "e\n";
  os << "# particle momentum is given in GeV\n";
  os << "# tabulated energy loss is given in MeV\n";
  os << "# delta is the total energy loss\n";
  os << "# delta_ion is the energy loss due to ionisation and excitation\n";
  os << "# delta_rad is the energy loss due to radiative effects\n";
  os << "# sigma is the width of the enery loss distribution\n";
  // column names
  os << std::left;
  os << std::setw(width) << "momentum" << separator;
  os << std::setw(width) << "beta" << separator;
  os << std::setw(width) << "beta_gamma" << separator;
  os << std::setw(width) << "delta" << separator;
  os << std::setw(width) << "delta_ion" << separator;
  os << std::setw(width) << "delta_rad" << separator;
  os << std::setw(width) << "sigma" << '\n';
}

static void printLine(std::ostream& os, float mass, float momentum, float delta,
                      float deltaIon, float deltaRad, float sigma) {
  const auto energy = std::sqrt(mass * mass + momentum * momentum);
  const auto beta = momentum / energy;
  const auto betaGamma = momentum / mass;
  os << std::right << std::fixed << std::setprecision(precision);
  os << std::setw(width) << momentum / 1_GeV << separator;
  os << std::setw(width) << beta << separator;
  os << std::setw(width) << betaGamma << separator;
  os << std::setw(width) << delta / 1_MeV << separator;
  os << std::setw(width) << deltaIon / 1_MeV << separator;
  os << std::setw(width) << deltaRad / 1_MeV << separator;
  os << std::setw(width) << sigma / 1_MeV << '\n';
}

int main(int argc, char const* argv[]) {
  // handle input arguments
  if (argc != 6) {
    std::cerr << "usage: " << argv[0] << " thickness pdg p_min p_max n\n";
    std::cerr << "\n";
    std::cerr << "tabulate energy loss over a configurable momentum range.\n";
    std::cerr << "\n";
    std::cerr << "parameters:\n";
    std::cerr << "  thickness: material thickness in mm\n";
    std::cerr << "  pdg: PDG particle type identifier\n";
    std::cerr << "  p_{min/max}: momentum range in GeV\n";
    std::cerr << "  n: number of points in the momentum range\n";
    return EXIT_FAILURE;
  }
  const auto thickness = atof(argv[1]) * 1_mm;
  const auto pdg = atoi(argv[2]);
  const auto mass = findMass(pdg);
  const auto charge = findCharge(pdg);
  const auto pmin = atof(argv[3]) * 1_GeV;
  const auto pmax = atof(argv[4]) * 1_GeV;
  const auto deltap = (pmax - pmin) / atoi(argv[5]);

  // use fixed material for now
  // TODO make material configurable by command line
  const auto& material = Acts::Test::makeBeryllium();
  const Acts::MaterialProperties slab(material, thickness);

  printHeader(std::cout, slab, pdg, mass, charge);
  // scan momentum
  for (auto p = pmin; p < pmax; p += deltap) {
    const auto qOverP = charge / p;

    // TODO make mean/mode configurable by command line
    const auto delta = computeEnergyLossMean(slab, pdg, mass, qOverP, charge);
    const auto deltaIon =
        Acts::computeEnergyLossBethe(slab, pdg, mass, qOverP, charge);
    const auto deltaRad =
        computeEnergyLossRadiative(slab, pdg, mass, qOverP, charge);
    const auto sigma =
        Acts::computeEnergyLossLandauSigma(slab, pdg, mass, qOverP, charge);

    printLine(std::cout, mass, p, delta, deltaIon, deltaRad, sigma);
  }

  return EXIT_SUCCESS;
}
