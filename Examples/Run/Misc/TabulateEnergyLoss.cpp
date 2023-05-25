// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @brief compute energy loss tables using the Acts implementation

#include "Acts/Definitions/Units.hpp"
#include "Acts/Material/Interactions.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/PdgParticle.hpp"
#include "ActsFatras/Utilities/ParticleData.hpp"

#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>

using namespace Acts::UnitLiterals;

static constexpr int width = 11;
static constexpr int precision = 3;
static constexpr char separator = ' ';

static void printHeader(std::ostream& os, const Acts::MaterialSlab& slab,
                        Acts::PdgParticle pdg, float mass, float charge) {
  os << "# material: " << slab << '\n';
  os << "# particle pdg id: " << pdg << '\n';
  os << "# particle mass: " << mass / 1_MeV << "MeV\n";
  os << "# particle charge: " << charge / 1_e << "e\n";
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
  const auto energy = std::hypot(mass, momentum);
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
  const auto pdg = static_cast<Acts::PdgParticle>(atoi(argv[2]));
  const auto mass = ActsFatras::findMass(pdg);
  const auto charge = ActsFatras::findCharge(pdg);
  const auto pmin = atof(argv[3]) * 1_GeV;
  const auto pmax = atof(argv[4]) * 1_GeV;
  const auto deltap = (pmax - pmin) / atoi(argv[5]);

  // use fixed material (beryllium) for now
  // TODO make material configurable by command line
  const auto material = Acts::Material::fromMassDensity(
      35.28_cm, 42.10_cm, 9.012, 4, 1.848_g / 1_cm3);
  const Acts::MaterialSlab slab(material, thickness);

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
