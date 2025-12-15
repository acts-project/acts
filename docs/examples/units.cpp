// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

// This file contains code snippets for the Units documentation

void basicUnits() {
  //! [Using Unit Constants]
  // Define input values using unit constants
  double length = 1 * Acts::UnitConstants::m;        // length == 1000.0
  double momentum = 100 * Acts::UnitConstants::MeV;  // momentum == 0.1

  // Convert output values from native units to specific units
  double lengthInMm = length / Acts::UnitConstants::mm;  // == 1000.0
  double lengthInM = length / Acts::UnitConstants::m;    // == 1.0
  //! [Using Unit Constants]
  static_cast<void>(length);
  static_cast<void>(momentum);
  static_cast<void>(lengthInMm);
  static_cast<void>(lengthInM);
}

void userLiterals() {
  //! [Using Unit Literals]
  using namespace Acts::UnitLiterals;

  // Define input values using user-defined literals
  double length = 1_m;         // == 1000.0 (native units: mm)
  double momentum = 1.25_TeV;  // == 1250.0 (native units: GeV)

  // Convert output values to specific units
  double lengthInUm = length / 1_um;        // == 1000000.0
  double momentumInMeV = momentum / 1_MeV;  // == 1250000.0
  //! [Using Unit Literals]
  static_cast<void>(length);
  static_cast<void>(momentum);
  static_cast<void>(lengthInUm);
  static_cast<void>(momentumInMeV);
}

void comprehensiveExample() {
  //! [Comprehensive Units Example]
  using namespace Acts::UnitLiterals;

  // Define input values with units (via unit constants)
  double width = 12 * Acts::UnitConstants::mm;
  double mmuon = 105.7 * Acts::UnitConstants::MeV;

  // Define input values with units (via user literals)
  double length = 23_cm;
  double time = 1214.2_ns;
  double angle = 123_degree;
  double particleMomentum = 2.5_TeV;
  double mass = 511_keV;
  double velocity = 345_m / 1_s;
  double bfield = 3.9_T;

  // All values are now in native units and can be used in calculations
  // Native units: mm, GeV, radian, GeV/(e*mm), etc.
  //! [Comprehensive Units Example]
  static_cast<void>(width);
  static_cast<void>(mmuon);
  static_cast<void>(length);
  static_cast<void>(time);
  static_cast<void>(angle);
  static_cast<void>(particleMomentum);
  static_cast<void>(mass);
  static_cast<void>(velocity);
  static_cast<void>(bfield);
}

void outputConversion() {
  Acts::GeometryContext gctx;
  Acts::BoundTrackParameters trackPars =
      Acts::BoundTrackParameters::createCurvilinear(
          Acts::Vector4::Zero(), Acts::Vector3::UnitX(), 1, std::nullopt,
          Acts::ParticleHypothesis::pion());
  //! [Converting Output Values]
  using namespace Acts::UnitLiterals;

  // Assume trackPars is a track parameter object with methods
  // that return values in native units

  // Convert output values (via unit constants)
  double t_in_ns = trackPars.time() / Acts::UnitConstants::ns;

  // Convert output values (via user literals)
  double x_in_mm = trackPars.position(gctx)[Acts::eBoundLoc0] / 1_mm;
  double p_in_TeV = trackPars.absoluteMomentum() / 1_TeV;
  //! [Converting Output Values]
  static_cast<void>(t_in_ns);
  static_cast<void>(x_in_mm);
  static_cast<void>(p_in_TeV);
}

void physicalConstants() {
  //! [Physical Constants]
  // Speed of light in vacuum (dimensionless in native units)
  constexpr double c = Acts::PhysicalConstants::c;  // == 1.0

  // Reduced Planck constant h/(2*pi)
  // In native units: GeV*mm
  constexpr double hbar = Acts::PhysicalConstants::hbar;
  //! [Physical Constants]
  static_cast<void>(c);
  static_cast<void>(hbar);
}

void goodPractice() {
  //! [Unit Best Practices]
  using namespace Acts::UnitLiterals;

  // GOOD: All input values explicitly specify units
  double momentum = 100 * Acts::UnitConstants::GeV;  // or 100_GeV
  double distance = 5 * Acts::UnitConstants::m;      // or 5_m

  // GOOD: Variable names indicate non-native units
  double momentumInMeV = 10.0;  // would be 0.01 in native units (GeV)

  // GOOD: Unqualified values are in native units
  double energyGeV = 100.0;  // native unit GeV, clear from context

  // Output conversion
  double output_in_meters = distance / 1_m;
  //! [Unit Best Practices]
  static_cast<void>(momentum);
  static_cast<void>(momentumInMeV);
  static_cast<void>(energyGeV);
  static_cast<void>(output_in_meters);
}
