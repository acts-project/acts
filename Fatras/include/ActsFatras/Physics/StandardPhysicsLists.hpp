// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsFatras/Kernel/PhysicsList.hpp"
#include "ActsFatras/Kernel/Process.hpp"
#include "ActsFatras/Physics/EnergyLoss/BetheBloch.hpp"
#include "ActsFatras/Physics/EnergyLoss/BetheHeitler.hpp"
#include "ActsFatras/Physics/Scattering/Highland.hpp"
#include "ActsFatras/Selectors/KinematicCasts.hpp"
#include "ActsFatras/Selectors/PdgSelectors.hpp"
#include "ActsFatras/Selectors/SelectorHelpers.hpp"

namespace ActsFatras {
namespace detail {
/// Select electrons and positrons only.
using SelectElectronLike = AbsPdgSelector<Acts::PdgParticle::eElectron>;
/// Select particles above a minimum absolute momentum.
using SelectPMin = Min<Casts::P>;

/// Highland multiple scattering that applies everywhere.
using StandardScattering =
    Process<HighlandScattering, EveryInput, EveryParticle, EveryParticle>;
/// Ionisation/excitation energy loss with a lower p cut on output particles.
///
/// Bethe-Bloch generates no particles and the child selector has no effect.
using StandardBetheBloch =
    Process<BetheBloch, EveryInput, SelectPMin, EveryParticle>;
/// Electron Bremsstrahlung energy loss with a lower p cut on output particles.
///
/// Only applies to electrons and positions. Bethe-Heitler generates no
/// particles and the child selector has no effect.
using StandardBetheHeitler =
    Process<BetheHeitler, AsInputSelector<SelectElectronLike>, SelectPMin,
            EveryParticle>;
}  // namespace detail

/// Electro-magnetic interactions for charged particles.
///
/// Scattering must come first so it is computed with the unmodified initial
/// energy before energy loss is applied.
///
/// @warning The list has no cuts on input particle charge or kinematics, i.e.
///          it relies on the simulator to preselect relevant charged particles
///          before application.
/// @todo Bethe-Bloch does not describe electrons; add correct ionisation loss
///       descriptions for electrons.
/// @todo Bethe-Heitler is applied after energy loss and thus sees the wrong
///       input energy.
using ChargedElectroMagneticPhysicsList =
    PhysicsList<detail::StandardScattering, detail::StandardBetheBloch,
                detail::StandardBetheHeitler>;

/// Construct the standard electro-magnetic physics list for charged particles.
///
/// @param minimumAbsMomentum lower p cut on output particles
ChargedElectroMagneticPhysicsList makeChargedElectroMagneticPhysicsList(
    double minimumAbsMomentum);

}  // namespace ActsFatras
