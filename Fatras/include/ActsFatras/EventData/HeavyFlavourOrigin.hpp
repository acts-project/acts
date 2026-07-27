// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>
#include <iosfwd>

namespace ActsFatras {

/// Simulation outcome identifier.
///
/// Encodes hadron heavy-flavour origin classification
/// the enum values correspond to the PDG codes for the quarks
/// Bottom: decay product of a bottomed hadron (it can also be B->D->X),
///         or particle originate in the fragmentation of a bottom quark
///         depending on the searchUpToHfQuark flag
/// Charm: decay product of a charmed hadron, or particle originate
///        in the fragmentation of a charm quark depending
///        on the searchUpToHfQuark flag
/// None: no heavy-flavour origin, depending on the criterion regulated
///       by the searchUpToHfQuark flag
enum class HeavyFlavorOrigin : std::uint8_t { None = 0, Charm = 4, Bottom = 5 };

/// Print simulation outcome to output stream
/// @param os Output stream
/// @param outcome Simulation outcome to print
/// @return Output stream
std::ostream &operator<<(std::ostream &os, HeavyFlavorOrigin outcome);

}  // namespace ActsFatras
