// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/IMultiAxis.hpp"

#include <memory>

namespace Acts::detail {

/// Build a 3-dimensional @c MultiAxis given an already-resolved concrete first
/// axis, dispatching the remaining two axes internally.
///
/// The three-axis factory instantiates @c MultiAxis for every combination of
/// concrete axis types (2 binning types x 3 boundary types per axis, i.e. 216
/// combinations). Instantiating all of them in a single translation unit is
/// very expensive in compiler memory. This helper is explicitly instantiated
/// per concrete first-axis type in dedicated translation units
/// (@c CreateMultiAxis3*.cpp) so that the combinatorial cost is spread across
/// several TUs instead of concentrated in one. The generated code and runtime
/// behaviour are identical to a single-TU dispatch.
///
/// @tparam FirstAxis the concrete type of the first axis
/// @param first the already-resolved first axis
/// @param second the second axis (resolved internally)
/// @param third the third axis (resolved internally)
/// @return unique pointer to the created 3D multi-axis
template <typename FirstAxis>
std::unique_ptr<IMultiAxis3D> createMultiAxis3(const FirstAxis& first,
                                               const IAxis& second,
                                               const IAxis& third);

}  // namespace Acts::detail
