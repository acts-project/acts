// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// @file Regression for qualified out-of-line template definitions.
#pragma once
namespace ActsAlignment {
/// Alignment mask.
enum class Mask { All };
/// Optional hierarchy.
struct Hierarchy {};
/// Alignment engine.
template <typename Fitter>
struct Alignment {
  /// Evaluate a track.
  void evaluate(const Mask& mask) const;
};
}  // namespace ActsAlignment
template <typename Fitter>
void ActsAlignment::Alignment<Fitter>::evaluate(
    const ActsAlignment::Mask& mask) const {}
