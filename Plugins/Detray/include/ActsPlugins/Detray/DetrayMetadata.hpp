// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <detray/definitions/algebra.hpp>
#include <detray/detectors/default_metadata.hpp>
#include <detray/detectors/odd_metadata.hpp>

namespace ActsPlugins::DetrayMetadata {

/// @addtogroup detray_plugin
/// @{

/// Metadata for the Open Data Detector
using Odd = detray::odd_metadata<detray::array<float>>;

/// Detray's generic default metadata
using Default = detray::default_metadata<detray::array<float>>;

/// @}

}  // namespace ActsPlugins::DetrayMetadata

/// X-macro over the closed set of supported detray metadata types.
///
/// Invoke with a single-argument macro; it is expanded once per metadata. This
/// is the single source of truth for the closed set and is used to generate the
/// `extern template` declarations of @c DetrayGeometryConverter::convert and the
/// corresponding explicit instantiations (one translation unit per metadata).
///
/// To add a metadata type, add the alias above and one line here, then create
/// the matching instantiation translation unit.
#define ACTS_DETRAY_METADATA_FOR_EACH(MACRO)  \
  MACRO(::ActsPlugins::DetrayMetadata::Odd)    \
  MACRO(::ActsPlugins::DetrayMetadata::Default)
