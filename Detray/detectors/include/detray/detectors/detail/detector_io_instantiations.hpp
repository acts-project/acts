// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// X-macro over the closed set of shipped detray metadata types for which the
/// heavy IO operations (write / read / consistency check / payload assembly)
/// are explicitly instantiated.
///
/// Invoke with a single-argument macro; it is expanded once per metadata name.
/// This is the single source of truth for the closed set and is used both to
/// generate the `extern template` declarations (see @c detector_io_array.hpp)
/// and the matching explicit instantiations (one translation unit per metadata,
/// see detectors/CMakeLists.txt). The algebra plugin is supplied by the
/// including header, so the same list serves every algebra variant.
///
/// To add a metadata type: add its header and one line here, then add the same
/// name to the CMake instantiation list.
#define DETRAY_IO_METADATA_FOR_EACH(MACRO) \
  MACRO(default_metadata)                  \
  MACRO(odd_metadata)                      \
  MACRO(toy_metadata)                      \
  MACRO(telescope_metadata)                \
  MACRO(wire_chamber_metadata)
