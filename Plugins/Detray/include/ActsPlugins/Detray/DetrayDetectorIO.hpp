// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsPlugins/Detray/DetrayMetadata.hpp"

#include <string>
#include <utility>
#include <vector>

#include <detray/core/detector.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace ActsPlugins::detail {

/// @ingroup detray_plugin
/// @brief Per-metadata detray detector operations (internal)
///
/// These are thin, internal wrappers around the heavy detray detector
/// algorithms (consistency checking and JSON I/O). They are function templates
/// over the detray metadata type and are explicitly instantiated (and declared
/// `extern template`) for the closed set of metadata types, so that the
/// expensive detray template trees are compiled only once per metadata in a
/// dedicated translation unit instead of in every consumer (Python bindings,
/// tests, …). The definitions live in @c DetrayDetectorIO.ipp; experiment code
/// may instantiate them for a custom metadata by including that header.

/// Check the consistency of a detray detector (throws on inconsistency).
/// @param detector the detray detector to check
template <typename metadata_t>
void checkDetrayConsistency(const detray::detector<metadata_t>& detector);

/// Write a detray detector to JSON file(s).
/// @param detector the detray detector to write
/// @param names the detray volume/surface name map
/// @param path the output path for the JSON file(s)
template <typename metadata_t>
void writeDetrayJson(const detray::detector<metadata_t>& detector,
                     const detray::name_map& names, const std::string& path);

/// Read a detray detector and its name map from JSON file(s).
/// @param mr the memory resource to build the detector with
/// @param files the JSON file(s) to read
/// @return the built detector together with its name map
template <typename metadata_t>
std::pair<detray::detector<metadata_t>, detray::name_map> readDetrayDetector(
    vecmem::memory_resource& mr, const std::vector<std::string>& files);

// Suppress implicit instantiation of the operations for the closed set; the
// matching definitions are emitted in generated translation units (see
// CMakeLists.txt).
#define ACTS_DETRAY_EXTERN_IO(META)                                         \
  extern template void checkDetrayConsistency<META>(                        \
      const detray::detector<META>&);                                       \
  extern template void writeDetrayJson<META>(const detray::detector<META>&, \
                                             const detray::name_map&,       \
                                             const std::string&);           \
  extern template std::pair<detray::detector<META>, detray::name_map>       \
  readDetrayDetector<META>(vecmem::memory_resource&,                        \
                           const std::vector<std::string>&);
ACTS_DETRAY_METADATA_FOR_EACH(ACTS_DETRAY_EXTERN_IO)
#undef ACTS_DETRAY_EXTERN_IO

}  // namespace ActsPlugins::detail
