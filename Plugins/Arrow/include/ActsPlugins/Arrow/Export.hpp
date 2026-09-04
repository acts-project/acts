// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Visibility annotation for the public API of libActsPluginArrow (and the
// Examples/Io/{Parquet,Arrow} sources donated into it as OBJECT libraries).
//
// The arrow-hosting library is built with -fvisibility=hidden so that arrow's
// thousands of symbols (statically linked in) do not leak out and collide
// with pyarrow's bundled libarrow under any dlopen mode. Types, methods, and
// free functions that must be callable from other .so files (pybind module,
// downstream code) are marked ACTS_ARROW_EXPORT to override the hidden
// default.
//
// Likely to be elevated to a library-level Acts utility macro once other
// subsystems adopt the same isolation pattern.
#define ACTS_ARROW_EXPORT __attribute__((visibility("default")))
