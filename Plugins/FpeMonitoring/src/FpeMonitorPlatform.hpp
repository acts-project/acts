// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsPlugins/FpeMonitoring/FpeMonitor.hpp"

#include <cstddef>
#include <csignal>
#include <optional>

namespace ActsPlugins::detail {

bool isRuntimeSupported();

std::optional<FpeType> decodeFpeType(int signal, siginfo_t* si, void* ctx);

void clearPendingExceptions(int excepts);
void enableExceptions(int excepts);
void disableExceptions(int excepts);
void maskTrapsInSignalContext(void* ctx, FpeType type);

std::size_t captureStackFromSignalContext(void* ctx, void* buffer,
                                          std::size_t bufferBytes);

std::size_t safeDumpSkipFrames();

bool shouldFailFastOnUnknownSignal();

void installSignalHandlers(void (*handler)(int, siginfo_t*, void*));

}  // namespace ActsPlugins::detail
