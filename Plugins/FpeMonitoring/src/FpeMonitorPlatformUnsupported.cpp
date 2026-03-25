// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "FpeMonitorPlatform.hpp"

namespace ActsPlugins::detail {

// This backend intentionally provides a complete no-op implementation so the
// public FpeMonitor API remains available even when platform trapping support
// is missing.
bool isRuntimeSupported() {
  return false;
}

std::optional<FpeType> decodeFpeType(int signal, const siginfo_t* si, void* ctx) {
  // No signal decoding support on unsupported platforms.
  static_cast<void>(signal);
  static_cast<void>(si);
  static_cast<void>(ctx);
  return std::nullopt;
}

// Trap-control hooks are intentionally inert.
void clearPendingExceptions(int excepts) {
  static_cast<void>(excepts);
}

void enableExceptions(int excepts) {
  static_cast<void>(excepts);
}

void disableExceptions(int excepts) {
  static_cast<void>(excepts);
}

void maskTrapsInSignalContext(void* ctx, FpeType type) {
  // No context mutation possible without platform-specific register layout.
  static_cast<void>(ctx);
  static_cast<void>(type);
}

std::size_t captureStackFromSignalContext(void* ctx, void* buffer,
                                          std::size_t bufferBytes) {
  // No signal-context stack unwinding backend on unsupported platforms.
  static_cast<void>(ctx);
  static_cast<void>(buffer);
  static_cast<void>(bufferBytes);
  return 0;
}

// Keep defaults aligned with safe_dump fallback behavior.
std::size_t safeDumpSkipFrames() {
  return 2;
}

bool shouldFailFastOnUnknownSignal() {
  return false;
}

void installSignalHandlers(void (*handler)(int, siginfo_t*, void*)) {
  // Signal handler installation is intentionally disabled.
  static_cast<void>(handler);
}

}  // namespace ActsPlugins::detail
