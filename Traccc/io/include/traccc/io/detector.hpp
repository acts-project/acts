/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Turn off exessive detector building debug logs (only works with gcc!)
// clang-format off
#include <detray/utils/quiet_log_start.hpp>
// clang-format on

// Detray include(s).
#include <detray/io/frontend/detector_reader.hpp>
#include <detray/io/frontend/detector_reader_config.hpp>
#include <detray/io/frontend/detector_writer.hpp>
#include <detray/io/frontend/impl/json_readers.hpp>

// clang-format off
#include <detray/utils/quiet_log_end.hpp>
// clang-format on
