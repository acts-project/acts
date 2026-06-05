/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <cassert>

#include "traccc/utils/logging.hpp"

namespace traccc {
class messaging {
    public:
    explicit messaging(std::unique_ptr<const Logger> ilogger)
        : m_logger(std::move(ilogger)) {}

    messaging() = delete;

    const Logger& logger() const {
        assert(m_logger.get() != nullptr);
        return *m_logger;
    }

    private:
    std::unique_ptr<const Logger> m_logger;
};
}  // namespace traccc
