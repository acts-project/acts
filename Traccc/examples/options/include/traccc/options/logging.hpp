/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <string>

#include "traccc/options/details/config_provider.hpp"
#include "traccc/options/details/interface.hpp"
#include "traccc/utils/logging.hpp"

namespace traccc::opts {

/// Options for logging
class logging : public interface,
                public config_provider<traccc::Logging::Level> {
    public:
    logging();

    virtual operator traccc::Logging::Level() const override;

    std::unique_ptr<configuration_printable> as_printable() const override;

    private:
    int m_verbosity_incr = 0;
    int m_verbosity_decr = 0;
};  // class output_data

}  // namespace traccc::opts
