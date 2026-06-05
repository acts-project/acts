/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/options/details/interface.hpp"
#include "traccc/utils/logging.hpp"

// Boost include(s).
#include <boost/program_options.hpp>

// System include(s).
#include <functional>
#include <string_view>
#include <vector>

namespace traccc::opts {

/// Top-level propgram options for an executable
class program_options {

    public:
    /// Constructor
    program_options(
        std::string_view description,
        const std::vector<std::reference_wrapper<interface> >& options,
        int argc, char* argv[],
        std::unique_ptr<const traccc::Logger> ilogger =
            traccc::getDummyLogger().clone());

    private:
    /// Description of all program options
    boost::program_options::options_description m_desc;

};  // class program_options

}  // namespace traccc::opts
