/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <boost/program_options.hpp>

#include "traccc/definitions/common.hpp"
#include "traccc/options/details/config_provider.hpp"
#include "traccc/options/details/interface.hpp"
#include "traccc/utils/seed_matching_config.hpp"

namespace traccc::opts {
class seed_matching : public interface,
                      public config_provider<seed_matching_config> {

    public:
    float m_matching_ratio = 0.5f;

    seed_matching();

    operator seed_matching_config() const override;

    std::unique_ptr<configuration_printable> as_printable() const override;
};
}  // namespace traccc::opts
