/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc::opts {
/**
 * @brief Mixin type to indicate that some set of program options can be
 * converted to some configuration type.
 *
 * @tparam Config The config type to which this can be converted
 */
template <typename Config>
class config_provider {
    public:
    using config_type = Config;

    virtual operator config_type() const = 0;
};
}  // namespace traccc::opts
