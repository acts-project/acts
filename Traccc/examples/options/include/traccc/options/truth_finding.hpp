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
#include "traccc/utils/truth_matching_config.hpp"

namespace traccc::opts {
class truth_finding : public interface,
                      public config_provider<truth_matching_config> {

    public:
    float m_pT_min = 0.5f * unit<float>::GeV;
    float m_z_min = -500.f * unit<float>::mm;
    float m_z_max = 500.f * unit<float>::mm;
    float m_r_max = 200.f * unit<float>::mm;
    float m_eta_max = 3.f;
    int m_process_id = 0;
    unsigned int m_min_track_candidates = 3;

    truth_finding();

    operator truth_matching_config() const override;

    void read(const boost::program_options::variables_map &) override;

    std::unique_ptr<configuration_printable> as_printable() const override;
};
}  // namespace traccc::opts
