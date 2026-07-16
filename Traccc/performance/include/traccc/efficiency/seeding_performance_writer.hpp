/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/utils/helpers.hpp"
#include "traccc/utils/messaging.hpp"
#include "traccc/utils/seed_matching_config.hpp"
#include "traccc/utils/truth_matching_config.hpp"

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/seed_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/utils/event_data.hpp"

// System include(s).
#include <map>
#include <memory>
#include <string>
#include <string_view>

namespace traccc {
namespace details {

/// Data members that should not pollute the API of
/// @c traccc::seeding_performance_writer
struct seeding_performance_writer_data;

}  // namespace details

class seeding_performance_writer : public messaging {

    public:
    /// Configuration for the tool
    struct config {

        /// Output filename.
        std::string file_path = "performance_track_seeding.root";
        /// Output file mode
        std::string file_mode = "RECREATE";

        /// Plot tool configurations.
        std::map<std::string, plot_helpers::binning> var_binning = {
            {"Eta", plot_helpers::binning("#eta", 40, -4.f, 4.f)},
            {"Phi", plot_helpers::binning("#phi", 100, -3.15f, 3.15f)},
            {"Pt", plot_helpers::binning("p_{T} [GeV/c]", 40, 0.f, 100.f)},
            {"Num", plot_helpers::binning("N", 30, -0.5f, 29.5f)}};

        /// Cut values
        truth_matching_config truth_config;
        seed_matching_config seed_truth_config;
    };

    /// Construct from configuration and log level.
    /// @param cfg The configuration
    ///
    seeding_performance_writer(const config& cfg,
                               std::unique_ptr<const traccc::Logger> logger);

    /// Destructor
    ~seeding_performance_writer();

    void write(const edm::seed_collection::const_view& seeds_view,
               const edm::spacepoint_collection::const_view& spacepoints_view,
               const edm::measurement_collection::const_view& measurements_view,
               const event_data& evt_data);

    void finalize();

    private:
    /// Configuration for the tool
    config m_cfg;

    /// Opaque data members for the class
    std::unique_ptr<details::seeding_performance_writer_data> m_data;

};  // class seeding_performance_writer

}  // namespace traccc
