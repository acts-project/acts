/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "traccc/options/throughput.hpp"

#include "traccc/examples/utils/printable.hpp"

// System include(s).
#include <format>

namespace traccc::opts {

/// Convenience namespace shorthand
namespace po = boost::program_options;

/// Type alias for the reconstruction stage enumeration
using stage_type = std::string;
/// Name of the reconstruction stage option
static const char* stage_option = "reco-stage";

throughput::throughput() : interface("Throughput Measurement Options") {

    m_desc.add_options()(
        stage_option, po::value<stage_type>()->default_value("full"),
        "Reconstruction stage to run (\"seeding\" or \"full\")");
    m_desc.add_options()(
        "processed-events",
        po::value(&processed_events)->default_value(processed_events),
        "Number of events to process");
    m_desc.add_options()(
        "cold-run-events",
        po::value(&cold_run_events)->default_value(cold_run_events),
        "Number of events to run 'cold'");
    m_desc.add_options()("deterministic",
                         po::value<bool>(&deterministic_event_order)
                             ->default_value(deterministic_event_order),
                         "Process events in deterministic order");
    m_desc.add_options()("random-seed",
                         po::value(&random_seed)->default_value(random_seed),
                         "Seed for event randomization (0 to use time)");
    m_desc.add_options()(
        "log-file", po::value(&log_file),
        "File where result logs will be printed (in append mode).");
}

void throughput::read(const po::variables_map& vm) {

    // Decode the input data format.
    if (vm.count(stage_option)) {
        const std::string stage_string = vm[stage_option].as<stage_type>();
        if (stage_string == "full") {
            reco_stage = stage::full;
        } else if (stage_string == "seeding") {
            reco_stage = stage::seeding;
        } else {
            throw std::invalid_argument("Unknown reconstruction stage");
        }
    }
}

std::unique_ptr<configuration_printable> throughput::as_printable() const {
    auto cat = std::make_unique<configuration_category>(m_description);

    std::string reco_stage_string;
    switch (reco_stage) {
        case stage::seeding:
            reco_stage_string = "seeding";
            break;
        case stage::full:
            reco_stage_string = "full";
            break;
        default:
            reco_stage_string = "unknown";
            break;
    }
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Reconstruction stage", reco_stage_string));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Cold run events", std::to_string(cold_run_events)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Processed events", std::to_string(processed_events)));
    cat->add_child(
        std::make_unique<configuration_kv_pair>("Log file", log_file));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Deterministic ordering",
        std::format("{}", deterministic_event_order)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Random seed",
        random_seed == 0 ? "time-based" : std::to_string(random_seed)));

    return cat;
}

}  // namespace traccc::opts
