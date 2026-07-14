/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <memory>
#include <optional>
#include <traccc/efficiency/nseed_performance_writer.hpp>
#include <traccc/efficiency/track_filter.hpp>
#include <traccc/efficiency/track_matcher.hpp>

namespace traccc {
nseed_performance_writer::nseed_performance_writer(
    const std::string& prefix, std::unique_ptr<track_filter>&& filter,
    std::unique_ptr<track_matcher>&& matcher)
    : _prefix(prefix),
      _filter(std::move(filter)),
      _matcher(std::move(matcher)) {}

void nseed_performance_writer::initialize() {
    output_seed_file.open(_prefix + "seeds.csv");
    output_track_file.open(_prefix + "tracks.csv");

    /*
     * This is a while away from proper file handling in C++, but it'll do.
     */
    assert(output_seed_file.good());
    assert(output_track_file.good());

    write_seed_header();
    write_track_header();
}

void nseed_performance_writer::finalize() {
    if (output_seed_file.is_open()) {
        output_seed_file.close();
    }

    if (output_track_file.is_open()) {
        output_track_file.close();
    }
}

void nseed_performance_writer::write_seed_header() {
    if (output_seed_file.good()) {
        output_seed_file << "event_id" << sep << "seed_id" << sep << "length"
                         << sep << "particle_id" << std::endl;
    }
}

void nseed_performance_writer::write_seed_row(
    std::size_t evt_id, std::size_t seed_id, std::size_t length,
    std::optional<std::size_t> part_id) {
    if (output_seed_file.good()) {
        output_seed_file << evt_id << sep << seed_id << sep << length << sep
                         << (part_id ? std::to_string(*part_id) : "-1")
                         << std::endl;
    }
}

void nseed_performance_writer::write_track_header() {
    if (output_track_file.good()) {
        output_track_file << "event_id" << sep << "particle_id" << sep
                          << "pass_cuts" << sep << "q" << sep << "eta" << sep
                          << "phi" << sep << "pt" << std::endl;
    }
}

void nseed_performance_writer::write_track_row(std::size_t evt_id,
                                               std::size_t part_id,
                                               bool pass_cuts, scalar q,
                                               scalar eta, scalar phi,
                                               scalar pt) {
    if (output_track_file.good()) {
        output_track_file << evt_id << sep << part_id << sep
                          << (pass_cuts ? "true" : "false") << sep << q << sep
                          << eta << sep << phi << sep << pt << std::endl;
    }
}

std::string nseed_performance_writer::generate_report_str() const {
    /*
     * Construct a stonking great string. Constructing a string of fixed format
     * like this isn't optimal, but hey it's an R&D project.
     */
    char buffer[512];
    std::string result;

    /*
     * Print a cute little header.
     */
    snprintf(buffer, 512, "==> Seed finding efficiency ...\n");
    result.append(buffer);

    /*
     * Print descriptors of the filter and matcher.
     */
    snprintf(buffer, 512, "- %-20s : %s\n", "Particle filter",
             _filter->get_name().c_str());
    result.append(buffer);
    snprintf(buffer, 512, "- %-20s : %s\n", "Particle matcher",
             _matcher->get_name().c_str());
    result.append(buffer);

    /*
     * Print some discrete statistics...
     */
    snprintf(buffer, 512, "- %-20s : %7ld\n", "Total seeds",
             (_stats.true_seeds + _stats.false_seeds));
    result.append(buffer);
    snprintf(buffer, 512, "- %-20s : %7ld\n", "True seeds", _stats.true_seeds);
    result.append(buffer);
    snprintf(buffer, 512, "- %-20s : %7ld\n", "False seeds",
             _stats.false_seeds);
    result.append(buffer);
    snprintf(buffer, 512, "- %-20s : %7ld\n", "Total tracks",
             (_stats.unmatched_tracks + _stats.matched_tracks));
    result.append(buffer);
    snprintf(buffer, 512, "- %-20s : %7ld\n", "Matched tracks",
             _stats.matched_tracks);
    result.append(buffer);
    snprintf(buffer, 512, "- %-20s : %7ld\n", "Unmatched tracks",
             _stats.unmatched_tracks);
    result.append(buffer);

    /*
     * ...and some compound statistics.
     */
    snprintf(buffer, 512, "- %-20s : %6.2f%%\n", "Precision",
             (100. * static_cast<double>(_stats.true_seeds) /
              static_cast<double>(_stats.true_seeds + _stats.false_seeds)));
    result.append(buffer);
    snprintf(buffer, 512, "- %-20s : %6.2f%%\n", "Fake rate",
             (100. * static_cast<double>(_stats.false_seeds) /
              static_cast<double>(_stats.true_seeds + _stats.false_seeds)));
    result.append(buffer);
    snprintf(
        buffer, 512, "- %-20s : %6.2f%%\n", "Recall/Efficiency",
        (100. * static_cast<double>(_stats.matched_tracks) /
         static_cast<double>(_stats.unmatched_tracks + _stats.matched_tracks)));
    result.append(buffer);

    return result;
}
}  // namespace traccc
