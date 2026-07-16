/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/examples/print_fitted_tracks_statistics.hpp"

namespace traccc::details {

void print_fitted_tracks_statistics(
    const edm::track_container<default_algebra>::host& tracks,
    const Logger& log) {

    std::size_t success = 0;
    std::size_t forward = 0;
    std::size_t backward = 0;
    std::size_t propagation_forward = 0;
    std::size_t propagation_backward = 0;
    std::size_t non_positive_ndf = 0;
    std::size_t not_all_fitted = 0;
    std::size_t not_all_smoothed = 0;
    std::size_t unknown = 0;

    for (track_fit_outcome outcome : tracks.tracks.fit_outcome()) {
        if (outcome == track_fit_outcome::SUCCESS) {
            ++success;
        } else if (outcome == track_fit_outcome::FAILURE_FITTER) {
            ++forward;
        } else if (outcome == track_fit_outcome::FAILURE_SMOOTHER) {
            ++backward;
        } else if (outcome == track_fit_outcome::FAILURE_FORWARD_PROPAGATION) {
            ++propagation_forward;
        } else if (outcome == track_fit_outcome::FAILURE_BACKWARD_PROPAGATION) {
            ++propagation_backward;
        } else if (outcome == track_fit_outcome::FAILURE_NON_POSITIVE_NDF) {
            ++non_positive_ndf;
        } else if (outcome == track_fit_outcome::FAILURE_NOT_ALL_FITTED) {
            ++not_all_fitted;
        } else if (outcome == track_fit_outcome::FAILURE_NOT_ALL_SMOOTHED) {
            ++not_all_smoothed;
        } else if (outcome == track_fit_outcome::UNKNOWN) {
            ++unknown;
        }
    }

    auto logger = [&log]() -> const Logger& { return log; };
    TRACCC_INFO("Success: "
                << success << "/" << tracks.tracks.size() << " ("
                << 100. * (static_cast<double>(success) /
                           static_cast<double>(tracks.tracks.size()))
                << "%)"
                << "\n  Fit failure during forward pass: " << forward
                << "\n  Fit failure during backward pass: " << backward
                << "\n  Failure during track propagation (forward): "
                << propagation_forward
                << "\n  Failure during track propagation (backward): "
                << propagation_backward
                << "\n  Non positive NDF: " << non_positive_ndf
                << "\n  Skipped track state(s) in fit: " << not_all_fitted
                << "\n  Skipped track state(s) in smoother: "
                << not_all_smoothed << "\n  Unknown: " << unknown
                << "\n  Total errors: "
                << (forward + backward + propagation_forward +
                    propagation_backward + non_positive_ndf + not_all_fitted +
                    not_all_smoothed + unknown));
}

}  // namespace traccc::details
