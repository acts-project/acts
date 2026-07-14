/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/cuda/utils/algorithm_base.hpp"

// Project include(s).
#include "traccc/cuda/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/finding/device/combinatorial_kalman_filter_algorithm.hpp"

namespace traccc::cuda {

/// CKF track finding algorithm using CUDA
class combinatorial_kalman_filter_algorithm
    : public device::combinatorial_kalman_filter_algorithm,
      public cuda::algorithm_base {

    public:
    /// Constructor with the algorithm's configuration
    combinatorial_kalman_filter_algorithm(
        const finding_config& config, const traccc::memory_resource& mr,
        const vecmem::copy& copy, const stream_wrapper& str,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone(),
        std::unique_ptr<traccc::cuda::kalman_fitting_algorithm> kf_fitter =
            nullptr);

    private:
    /// @name Function(s) inherited from
    ///       @c traccc::device::combinatorial_kalman_filter_algorithm
    /// @{

    /// Function meant to perform sanity checks on the input data
    ///
    /// @param measurements All measurements in an event
    /// @return @c true if the input data is valid, @c false otherwise
    ///
    bool input_is_valid(const edm::measurement_collection::const_view&
                            measurements) const override;

    /// Function building the measurement ranges buffer
    ///
    /// @param det The detector object
    /// @param n_measurements The number of measurements in the event
    /// @param measurements All measurements in an event
    /// @return The measurement ranges buffer
    ///
    vecmem::data::vector_buffer<
        edm::measurement_collection::const_view::size_type>
    build_measurement_ranges_buffer(
        const detector_buffer& det,
        const edm::measurement_collection::const_view::size_type n_measurements,
        const edm::measurement_collection::const_view& measurements)
        const override;

    /// Launch the @c progressive_kalman_filter kernel
    ///
    /// @param n_threads The number of threads to launch the kernel with
    /// @param config The track finding configuration
    /// @param det The detector object
    /// @param bfield The magnetic field object
    /// @param payload The payload for the kernel
    ///
    /// @returns the type-erased surface sequence buffer
    ///
    void progressive_kalman_filter_kernel(
        unsigned int n_threads, const finding_config& config,
        const detector_buffer& det, const magnetic_field& bfield,
        const device::progressive_kalman_filter_payload& payload,
        const device::kalman_fitting_algorithm::fit_payload& smoothing_payload)
        const override;

    /// Track finding kernel launcher
    ///
    /// @param n_threads The number of threads to launch the kernel with
    /// @param config The track finding configuration
    /// @param det The detector object
    /// @param payload The payload for the kernel
    ///
    void find_tracks_kernel(
        unsigned int n_threads, const finding_config& config,
        const detector_buffer& det,
        const device::find_tracks_payload& payload) const override;

    /// @brief Track condensing kernel launcher
    ///
    /// @param n_threads The number of threads to launch the kernel with
    /// @param out_params_per_in_param Vector of output parameters per input
    ///                                parameter, coming from @c find_tracks
    /// @param params_index Vector to hold the starting index of the output
    ///                     parameters per input parameter. To be used by
    ///                     @c condense_tracks.
    /// @param payload  The payload for the kernel
    ///
    void condense_tracks_kernel(
        unsigned int n_threads,
        const vecmem::data::vector_view<const unsigned int>&
            out_params_per_in_param,
        vecmem::data::vector_view<unsigned int>& params_index,
        const device::condense_tracks_payload& payload) const override;

    /// Launch the kernel sorting the parameters before duplicate removal
    ///
    /// @param n_threads The number of threads to launch the kernel with
    /// @param payload The payload for the kernel
    ///
    void fill_finding_duplicate_removal_sort_keys_kernel(
        unsigned int n_threads,
        const device::fill_finding_duplicate_removal_sort_keys_payload& payload)
        const override;

    /// Sort the parameter IDs according to the last measurement index
    ///
    /// @param link_last_measurement The last measurement index per link
    /// @param param_ids The parameter IDs to sort
    ///
    void sort_param_ids_by_last_measurement(
        vecmem::data::vector_view<unsigned int>& link_last_measurement,
        vecmem::data::vector_view<unsigned int>& param_ids) const override;

    /// Duplicate removal kernel launcher
    ///
    /// @param n_threads The number of threads to launch the kernel with
    /// @param config The track finding configuration
    /// @param payload The payload for the kernel
    ///
    void remove_duplicates_kernel(
        unsigned int n_threads, const finding_config& config,
        const device::remove_duplicates_payload& payload) const override;

    /// Launch the @c fill_finding_propagation_sort_keys kernel
    ///
    /// @param n_threads The number of threads to launch the kernel with
    /// @param payload The payload for the kernel
    ///
    void fill_finding_propagation_sort_keys_kernel(
        unsigned int n_threads,
        const device::fill_finding_propagation_sort_keys_payload& payload)
        const override;

    /// Sort the parameter IDs according to a custom set of keys
    ///
    /// @param keys The sort keys
    /// @param param_ids The parameter IDs to sort
    ///
    void sort_param_ids_by_keys(
        vecmem::data::vector_view<device::sort_key>& keys,
        vecmem::data::vector_view<unsigned int>& param_ids) const override;

    /// Launch the @c propagate_to_next_surface kernel
    ///
    /// @param n_threads The number of threads to launch the kernel with
    /// @param config The track finding configuration
    /// @param det The detector object
    /// @param bfield The magnetic field object
    /// @param payload The payload for the kernel
    ///
    void propagate_to_next_surface_kernel(
        unsigned int n_threads, const finding_config& config,
        const detector_buffer& det, const magnetic_field& bfield,
        const device::propagate_to_next_surface_payload& payload)
        const override;

    /// Launch the @c gather_best_tips_per_measurement kernel
    ///
    /// @param n_threads The number of threads to launch the kernel with
    /// @param payload The payload for the kernel
    ///
    void gather_best_tips_per_measurement_kernel(
        unsigned int n_threads,
        const device::gather_best_tips_per_measurement_payload<default_algebra>&
            payload) const override;

    /// Launch the @c gather_measurement_votes kernel
    ///
    /// @param n_threads The number of threads to launch the kernel with
    /// @param payload The payload for the kernel
    ///
    void gather_measurement_votes_kernel(
        unsigned int n_threads,
        const device::gather_measurement_votes_payload& payload) const override;

    /// Launch the @c update_tip_length_buffer kernel
    ///
    /// @param n_threads The number of threads to launch the kernel with
    /// @param payload The payload for the kernel
    ///
    void update_tip_length_buffer_kernel(
        unsigned int n_threads,
        const device::update_tip_length_buffer_payload& payload) const override;

    /// Launch the @c build_tracks kernel
    ///
    /// @param n_threads The number of threads to launch the kernel with
    /// @param run_mbf_smoother Whether the MBF smoother was run
    /// @param calib_cfg The measurement selector calibration configuration
    /// @param payload The payload for the kernel
    ///
    void build_tracks_kernel(
        unsigned int n_threads, bool run_mbf_smoother,
        const measurement_selector::config& calib_cfg,
        const device::build_tracks_payload& payload) const override;

    /// @}

};  // class combinatorial_kalman_filter_algorithm

}  // namespace traccc::cuda
