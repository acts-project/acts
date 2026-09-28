/* TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// core
#include "traccc/geometry/detector.hpp"

// io
#include "traccc/io/read_cells.hpp"
#include "traccc/io/read_detector.hpp"
#include "traccc/io/read_detector_description.hpp"
#include "traccc/io/utils.hpp"

// algorithms
#include "traccc/device/container_d2h_copy_alg.hpp"
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"
#include "traccc/sycl/clusterization/clusterization_algorithm.hpp"
#include "traccc/sycl/clusterization/measurement_sorting_algorithm.hpp"
#include "traccc/sycl/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/sycl/gbts_seeding/gbts_seeding_algorithm.hpp"
#include "traccc/sycl/seeding/seed_parameter_estimation_algorithm.hpp"
#include "traccc/sycl/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"
#include "traccc/sycl/seeding/triplet_seeding_algorithm.hpp"
#include "traccc/sycl/utils/make_magnetic_field.hpp"

// performance
#include "traccc/efficiency/seeding_performance_writer.hpp"
#include "traccc/performance/timer.hpp"

// options
#include "traccc/options/accelerator.hpp"
#include "traccc/options/clusterization.hpp"
#include "traccc/options/detector.hpp"
#include "traccc/options/input_data.hpp"
#include "traccc/options/magnetic_field.hpp"
#include "traccc/options/performance.hpp"
#include "traccc/options/program_options.hpp"
#include "traccc/options/track_finding.hpp"
#include "traccc/options/track_fitting.hpp"
#include "traccc/options/track_gbts_seeding.hpp"
#include "traccc/options/track_propagation.hpp"
#include "traccc/options/track_seeding.hpp"

// examples
#include "traccc/examples/make_magnetic_field.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/sycl/device_memory_resource.hpp>
#include <vecmem/memory/sycl/host_memory_resource.hpp>
#include <vecmem/memory/sycl/shared_memory_resource.hpp>
#include <vecmem/utils/sycl/async_copy.hpp>

// Project include(s).
#include "traccc/utils/memory_resource.hpp"

// System include(s).
#include <exception>
#include <iomanip>
#include <iostream>
#include <memory>

int seq_run(const traccc::opts::detector& detector_opts,
            const traccc::opts::magnetic_field& bfield_opts,
            const traccc::opts::input_data& input_opts,
            const traccc::opts::clusterization& clusterization_opts,
            const traccc::opts::track_seeding& seeding_opts,
            const traccc::opts::track_gbts_seeding& seeding_gbts_opts,
            const traccc::opts::track_finding& finding_opts,
            const traccc::opts::track_propagation& propagation_opts,
            const traccc::opts::track_fitting& /*fitting_opts*/,
            const traccc::opts::performance& performance_opts,
            [[maybe_unused]] const traccc::opts::accelerator& accelerator_opts,
            std::unique_ptr<const traccc::Logger> ilogger, bool usingGBTS) {
  TRACCC_LOCAL_LOGGER(std::move(ilogger));

  // Creating SYCL queue object
  vecmem::sycl::queue_wrapper vecmem_queue;
  traccc::sycl::queue_wrapper traccc_queue{vecmem_queue.queue()};
  TRACCC_INFO("Running on device: " << vecmem_queue.device_name());

  // Memory resources used by the application.
  vecmem::host_memory_resource host_mr;
  vecmem::sycl::host_memory_resource sycl_host_mr{vecmem_queue};
  vecmem::sycl::device_memory_resource device_mr{vecmem_queue};
  traccc::memory_resource mr{device_mr, &sycl_host_mr};

  // Copy object for asynchronous data transfers.
  vecmem::sycl::async_copy copy{vecmem_queue};

  // Construct the detector description object.
  traccc::detector_design_description::host host_det_descr{host_mr};
  traccc::detector_conditions_description::host host_det_cond{host_mr};
  traccc::io::read_detector_description(
      host_det_descr, host_det_cond, detector_opts.detector_file,
      detector_opts.digitization_file, detector_opts.conditions_file,
      traccc::data_format::json);
  traccc::detector_design_description::data host_det_descr_data{
      vecmem::get_data(host_det_descr)};
  traccc::detector_conditions_description::data host_det_cond_data{
      vecmem::get_data(host_det_cond)};
  traccc::detector_design_description::buffer device_det_descr{
      [&]() {
        // number of elements in the detector design description
        std::vector<unsigned int> sizes(host_det_descr.size());
        for (std::size_t i = 0; i < host_det_descr.size(); ++i) {
          auto this_design = host_det_descr.at(i);
          // now for each element, set the size to the largest size of
          // that element across all modules
          sizes[i] = std::max(
              static_cast<unsigned int>(((this_design.bin_edges_x()).size())),
              static_cast<unsigned int>(((this_design.bin_edges_y()).size())));
        }
        return sizes;
      }(),
      device_mr, &host_mr, vecmem::data::buffer_type::resizable};
  copy.setup(device_det_descr)->wait();
  copy(host_det_descr_data, device_det_descr)->wait();

  traccc::detector_conditions_description::buffer device_det_cond{
      static_cast<traccc::detector_conditions_description::buffer::size_type>(
          host_det_cond.size()),
      device_mr};
  copy.setup(device_det_cond)->wait();
  copy(host_det_cond_data, device_det_cond)->wait();

  // Construct a Detray detector object, if supported by the configuration.
  traccc::host_detector host_det;
  traccc::io::read_detector(host_det, host_mr, detector_opts.detector_file,
                            detector_opts.material_file,
                            detector_opts.grid_file);

  const traccc::detector_buffer detector_buffer =
      traccc::buffer_from_host_detector(host_det, device_mr, copy);

  // Output stats
  std::uint64_t n_cells = 0;
  // std::uint64_t n_clusters = 0;
  std::uint64_t n_spacepoints_sycl = 0;
  std::uint64_t n_seeds_sycl = 0;
  std::uint64_t n_found_tracks_sycl = 0;

  // Constant B field for the track finding and fitting
  const auto host_field = traccc::details::make_magnetic_field(bfield_opts);
  const auto device_field =
      traccc::sycl::make_magnetic_field(host_field, traccc_queue);

  // Algorithm configuration(s).
  const traccc::seedfinder_config seedfinder_config(seeding_opts);
  const traccc::seedfilter_config seedfilter_config(seeding_opts);
  const traccc::spacepoint_grid_config spacepoint_grid_config(seeding_opts);

  // GBTS seeding configuration
  const traccc::gbts_seedfinder_config gbts_config(seeding_gbts_opts);

  detray::propagation::config propagation_config(propagation_opts);

  traccc::finding_config finding_cfg(finding_opts);
  finding_cfg.propagation = propagation_config;

  // Algorithms.
  traccc::track_params_estimation_config track_params_estimation_config;

  traccc::sycl::clusterization_algorithm ca_sycl(
      mr, copy, traccc_queue, clusterization_opts,
      logger().clone("SyclClusteringAlg"));
  traccc::sycl::measurement_sorting_algorithm ms_sycl(
      mr, copy, traccc_queue, logger().clone("SyclMeasSortingAlg"));
  traccc::sycl::silicon_pixel_spacepoint_formation_algorithm sf_sycl(
      mr, copy, traccc_queue, logger().clone("SyclSpFormationAlg"));
  traccc::sycl::triplet_seeding_algorithm sa_sycl(
      seedfinder_config, spacepoint_grid_config, seedfilter_config, mr, copy,
      traccc_queue, logger().clone("SyclSeedingAlg"));
  traccc::sycl::gbts_seeding_algorithm gbts_sa_sycl(
      gbts_config, mr, copy, traccc_queue,
      logger().clone("SyclGbtsSeedingAlg"));
  traccc::sycl::seed_parameter_estimation_algorithm tp_sycl(
      track_params_estimation_config, mr, copy, traccc_queue,
      logger().clone("SyclTrackParEstAlg"));
  traccc::sycl::combinatorial_kalman_filter_algorithm finding_alg_sycl{
      finding_cfg, mr, copy, traccc_queue, logger().clone("SyclFindingAlg")};

  // performance writer
  traccc::seeding_performance_writer sd_performance_writer(
      traccc::seeding_performance_writer::config{},
      logger().clone("SeedingPerformanceWriter"));

  traccc::performance::timing_info elapsedTimes;

  // Loop over events
  for (std::size_t event = input_opts.skip;
       event < input_opts.events + input_opts.skip; ++event) {
    // Instantiate SYCL containers/collections
    traccc::edm::measurement_collection::buffer measurements_sycl_buffer;
    traccc::sycl::silicon_pixel_spacepoint_formation_algorithm::output_type
        spacepoints_sycl_buffer;
    traccc::sycl::triplet_seeding_algorithm::output_type seeds_sycl_buffer;
    traccc::sycl::seed_parameter_estimation_algorithm::output_type
        params_sycl_buffer(0, *mr.host);
    traccc::sycl::combinatorial_kalman_filter_algorithm::output_type
        track_candidates_sycl_buffer;

    {
      traccc::performance::timer wall_t("Wall time", elapsedTimes);

      traccc::edm::silicon_cell_collection::host cells_per_event{host_mr};

      {
        traccc::performance::timer t("File reading  (cpu)", elapsedTimes);
        // Read the cells from the relevant event file into host memory.
        static constexpr bool DEDUPLICATE = true;
        traccc::io::read_cells(cells_per_event, event, input_opts.directory,
                               logger().clone(), &host_det_cond,
                               input_opts.format, DEDUPLICATE,
                               input_opts.use_acts_geom_source);
      }  // stop measuring file reading timer

      n_cells += cells_per_event.size();

      // Create device copy of input collections
      traccc::edm::silicon_cell_collection::buffer cells_buffer(
          static_cast<unsigned int>(cells_per_event.size()), mr.main);
      copy(vecmem::get_data(cells_per_event), cells_buffer)->wait();

      // SYCL
      {
        traccc::performance::timer t("Clusterization (sycl)", elapsedTimes);
        // Reconstruct it into spacepoints on the device.
        auto unsorted_measurements =
            ca_sycl(cells_buffer, device_det_descr, device_det_cond);
        measurements_sycl_buffer = ms_sycl(unsorted_measurements);
        vecmem_queue.synchronize();
      }  // stop measuring clusterization sycl timer

      // Perform seeding, track finding and fitting only when using a
      // Detray geometry.
      // SYCL
      {
        traccc::performance::timer t("Spacepoint formation (sycl)",
                                     elapsedTimes);
        // Reconstruct it into spacepoints on the device.
        spacepoints_sycl_buffer =
            sf_sycl(detector_buffer, measurements_sycl_buffer);
        vecmem_queue.synchronize();
      }  // stop measuring clusterization sycl timer

      // SYCL
      {
        traccc::performance::timer t("Seeding (sycl)", elapsedTimes);
        if (usingGBTS) {
          seeds_sycl_buffer =
              gbts_sa_sycl(spacepoints_sycl_buffer, measurements_sycl_buffer);
        } else {
          seeds_sycl_buffer = sa_sycl(spacepoints_sycl_buffer);
        }
        vecmem_queue.synchronize();
      }  // stop measuring seeding sycl timer

      // SYCL
      {
        traccc::performance::timer t("Track params (sycl)", elapsedTimes);
        params_sycl_buffer =
            tp_sycl(device_field, measurements_sycl_buffer,
                    spacepoints_sycl_buffer, seeds_sycl_buffer);
        vecmem_queue.synchronize();
      }  // stop measuring track params timer

      // SYCL
      {
        traccc::performance::timer timer{"Track finding (sycl)", elapsedTimes};
        track_candidates_sycl_buffer =
            finding_alg_sycl(detector_buffer, device_field,
                             measurements_sycl_buffer, params_sycl_buffer);
        vecmem_queue.synchronize();
      }
    }  // stop measuring wall time

    traccc::edm::measurement_collection::host measurements_per_event_sycl{
        host_mr};
    traccc::edm::spacepoint_collection::host spacepoints_per_event_sycl{
        host_mr};
    traccc::edm::seed_collection::host seeds_sycl{host_mr};
    traccc::edm::track_collection<traccc::default_algebra>::host
        track_candidates_sycl{host_mr};

    copy(measurements_sycl_buffer, measurements_per_event_sycl)->wait();
    copy(spacepoints_sycl_buffer, spacepoints_per_event_sycl)->wait();
    copy(seeds_sycl_buffer, seeds_sycl)->wait();
    copy(track_candidates_sycl_buffer.tracks, track_candidates_sycl,
         vecmem::copy::type::device_to_host)
        ->wait();

    /// Statistics
    n_spacepoints_sycl += spacepoints_per_event_sycl.size();
    n_seeds_sycl += seeds_sycl.size();
    n_found_tracks_sycl += track_candidates_sycl.size();

    if (performance_opts.run) {
      traccc::event_data evt_data(input_opts.directory, event, host_mr,
                                  input_opts.use_acts_geom_source, &host_det,
                                  input_opts.format, true);

      sd_performance_writer.write(vecmem::get_data(seeds_sycl),
                                  vecmem::get_data(spacepoints_per_event_sycl),
                                  vecmem::get_data(measurements_per_event_sycl),
                                  evt_data);
    }
  }

  if (performance_opts.run) {
    sd_performance_writer.finalize();
  }

  TRACCC_INFO("==> Statistics ... ");
  TRACCC_INFO("- read    " << n_cells << " cells");
  TRACCC_INFO("- created (sycl) " << n_spacepoints_sycl << " spacepoints     ");

  TRACCC_INFO("- created (sycl) " << n_seeds_sycl << " seeds");
  TRACCC_INFO("- found (sycl)   " << n_found_tracks_sycl << " tracks");
  TRACCC_INFO("==>Elapsed times...\n" << elapsedTimes);

  return 0;
}
//
// The main routine
//
int main(int argc, char* argv[]) {
  std::unique_ptr<const traccc::Logger> logger = traccc::getDefaultLogger(
      "TracccExampleSeqSycl", traccc::Logging::Level::INFO);

  // Program options.
  traccc::opts::detector detector_opts;
  traccc::opts::magnetic_field bfield_opts;
  traccc::opts::input_data input_opts;
  traccc::opts::clusterization clusterization_opts;
  traccc::opts::track_seeding seeding_opts;
  traccc::opts::track_gbts_seeding seeding_gbts_opts;
  traccc::opts::track_finding finding_opts;
  traccc::opts::track_propagation propagation_opts;
  traccc::opts::track_fitting fitting_opts;
  traccc::opts::performance performance_opts;
  traccc::opts::accelerator accelerator_opts;
  traccc::opts::program_options program_opts{
      "Full Tracking Chain Using SYCL",
      {detector_opts, bfield_opts, input_opts, clusterization_opts,
       seeding_opts, seeding_gbts_opts, finding_opts, propagation_opts,
       fitting_opts, performance_opts, accelerator_opts},
      argc,
      argv,
      logger->cloneWithSuffix("Options")};

  // Run the application.
  return seq_run(detector_opts, bfield_opts, input_opts, clusterization_opts,
                 seeding_opts, seeding_gbts_opts, finding_opts,
                 propagation_opts, fitting_opts, performance_opts,
                 accelerator_opts, logger->clone(), seeding_gbts_opts.useGBTS);
}
