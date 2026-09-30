/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/cuda/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"
#include "traccc/cuda/clusterization/clusterization_algorithm.hpp"
#include "traccc/cuda/clusterization/measurement_sorting_algorithm.hpp"
#include "traccc/cuda/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/cuda/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/cuda/gbts_seeding/gbts_seeding_algorithm.hpp"
#include "traccc/cuda/seeding/seed_parameter_estimation_algorithm.hpp"
#include "traccc/cuda/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"
#include "traccc/cuda/seeding/triplet_seeding_algorithm.hpp"
#include "traccc/cuda/utils/make_magnetic_field.hpp"
#include "traccc/cuda/utils/stream_wrapper.hpp"
#include "traccc/device/container_d2h_copy_alg.hpp"
#include "traccc/efficiency/seeding_performance_writer.hpp"
#include "traccc/examples/make_magnetic_field.hpp"
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/detector_buffer.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/io/read_cells.hpp"
#include "traccc/io/read_detector.hpp"
#include "traccc/io/read_detector_description.hpp"
#include "traccc/io/utils.hpp"
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
#include "traccc/options/track_resolution.hpp"
#include "traccc/options/track_seeding.hpp"
#include "traccc/performance/container_comparator.hpp"
#include "traccc/performance/timer.hpp"
#include "traccc/utils/propagation.hpp"

// VecMem include(s).
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/host_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>
#include <vecmem/utils/cuda/stream_wrapper.hpp>

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
            const traccc::opts::track_resolution& resolution_opts,
            const traccc::opts::track_fitting& fitting_opts,
            const traccc::opts::performance& performance_opts,
            const traccc::opts::accelerator& accelerator_opts,
            std::unique_ptr<const traccc::Logger> ilogger, bool usingGBTS) {
  TRACCC_LOCAL_LOGGER(std::move(ilogger));

  // Memory resources used by the application.
  vecmem::host_memory_resource host_mr;
  vecmem::cuda::host_memory_resource cuda_host_mr;
  vecmem::cuda::device_memory_resource device_mr;
  traccc::memory_resource mr{device_mr, &cuda_host_mr};

  // CUDA types used.
  vecmem::cuda::stream_wrapper vecmem_stream;
  traccc::cuda::stream_wrapper stream{vecmem_stream.stream()};
  vecmem::cuda::async_copy copy{stream.cudaStream()};

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
  traccc::host_detector host_detector;
  traccc::io::read_detector(host_detector, host_mr, detector_opts.detector_file,
                            detector_opts.material_file,
                            detector_opts.grid_file);
  const traccc::detector_buffer device_detector =
      traccc::buffer_from_host_detector(host_detector, device_mr, copy);
  stream.synchronize();

  // Output stats
  std::uint64_t n_cells = 0;
  std::uint64_t n_measurements_cuda = 0;
  std::uint64_t n_spacepoints_cuda = 0;
  std::uint64_t n_seeds_cuda = 0;
  std::uint64_t n_found_tracks_cuda = 0;
  std::uint64_t n_ambiguity_free_tracks_cuda = 0;
  std::uint64_t n_fitted_tracks_cuda = 0;

  // Type definitions
  using device_spacepoint_formation_algorithm =
      traccc::cuda::silicon_pixel_spacepoint_formation_algorithm;

  using device_finding_algorithm =
      traccc::cuda::combinatorial_kalman_filter_algorithm;

  using device_fitting_algorithm = traccc::cuda::kalman_fitting_algorithm;

  // Algorithm configuration(s).
  detray::propagation::config propagation_config(propagation_opts);

  const traccc::seedfinder_config seedfinder_config(seeding_opts);
  const traccc::seedfilter_config seedfilter_config(seeding_opts);
  const traccc::spacepoint_grid_config spacepoint_grid_config(seeding_opts);

  // GBTS seeding configuration
  traccc::gbts_seedfinder_config gbts_config(seeding_gbts_opts);

  traccc::finding_config finding_cfg(finding_opts);
  finding_cfg.propagation = propagation_config;

  traccc::cuda::greedy_ambiguity_resolution_algorithm::config_type
      resolution_config(resolution_opts);

  traccc::fitting_config fitting_cfg(fitting_opts);
  fitting_cfg.propagation = propagation_config;

  // Constant B field for the track finding and fitting
  const auto host_field = traccc::details::make_magnetic_field(bfield_opts);
  const auto device_field = traccc::cuda::make_magnetic_field(
      host_field, (accelerator_opts.use_gpu_texture_memory
                       ? traccc::cuda::magnetic_field_storage::texture_memory
                       : traccc::cuda::magnetic_field_storage::global_memory));

  traccc::track_params_estimation_config track_params_estimation_config;

  traccc::cuda::clusterization_algorithm ca_cuda(
      mr, copy, stream, clusterization_opts,
      logger().clone("CudaClusteringAlg"));
  traccc::cuda::measurement_sorting_algorithm ms_cuda(
      mr, copy, stream, logger().clone("CudaMeasSortingAlg"));
  device_spacepoint_formation_algorithm sf_cuda(
      mr, copy, stream, logger().clone("CudaSpFormationAlg"));
  traccc::cuda::triplet_seeding_algorithm sa_cuda(
      seedfinder_config, spacepoint_grid_config, seedfilter_config, mr, copy,
      stream, logger().clone("CudaSeedingAlg"));
  traccc::cuda::gbts_seeding_algorithm gbts_sa_cuda(
      gbts_config, mr, copy, stream, logger().clone("CudaGbtsSeedingAlg"));
  traccc::cuda::seed_parameter_estimation_algorithm tp_cuda(
      track_params_estimation_config, mr, copy, stream,
      logger().clone("CudaTrackParEstAlg"));
  device_finding_algorithm finding_alg_cuda(finding_cfg, mr, copy, stream,
                                            logger().clone("CudaFindingAlg"));
  traccc::cuda::greedy_ambiguity_resolution_algorithm resolution_alg_cuda(
      resolution_config, mr, copy, stream,
      logger().clone("CudaAmbiguityResolutionAlg"));
  device_fitting_algorithm fitting_alg_cuda(fitting_cfg, mr, copy, stream,
                                            logger().clone("CudaFittingAlg"));

  // performance writer
  traccc::seeding_performance_writer sd_performance_writer(
      traccc::seeding_performance_writer::config{},
      logger().clone("SeedingPerformanceWriter"));

  traccc::performance::timing_info elapsedTimes;

  // Loop over events
  for (std::size_t event = input_opts.skip;
       event < input_opts.events + input_opts.skip; ++event) {
    // Instantiate cuda containers/collections
    traccc::edm::measurement_collection::buffer measurements_cuda_buffer;
    traccc::edm::spacepoint_collection::buffer spacepoints_cuda_buffer;
    traccc::edm::seed_collection::buffer seeds_cuda_buffer;
    traccc::bound_track_parameters_collection_types::buffer params_cuda_buffer(
        0, *mr.host);
    traccc::edm::track_container<traccc::default_algebra>::buffer
        track_candidates_buffer;
    traccc::edm::track_container<traccc::default_algebra>::buffer
        res_track_candidates_buffer;
    traccc::edm::track_container<traccc::default_algebra>::buffer
        track_states_buffer;

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
      copy.setup(cells_buffer)->wait();
      copy(vecmem::get_data(cells_per_event), cells_buffer)->wait();

      // CUDA
      {
        traccc::performance::timer t("Clusterization (cuda)", elapsedTimes);
        // Reconstruct it into spacepoints on the device.
        auto unsorted_measurements =
            ca_cuda(cells_buffer, device_det_descr, device_det_cond);
        measurements_cuda_buffer = ms_cuda(unsorted_measurements);
        stream.synchronize();
      }  // stop measuring clusterization cuda timer

      // Perform seeding, track finding and fitting only when using a
      // Detray geometry.
      // CUDA
      {
        traccc::performance::timer t("Spacepoint formation (cuda)",
                                     elapsedTimes);
        spacepoints_cuda_buffer =
            sf_cuda(device_detector, measurements_cuda_buffer);
        stream.synchronize();
      }  // stop measuring spacepoint formation cuda timer

      // CUDA
      {
        traccc::performance::timer t("Seeding (cuda)", elapsedTimes);
        if (usingGBTS) {
          seeds_cuda_buffer =
              gbts_sa_cuda(spacepoints_cuda_buffer, measurements_cuda_buffer);
        } else {
          seeds_cuda_buffer = sa_cuda(spacepoints_cuda_buffer);
        }
        stream.synchronize();
      }  // stop measuring seeding cuda timer

      // CUDA
      {
        traccc::performance::timer t("Track params (cuda)", elapsedTimes);
        params_cuda_buffer =
            tp_cuda(device_field, measurements_cuda_buffer,
                    spacepoints_cuda_buffer, seeds_cuda_buffer);
        stream.synchronize();
      }  // stop measuring track params timer

      // CUDA
      {
        traccc::performance::timer timer{"Track finding (cuda)", elapsedTimes};
        track_candidates_buffer =
            finding_alg_cuda(device_detector, device_field,
                             measurements_cuda_buffer, params_cuda_buffer);
      }

      // CUDA
      {
        traccc::performance::timer timer{"Ambiguity resolution (cuda)",
                                         elapsedTimes};
        res_track_candidates_buffer =
            resolution_alg_cuda(track_candidates_buffer);
      }

      // CUDA
      {
        traccc::performance::timer timer{"Track fitting (cuda)", elapsedTimes};
        track_states_buffer = fitting_alg_cuda(device_detector, device_field,
                                               res_track_candidates_buffer);
      }

    }  // Stop measuring wall time

    traccc::edm::measurement_collection::host measurements_per_event_cuda{
        host_mr};
    traccc::edm::spacepoint_collection::host spacepoints_per_event_cuda{
        host_mr};
    traccc::edm::seed_collection::host seeds_cuda{host_mr};
    traccc::edm::track_collection<traccc::default_algebra>::host
        track_candidates_cuda{host_mr};
    traccc::edm::track_collection<traccc::default_algebra>::host
        res_track_candidates_cuda{host_mr};
    traccc::edm::track_container<traccc::default_algebra>::host
        track_states_cuda{host_mr};

    copy(measurements_cuda_buffer, measurements_per_event_cuda)->wait();
    copy(spacepoints_cuda_buffer, spacepoints_per_event_cuda)->wait();
    copy(seeds_cuda_buffer, seeds_cuda)->wait();
    copy(track_candidates_buffer.tracks, track_candidates_cuda,
         vecmem::copy::type::device_to_host)
        ->wait();
    copy(res_track_candidates_buffer.tracks, res_track_candidates_cuda,
         vecmem::copy::type::device_to_host)
        ->wait();
    copy(track_states_buffer.tracks, track_states_cuda.tracks,
         vecmem::copy::type::device_to_host)
        ->wait();
    copy(track_states_buffer.states, track_states_cuda.states,
         vecmem::copy::type::device_to_host)
        ->wait();
    track_states_cuda.measurements =
        vecmem::get_data(measurements_per_event_cuda);
    stream.synchronize();

    /// Statistics
    n_measurements_cuda += measurements_per_event_cuda.size();
    n_spacepoints_cuda += spacepoints_per_event_cuda.size();
    n_seeds_cuda += seeds_cuda.size();
    n_found_tracks_cuda += track_candidates_cuda.size();
    n_ambiguity_free_tracks_cuda += res_track_candidates_cuda.size();
    n_fitted_tracks_cuda += track_states_cuda.tracks.size();

    if (performance_opts.run) {
      // TODO: Do evt_data.fill_cca_result(...) with cuda clusters and
      // measurements
    }
  }

  if (performance_opts.run) {
    sd_performance_writer.finalize();
  }

  TRACCC_INFO("==> Statistics ... ");
  TRACCC_INFO("- read    " << n_cells << " cells");
  TRACCC_INFO("- created (cuda)  " << n_measurements_cuda
                                   << " measurements     ");
  TRACCC_INFO("- created (cuda) " << n_spacepoints_cuda << " spacepoints     ");

  TRACCC_INFO("- created (cuda) " << n_seeds_cuda << " seeds");
  TRACCC_INFO("- found (cuda)   " << n_found_tracks_cuda << " tracks");
  TRACCC_INFO("- resolved (cuda)   " << n_ambiguity_free_tracks_cuda
                                     << " tracks");
  TRACCC_INFO("- fitted (cuda)  " << n_fitted_tracks_cuda << " tracks");
  TRACCC_INFO("==>Elapsed times... " << elapsedTimes);

  return 0;
}

// The main routine
//
int main(int argc, char* argv[]) {
  std::unique_ptr<const traccc::Logger> logger =
      traccc::getDefaultLogger("CudaSeqExample", traccc::Logging::Level::INFO);

  // Program options.
  traccc::opts::detector detector_opts;
  traccc::opts::magnetic_field bfield_opts;
  traccc::opts::input_data input_opts;
  traccc::opts::clusterization clusterization_opts;
  traccc::opts::track_seeding seeding_opts;
  traccc::opts::track_gbts_seeding seeding_gbts_opts;
  traccc::opts::track_finding finding_opts;
  traccc::opts::track_propagation propagation_opts;
  traccc::opts::track_resolution resolution_opts;
  traccc::opts::track_fitting fitting_opts;
  traccc::opts::performance performance_opts;
  traccc::opts::accelerator accelerator_opts;
  traccc::opts::program_options program_opts{
      "Full Tracking Chain Using CUDA",
      {detector_opts, bfield_opts, input_opts, clusterization_opts,
       seeding_opts, seeding_gbts_opts, finding_opts, propagation_opts,
       resolution_opts, performance_opts, fitting_opts, accelerator_opts},
      argc,
      argv,
      logger->cloneWithSuffix("Options")};

  // Run the application.
  return seq_run(detector_opts, bfield_opts, input_opts, clusterization_opts,
                 seeding_opts, seeding_gbts_opts, finding_opts,
                 propagation_opts, resolution_opts, fitting_opts,
                 performance_opts, accelerator_opts, logger->clone(),
                 seeding_gbts_opts.useGBTS);
}
