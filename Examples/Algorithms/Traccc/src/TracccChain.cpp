// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/TracccChain.hpp"

// ================================================================
// Backend selection: ACTS_ENABLE_CUDA is set by the ACTS build
// system when -DACTS_ENABLE_CUDA=ON. No definition means CPU.
// ================================================================

#if defined(ACTS_ENABLE_CUDA)

#include "traccc/cuda/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"
#include "traccc/cuda/clusterization/clusterization_algorithm.hpp"
#include "traccc/cuda/clusterization/measurement_sorting_algorithm.hpp"
#include "traccc/cuda/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/cuda/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/cuda/seeding/seed_parameter_estimation_algorithm.hpp"
#include "traccc/cuda/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"
#include "traccc/cuda/seeding/triplet_seeding_algorithm.hpp"
#include "traccc/cuda/utils/make_magnetic_field.hpp"
#include "traccc/cuda/utils/stream.hpp"
#include "vecmem/memory/cuda/device_memory_resource.hpp"
#include "vecmem/memory/cuda/host_memory_resource.hpp"
#include "vecmem/memory/host_memory_resource.hpp"
#include "vecmem/utils/cuda/async_copy.hpp"

#else  // CPU

#include "traccc/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"
#include "traccc/clusterization/clusterization_algorithm.hpp"
#include "traccc/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"
#include "vecmem/memory/host_memory_resource.hpp"
#include "vecmem/utils/copy.hpp"  // FIX: was vecmem/utils/cuda/async_copy.hpp — wrong for CPU build

#endif

// traccc EDM
#include "traccc/edm/silicon_cell_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/edm/track_collection.hpp"

// traccc geometry + IO
#include "traccc/geometry/detector_buffer.hpp"
#include "traccc/geometry/detector_conditions_description.hpp"
#include "traccc/geometry/detector_design_description.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/io/data_format.hpp"
#include "traccc/io/read_detector.hpp"
#include "traccc/io/read_detector_description.hpp"
#include "traccc/io/read_magnetic_field.hpp"

// traccc magnetic field
#include "Acts/Utilities/Logger.hpp"

#include "traccc/bfield/magnetic_field.hpp"

namespace ActsExamples::Traccc {

struct TracccChain::Impl {
  std::unique_ptr<const Acts::Logger> m_logger;
  const Acts::Logger& logger() const { return *m_logger; }
  TracccChainConfig m_cfg;

  // ---- memory resources ----
  struct BaseMemory {
    vecmem::copy* copy = nullptr;
    vecmem::memory_resource* mr = nullptr;  // host-accessible resource for IO
  };
  struct HostMemory : public BaseMemory {
    vecmem::host_memory_resource host_mr;
    vecmem::copy host_copy;
    HostMemory() : host_mr(), host_copy() {
      copy = &host_copy;
      mr = &host_mr;
    }
  };
#if defined(ACTS_ENABLE_CUDA)
  struct CudaMemory : public BaseMemory {
    vecmem::cuda::host_memory_resource cuda_host_mr;
    vecmem::cuda::device_memory_resource device_mr;
    traccc::memory_resource mr_device{device_mr, &cuda_host_mr};
    traccc::cuda::stream stream;
    vecmem::cuda::async_copy cuda_copy;
    CudaMemory()
        : cuda_host_mr(),
          device_mr(),
          stream(),
          cuda_copy(stream.cudaStream()) {
      copy = &cuda_copy;
      mr = &cuda_host_mr;  // host-pinned; usable by IO and host-side allocs
    }
  };
#endif

  std::unique_ptr<BaseMemory> m_memory;

  // ---- detector description ----
  // unique_ptr: constructed in body once memory resource is known
  std::unique_ptr<traccc::detector_design_description::host> m_det_descr;
  std::unique_ptr<traccc::detector_conditions_description::host> m_det_cond;
  traccc::host_detector m_host_detector;

  // Device-side copies — only populated for the CUDA backend
  // you can imagine similar setup for HIP, Alpaka etc.
#if defined(ACTS_ENABLE_CUDA)
  std::unique_ptr<traccc::detector_buffer> m_device_detector;
  std::unique_ptr<traccc::detector_design_description::buffer>
      m_device_det_descr;
  std::unique_ptr<traccc::detector_conditions_description::buffer>
      m_device_det_cond;
#endif

  // ---- magnetic field ----
  traccc::magnetic_field m_host_field;
  traccc::vector3
      m_field_vec;  // uniform approximation for CPU params estimation
#if defined(ACTS_ENABLE_CUDA)
  traccc::magnetic_field m_device_field;
#endif

  // ---- algorithm objects (host / CPU) ----
  std::unique_ptr<traccc::host::clusterization_algorithm> m_host_clusterization;
  std::unique_ptr<traccc::host::silicon_pixel_spacepoint_formation_algorithm>
      m_host_spacepoint_formation;
  std::unique_ptr<traccc::host::seeding_algorithm> m_host_seeding;
  std::unique_ptr<traccc::host::track_params_estimation> m_host_params_est;
  std::unique_ptr<traccc::host::combinatorial_kalman_filter_algorithm>
      m_host_finding;
  std::unique_ptr<traccc::host::greedy_ambiguity_resolution_algorithm>
      m_host_ambiguity_resolution;
  std::unique_ptr<traccc::host::kalman_fitting_algorithm> m_host_fitting;

  // ---- algorithm objects (device / CUDA) ----
#if defined(ACTS_ENABLE_CUDA)
  std::unique_ptr<traccc::device::clusterization_algorithm>
      m_device_clusterization;
  std::unique_ptr<traccc::device::silicon_pixel_spacepoint_formation_algorithm>
      m_device_spacepoint_formation;
  std::unique_ptr<traccc::device::triplet_seeding_algorithm> m_device_seeding;
  std::unique_ptr<
      traccc::algorithm<traccc::bound_track_parameters_collection_types::buffer(
          const traccc::magnetic_field&,
          const traccc::edm::measurement_collection::const_view&,
          const traccc::edm::spacepoint_collection::const_view&,
          const traccc::edm::seed_collection::const_view&)>>
      m_device_params_est;
  std::unique_ptr<traccc::algorithm<
      traccc::edm::track_container<traccc::default_algebra>::buffer(
          const traccc::detector_buffer&, const traccc::magnetic_field&,
          const traccc::edm::measurement_collection::const_view&,
          const traccc::bound_track_parameters_collection_types::const_view&)>>
      m_device_finding;
  std::unique_ptr<traccc::algorithm<
      traccc::edm::track_container<traccc::default_algebra>::buffer(
          const traccc::detector_buffer&, const traccc::magnetic_field&,
          const traccc::edm::track_container<
              traccc::default_algebra>::const_view&)>>
      m_device_fitting;
  // ambiguity resolution is not yet implemented on device
#endif

  // --------------------------------------------------------------
  explicit Impl(const TracccChainConfig& cfg, Acts::Logging::Level logLevel)
      : m_logger(Acts::getDefaultLogger("TracccChain", logLevel)), m_cfg(cfg) {
    if (cfg.backend == TracccChainConfig::Backend::CUDA) {
#if !defined(ACTS_ENABLE_CUDA)
      throw std::runtime_error(
          "TracccChainConfig requests CUDA backend but ACTS_ENABLE_CUDA is not "
          "set");
#else
      ACTS_INFO("TracccChain configured to use CUDA backend");
      m_memory = std::make_unique<CudaMemory>();
      auto* cuda_mem = static_cast<CudaMemory*>(m_memory.get());

      m_det_descr = std::make_unique<traccc::detector_design_description::host>(
          cuda_mem->cuda_host_mr);
      m_det_cond =
          std::make_unique<traccc::detector_conditions_description::host>(
              cuda_mem->cuda_host_mr);

      m_device_clusterization =
          std::make_unique<traccc::cuda::clusterization_algorithm>(
              cuda_mem->mr_device, cuda_mem->cuda_copy, cuda_mem->stream,
              cfg.clusteringConfig,
              Acts::getDefaultLogger("TracccClusterization", logLevel));
      m_device_spacepoint_formation = std::make_unique<
          traccc::cuda::silicon_pixel_spacepoint_formation_algorithm>(
          cuda_mem->mr_device, cuda_mem->cuda_copy, cuda_mem->stream,
          Acts::getDefaultLogger("TracccSpFormation", logLevel));
      m_device_seeding =
          std::make_unique<traccc::cuda::triplet_seeding_algorithm>(
              cfg.seedfinderConfig,
              traccc::spacepoint_grid_config{cfg.seedfinderConfig},
              cfg.seedfilterConfig, cuda_mem->mr_device, cuda_mem->cuda_copy,
              cuda_mem->stream,
              Acts::getDefaultLogger("TracccSeeding", logLevel));
      m_device_params_est =
          std::make_unique<traccc::cuda::seed_parameter_estimation_algorithm>(
              cfg.trackParamsEstConfig, cuda_mem->mr_device,
              cuda_mem->cuda_copy, cuda_mem->stream,
              Acts::getDefaultLogger("TracccParamsEst", logLevel));
      m_device_finding =
          std::make_unique<traccc::cuda::combinatorial_kalman_filter_algorithm>(
              cfg.findingConfig, cuda_mem->mr_device, cuda_mem->cuda_copy,
              cuda_mem->stream,
              Acts::getDefaultLogger("TracccFinding", logLevel));
      m_device_fitting =
          std::make_unique<traccc::cuda::kalman_fitting_algorithm>(
              cfg.fittingConfig, cuda_mem->mr_device, cuda_mem->cuda_copy,
              cuda_mem->stream,
              Acts::getDefaultLogger("TracccFitting", logLevel));
#endif
    } else {
      ACTS_INFO("TracccChain configured to use CPU backend");
      m_memory = std::make_unique<HostMemory>();
      auto* host_mem = static_cast<HostMemory*>(m_memory.get());

      m_det_descr = std::make_unique<traccc::detector_design_description::host>(
          host_mem->host_mr);
      m_det_cond =
          std::make_unique<traccc::detector_conditions_description::host>(
              host_mem->host_mr);

      m_host_clusterization =
          std::make_unique<traccc::host::clusterization_algorithm>(
              host_mem->host_mr,
              Acts::getDefaultLogger("TracccClusterization", logLevel));
      m_host_spacepoint_formation = std::make_unique<
          traccc::host::silicon_pixel_spacepoint_formation_algorithm>(
          host_mem->host_mr,
          Acts::getDefaultLogger("TracccSpFormation", logLevel));
      m_host_seeding = std::make_unique<traccc::host::seeding_algorithm>(
          cfg.seedfinderConfig,
          traccc::spacepoint_grid_config{cfg.seedfinderConfig},
          cfg.seedfilterConfig, host_mem->host_mr,
          Acts::getDefaultLogger("TracccSeeding", logLevel));
      m_host_params_est =
          std::make_unique<traccc::host::track_params_estimation>(
              cfg.trackParamsEstConfig, host_mem->host_mr,
              Acts::getDefaultLogger("TracccParamsEst", logLevel));
      m_host_finding =
          std::make_unique<traccc::host::combinatorial_kalman_filter_algorithm>(
              cfg.findingConfig, host_mem->host_mr,
              Acts::getDefaultLogger("TracccFinding", logLevel));
      m_host_ambiguity_resolution =
          std::make_unique<traccc::host::greedy_ambiguity_resolution_algorithm>(
              cfg.resolutionConfig, host_mem->host_mr,
              Acts::getDefaultLogger("TracccAmbiguity", logLevel));
      m_host_fitting = std::make_unique<traccc::host::kalman_fitting_algorithm>(
          cfg.fittingConfig, host_mem->host_mr, *m_memory->copy,
          Acts::getDefaultLogger("TracccFitting", logLevel));
    }

    // Load detector description
    if (!cfg.digitizationFile.empty()) {
      traccc::io::read_detector_description(
          *m_det_descr, *m_det_cond, cfg.detectorFile, cfg.digitizationFile,
          cfg.conditionsFile, traccc::data_format::json);
    }

    // Load detector geometry — requires a host memory resource
    traccc::io::read_detector(m_host_detector, *m_memory->mr, cfg.detectorFile,
                              cfg.materialFile, cfg.gridFile);

    // Load magnetic field
    traccc::io::read_magnetic_field(m_host_field, cfg.magneticFieldFile);
    m_field_vec = {0.f, 0.f, cfg.bFieldInZ};

    if (cfg.backend == TracccChainConfig::Backend::CUDA) {
#if defined(ACTS_ENABLE_CUDA)
      auto* cuda_mem = static_cast<CudaMemory*>(m_memory.get());

      m_device_field = traccc::cuda::make_magnetic_field(
          m_host_field, traccc::cuda::magnetic_field_storage::global_memory);

      m_device_detector = std::make_unique<traccc::detector_buffer>(
          traccc::buffer_from_host_detector(
              m_host_detector, cuda_mem->device_mr, cuda_mem->cuda_copy));

      std::vector<unsigned int> descr_sizes(m_det_descr->size());
      for (std::size_t i = 0; i < m_det_descr->size(); ++i) {
        auto design = m_det_descr->at(i);
        descr_sizes[i] =
            std::max(static_cast<unsigned int>(design.bin_edges_x().size()),
                     static_cast<unsigned int>(design.bin_edges_y().size()));
      }
      m_device_det_descr =
          std::make_unique<traccc::detector_design_description::buffer>(
              descr_sizes, cuda_mem->device_mr, &cuda_mem->cuda_host_mr,
              vecmem::data::buffer_type::resizable);

      cuda_mem->copy->setup(*m_device_det_descr)->wait();
      cuda_mem->copy
          ->operator()(vecmem::get_data(*m_det_descr), *m_device_det_descr)
          ->wait();

      m_device_det_cond =
          std::make_unique<traccc::detector_conditions_description::buffer>(
              static_cast<
                  traccc::detector_conditions_description::buffer::size_type>(
                  m_det_cond->size()),
              cuda_mem->device_mr);
      cuda_mem->copy->setup(*m_device_det_cond)->wait();
      cuda_mem->copy
          ->operator()(vecmem::get_data(*m_det_cond), *m_device_det_cond)
          ->wait();

      cuda_mem->stream.synchronize();
#endif
    }

    ACTS_INFO(
        "TracccChain initialised ("
        << (cfg.backend == TracccChainConfig::Backend::CPU ? "CPU" : "CUDA")
        << "\n"
        << "), detector: " << cfg.detectorFile << "\n"
        << ", digitization: " << cfg.digitizationFile << "\n"
        << ", conditions: " << cfg.conditionsFile << "\n"
        << ", material: " << cfg.materialFile << "\n"
        << ", grid: " << cfg.gridFile << "\n"
        << ", bfield: " << cfg.magneticFieldFile);
  }

  // --------------------------------------------------------------
  EventResult run(traccc::edm::measurement_collection::host& meas_host,
                  traccc::edm::spacepoint_collection::host& sp_host) const {
    EventResult result;
    result.n_measurements = meas_host.size();
    result.n_spacepoints = sp_host.size();

    ACTS_DEBUG("Running traccc chain with "
               << meas_host.size() << " measurements and " << sp_host.size()
               << " spacepoints.");

    if (m_cfg.backend != TracccChainConfig::Backend::CPU) {
#if defined(ACTS_ENABLE_CUDA)
      auto* cuda_mem = static_cast<CudaMemory*>(m_memory.get());

      traccc::edm::measurement_collection::buffer meas_buf(
          static_cast<unsigned int>(meas_host.size()), cuda_mem->device_mr);
      traccc::edm::spacepoint_collection::buffer sp_buf(
          static_cast<unsigned int>(sp_host.size()), cuda_mem->device_mr);

      cuda_mem->copy->setup(meas_buf)->wait();
      cuda_mem->copy->setup(sp_buf)->wait();
      cuda_mem->copy->operator()(vecmem::get_data(meas_host), meas_buf)->wait();
      cuda_mem->copy->operator()(vecmem::get_data(sp_host), sp_buf)->wait();

      auto seeds_buf = (*m_device_seeding)(sp_buf);
      result.n_seeds = cuda_mem->copy->get_size(seeds_buf);

      auto params_buf =
          (*m_device_params_est)(m_device_field, meas_buf, sp_buf, seeds_buf);

      auto track_candidates_buf = (*m_device_finding)(
          *m_device_detector, m_device_field, meas_buf, params_buf);
      result.n_track_candidates =
          cuda_mem->copy->get_size(track_candidates_buf.tracks);

      // No device ambiguity resolution yet — pass candidates directly to fitter
      // if needed Fitting happened through MBF track finding

      //   auto fitted_tracks_buf = (*m_device_fitting)(
      //       *m_device_detector, m_device_field,
      //       traccc::edm::track_container<traccc::default_algebra>::const_view(
      //           track_candidates_buf));

      //   traccc::edm::track_container<traccc::default_algebra>::host
      //   found_tracks_host;
      //   cuda_mem->copy->operator()(track_candidates_buf.tracks,
      //   found_tracks_host.tracks,
      //                              vecmem::copy::type::device_to_host)->wait();
      //   cuda_mem->copy->operator()(track_candidates_buf.states,
      //   found_tracks_host.states,
      //                              vecmem::copy::type::device_to_host)->wait();
      //   cuda_mem->stream.synchronize();

      result.n_fitted_tracks =
          cuda_mem->copy->get_size(track_candidates_buf.tracks);
#endif
    } else {
      // CPU execution
      auto seeds_host = (*m_host_seeding)(vecmem::get_data(sp_host));
      result.n_seeds = seeds_host.size();

      auto params_host = (*m_host_params_est)(
          vecmem::get_data(meas_host), vecmem::get_data(sp_host),
          vecmem::get_data(seeds_host), m_field_vec);

      auto track_candidates_host = (*m_host_finding)(
          m_host_detector, m_host_field, vecmem::get_data(meas_host),
          vecmem::get_data(params_host));
      result.n_track_candidates = track_candidates_host.tracks.size();

      auto fitted_tracks_host = (*m_host_fitting)(
          m_host_detector, m_host_field,
          traccc::edm::track_container<traccc::default_algebra>::const_data(
              track_candidates_host));
      result.n_fitted_tracks = fitted_tracks_host.tracks.size();
    }

    return result;
  }
};

TracccChain::TracccChain(const TracccChainConfig& cfg,
                         Acts::Logging::Level logLevel)
    : m_impl(std::make_unique<Impl>(cfg, logLevel)) {}

TracccChain::~TracccChain() = default;
TracccChain::TracccChain(TracccChain&&) noexcept = default;
TracccChain& TracccChain::operator=(TracccChain&&) noexcept = default;

EventResult TracccChain::operator()(
    traccc::edm::measurement_collection::host& measurements,
    traccc::edm::spacepoint_collection::host& spacepoints) const {
  return m_impl->run(measurements, spacepoints);
}

}  // namespace ActsExamples::Traccc
