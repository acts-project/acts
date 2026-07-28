/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/bfield/construct_const_bfield.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/io/detector.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/performance/kalman_filter_comparison.hpp"
#include "traccc/simulation/event_generators.hpp"
#include "traccc/simulation/simulator.hpp"
#include "traccc/utils/fill_track_container.hpp"

// Test include(s).
#include "tests/test_detectors.hpp"
#include "tests/toy_detector_fixture.hpp"

// GTest include(s).
#include <gtest/gtest.h>

// System include(s)
#include <algorithm>
#include <filesystem>
#include <optional>
#include <string>

using namespace traccc;

constexpr std::size_t n_events{1u};
constexpr std::size_t n_tracks{1000u};

constexpr detray::pdg_particle ptc_type{detray::muon<scalar>()};

// Truth particle cuts
constexpr scalar min_p{50.f * traccc::unit<traccc::scalar>::MeV};
constexpr scalar max_r{75.f * traccc::unit<traccc::scalar>::mm};

/// Test suite for navigation tests for the CKF
class KF_integration_test_toy_detector
    : public ToyDetectorFixture,
      public ::testing::WithParamInterface<
          std::tuple<float, float, bool, bool, bool, float, float>> {};

/// Test the detray navigation on simulated tracks
TEST_P(KF_integration_test_toy_detector, toy_detector) {
  // TODO: Enable these tests at a later date.
  GTEST_SKIP();

  using detector_t = traccc::default_detector::host;
  using algebra_t = typename detector_t::algebra_type;
  using b_field_t = covfie::field<traccc::const_bfield_backend_t<scalar>>;
  using track_t = traccc::free_track_parameters<algebra_t>;

  using sf_candidate_t =
      traccc::propagation_validator::candidate_type<detector_t>;
  using generator_t = detray::random_track_generator<track_t>;
  using writer_t = smearing_writer<measurement_smearer<algebra_t>>;
  using simulator_t = simulator<detector_t, b_field_t, generator_t, writer_t>;

  vecmem::host_memory_resource host_mr;

  std::unique_ptr<const traccc::Logger> logger =
      traccc::getDefaultLogger("KF_integration_test_toy_detector_toy_detector",
                               traccc::Logging::Level::INFO);

  const ::testing::TestInfo& test_info =
      *::testing::UnitTest::GetInstance()->current_test_info();
  std::string name =
      std::string{test_info.test_suite_name()} + "_" + test_info.name();
  std::replace(name.begin(), name.end(), '/', '_');
  const std::filesystem::path det_dir{name};

  WriteDetector(true, name);

  detray::io::detector_reader_config reader_cfg{};
  reader_cfg.add_file((det_dir / "toy_detector_geometry.json").native())
      .add_file((det_dir / "toy_detector_surface_grids.json").native())
      .do_check(true);
  if (std::get<2>(GetParam())) {
    reader_cfg.add_file((det_dir / "toy_detector_material_maps.json").native());
  }

  auto [io_det, names] =
      detray::io::read_detector<traccc::default_detector::host>(host_mr,
                                                                reader_cfg);
  traccc::host_detector host_det{};
  host_det.template set<detector_traits<typename detector_t::metadata>>(
      std::move(io_det));
  const auto& det = host_det.template as<
      traccc::detector_traits<typename detector_t::metadata>>();

  // Create B field
  traccc::magnetic_field field = traccc::construct_const_bfield(B);

  // Create track generator
  const scalar pT{std::get<0>(GetParam())};
  generator_t::configuration gen_cfg{};
  gen_cfg.n_tracks(n_tracks).eta_range(-3, 3).p_T(pT).randomize_charge(true);
  // Choose different random seed than detray for more test coverage
  gen_cfg.seed(135346);

  // Create data directory
  std::filesystem::path data_dir{traccc::io::data_directory()};
  std::filesystem::path outdir{std::filesystem::path{"fast_track_simulation"} /
                               name};

  std::filesystem::path full_path = data_dir / outdir;
  if (!std::filesystem::exists(full_path)) {
    if (std::error_code err;
        !std::filesystem::create_directories(full_path, err)) {
      throw std::runtime_error(err.message());
    }
  }

  // Create measurement smearer
  measurement_smearer<algebra_t> smearer(smearing[0], smearing[1]);
  auto sim = simulator_t(
      ptc_type, n_events, det,
      field.template as_field<traccc::const_bfield_backend_t<scalar>>(),
      generator_t{gen_cfg}, writer_t::config{smearer}, full_path.string());

  // Propagation config for the simulation
  detray::propagation::config prop_cfg{};
  prop_cfg.navigation.search_window = search_window;  //< toy detector grids

  sim.get_config().propagation = prop_cfg;
  sim.get_config().propagation.stepping.step_constraint = step_constraint;
  sim.get_config().do_multiple_scattering = std::get<3>(GetParam());
  sim.get_config().do_energy_loss = std::get<4>(GetParam());
  sim.get_config().min_pT(10.f * traccc::unit<scalar>::MeV);

  // Run the simulation: Produces data files
  sim.run();

  // Configure the test
  detray::propagation_validation_config test_cfg{};
  traccc::seed_generator<detector_t>::config smearing_cfg{};

  // General config
  test_cfg.display_svg = false;

  // Specific config for the navigation test
  test_cfg.propagation = prop_cfg;
  test_cfg.propagation.navigation.intersection.mask_tolerance_scalor = 1.5f;
  test_cfg.propagation.navigation.intersection.min_mask_tolerance =
      std::get<1>(GetParam());
  test_cfg.propagation.navigation.intersection.max_mask_tolerance =
      5.f * traccc::unit<float>::mm;
  test_cfg.propagation.navigation.estimate_scattering_noise = false;
  test_cfg.propagation.stepping.path_limit = 2.f * traccc::unit<float>::m;

  // Configure the material interaction
  test_cfg.particle = ptc_type;
  test_cfg.do_multiple_scattering = std::get<3>(GetParam());
  test_cfg.do_energy_loss = std::get<4>(GetParam());

  // Prepare the data for the test

  // Initial track parameters from truth particle
  std::vector<traccc::free_track_parameters<algebra_t>> tracks{};
  // The traces of truth hits forward and in reverse order
  std::vector<vecmem::vector<sf_candidate_t>> truth_traces_fw{};
  // Measurements
  typename edm::measurement_collection::host measurements{host_mr};
  // Collection for bound track parameters and track state/measurement links
  // Input/results of the KF
  traccc::edm::track_container<algebra_t>::host track_container{host_mr};

  fill_track_containers(logger->clone(), &host_det, outdir, n_events, false,
                        min_p, max_r, tracks, truth_traces_fw, measurements,
                        track_container);

  // Save the original truth traces before dummy records are inserted for
  // missing sensitive surfaces
  auto truth_traces_fw_KF = truth_traces_fw;

  // Run the test without KF
  test_cfg.max_percent_missed = 25.f;
  test_cfg.max_percent_additional = 60.f;

  std::optional<b_field_t::view_t> field_view{
      field.template as_field<traccc::const_bfield_backend_t<scalar>>()};

  bool success = detray::propagation_validation(
      det, names, field_view, test_cfg, tracks, truth_traces_fw);

  ASSERT_TRUE(success);

  // Run the test with KF
  test_cfg.max_percent_missed = std::get<5>(GetParam());
  test_cfg.max_percent_additional = std::get<6>(GetParam());

  // Set up the containers for the Kalman actor
  edm::measurement_collection::const_device device_measurements{
      vecmem::get_data(measurements)};

  success = kalman_filter_comparison(
      det, names, field, test_cfg, smearing_cfg, logger->clone(), tracks,
      truth_traces_fw_KF, device_measurements, track_container);

  ASSERT_TRUE(success);
}

// Parameters:
// 1: p_T
// 2: min mask tolerance
// 3: Build detector with material
// 4: Do multiple scattering
// 5: Do energy loss
// 6: max % missing surfaces
// 7: max % additional surfaces

// No material - navigation should work
INSTANTIATE_TEST_SUITE_P(
    pT_05GeV_no_mat, KF_integration_test_toy_detector,
    ::testing::Values(std::make_tuple(0.5f * traccc::unit<scalar>::GeV,
                                      1e-5f * traccc::unit<float>::mm, false,
                                      false, false, 0.8f, 1.f)));

// No scattering - navigation should work (material interactor models e-loss)
INSTANTIATE_TEST_SUITE_P(
    pT_05GeV_only_eloss, KF_integration_test_toy_detector,
    ::testing::Values(std::make_tuple(0.5f * traccc::unit<scalar>::GeV,
                                      0.1f * traccc::unit<float>::mm, true,
                                      false, true, 0.5f, 7.5f)));

// Nominal (e-loss + scattering)
INSTANTIATE_TEST_SUITE_P(
    pT_05GeV_nominal, KF_integration_test_toy_detector,
    ::testing::Values(std::make_tuple(0.5f * traccc::unit<scalar>::GeV,
                                      1.f * traccc::unit<float>::mm, true, true,
                                      true, 1.7f, 54.4f)));
INSTANTIATE_TEST_SUITE_P(
    pT_100GeV_nominal, KF_integration_test_toy_detector,
    ::testing::Values(std::make_tuple(100.f * traccc::unit<scalar>::GeV,
                                      0.1f * traccc::unit<float>::mm, true,
                                      true, true, 0.2f, 8.5f)));
INSTANTIATE_TEST_SUITE_P(
    pT_5GeV_nominal, KF_integration_test_toy_detector,
    ::testing::Values(std::make_tuple(5.f * traccc::unit<scalar>::GeV,
                                      0.1f * traccc::unit<float>::mm, true,
                                      true, true, 0.1f, 8.3f)));
INSTANTIATE_TEST_SUITE_P(
    pT_1GeV_nominal, KF_integration_test_toy_detector,
    ::testing::Values(std::make_tuple(1.f * traccc::unit<scalar>::GeV,
                                      0.15f * traccc::unit<float>::mm, true,
                                      true, true, 1.7f, 11.8f)));
