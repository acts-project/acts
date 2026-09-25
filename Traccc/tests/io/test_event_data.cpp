/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/io/data_format.hpp"
#include "traccc/io/detector.hpp"
#include "traccc/io/read_digitization_config.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/utils/event_data.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

TEST(event_data, acts_odd) {
  /// Type declarations
  using host_detector_type = traccc::default_detector_traits::host;

  vecmem::host_memory_resource resource;

  const std::string path = "odd/geant4_1muon_100GeV";
  const std::string det_file = "geometries/odd/odd-detray_geometry_detray.json";

  // Read detector file
  detray::io::detector_reader_config reader_cfg{};
  reader_cfg.add_file(traccc::io::data_directory() + det_file);

  auto [host_det, names] =
      detray::io::read_detector<host_detector_type>(resource, reader_cfg);

  traccc::host_detector polymorphic_detector;
  polymorphic_detector.set<traccc::default_detector>(std::move(host_det));

  {
    // without cell
    traccc::event_data evt_data(path, 0u, resource, true, &polymorphic_detector,
                                traccc::data_format::csv, false);
    EXPECT_EQ(evt_data.m_particle_map.size(), 4515u);
    EXPECT_EQ(evt_data.m_meas_to_ptc_map.size(), 58u);
    EXPECT_EQ(evt_data.m_meas_to_param_map.size(), 58u);
  }
  {
    // with cell
    traccc::event_data evt_data(path, 0u, resource, true, &polymorphic_detector,
                                traccc::data_format::csv, true);
    EXPECT_EQ(evt_data.m_particle_map.size(), 4515u);
    EXPECT_EQ(evt_data.m_meas_to_ptc_map.size(), 58u);
    EXPECT_EQ(evt_data.m_meas_to_param_map.size(), 58u);
  }
}

TEST(event_data, mock_data) {
  /***
   * Mock data test
   *
   * Mock data consists of three particles each of which has one hit
   *
   * first particle: one hit, three cells
   * second particle: one hit, four cells
   * third particle: one hit, three cells
   *
   * [ ] [1] [ ] [ ] [ ] [ ] [ ] [ ]
   * [1][1,2][2] [ ] [ ] [ ] [ ] [ ]
   * [ ] [2] [2] [ ] [ ] [ ] [ ] [ ]
   * [ ] [ ] [ ] [ ] [ ] [ ] [ ] [ ]
   * [ ] [ ] [ ] [ ] [ ] [ ] [ ] [ ]
   * [ ] [ ] [ ] [ ] [ ] [3] [3] [3]
   * [ ] [ ] [ ] [ ] [ ] [ ] [ ] [ ]
   * [ ] [ ] [ ] [ ] [ ] [ ] [ ] [ ]
   *
   */

  /// Type declarations
  using host_detector_type = traccc::default_detector_traits::host;

  vecmem::host_memory_resource resource;

  const std::string path = TRACCC_TEST_IO_MOCK_DATA_DIR;

  // Dummy detector file
  const std::string det_file = "geometries/odd/odd-detray_geometry_detray.json";

  // Read detector file
  detray::io::detector_reader_config reader_cfg{};
  reader_cfg.add_file(traccc::io::data_directory() + det_file);

  auto [host_det, names] =
      detray::io::read_detector<host_detector_type>(resource, reader_cfg);

  traccc::host_detector polymorphic_detector;
  polymorphic_detector.set<traccc::default_detector>(std::move(host_det));

  traccc::event_data evt_data(path, 0u, resource, true, &polymorphic_detector,
                              traccc::data_format::csv, true);

  // There are three measurements
  EXPECT_EQ(evt_data.m_meas_to_ptc_map.size(), 3u);
  EXPECT_EQ(evt_data.m_meas_to_param_map.size(), 3u);
  // There are three particles
  EXPECT_EQ(evt_data.m_particle_map.size(), 3u);
  EXPECT_EQ(evt_data.m_ptc_to_meas_map.size(), 3u);

  for (auto const& [meas, ptcs] : evt_data.m_meas_to_ptc_map) {
    EXPECT_EQ(ptcs.size(), 1);

    for (auto const& [ptc, count] : ptcs) {
      if (ptc.particle_id == 4503599644147712) {
        // number of cells from 1st particle
        EXPECT_EQ(count, 3);
      } else if (ptc.particle_id == 4503599660924928) {
        // number of cells from 2nd particle
        EXPECT_EQ(count, 4);
      } else if (ptc.particle_id == 4503599744811008) {
        // number of cells from 3rd particle
        EXPECT_EQ(count, 3);
      }
    }
  }
}
