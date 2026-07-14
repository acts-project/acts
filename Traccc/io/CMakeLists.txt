# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Project include(s).
include( traccc-compiler-options-cpp )

# Look for OpenMP.
find_package( OpenMP COMPONENTS CXX )

# Set up the "build" of the traccc::io library.
traccc_add_library( traccc_io io TYPE SHARED
  # Public headers
  "include/traccc/io/detector.hpp"
  "include/traccc/io/digitization_config.hpp"
  "include/traccc/io/read_cells.hpp"
  "include/traccc/io/read_detector.hpp"
  "include/traccc/io/read_detector_description.hpp"
  "include/traccc/io/read_digitization_config.hpp"
  "include/traccc/io/read_conditions_config.hpp"
  "include/traccc/io/read_magnetic_field.hpp"
  "include/traccc/io/read_measurements.hpp"
  "include/traccc/io/read_particles.hpp"
  "include/traccc/io/read_spacepoints.hpp"
  "include/traccc/io/data_format.hpp"
  "include/traccc/io/write.hpp"
  "include/traccc/io/utils.hpp"
  "include/traccc/io/csv/cell.hpp"
  "include/traccc/io/csv/hit.hpp"
  "include/traccc/io/csv/measurement_hit_id.hpp"
  "include/traccc/io/csv/measurement.hpp"
  "include/traccc/io/csv/particle.hpp"
  "include/traccc/io/csv/surface.hpp"
  "include/traccc/io/csv/make_cell_reader.hpp"
  "include/traccc/io/csv/make_hit_reader.hpp"
  "include/traccc/io/csv/make_measurement_edm.hpp"
  "include/traccc/io/csv/make_measurement_hit_id_reader.hpp"
  "include/traccc/io/csv/make_particle_reader.hpp"
  "include/traccc/io/csv/make_surface_reader.hpp"
  # Implementation
  "src/data_format.cpp"
  "src/read_cells.cpp"
  "src/read_detector.cpp"
  "src/read_detector_description.cpp"
  "src/read_digitization_config.cpp"
  "src/read_conditions_config.cpp"
  "src/read_magnetic_field.cpp"
  "src/read_measurements.cpp"
  "src/read_particles.cpp"
  "src/read_spacepoints.cpp"
  "src/write.cpp"
  "src/utils.cpp"
  "src/read_binary.hpp"
  "src/write_binary.hpp"
  "src/csv/make_surface_reader.cpp"
  "src/csv/make_cell_reader.cpp"
  "src/csv/read_cells.hpp"
  "src/csv/read_cells.cpp"
  "src/csv/write_cells.hpp"
  "src/csv/write_cells.cpp"
  "src/csv/read_spacepoints.hpp"
  "src/csv/read_spacepoints.cpp"
  "src/csv/make_measurement_reader.cpp"
  "src/csv/read_measurements.hpp"
  "src/csv/read_measurements.cpp"
  "src/csv/make_measurement_edm.cpp"
  "src/csv/make_measurement_hit_id_reader.cpp"
  "src/csv/make_hit_reader.cpp"
  "src/csv/make_particle_reader.cpp"
  "src/csv/read_particles.hpp"
  "src/csv/read_particles.cpp"
  "src/json/read_conditions_config.hpp"
  "src/json/read_conditions_config.cpp"
  "src/json/read_digitization_config.hpp"
  "src/json/read_digitization_config.cpp"
  "src/json/write_digitization_config.hpp"
  "src/json/write_digitization_config.cpp"
  "src/obj/write_seeds.hpp"
  "src/obj/write_seeds.cpp"
  "src/obj/write_spacepoints.hpp"
  "src/obj/write_spacepoints.cpp"
  "src/obj/write_tracks.hpp"
  "src/obj/write_tracks.cpp"
  )
target_link_libraries( traccc_io
  PUBLIC vecmem::core traccc::core Acts::Core
  PRIVATE detray::core detray::io Acts::PluginJson )
target_compile_definitions( traccc_io
  PRIVATE TRACCC_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/data" )
if( OpenMP_CXX_FOUND )
  target_link_libraries( traccc_io PRIVATE OpenMP::OpenMP_CXX )
endif()
