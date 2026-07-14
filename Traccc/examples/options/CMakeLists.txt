# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2022-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the "build" of the traccc::option library.
traccc_add_library( traccc_options options TYPE SHARED
  # header files
  "include/traccc/options/details/interface.hpp"
  "include/traccc/options/details/value_array.hpp"
  "include/traccc/options/details/value_array.ipp"
  "include/traccc/options/accelerator.hpp"
  "include/traccc/options/clusterization.hpp"
  "include/traccc/options/detector.hpp"
  "include/traccc/options/generation.hpp"
  "include/traccc/options/input_data.hpp"
  "include/traccc/options/magnetic_field.hpp"
  "include/traccc/options/output_data.hpp"
  "include/traccc/options/performance.hpp"
  "include/traccc/options/program_options.hpp"
  "include/traccc/options/telescope_detector.hpp"
  "include/traccc/options/threading.hpp"
  "include/traccc/options/throughput.hpp"
  "include/traccc/options/track_finding.hpp"
  "include/traccc/options/track_fitting.hpp"
  "include/traccc/options/track_propagation.hpp"
  "include/traccc/options/track_resolution.hpp"
  "include/traccc/options/track_seeding.hpp"
  "include/traccc/options/track_gbts_seeding.hpp"
  # source files
  "src/details/interface.cpp"
  "src/accelerator.cpp"
  "src/clusterization.cpp"
  "src/detector.cpp"
  "src/generation.cpp"
  "src/input_data.cpp"
  "src/logging.cpp"
  "src/magnetic_field.cpp"
  "src/output_data.cpp"
  "src/performance.cpp"
  "src/program_options.cpp"
  "src/seed_matching.cpp"
  "src/telescope_detector.cpp"
  "src/threading.cpp"
  "src/throughput.cpp"
  "src/track_finding.cpp"
  "src/track_fitting.cpp"
  "src/track_matching.cpp"
  "src/track_propagation.cpp"
  "src/track_resolution.cpp"
  "src/track_seeding.cpp"
  "src/track_gbts_seeding.cpp"
  "src/truth_finding.cpp"
  )
target_link_libraries( traccc_options
   PUBLIC
      Boost::program_options
      traccc::io
      traccc::performance
      traccc::utils
)
