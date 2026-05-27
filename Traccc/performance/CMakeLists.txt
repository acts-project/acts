# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2022-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Project include(s).
include( traccc-compiler-options-cpp )

# Set up the build of the traccc::performance library.
traccc_add_library( traccc_performance performance TYPE SHARED
   # Efficiency calculation code.
   "src/efficiency/duplication_plot_tool.hpp"
   "src/efficiency/eff_plot_tool.hpp"
   "include/traccc/efficiency/seeding_performance_writer.hpp"
   "src/efficiency/seeding_performance_writer.cpp"
   "include/traccc/efficiency/finding_performance_writer.hpp"
   "src/efficiency/finding_performance_writer.cpp"
   "src/efficiency/track_classification.hpp"
   "include/traccc/efficiency/nseed_performance_writer.hpp"
   "src/efficiency/nseed_performance_writer.cpp"
   "src/efficiency/track_filter.cpp"
   "src/efficiency/track_matcher.cpp"
   # Resolution calculation code.
   "include/traccc/resolution/fitting_performance_writer.hpp"
   "include/traccc/resolution/res_plot_tool_config.hpp"
   "include/traccc/resolution/stat_plot_tool_config.hpp"
   "src/resolution/fitting_performance_writer.cpp"
   "src/resolution/res_plot_tool.hpp"
   "src/resolution/res_plot_tool.cpp"
   "src/resolution/stat_plot_tool.hpp"
   "src/resolution/stat_plot_tool.cpp"
   # Utils
   "include/traccc/utils/event_data.hpp"
   "include/traccc/utils/fill_track_container.hpp"
   "src/utils/event_data.cpp"
   "include/traccc/utils/helpers.hpp"
   "include/traccc/utils/ranges.hpp"
   "include/traccc/utils/transcribe_to_trace.hpp"
   "src/utils/helpers.hpp"
   "src/utils/helpers.cpp"
   # Value/object comparison code.
   "include/traccc/performance/details/is_same_angle.hpp"
   "src/performance/details/is_same_angle.cpp"
   "include/traccc/performance/details/is_same_scalar.hpp"
   "src/performance/details/is_same_scalar.cpp"
   "include/traccc/performance/details/is_same_object.hpp"
   "include/traccc/performance/impl/is_same_object.ipp"
   "include/traccc/performance/impl/is_same_seed.ipp"
   "include/traccc/performance/impl/is_same_measurement.ipp"
   "include/traccc/performance/impl/is_same_spacepoint.ipp"
   "include/traccc/performance/impl/is_same_track_parameters.ipp"
   "include/traccc/performance/impl/is_same_track_state.ipp"
   "include/traccc/performance/impl/is_same_track.ipp"
   "src/performance/details/is_same_object.cpp"
   "include/traccc/performance/details/comparator_factory.hpp"
   "include/traccc/performance/impl/comparator_factory.ipp"
   "include/traccc/performance/impl/seed_comparator_factory.ipp"
   "include/traccc/performance/impl/track_comparator_factory.ipp"
   # Collection/container comparison code.
   "include/traccc/performance/collection_comparator.hpp"
   "include/traccc/performance/impl/collection_comparator.ipp"
   "include/traccc/performance/container_comparator.hpp"
   "include/traccc/performance/impl/container_comparator.ipp"
   "include/traccc/performance/soa_comparator.hpp"
   "include/traccc/performance/impl/soa_comparator.ipp"
   # Truth propagation comparison
   "include/traccc/performance/kalman_filter_comparison.hpp"
   "src/performance/kalman_filter_comparison.cpp"
   # Performance time measurement code.
   "include/traccc/performance/timer.hpp"
   "src/performance/timer.cpp"
   "include/traccc/performance/timing_info.hpp"
   "src/performance/timing_info.cpp"
   "include/traccc/performance/throughput.hpp"
   "src/performance/throughput.cpp" )
target_link_libraries( traccc_performance
   PUBLIC traccc::core traccc::io covfie::core detray::test_utils detray::validation_utils)

# Use ROOT in traccc::performance, if requested.
if( TRACCC_USE_ROOT )
   find_package( ROOT COMPONENTS Core RIO Hist REQUIRED )
   target_link_libraries( traccc_performance
      PRIVATE ROOT::Core ROOT::RIO ROOT::Hist )
   target_compile_definitions( traccc_performance
      PRIVATE TRACCC_HAVE_ROOT )
endif()

# Use CUDA/NVTX in traccc::performance, if it's available.
find_package( CUDAToolkit )
if( CUDAToolkit_FOUND )
   target_link_libraries( traccc_performance
      PRIVATE CUDA::cudart ${CMAKE_DL_LIBS} )
   target_compile_definitions( traccc_performance
      PRIVATE TRACCC_HAVE_NVTX )
endif()
