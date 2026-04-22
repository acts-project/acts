# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building detray host tutorials")

# Set up the rest of the tutorials
set(source_files
    "detector/build_predefined_detectors.cpp"
    "detector/detector_ray_scan.cpp"
    "detector/detector_to_dot.cpp"
    "detector/do_it_yourself_detector.cpp"
    "detector/read_detector_from_file.cpp"
    "detector/track_geometry_changes.cpp"
    "detector/write_detector_to_file.cpp"
    "propagation/navigation_inspection.cpp"
    "propagation/propagation.cpp"
    "propagation/define_an_actor.cpp"
)

foreach(file ${source_files})
    # Extract name of the tutorial
    get_filename_component(file_name ${file} NAME_WE)

    detray_add_tutorial(${file_name}
      ${file}
      LINK_LIBRARIES detray::core_array covfie::core vecmem::core detray::io
                     detray::tutorial detray::test_common detray::test_utils detray::validation_utils
    )
endforeach()
