add_subdirectory(Algorithms)
add_subdirectory(Detectors)
add_subdirectory_if(Nodd ACTS_BUILD_EXAMPLES_DD4HEP)
add_subdirectory(Framework)
add_subdirectory(Io)
add_subdirectory_if(Scripts ACTS_BUILD_ANALYSIS_APPS)
add_subdirectory_if(
    CudaMuonSegmentFitting
    ACTS_BUILD_EXAMPLES_CUDA_MUON_SEGMENT_FITTING
)
