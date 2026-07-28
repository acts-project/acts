add_subdirectory(Algorithms)
add_subdirectory(EventData)
add_subdirectory(Io)
add_subdirectory(Framework)
add_subdirectory_if(
    CudaMuonSegmentFitting
    ACTS_BUILD_EXAMPLES_CUDA_MUON_SEGMENT_FITTING
)
