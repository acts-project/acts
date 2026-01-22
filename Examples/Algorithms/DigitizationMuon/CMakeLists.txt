acts_add_library(
    ExamplesDigitizationMuon
    src/MuonSpacePointDigitizer.cpp
)

target_include_directories(
    ActsExamplesDigitizationMuon
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesDigitizationMuon
    PUBLIC
    ActsExamplesFramework
    ActsExamplesDigitization
    PRIVATE
    ROOT::Core
    ROOT::Geom
    ROOT::Graf
    ROOT::Hist
    ROOT::Gpad
)

acts_add_headers(ActsExamplesDigitizationMuon include/**/*.hpp)
