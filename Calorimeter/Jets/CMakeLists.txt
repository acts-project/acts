acts_add_library(
    CalorimeterJets
    src/CalorimeterJetReconstruction.cpp
    ACTS_INCLUDE_FOLDER include/ActsCalorimeter
)
target_include_directories(
    ActsCalorimeterJets
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(
    ActsCalorimeterJets
    PUBLIC Acts::Calorimeter
    PRIVATE FastJet::FastJet
)
acts_compile_headers(CalorimeterJets GLOB "include/**/*.hpp")

add_subdirectory(Examples)
