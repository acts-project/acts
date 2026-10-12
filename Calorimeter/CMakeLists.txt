acts_add_library(
    Calorimeter
    src/Clustering/CalorimeterClusterer.cpp
    src/Digitization/CalorimeterResponse.cpp
    ACTS_INCLUDE_FOLDER include/ActsCalorimeter
)
target_include_directories(
    ActsCalorimeter
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(ActsCalorimeter PUBLIC Acts::Core)
acts_compile_headers(Calorimeter GLOB "include/**/*.hpp")

add_subdirectory(Examples)
