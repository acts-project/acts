acts_add_library(
    ExamplesCalorimeter
    src/CalorimeterClusteringAlgorithm.cpp
    src/CalorimeterDigitizationAlgorithm.cpp
    ACTS_INCLUDE_FOLDER include/ActsExamples
)
target_include_directories(
    ActsExamplesCalorimeter
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(
    ActsExamplesCalorimeter
    PUBLIC Acts::Calorimeter Acts::ExamplesFramework
)
acts_compile_headers(ExamplesCalorimeter GLOB "include/**/*.hpp")

add_executable(ActsExampleCalorimeter CalorimeterExample.cpp)
target_link_libraries(ActsExampleCalorimeter PRIVATE Acts::ExamplesCalorimeter)
install(
    TARGETS ActsExampleCalorimeter
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
