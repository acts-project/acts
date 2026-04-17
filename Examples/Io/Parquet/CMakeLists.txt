acts_add_library(
    ExamplesIoParquet
    src/ParquetReader.cpp
    src/ParquetWriter.cpp
    src/ArrowInputConverter.cpp
)

target_include_directories(
    ActsExamplesIoParquet
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesIoParquet
    PUBLIC
        Acts::ExamplesFramework
        Acts::PluginArrow
        Arrow::arrow_shared
        Parquet::parquet_shared
)

acts_compile_headers(ExamplesIoParquet GLOB include/ActsExamples/**/*.hpp)
