acts_add_library(
    FatrasJson
    SHARED
    src/DataDirectory.cpp
    src/DetectorDescriptionJsonConverter.cpp
    src/EventConfigJsonConverter.cpp
    ACTS_INCLUDE_FOLDER include/ActsFatras
)
target_include_directories(
    ActsFatrasJson
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(
    ActsFatrasJson
    PUBLIC Acts::Fatras nlohmann_json::nlohmann_json
)

# An enum written as the integer it happens to be is a file no one can read and
# a value that moves when the enum does.
target_compile_definitions(
    ActsFatrasJson
    PRIVATE JSON_DISABLE_ENUM_SERIALIZATION=1
)

# The shipped descriptions, their material and the configurations fitted to
# them. `dataPath` resolves a name against the source tree in a build and
# against the installed copy otherwise.
target_compile_definitions(
    ActsFatrasJson
    PRIVATE
        "ACTS_FATRAS_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/../data\""
        "ACTS_FATRAS_INSTALLED_DATA_DIR=\"${CMAKE_INSTALL_FULL_DATADIR}/Acts/Fatras\""
)
install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../data/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/Acts/Fatras
    FILES_MATCHING
    PATTERN "*.json"
)

acts_compile_headers(FatrasJson GLOB include/**/*.hpp)
