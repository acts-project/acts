# Create relocatable CMake package config files

include(CMakePackageConfigHelpers)

# use path suggested by
# https://cmake.org/cmake/help/v3.18/manual/cmake-packages.7.html
set(install_package_config_dir "${CMAKE_INSTALL_LIBDIR}/cmake/Acts")

# version is taken automatically from PROJECT_VERSION; no need to specify
write_basic_package_version_file(
    ${PROJECT_BINARY_DIR}/ActsConfigVersion.cmake
    COMPATIBILITY SameMajorVersion
)
configure_package_config_file(
    ${CMAKE_CURRENT_LIST_DIR}/ActsConfig.cmake.in
    ${PROJECT_BINARY_DIR}/ActsConfig.cmake
    INSTALL_DESTINATION ${install_package_config_dir}
    PATH_VARS CMAKE_INSTALL_BINDIR CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_LIBDIR
)

# install cmake package configs
install(
    FILES
        ${PROJECT_BINARY_DIR}/ActsConfigVersion.cmake
        ${PROJECT_BINARY_DIR}/ActsConfig.cmake
    DESTINATION ${install_package_config_dir}
)

# For git repositories, patch ActsConfig.cmake at install time with fresh git hash
if(_acts_is_git_repo)
    install(
        CODE
            "
        # Get git hash at install time
        execute_process(
            COMMAND git rev-parse HEAD
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE _install_git_hash
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        execute_process(
            COMMAND git rev-parse --short HEAD
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE _install_git_hash_short
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        execute_process(
            COMMAND git diff-index --quiet HEAD --
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            RESULT_VARIABLE _git_is_dirty
            ERROR_QUIET
        )

        if(NOT _git_is_dirty EQUAL 0)
            set(_install_git_hash \"\${_install_git_hash}-dirty\")
            set(_install_git_hash_short \"\${_install_git_hash_short}-dirty\")
        endif()

        # Read the installed ActsConfig.cmake
        set(_config_file \"${install_package_config_dir}/ActsConfig.cmake\")
        file(READ \"\${_config_file}\" _config_content)

        # Replace the commit hash lines
        string(REGEX REPLACE
            \"set\\\\(Acts_COMMIT_HASH \\\"[^\\\"]*\\\"\\\\)\"
            \"set(Acts_COMMIT_HASH \\\"\${_install_git_hash}\\\")\"
            _config_content \"\${_config_content}\")
        string(REGEX REPLACE
            \"set\\\\(Acts_COMMIT_HASH_SHORT \\\"[^\\\"]*\\\"\\\\)\"
            \"set(Acts_COMMIT_HASH_SHORT \\\"\${_install_git_hash_short}\\\")\"
            _config_content \"\${_config_content}\")

        # Write back the patched content
        file(WRITE \"\${_config_file}\" \"\${_config_content}\")
        "
    )
endif()

# install third party FindXXX.cmake files
file(GLOB_RECURSE _pckg_find_files "${CMAKE_CURRENT_LIST_DIR}/Find*.cmake")
install(
    FILES ${_pckg_find_files}
    DESTINATION ${install_package_config_dir}/Modules
)

# install target configs for all available components
foreach(_component ${_components})
    install(
        EXPORT Acts${_component}Targets
        DESTINATION ${install_package_config_dir}
    )
endforeach()
