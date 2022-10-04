# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

# originally from https://github.com/vector-of-bool/CMakeCM/blob/master/modules/FindFilesystem.cmake

if(TARGET std::filesystem)
    # This module has already been processed. Don't do it again.
    return()
endif()

cmake_minimum_required(VERSION 3.10)

include(CMakePushCheckState)

# If we're not cross-compiling, try to run test executables.
# Otherwise, assume that compile + link is a sufficient check.
if(CMAKE_CROSSCOMPILING)
    include(CheckCXXSourceCompiles)
    macro(_cmcm_check_cxx_source code var)
        check_cxx_source_compiles("${code}" ${var})
    endmacro()
else()
    include(CheckCXXSourceRuns)
    macro(_cmcm_check_cxx_source code var)
        check_cxx_source_runs("${code}" ${var})
    endmacro()
endif()

cmake_push_check_state()

set(CMAKE_REQUIRED_QUIET ${Filesystem_FIND_QUIETLY})

# All of our tests required C++17 or later
set(CMAKE_CXX_STANDARD 17)

set(CXX_FILESYSTEM_HEADER "filesystem" CACHE STRING "The header that should be included to obtain the filesystem APIs")
set(CXX_FILESYSTEM_NAMESPACE "std::filesystem" CACHE STRING "The C++ namespace that contains the filesystem APIs")

set(_found FALSE)

macro(_check_filesystem)
    # We have some filesystem library available. Do link checks
    string(CONFIGURE [[
        #include <cstdlib>
        #include <@CXX_FILESYSTEM_HEADER@>

        int main() {
            auto cwd = @CXX_FILESYSTEM_NAMESPACE@::current_path();
            printf("%s", cwd.c_str());
            return EXIT_SUCCESS;
        }
    ]] _code @ONLY)
    message("code ${_code}")

    set(prev_libraries ${CMAKE_REQUIRED_LIBRARIES})

    if(NOT _found)
        _cmcm_check_cxx_source("${_code}" CXX_FILESYSTEM_NO_LINK_NEEDED)
        set(_found ${CXX_FILESYSTEM_NO_LINK_NEEDED})
    endif()
    if(NOT _found)
        # Add the libstdc++ flag
        set(CMAKE_REQUIRED_LIBRARIES ${prev_libraries} -lstdc++fs)
        _cmcm_check_cxx_source("${_code}" CXX_FILESYSTEM_STDCPPFS_NEEDED)
        set(_found ${CXX_FILESYSTEM_STDCPPFS_NEEDED})
    endif()
    if(NOT _found)
        # Try the libc++ flag
        set(CMAKE_REQUIRED_LIBRARIES ${prev_libraries} -lc++fs)
        _cmcm_check_cxx_source("${_code}" CXX_FILESYSTEM_STDCPPFS_NEEDED)
        set(_found ${CXX_FILESYSTEM_STDCPPFS_NEEDED})
    endif()
endmacro()

_check_filesystem()
if(NOT _found)
    set(_cxx17_flag "-std=c++17")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_cxx17_flag}")
    _check_filesystem()
endif()

if(_found)
    add_library(std::filesystem INTERFACE IMPORTED)
    set_property(TARGET std::filesystem APPEND PROPERTY INTERFACE_COMPILE_FEATURES cxx_std_17)

    if(CXX_FILESYSTEM_NO_LINK_NEEDED)
        # Nothing to add...
    elseif(CXX_FILESYSTEM_STDCPPFS_NEEDED)
        set_property(TARGET std::filesystem APPEND PROPERTY INTERFACE_LINK_LIBRARIES -lstdc++fs)
    elseif(CXX_FILESYSTEM_CPPFS_NEEDED)
        set_property(TARGET std::filesystem APPEND PROPERTY INTERFACE_LINK_LIBRARIES -lc++fs)
    endif()
endif()

cmake_pop_check_state()

set(Filesystem_FOUND ${_found} CACHE BOOL "TRUE if we can run a program using std::filesystem" FORCE)

if(Filesystem_FIND_REQUIRED AND NOT Filesystem_FOUND)
    message(FATAL_ERROR "Cannot run simple program using std::filesystem")
endif()
