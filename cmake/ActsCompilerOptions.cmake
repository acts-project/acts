# Acts compiler flags
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        RelWithDebInfo
        CACHE STRING
        "Build type configuration"
        FORCE
    )
    message(STATUS "Setting default build type: ${CMAKE_BUILD_TYPE}")
endif()

# Warning flags are *not* set here: they are applied per-directory by
# acts_apply_warning_flags() from ActsWarningFlags.cmake, so that they only
# affect Acts' own code and not the dependencies built as part of this project.
set(cxx_flags "")

# Add assertions to standard libraries
if(ACTS_FORCE_ASSERTIONS)
    set(cxx_flags "${cxx_flags} -D_GLIBCXX_ASSERTIONS -D_LIBCPP_DEBUG")
endif()

set(ACTS_CXX_STANDARD 20)
set(ACTS_CXX_STANDARD_FEATURE cxx_std_20)
if(DEFINED CMAKE_CXX_STANDARD)
    if(${CMAKE_CXX_STANDARD} GREATER_EQUAL 20)
        set(ACTS_CXX_STANDARD ${CMAKE_CXX_STANDARD})
        set(ACTS_CXX_STANDARD_FEATURE "cxx_std_${CMAKE_CXX_STANDARD}")
    else()
        message(
            SEND_ERROR
            "CMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}, but ACTS requires C++ >=20"
        )
    endif()
endif()

message(STATUS "C++ standard: ${ACTS_CXX_STANDARD}")

if(ACTS_ENABLE_CPU_PROFILING OR ACTS_ENABLE_MEMORY_PROFILING)
    message(STATUS "Added debug symbol compile flag")
    set(cxx_flags "${cxx_flags} ${CMAKE_CXX_FLAGS_DEBUG_INIT}")
endif()

# assign to global CXX flags
set(CMAKE_CXX_FLAGS "${cxx_flags} ${CMAKE_CXX_FLAGS}")
message(STATUS "Using compiler flags: ${CMAKE_CXX_FLAGS}")

# do not scan for C++ modules
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

# silence warning about missing RPATH on Mac OSX
set(CMAKE_MACOSX_RPATH 1)

# bake where we found external dependencies, if they
# were not in the default library directories
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

message(CHECK_START "Checking C++20 std::format support")
try_compile(
    HAVE_STDFORMAT
    SOURCES ${PROJECT_SOURCE_DIR}/cmake/src/format.cpp
    CXX_STANDARD 20
)

if(NOT HAVE_STDFORMAT)
    message(CHECK_FAIL "missing")
    message(
        SEND_ERROR
        "C++20 std::format support is required, "
        "but not available in this compiler"
    )
else()
    message(CHECK_PASS "ok")
endif()
