set(_fpe_platform_source src/FpeMonitorPlatformUnsupported.cpp)
string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" _acts_fpe_arch)

if(APPLE)
    if(_acts_fpe_arch MATCHES "^(x86_64|amd64)$")
        set(_fpe_platform_source src/FpeMonitorPlatformDarwinX86_64.cpp)
    elseif(_acts_fpe_arch MATCHES "^(arm64|aarch64)$")
        set(_fpe_platform_source src/FpeMonitorPlatformDarwinArm64.cpp)
    endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    if(_acts_fpe_arch MATCHES "^(x86_64|amd64)$")
        set(_fpe_platform_source src/FpeMonitorPlatformLinuxX86_64.cpp)
    endif()
endif()

acts_add_library(
    PluginFpeMonitoring
    src/FpeMonitor.cpp
    ${_fpe_platform_source}
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

target_include_directories(
    ActsPluginFpeMonitoring
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(ActsPluginFpeMonitoring PUBLIC Acts::Core)

acts_compile_headers(PluginFpeMonitoring GLOB include/**/*.hpp)

# FPE compile definitions
set(_fpe_definitions "")

# In case we're going to run clang-tidy, let's unconditionally turn off backtrace
if(ACTS_RUN_CLANG_TIDY)
    list(APPEND _fpe_definitions BOOST_STACKTRACE_USE_NOOP=1)
else()
    find_package(Backtrace)
    find_program(addr2line_EXECUTABLE addr2line)
    if(APPLE)
        list(APPEND _fpe_definitions _GNU_SOURCE)
    elseif(CMAKE_DL_LIBS)
        target_link_libraries(ActsPluginFpeMonitoring PUBLIC ${CMAKE_DL_LIBS})

        set(_backtrace_setup_complete FALSE)
        set(_backtrace_src
            "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/backtrace.cpp"
        )
        file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp")

        find_path(
            boost_stacktrace_include
            NAMES "boost/stacktrace.hpp"
            REQUIRED
            PATHS ${Boost_INCLUDE_DIRS}
        )

        file(
            WRITE "${_backtrace_src}"
            "#include <boost/stacktrace.hpp>\n"
            "#include <iostream>\n"
            "int main() { std::cout << boost::stacktrace::stacktrace(); }\n"
        )

        message(
            CHECK_START
            "Does backtrace work without libbacktrace linker flag"
        )
        try_compile(
            _backtrace_nolink
            "${CMAKE_CURRENT_BINARY_DIR}"
            "${_backtrace_src}"
            LINK_LIBRARIES ${CMAKE_DL_LIBS}
            CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${boost_stacktrace_include}"
            COMPILE_DEFINITIONS BOOST_STACKTRACE_USE_BACKTRACE
            OUTPUT_VARIABLE __OUTPUT
        )

        if(_backtrace_nolink)
            message(CHECK_PASS "yes")
            list(APPEND _fpe_definitions BOOST_STACKTRACE_USE_BACKTRACE)
            set(_backtrace_setup_complete TRUE)
        else()
            message(CHECK_FAIL "no")

            file(GLOB hints "/usr/lib/gcc/*/*/include")
            find_file(backtrace_header NAMES "backtrace.h" HINTS ${hints})

            if(backtrace_header STREQUAL "backtrace_header-NOTFOUND")
                message(STATUS "Could not find backtrace header file")
            else()
                message(CHECK_START "Does backtrace work with explicit include")
                file(
                    WRITE "${_backtrace_src}"
                    "#include <boost/stacktrace.hpp>\n"
                    "#include <iostream>\n"
                    "int main() { std::cout << boost::stacktrace::stacktrace(); }\n"
                )
                try_compile(
                    _backtrace_explicit
                    "${CMAKE_CURRENT_BINARY_DIR}"
                    "${_backtrace_src}"
                    LINK_LIBRARIES ${CMAKE_DL_LIBS}
                    CMAKE_FLAGS
                        "-DINCLUDE_DIRECTORIES=${boost_stacktrace_include}"
                    COMPILE_DEFINITIONS BOOST_STACKTRACE_USE_BACKTRACE
                    OUTPUT_VARIABLE __OUTPUT
                )

                if(_backtrace_explicit)
                    message(CHECK_PASS "yes")
                    list(APPEND _fpe_definitions BOOST_STACKTRACE_USE_BACKTRACE)
                    set(_backtrace_setup_complete TRUE)
                else()
                    message(CHECK_FAIL "no")

                    set(_backtrace_include_def
                        "BOOST_STACKTRACE_BACKTRACE_INCLUDE_FILE=\"${backtrace_header}\""
                    )

                    file(
                        WRITE "${_backtrace_src}"
                        "#include <boost/stacktrace.hpp>\n"
                        "int main() {}\n"
                    )

                    message(
                        CHECK_START
                        "Does backtrace work with explicit include"
                    )
                    try_compile(
                        _backtrace_explicit_header
                        "${CMAKE_CURRENT_BINARY_DIR}"
                        "${_backtrace_src}"
                        LINK_LIBRARIES ${CMAKE_DL_LIBS}
                        CMAKE_FLAGS
                            "-DINCLUDE_DIRECTORIES=${boost_stacktrace_include}"
                        COMPILE_DEFINITIONS
                            BOOST_STACKTRACE_USE_BACKTRACE
                            ${_backtrace_include_def}
                        OUTPUT_VARIABLE __OUTPUT
                    )

                    if(_backtrace_explicit_header)
                        message(CHECK_PASS "yes")
                        list(
                            APPEND _fpe_definitions
                            BOOST_STACKTRACE_USE_BACKTRACE
                            ${_backtrace_include_def}
                        )
                        set(_backtrace_setup_complete TRUE)
                    else()
                        message(CHECK_FAIL "no")

                        file(
                            WRITE "${_backtrace_src}"
                            "#include <boost/stacktrace.hpp>\n"
                            "#include <iostream>\n"
                            "int main() { std::cout << boost::stacktrace::stacktrace(); }\n"
                        )

                        message(
                            CHECK_START
                            "Does backtrace work without linker flag"
                        )
                        try_compile(
                            _backtrace_nolink
                            "${CMAKE_CURRENT_BINARY_DIR}"
                            "${_backtrace_src}"
                            LINK_LIBRARIES ${CMAKE_DL_LIBS}
                            CMAKE_FLAGS
                                "-DINCLUDE_DIRECTORIES=${boost_stacktrace_include}"
                            COMPILE_DEFINITIONS
                                BOOST_STACKTRACE_USE_BACKTRACE
                                ${_backtrace_include_def}
                            OUTPUT_VARIABLE __OUTPUT
                        )

                        if(_backtrace_nolink)
                            message(CHECK_PASS "yes")
                            list(
                                APPEND _fpe_definitions
                                BOOST_STACKTRACE_USE_BACKTRACE
                                ${_backtrace_include_def}
                            )
                            set(_backtrace_setup_complete TRUE)
                        else()
                            message(CHECK_FAIL "no")

                            file(
                                WRITE "${_backtrace_src}"
                                "#include <boost/stacktrace.hpp>\n"
                                "#include <iostream>\n"
                                "int main() { std::cout << boost::stacktrace::stacktrace(); }\n"
                            )

                            message(
                                CHECK_START
                                "Does backtrace work with linker flag"
                            )
                            try_compile(
                                _backtrace_link
                                "${CMAKE_CURRENT_BINARY_DIR}"
                                "${_backtrace_src}"
                                LINK_LIBRARIES backtrace ${CMAKE_DL_LIBS}
                                CMAKE_FLAGS
                                    "-DINCLUDE_DIRECTORIES=${boost_stacktrace_include}"
                                COMPILE_DEFINITIONS
                                    BOOST_STACKTRACE_USE_BACKTRACE
                                    ${_backtrace_include_def}
                                OUTPUT_VARIABLE __OUTPUT
                            )

                            if(_backtrace_link)
                                message(CHECK_PASS "yes")
                                list(
                                    APPEND _fpe_definitions
                                    BOOST_STACKTRACE_USE_BACKTRACE
                                    ${_backtrace_include_def}
                                )
                                target_link_libraries(
                                    ActsPluginFpeMonitoring
                                    PUBLIC backtrace
                                )

                                set(_backtrace_setup_complete TRUE)
                            else()
                                message(CHECK_FAIL "no")
                            endif()
                        endif()
                    endif()
                endif()
            endif()
        endif()

        if(NOT _backtrace_setup_complete)
            message(
                STATUS
                "libdl available but backtrace backend probe failed: use boost basic fallback"
            )
        endif()
    endif()
endif()

if(_fpe_definitions)
    target_compile_definitions(
        ActsPluginFpeMonitoring
        PUBLIC ${_fpe_definitions}
    )
endif()
