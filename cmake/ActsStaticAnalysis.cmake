if(ACTS_RUN_CLANG_TIDY)
    find_program(CLANG_TIDY_COMMAND NAMES clang-tidy)

    if(NOT CLANG_TIDY_COMMAND)
        message(
            WARNING
            "ACTS_RUN_CLANG_TIDY is ON but clang-tidy is not found!"
        )
        set(ACTS_CXX_CLANG_TIDY "" CACHE STRING "" FORCE)
    else()
        message(STATUS "Setting up clang-tidy run")

        set(CLANG_TIDY_HEADER_FILTER ".*")

        set(ACTS_CXX_CLANG_TIDY
            "${CLANG_TIDY_COMMAND};-header-filter=${CLANG_TIDY_HEADER_FILTER};-config-file=${CMAKE_SOURCE_DIR}/.clang-tidy"
        )
    endif()
endif()

if(ACTS_RUN_CPPCHECK)
    find_program(CPPCHECK_COMMAND NAMES cppcheck)

    if(NOT CPPCHECK_COMMAND)
        message(
            WARNING
            "ACTS_RUN_CPPCHECK is ON but cppcheck is not found!"
        )
        set(ACTS_CXX_CPPCHECK "" CACHE STRING "" FORCE)
    else()
        message(STATUS "Setting up cppcheck run")

        set(ACTS_CXX_CPPCHECK
            "${CPPCHECK_COMMAND};--enable=warning,performance,portability;--std=c++20;--suppress=missingIncludeSystem;--error-exitcode=1"
        )
    endif()
endif()

function(acts_enable_static_analysis)
    set(CMAKE_CXX_CLANG_TIDY "${ACTS_CXX_CLANG_TIDY}" CACHE STRING "" FORCE)
    set(CMAKE_CXX_CPPCHECK "${ACTS_CXX_CPPCHECK}" CACHE STRING "" FORCE)
endfunction()

function(acts_disable_static_analysis)
    set(CMAKE_CXX_CLANG_TIDY "" CACHE STRING "" FORCE)
    set(CMAKE_CXX_CPPCHECK "" CACHE STRING "" FORCE)
endfunction()
