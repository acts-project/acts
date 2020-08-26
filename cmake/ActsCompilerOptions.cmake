# Acts compiler flags
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type configuration" FORCE)
  message(STATUS "Setting default build type: ${CMAKE_BUILD_TYPE}")
endif() 

set(ACTS_CXX_FLAGS "-Wall -Wextra -Wpedantic -Wshadow -Wunused-local-typedefs")
set(ACTS_CXX_FLAGS_DEBUG "--coverage")
set(ACTS_CXX_FLAGS_MINSIZEREL "")
set(ACTS_CXX_FLAGS_RELEASE "")
set(ACTS_CXX_FLAGS_RELWITHDEBINFO "")

# Acts linker flags
set(ACTS_EXE_LINKER_FLAGS_DEBUG "--coverage")
set(ACTS_SHARED_LINKER_FLAGS_DEBUG "--coverage ")

# assign to global CXX flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ACTS_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${ACTS_CXX_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} ${ACTS_CXX_FLAGS_MINSIZEREL}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${ACTS_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${ACTS_CXX_FLAGS_RELWITHDEBINFO}")

# assign to global linker flags
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${ACTS_EXE_LINKER_FLAGS_DEBUG}")
set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} ${ACTS_SHARED_LINKER_FLAGS_DEBUG}")

# silence warning about missing RPATH on Mac OSX
set(CMAKE_MACOSX_RPATH 1)

if(ACTS_RUN_CLANG_TIDY)
  find_program(CLANG_TIDY_COMMAND NAMES clang-tidy)
  if(NOT CLANG_TIDY_COMMAND)
    message(WARNING "CMake_RUN_CLANG_TIDY is ON but clang-tidy is not found!")
    set(CMAKE_CXX_CLANG_TIDY "" CACHE STRING "" FORCE)
  else()
    message(STATUS "Setting up clang-tidy run")
    set(CLANG_TIDY_CHECKS "-*,readability-*,-readability-redundant-member-init,misc-*,-misc-unused-parameters,bugprone-*,performance-*,modernize-*,-modernize-use-auto,clang-analyzer-deadcode.*,clang-analyzer-*,-clang-analyzer-osx.*,-clang-analyzer-unix.*,cppcoreguidelines-*,-cppcoreguidelines-pro-type-vararg,-cppcoreguidelines-owning-memory,-cppcoreguidelines-pro-bounds-constant-array")
    # performance-move-const-arg
    set(CLANG_TIDY_ERRORS "readability-inconsistent-declaration-parameter-name,readability-named-parameter,readability-container-size-empty,modernize-use-using,readability-braces-around-statements,modernize-use-override,modernize-use-equals-default,readability-implicit-bool-cast,modernize-use-default-member-init,performance-unnecessary-value-param,performance-move-const-arg,modernize-use-equals-default,modernize-use-nullptr")
    set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND};-checks=${CLANG_TIDY_CHECKS};-header-filter='.*';-warnings-as-errors=${CLANG_TIDY_ERRORS}")
  endif()
endif()
