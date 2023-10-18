# Acts compiler flags
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type configuration" FORCE)
  message(STATUS "Setting default build type: ${CMAKE_BUILD_TYPE}")
endif() 

set(ACTS_CXX_FLAGS "-Wall -Wextra -Wpedantic -Wshadow -Wno-unused-local-typedefs")
set(ACTS_CXX_FLAGS_DEBUG "--coverage")
set(ACTS_CXX_FLAGS_MINSIZEREL "")
set(ACTS_CXX_FLAGS_RELEASE "")
set(ACTS_CXX_FLAGS_RELWITHDEBINFO "")

set(ACTS_CXX_STANDARD 17)
set(ACTS_CXX_STANDARD_FEATURE cxx_std_17)
if(DEFINED CMAKE_CXX_STANDARD)
  if(${CMAKE_CXX_STANDARD} GREATER_EQUAL 17)
    set(ACTS_CXX_STANDARD ${CMAKE_CXX_STANDARD})
    set(ACTS_CXX_STANDARD_FEATURE "cxx_std_${CMAKE_CXX_STANDARD}")
  else()
    message(ERROR "CMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}, but ACTS requires C++ >=17")
  endif()
endif()

# This adds some useful conversion checks like float-to-bool, float-to-int, etc.
# However, at the moment this is only added to clang builds, since GCC's -Wfloat-conversion 
# is much more aggressive and also triggers on e.g., double-to-float
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang") 
 set(ACTS_CXX_FLAGS "${ACTS_CXX_FLAGS} -Wfloat-conversion")
endif()

if(ACTS_ENABLE_CPU_PROFILING OR ACTS_ENABLE_MEMORY_PROFILING)
  message(STATUS "Added -g compile flag")
  set(ACTS_CXX_FLAGS "${ACTS_CXX_FLAGS} -g")
endif()

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

# bake where we found external dependencies, if they were not in the default library directories
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# set relative library path for ACTS libraries
set(CMAKE_INSTALL_RPATH "\$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")

if(${ACTS_CXX_STANDARD} GREATER_EQUAL 20)
  file(WRITE
    "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/concepts.cpp"
    "#include <concepts>\n"
    "template<class T, class U>\n"
    "concept Derived = std::is_base_of<U, T>::value;\n"
    "struct A {}; struct B : public A {};"
    "int main() { static_assert(Derived<B, A>, \"works\");  }\n" )

  message(CHECK_START "Are C++20 concepts supported")
  try_compile(ACTS_CONCEPTS_SUPPORTED "${CMAKE_BINARY_DIR}"
      "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/concepts.cpp"
      CXX_STANDARD 20
      OUTPUT_VARIABLE __OUTPUT)

  if(ACTS_CONCEPTS_SUPPORTED)
    message(CHECK_PASS "yes")
  else()
    message(CHECK_FAIL "no")
  endif()
endif()
