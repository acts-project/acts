# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/home/aicha/Atlas/acts/acts_build/_deps/dfelibs-src")
  file(MAKE_DIRECTORY "/home/aicha/Atlas/acts/acts_build/_deps/dfelibs-src")
endif()
file(MAKE_DIRECTORY
  "/home/aicha/Atlas/acts/acts_build/_deps/dfelibs-build"
  "/home/aicha/Atlas/acts/acts_build/_deps/dfelibs-subbuild/dfelibs-populate-prefix"
  "/home/aicha/Atlas/acts/acts_build/_deps/dfelibs-subbuild/dfelibs-populate-prefix/tmp"
  "/home/aicha/Atlas/acts/acts_build/_deps/dfelibs-subbuild/dfelibs-populate-prefix/src/dfelibs-populate-stamp"
  "/home/aicha/Atlas/acts/acts_build/_deps/dfelibs-subbuild/dfelibs-populate-prefix/src"
  "/home/aicha/Atlas/acts/acts_build/_deps/dfelibs-subbuild/dfelibs-populate-prefix/src/dfelibs-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/aicha/Atlas/acts/acts_build/_deps/dfelibs-subbuild/dfelibs-populate-prefix/src/dfelibs-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/aicha/Atlas/acts/acts_build/_deps/dfelibs-subbuild/dfelibs-populate-prefix/src/dfelibs-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
