# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Originally taken from https://github.com/facebook/folly/blob/main/CMake/FindZstd.cmake

#
# - Try to find Facebook zstd library
# This will define
# ZSTD_FOUND
# ZSTD_INCLUDE_DIR
# ZSTD_LIBRARY
#

find_path(ZSTD_INCLUDE_DIR NAMES zstd.h)

find_library(ZSTD_LIBRARY_DEBUG NAMES zstdd zstd_staticd)
find_library(ZSTD_LIBRARY_RELEASE NAMES zstd zstd_static)

include(SelectLibraryConfigurations)
select_library_configurations(ZSTD)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    ZSTD
    DEFAULT_MSG
    ZSTD_LIBRARY
    ZSTD_INCLUDE_DIR
)

if(ZSTD_FOUND)
    add_library(ZSTD UNKNOWN IMPORTED)
    set_target_properties(ZSTD PROPERTIES IMPORTED_LOCATION ${ZSTD_LIBRARY})
    target_include_directories(ZSTD INTERFACE ${ZSTD_INCLUDE_DIR})
    add_library(ZSTD::ZSTD ALIAS ZSTD)
endif()

mark_as_advanced(ZSTD_INCLUDE_DIR ZSTD_LIBRARY)
