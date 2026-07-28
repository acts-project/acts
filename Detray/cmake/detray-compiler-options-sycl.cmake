# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 3.21)

# All SYCL warning flags come from the shared module. Note that CMake has no
# COMPILE_WARNING_AS_ERROR implementation for SYCL (it is a vecmem-provided
# pseudo-language), so the module emits -Werror for it by hand.
acts_apply_warning_flags(PROFILE detray LANGUAGES SYCL)
