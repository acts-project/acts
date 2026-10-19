# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

traccc_add_executable(convert_csv_bfield "convert_csv_bfield.cpp"
    LINK_LIBRARIES traccc::core covfie::core Boost::program_options
)

traccc_add_executable(generate_constant_bfield "generate_constant_bfield.cpp"
    LINK_LIBRARIES traccc::core covfie::core Boost::program_options
)
