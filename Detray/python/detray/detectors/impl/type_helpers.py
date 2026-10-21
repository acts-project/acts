# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# python includes
from collections import namedtuple

""" The type can be 'single' (index) or 'range' (index range) """
link = namedtuple(
    "link",
    "link_type data_type",
    defaults=["single", None],
)

""" A detray class with a name and its include path """
cpp_class = namedtuple("cpp_class", "specifier path param", defaults=["", "", {}])
