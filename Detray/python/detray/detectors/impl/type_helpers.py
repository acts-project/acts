# Detray library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

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
