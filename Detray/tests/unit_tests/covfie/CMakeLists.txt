# Detray library, part of the ACTS project (R&D line)
#
# (c) 2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the covfie tests.
detray_add_unit_test( covfie
   "constant_field.cpp"
   LINK_LIBRARIES GTest::gtest_main covfie::core
)
