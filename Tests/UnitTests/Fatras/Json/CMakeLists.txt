set(unittest_extra_libraries ActsFatrasJson)

add_unittest(
    FatrasSyntheticDescriptionJson
    DetectorDescriptionJsonConverterTests.cpp
)
add_unittest(FatrasSyntheticEventConfigJson EventConfigJsonConverterTests.cpp)
add_unittest(FatrasSyntheticShippedDetectors ShippedDetectorsTests.cpp)
