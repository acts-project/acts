set(unittest_extra_libraries ActsExamplesDigitization ActsExamplesIoJson)

add_unittest(JsonDigitizationConfig JsonDigitizationConfigTests.cpp)

add_unittest(
    TrackingGeometryMaterialJsonWriter
    TrackingGeometryMaterialJsonWriterTests.cpp
)
