set(unittest_extra_libraries ActsPluginDetray)
add_unittest(DetrayPayloadConverter DetrayPayloadConverterTests.cpp)
target_link_libraries(
    ActsUnitTestDetrayPayloadConverter
    PUBLIC detray::svgtools
)

add_unittest(DetrayMaterialConversion DetrayMaterialConversionTests.cpp)
