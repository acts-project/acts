set(unittest_extra_libraries
    ActsExamplesFramework
    ActsExamplesIoRoot
    ActsPluginRoot
)
add_unittest(DataHandle DataHandleTest.cpp)
add_unittest(HistogramFitInterface HistogramFitInterfaceTests.cpp)
