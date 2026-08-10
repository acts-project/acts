set(unittest_extra_libraries
    ActsExamplesFramework
    ActsExamplesIoRoot
    ActsPluginRoot
)
add_unittest(DataHandle DataHandleTest.cpp)
add_unittest(GaussianHistogramFit GaussianHistogramFitTests.cpp)
add_unittest(GaussianHistogramFitRootBaseline GaussianHistogramFitRootBaselineTests.cpp)
