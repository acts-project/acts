set(unittest_extra_libraries
    ActsExamplesFramework
    ActsExamplesIoRoot
    ActsPluginRoot
)
add_unittest(DataHandle DataHandleTest.cpp)
add_unittest(HistogramFitInterface HistogramFitInterfaceTests.cpp)
add_unittest(MeasurementResidual MeasurementResidualTests.cpp)
add_unittest(ResPlotTool ResPlotToolTests.cpp)
