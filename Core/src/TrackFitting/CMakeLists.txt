target_sources(
    ActsCore
    PRIVATE
        KalmanFitterError.cpp
        GainMatrixUpdater.cpp
        GainMatrixSmoother.cpp
        GlobalChiSquareFitterError.cpp
        GsfError.cpp
        GsfUtils.cpp
        BetheHeitlerApprox.cpp
        GsfMixtureReduction.cpp
        GlobalChiSquareFitter.cpp
        MbfSmoother.cpp
        detail/GsfComponentMerging.cpp
)

foreach(DIM RANGE 1 6)
    set(dim_file
        ${CMAKE_CURRENT_BINARY_DIR}/detail/GainMatrixUpdaterImpl${DIM}.cpp
    )
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/detail/GainMatrixUpdaterImpl.cpp.in
        ${dim_file}
        @ONLY
    )

    set_source_files_properties(
        ${dim_file}
        PROPERTIES
            COMPILE_DEFINITIONS ACTS_GAIN_MATRIX_UPDATER_INSTANTIATE=${DIM}
    )
    target_sources(ActsCore PRIVATE ${dim_file})
endforeach()

foreach(DIM RANGE 1 6)
    set(dim_file ${CMAKE_CURRENT_BINARY_DIR}/detail/MbfSmootherImpl${DIM}.cpp)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/detail/MbfSmootherImpl.cpp.in
        ${dim_file}
        @ONLY
    )

    set_source_files_properties(
        ${dim_file}
        PROPERTIES COMPILE_DEFINITIONS ACTS_MBF_SMOOTHER_INSTANTIATE=${DIM}
    )
    target_sources(ActsCore PRIVATE ${dim_file})
endforeach()
