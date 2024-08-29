set(ACTS_ACTSVG_SOURCE
    "URL;https://github.com/acts-project/actsvg/archive/refs/tags/v0.4.47.tar.gz;URL_HASH;SHA256=aeb3927f8db1c7c9b8afa76adaa720013aa03a01032c8d5f457f623e2b04a172"
    CACHE STRING
    "Source to take ACTSVG from"
)
mark_as_advanced(ACTS_ACTSVG_SOURCE)

set(ACTS_VECMEM_SOURCE
    "URL;https://github.com/acts-project/vecmem/archive/refs/tags/v1.4.0.tar.gz;URL_HASH;SHA256=545dfb4de4f9f3d773eef6a0e3297ebf981bb81950930d0991ad739e31ab16af"
    CACHE STRING
    "Source to take VECMEM from"
)
mark_as_advanced(ACTS_VECMEM_SOURCE)

set(ACTS_ALGEBRAPLUGINS_SOURCE
    "URL;https://github.com/acts-project/algebra-plugins/archive/refs/tags/v0.22.0.tar.gz;URL_HASH;SHA256=6fde02181c1b856c0a17a1925f0969798eecd5e3d6f2a87ea2eb365b6c948cc1"
    CACHE STRING
    "Source to take ALGEBRAPLUGINS from"
)
mark_as_advanced(ACTS_ALGEBRAPLUGINS_SOURCE)

set(ACTS_COVFIE_SOURCE
    "URL;https://github.com/acts-project/covfie/archive/refs/tags/v0.10.0.tar.gz;URL_HASH;SHA256=d44142b302ffc193ad2229f1d2cc6d8d720dd9da8c37989ada4f23018f86c964"
    CACHE STRING
    "Source to take COVFIE from"
)
mark_as_advanced(ACTS_COVFIE_SOURCE)

set(ACTS_DETRAY_SOURCE
    "URL;https://github.com/acts-project/detray/archive/refs/tags/v0.72.1.tar.gz;URL_HASH;SHA256=6cc8d34bc0d801338e9ab142c4a9884d19d9c02555dbb56972fab86b98d0f75b"
    CACHE STRING
    "Source to take DETRAY from"
)
mark_as_advanced(ACTS_DETRAY_SOURCE)

set(ACTS_TRACCC_SOURCE
    "URL;https://github.com/acts-project/traccc/archive/refs/tags/v0.15.0.tar.gz;URL_HASH;SHA256=1a9a1d0d2f6c4a7773eae3119b1e044fb52031ca49dfc88e6dc4ab8a11df167e"
    CACHE STRING
    "Source to take TRACCC from"
)
mark_as_advanced(ACTS_TRACCC_SOURCE)

set(ACTS_FRNN_SOURCE
    "GIT_REPOSITORY;https://github.com/hrzhao76/FRNN/;GIT_TAG;5f8a48b0022300cd2863119f5646a5f31373e0c8"
    CACHE STRING
    "Source to take FRNN from"
)
mark_as_advanced(ACTS_FRNN_SOURCE)

set(ACTS_NLOHMANNJSON_SOURCE
    "URL;https://github.com/nlohmann/json/archive/refs/tags/v3.10.5.tar.gz;URL_HASH;SHA256=5daca6ca216495edf89d167f808d1d03c4a4d929cef7da5e10f135ae1540c7e4"
    CACHE STRING
    "Source to take nlohmann_json from"
)
mark_as_advanced(ACTS_NLOHMANN_JSON_SOURCE)

set(ACTS_EIGEN3_SOURCE
    "URL;https://gitlab.com/libeigen/eigen/-/archive/${_acts_eigen3_version}/${_acts_eigen3_version}.tar.gz;URL_HASH;SHA256=ba6ef66ba2d319e0a871a267889411c550d4bdf5bc7c62f86c60276913f3f4ba"
    CACHE STRING
    "Source to take eigen3 from"
)
mark_as_advanced(ACTS_EIGEN3_SOURCE)

set(ACTS_PYBIND11_SOURCE
    "GIT_REPOSITORY;https://github.com/pybind/pybind11.git;GIT_TAG;v2.13.1"
    CACHE STRING
    "Source to take pybind11 from"
)
mark_as_advanced(ACTS_PYBIND11_SOURCE)
