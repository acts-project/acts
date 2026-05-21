set(ACTS_ACTSVG_SOURCE
    "URL;https://github.com/acts-project/actsvg/archive/refs/tags/v${_acts_actsvg_version}.tar.gz;URL_HASH;SHA256=75944cbc4444d3099764f19f6fa5fbf9f2b81e3ac776b5874746add74a7cd0f8"
    CACHE STRING
    "Source to take ACTSVG from"
)
mark_as_advanced(ACTS_ACTSVG_SOURCE)

set(ACTS_VECMEM_SOURCE
    "URL;https://github.com/acts-project/vecmem/archive/refs/tags/v${_acts_vecmem_version}.tar.gz;URL_HASH;SHA256=a1dd195e154ed23a0e50c52e22fb9f986fc65cd99860020fc47a292f597fa88d"
    CACHE STRING
    "Source to take VECMEM from"
)
mark_as_advanced(ACTS_VECMEM_SOURCE)

set(ACTS_COVFIE_SOURCE
    "URL;https://github.com/acts-project/covfie/archive/refs/tags/v${_acts_covfie_version}.tar.gz;URL_HASH;SHA256=6eff65e05118d3007c689e3529a62bb1674348ac1b0f0f32afd953c62d1b8890"
    CACHE STRING
    "Source to take COVFIE from"
)
mark_as_advanced(ACTS_COVFIE_SOURCE)

set(ACTS_TRACCC_SOURCE
    "URL;https://github.com/acts-project/traccc/archive/refs/tags/v${_acts_traccc_version}.tar.gz;URL_HASH;SHA256=9c7e772a088cadfbdcb63abe51e4064241c84f10fbfc27594cbd0cc78d140d64"
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
    "URL;https://github.com/nlohmann/json/archive/refs/tags/v${_acts_nlohmanjson_version}.tar.gz;URL_HASH;SHA256=0d8ef5af7f9794e3263480193c491549b2ba6cc74bb018906202ada498a79406"
    CACHE STRING
    "Source to take nlohmann_json from"
)
mark_as_advanced(ACTS_NLOHMANN_JSON_SOURCE)

# translate version string to the historical Mille release naming convention
string(REPLACE "." "-" _acts_mille_release_string ${_acts_mille_version})
set(ACTS_MILLE_SOURCE
    "URL;https://gitlab.desy.de/millepede/mille/-/archive/V${_acts_mille_release_string}/mille-V${_acts_mille_release_string}.tar.gz;URL_HASH;SHA256=ae4bf37de8d835aa8adc2960bb795a2080233a4c8af3d4b55adf395e20df0f3e"
    CACHE STRING
    "Source to take Mille from"
)
mark_as_advanced(ACTS_MILLE_SOURCE)

# translate version string to the historical GBL release naming convention
string(REPLACE "." "-" _acts_gbl_release_string ${_acts_gbl_version})
set(ACTS_GBL_SOURCE
    "URL;https://gitlab.desy.de/millepede/general-broken-lines/-/archive/V${_acts_gbl_release_string}/general-broken-lines-V${_acts_gbl_release_string}.tar.gz;URL_HASH;SHA256=e40401a77a828c81a9217d8df3201e7712ac86b4fd5058d526ae1e1f6664304f"
    CACHE STRING
    "Source to take General Broken Lines (GBL) from"
)
mark_as_advanced(ACTS_GBL_SOURCE)
# translate version string to the historical Millepede release naming convention
string(REPLACE "." "-" _acts_mp2_release_string ${_acts_mp2_version})
set(ACTS_MP2_SOURCE
    "URL;https://gitlab.desy.de/millepede/millepede-ii/-/archive/V${_acts_mp2_release_string}/millepede-ii-V${_acts_mp2_release_string}.tar.gz;URL_HASH;SHA256=b6a316e4b1ebf93cbf72ddd57a157e09f4446e4677352ef288748731ac2c0297"
    CACHE STRING
    "Source to take Millepede-II from"
)
mark_as_advanced(ACTS_MP2_SOURCE)

set(ACTS_EIGEN3_SOURCE
    "URL;https://gitlab.com/libeigen/eigen/-/archive/${_acts_eigen3_version}/${_acts_eigen3_version}.tar.gz;URL_HASH;SHA256=ba6ef66ba2d319e0a871a267889411c550d4bdf5bc7c62f86c60276913f3f4ba"
    CACHE STRING
    "Source to take eigen3 from"
)
mark_as_advanced(ACTS_EIGEN3_SOURCE)

set(ACTS_PYBIND11_SOURCE
    "GIT_REPOSITORY;https://github.com/pybind/pybind11.git;GIT_TAG;v${_acts_pybind11_version}"
    CACHE STRING
    "Source to take pybind11 from"
)
mark_as_advanced(ACTS_PYBIND11_SOURCE)

set(ACTS_ANNOY_SOURCE
    "URL;https://github.com/spotify/annoy/archive/refs/tags/v${_acts_annoy_version}.tar.gz;URL_HASH;SHA256=c121d38cacd98f5103b24ca4e94ca097f18179eed3037e9eb93ad70ec1e6356e"
    CACHE STRING
    "Source to take Annoy from"
)
mark_as_advanced(ACTS_ANNOY_SOURCE)

set(ACTS_ODD_SOURCE
    "GIT_REPOSITORY;https://github.com/OpenDataDetector/OpenDataDetector.git;GIT_TAG;v4.0.4"
    CACHE STRING
    "Source to take OpenDataDetector from"
)
mark_as_advanced(ACTS_ODD_SOURCE)

set(ACTS_MODULEMAPGRAPH_SOURCE
    "GIT_REPOSITORY;https://gitlab.cern.ch/gnn4itkteam/ModuleMapGraph;GIT_TAG;1.4.1"
    CACHE STRING
    "Source to take ModuleMapGraph from"
)
mark_as_advanced(ACTS_MODULEMAPGRAPH_SOURCE)
