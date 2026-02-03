set(ACTS_ACTSVG_SOURCE
    "URL;https://github.com/acts-project/actsvg/archive/refs/tags/v${_acts_actsvg_version}.tar.gz;URL_HASH;SHA256=75944cbc4444d3099764f19f6fa5fbf9f2b81e3ac776b5874746add74a7cd0f8"
    CACHE STRING
    "Source to take ACTSVG from"
)
mark_as_advanced(ACTS_ACTSVG_SOURCE)

set(ACTS_VECMEM_SOURCE
    "URL;https://github.com/acts-project/vecmem/archive/refs/tags/v${_acts_vecmem_version}.tar.gz;URL_HASH;SHA256=a6e30da7dcb9a792eb0f2c0ceb821f52cf867c9745a65468b2ce6f026b23a49d"
    CACHE STRING
    "Source to take VECMEM from"
)
mark_as_advanced(ACTS_VECMEM_SOURCE)

set(ACTS_ALGEBRAPLUGINS_SOURCE
    "URL;https://github.com/acts-project/algebra-plugins/archive/refs/tags/v${_acts_algebraplugins_version}.tar.gz;URL_HASH;SHA256=d798ba2129bf922f54627233ef947b8bb2345db9199e3868cc48bc1da86d5f15"
    CACHE STRING
    "Source to take ALGEBRAPLUGINS from"
)
mark_as_advanced(ACTS_ALGEBRAPLUGINS_SOURCE)

set(ACTS_COVFIE_SOURCE
    "URL;https://github.com/acts-project/covfie/archive/refs/tags/v${_acts_covfie_version}.tar.gz;URL_HASH;SHA256=6eff65e05118d3007c689e3529a62bb1674348ac1b0f0f32afd953c62d1b8890"
    CACHE STRING
    "Source to take COVFIE from"
)
mark_as_advanced(ACTS_COVFIE_SOURCE)

set(ACTS_DETRAY_SOURCE
    "URL;https://github.com/acts-project/detray/archive/refs/tags/v${_acts_detray_version}.tar.gz;URL_HASH;SHA256=7b3d3c94cf42be7450e9fe008b567a2f425e6f1986b61d8a3a66814383599043"
    CACHE STRING
    "Source to take DETRAY from"
)
mark_as_advanced(ACTS_DETRAY_SOURCE)

set(ACTS_TRACCC_SOURCE
    "URL;https://github.com/acts-project/traccc/archive/refs/tags/v${_acts_traccc_version}.tar.gz;URL_HASH;SHA256=d2dd856a83503ab452936f79ebb2496ce59e2f63984b0e9403d1c881546b8872"
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
    "GIT_REPOSITORY;https://gitlab.cern.ch/acts/OpenDataDetector.git;GIT_TAG;v4.0.4"
    CACHE STRING
    "Source to take OpenDataDetector from"
)
mark_as_advanced(ACTS_ODD_SOURCE)

set(ACTS_MODULEMAPGRAPH_SOURCE
    "GIT_REPOSITORY;https://gitlab.cern.ch/gnn4itkteam/ModuleMapGraph;GIT_TAG;1.1.16"
    CACHE STRING
    "Source to take ModuleMapGraph from"
)
mark_as_advanced(ACTS_MODULEMAPGRAPH_SOURCE)
