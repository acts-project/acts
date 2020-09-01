# plugins are manually ordered such that dependent plugins appear after their
# dependencies in the component list.
# independent plugins
add_component_if(Cuda CudaPlugin ACTS_BUILD_CUDA_PLUGIN)
add_component_if(Sycl SyclPlugin ACTS_BUILD_SYCL_PLUGIN)
add_component_if(Digitization DigitizationPlugin ACTS_BUILD_DIGITIZATION_PLUGIN)
add_component_if(Identification IdentificationPlugin ACTS_BUILD_IDENTIFICATION_PLUGIN)
add_component_if(Json JsonPlugin ACTS_BUILD_JSON_PLUGIN)
# legacy behaves like a plugin but is not named like one for historical reasons
add_component_if(Legacy Legacy ACTS_BUILD_LEGACY)
# dependent plugins. depend either on a independent plugins or on one another
add_component_if(TGeo PluginTGeo ACTS_BUILD_PLUGIN_TGEO)
add_component_if(DD4hep PluginDD4hep ACTS_BUILD_PLUGIN_DD4HEP)

propagate_components_to_parent()
