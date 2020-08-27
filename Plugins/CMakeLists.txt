# plugins are manually ordered such that dependent plugins appear after their
# dependencies in the component list.
# independent plugins
add_component_if(Cuda CudaPlugin ACTS_BUILD_PLUGIN_CUDA)
add_component_if(Digitization DigitizationPlugin ACTS_BUILD_PLUGIN_DIGITIZATION)
add_component_if(Identification IdentificationPlugin ACTS_BUILD_PLUGIN_IDENTIFICATION)
add_component_if(Json JsonPlugin ACTS_BUILD_PLUGIN_JSON)
# legacy behaves like a plugin but is not named like one for historical reasons
add_component_if(Legacy Legacy ACTS_BUILD_PLUGIN_LEGACY)
# dependent plugins. depend either on a independent plugins or on one another
add_component_if(TGeo TGeoPlugin ACTS_BUILD_PLUGIN_TGEO)
add_component_if(DD4hep DD4hepPlugin ACTS_BUILD_PLUGIN_DD4HEP)

propagate_components_to_parent()
