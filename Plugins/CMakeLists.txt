# plugins are manually ordered such that dependent plugins appear after their
# dependencies in the component list.
# independent plugins
add_component_if(ActSVG PluginActSVG ACTS_BUILD_PLUGIN_ACTSVG)
add_component_if(FpeMonitoring PluginFpeMonitoring ACTS_BUILD_PLUGIN_FPEMON)
add_component_if(Geant4 PluginGeant4 ACTS_BUILD_PLUGIN_GEANT4)
add_component_if(GeoModel PluginGeoModel ACTS_BUILD_PLUGIN_GEOMODEL)
add_component_if(Json PluginJson ACTS_BUILD_PLUGIN_JSON)
add_component_if(Mille PluginMille ACTS_BUILD_PLUGIN_MILLE)
add_component_if(Onnx PluginOnnx ACTS_BUILD_PLUGIN_ONNX)
add_component_if(Gnn PluginGnn ACTS_BUILD_PLUGIN_GNN)
add_component_if(Detray PluginDetray ACTS_BUILD_PLUGIN_TRACCC)
add_component_if(Covfie PluginCovfie ACTS_BUILD_PLUGIN_TRACCC)

# dependent plugins. depend either on a independent plugins or on one another
add_component_if(Root PluginRoot ACTS_BUILD_PLUGIN_ROOT)
add_component_if(DD4hep PluginDD4hep ACTS_BUILD_PLUGIN_DD4HEP)
add_component_if(EDM4hep PluginEDM4hep ACTS_BUILD_PLUGIN_EDM4HEP)

# Arrow must be registered last. When the examples are built, the Arrow plugin
# absorbs the Examples/Io/{Arrow,Parquet} converter object files and therefore
# links Acts::ExamplesFramework, which is dragged into the ActsPluginArrowTargets
# export set. ExamplesFramework publicly links other Acts components
# (Acts::Fatras, Acts::PluginFpeMonitoring, and Acts::PluginOnnx when enabled),
# so the generated ActsPluginArrowTargets.cmake guards that those targets already
# exist. Components are loaded in registration order by ActsConfig.cmake, so
# PluginArrow has to come after every component ExamplesFramework can reference;
# keeping it last guarantees that. See Examples/Io/Parquet/CMakeLists.txt.
add_component_if(Arrow PluginArrow ACTS_BUILD_PLUGIN_ARROW)

propagate_components_to_parent()
