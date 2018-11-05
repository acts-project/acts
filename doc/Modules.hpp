// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
// clang-format off

/// @defgroup Design Design and concept descriptions
/// @brief description of general concepts used in Acts

/// @defgroup Logging Debug output options
/// @ingroup Design
/// @brief description of debug output options
///
/// The Acts logging facility supports several severity levels which allow you
/// to control the amount of information displayed at run-time. Logger objects
/// can easily be created using the Acts::getDefaultLogger function which
/// should be sufficient to get you started. In case you need more customized
/// debug output, you can make use of the output decorators defined in
/// Acts::Logging or even write your own implementation of
/// Acts::Logging::OutputDecorator.
///
/// In order to add debug messages to your program, you should use the provided
/// macros for the different severity levels:
/// - #ACTS_VERBOSE
/// - #ACTS_DEBUG
/// - #ACTS_INFO
/// - #ACTS_WARNING
/// - #ACTS_ERROR
/// - #ACTS_FATAL
///
/// All of these macros require that a function <tt>logger()</tt> returning a
/// Acts::Logger object is available in the scope in which the macros are used.
/// Inside classes containing an Acts::Logger object as member variable, this
/// could be achieved by providing a private class method called <tt>logger()</tt>
/// (for an example see e.g. Acts::CylinderVolumeBuilder::logger). Inside free
/// functions or member methods with local logger objects, the same effect can
/// be achieved by using the macro #ACTS_LOCAL_LOGGER which is provided for your
/// convenience.
///
/// @par Code example illustrating the usage
/// @code{.cpp}
/// #include <fstream>
/// #include <memory>
/// #include "Acts/Utilities/Logger.hpp"
///
/// void myFunction() {
///    // open the logfile
///    std::ofstream logfile("log.txt");
///
///    // setup a logger instance for >= INFO messages, streaming into the log file
///    // make sure you do NOT call the variable 'logger'
///    std::unique_ptr<const Acts::Logger> myLogger
///        = Acts::getDefaultLogger("MyLogger", Acts::Logging::INFO, &logfile);
///
///    // make sure the Acts debug macros can work with your logger
///    ACTS_LOCAL_LOGGER(myLogger);
///
///    ACTS_VERBOSE("This message will not appear in the logfile.");
///    ACTS_INFO("But this one will: Hello World!");
///
///    // do not forget to close the logfile
///    logfile.close();
/// }
/// @endcode
///
/// In case you are using Acts in another framework which comes with its own
/// logging facility (e.g. Gaudi) you can pipe the logging output from Acts
/// tools and algorithms to your framework's logging system by supplying different
/// implementations of:
/// - Acts::Logging::OutputFilterPolicy (for mapping logging levels)
/// - Acts::Logging::OutputPrintPolicy (for passing the Acts output to your internal logging system)
///
/// Since Acts makes extensive use of Acts::getDefaultLogger to provide
/// sufficient information for debugging, you would need to provide a modified
/// implementation of this function (using your output filter and printing
/// policies) to also pipe this output to your framework.
///
/// Changing the implementation of an already defined function is a non-trivial
/// task in C++. We recommend the following approach using the possibility to inject
/// custom code by pre-loading shared libraries with <tt>LD_PRELOAD</tt>. You
/// need to provide an appropriate implementation for a function of the following
/// signature into a separate source file and compile it in a shared library
///
/// @code{.cpp}
/// namespace Acts {
///   std::unique_ptr<const Logger> getDefaultLogger(const std::string&,
///                                                  const Logging::Level&,
///                                                  std::ostream*);
/// }
/// @endcode
///
/// Then you can run your executable, which uses Acts tools and algorithms, in
/// the following way (tested under Unix)
///
/// @code{bash}
/// LD_PRELOAD=<YOUR_SHARED_LIBRARY> path/to/your/exectuable
/// @endcode
///
/// For an example have a look at CustomDefaultLogger.cpp which you can use as follows:
///
/// @code{bash}
/// cd <ACTS/INSTALL/DIRECTORY>
/// source bin/setup.sh
/// LD_PRELOAD=lib/libActsCustomLogger.so bin/Examples/ActsGenericDetector
/// @endcode

/// @defgroup Core Core classes
/// @brief Acts core classes

/// @defgroup Detector Tracking geometry
/// @ingroup Core
/// @brief Description of the tracking geometry

/// @defgroup EventData Event data model
/// @ingroup Core
/// @brief Event data model

/// @defgroup Extrapolation Track extrapolation
/// @ingroup Core
/// @brief Algorithms for extrapolation of track parameters

/// @defgroup Fitter Track fitters
/// @ingroup Core
/// @brief Algorithms for track fitting

/// @defgroup Layers Layers
/// @ingroup Core
/// @brief Description of detector layers

/// @defgroup MagneticField Magnetic field
/// @ingroup Core
/// @brief Description of magnetic field configurations
///
/// This module collects information about classes and typedefs useful for
/// describing different magnetic field configurations.
///
/// Acts is independent of the magnetic field implementation used.
/// Algorithms which need magnetic field information (e.g. Acts::RungeKuttaEngine,
/// Acts::AtlasStepper, Acts::EigenStepper) are templated on the magnetic field.
/// The requirements for the magnetic field implementation are the implementation of the following functions:
/// @todo implement concept for magnetic field implementation
/// @code
/// // retrieve field at given position
/// Acts::Vector3D getField(const Acts::Vector3D& position) const
/// // retrieve magnetic field cell at given position
/// // where Cache is the specific Cache struct defined by the magnetic field implementation
/// Acts::Vector3D getField(const Acts::Vector3D& position, Cache& cache) const
/// // retrieve gradient of magnetic field value
/// Acts::Vector3D getFieldGradient(const Acts::Vector3D& position, Acts::ActsMatrixD<3, 3>& derivative) const
/// // check whether given 3D position is inside this field cell
/// bool isInside(const Acts::Vector3D& position) const
/// @endcode
///
/// Each magnetic field implementation expects to be passed a reference to an 
/// implementation specific cache object. This can usually be achieved through
/// <tt>typename BField::Cache cache</tt>, where @c BField is a template parameter.
///
/// Acts comes with these implementations of this (implicit) interface:
///
/// <B>1. Constant magnetic field implementation</B>
///
/// Should be used to describe a constant magnetic field.
/// The Acts::ConstantBField returns a given constant magnetic field value at every point and can be set by the
/// user either at construction or with a set function.
///
/// <B>2. Interpolated magnetic field  implementation</B>
///
/// For more complex magnetic field implementations the Acts::InterpolatedBFieldMap can be used.
///
/// The Acts::InterpolatedBFieldMap internally uses a field mapper which follows
/// the <a href="http://www.boost.org/doc/libs/1_55_0/doc/html/boost_typeerasure/basic.html">boost concept</a>
/// Acts::concept::AnyFieldLookup. This allows users to provide their own field mapper implementation
/// using the Acts::InterpolatedBFieldMap interface.
///
/// Acts provides a default field mapper implementation: Acts::InterpolatedBFieldMap::FieldMapper, which maps global cartesian 3D positions
/// to magnetic field values. It uses an underlying grid which follows the <a href="http://www.boost.org/doc/libs/1_55_0/doc/html/boost_typeerasure/basic.html">boost concept</a>
/// Acts::concept::AnyNDimGrid which can be a grid of any dimension and allows users to provide their own
/// grid implementation. Furthermore users also need to provide two functions in order to use the Acts::InterpolatedBFieldMap::FieldMapper:
/// 1. a function mapping cartesian global 3D coordinates onto the grid coordinates (of dimension N)
/// 2. a function calculating cartesian global 3D coordinates of the magnetic field with the local N dimensional field and the global 3D position as an input
///
/// A default Acts::detail::Grid implementation is provided following the Acts::concept::AnyNDimGrid,
/// which is flexible (using template parameters) on the dimension of the grid and the value stored in the grid.
///
/// Two convenience functions in order ease the creation of an Acts::InterpolatedBFieldMap::FieldMapper e.g. when reading in a field map from a file,
/// are provided:
/// 1. Acts::InterpolatedBFieldMap::FieldMapper<2, 2> fieldMapperRZ()
/// 2. Acts::InterpolatedBFieldMap::FieldMapper<3, 3> fieldMapperXYZ()
///
/// <B>3. SharedBField</B>
/// 
/// Wraps another @c BField type, which it holds as a @c shared_ptr. The instance can then be copied without having to duplicate
/// the underlying field implementation. This is useful in the case of a large B-Field map.

/// @defgroup Material Material
/// @ingroup Core
/// @brief Description of material properties

/// @defgroup Surfaces Geometric surfaces
/// @ingroup Core
/// @brief Description of geometric surfaces

/// @defgroup Tools Tools
/// @ingroup Core
/// @brief Geometry building tools

/// @defgroup Utilities Helper classes
/// @ingroup Core
/// @brief Helper utilities

/// @defgroup Volumes Volumes
/// @ingroup Core
/// @brief Description of geometric volumes

/// @defgroup Examples Examples
/// @brief Acts Examples


/// @defgroup Plugins Plugins
/// @brief Acts extensions

/// @defgroup DD4hepPlugins DD4hepPlugins
/// @ingroup Plugins
/// @brief Build Acts tracking geometry from \a %DD4hep input.
///
/// The DD4hepPlugin allows building of the Acts TrackingGeometry from <a href="http://aidasoft.web.cern.ch/DD4hep">DD4hep</a> input.
/// \a %DD4hep uses <a href="https://root.cern.ch">ROOT</a> TGeo as an underlying geometry model.
///
/// <B>General</B>
///
/// The basic input for building the detector are a detector description in the
/// XML file format and a corresponding detector constructor written in C++. These
/// two components have to be changed accordingly. Detector constructors use these
/// XML files as an input and construct a detector in the \a %DD4hep geometry format.
/// The whole detector is segmented into different detector parts, e.g. barrels
/// and endcaps which describe different sub-detectors. Since these sub-detectors
/// are built differently, they need different detector constructors. In this way
/// one detector description in XML can use various detector constructors, needed
/// for the different kinds of detector parts.
///
/// In the detector description model of \a %DD4hep any detector is a tree of instances
/// of the so-called \a DetElement class. This \a DetElement class provides all needed
/// detector information, e.g. readout segmentation, geometrical information,
/// environmental conditions. This tree is parallel to the volume tree, which
/// provides the TGeoVolumes and their placements.
/// The relation between these two trees is one-directional, i.e. every volume can be
/// accessed via its corresponding \a DetElement, but not vice versa.
/// Not every supporting material will be declared as a detector element, hence, the
/// geometrical tree can have a deeper hierarchy structure. In order to access both,
/// detector specific and geometrical information, the conversion to the tracking
/// geometry navigates through the detector tree.
/// The \a DetElement can also be extended, to add specific features or to access
/// information. This extension mechanism is used during the translation process.
///
/// <B>ActsExtensions</B>
///
/// \a %DD4hep provides a special extension mechanism for the \a DetElement which allows to
/// add custom features. In Acts this functionality is used for the conversion from
/// \a %DD4hep into Acts.
/// The extensions are used to indicate certain volumes, e.g. if a \a DetElement is the
/// beam pipe or if a \a DetElement is a layer carrying the sensitive modules. In addition
/// the extensions are used in order to distinguish if a sub detector is a barrel (described
/// as a cylinder volume in Acts) or an endcap (which is described as a disc volume in Acts).
/// Furthermore the
/// extensions are used to hand over specific information needed for tracking, e.g.
/// paramters for material mapping.
/// Please find further information in Acts::ActsExtension.
///
/// <B>DD4hepDetectorElement</B>
///
/// In Acts the surfaces describing the sensitive modules of a detector are directly
/// linked to these of the initial geometry input. In the case of \a %DD4hep the
/// Acts::DD4hepDetectorElement was introduced which is the direct link of Acts to \a %DD4hep.
/// In the case for tracking relevant paramters in the \a %DD4hep geometry description
/// are changed (e.g. alignment) it will be automatically changed in Acts.
///
/// <B>Build</B>
///
/// The DD4hepPlugin is only build on demand. The DD4hepPlugin depends on the TGeoPlugin
/// therefore both Plugins need to be installed.
/// During the cmake configuration the flags \c ACTS_BUILD_DD4HEP_PLUGIN and \c ACTS_BUILD_TGEO_PLUGIN need to be set \a ON.
/// In addition \a ROOT and \a %DD4hep need to be added to the \c CMAKE_PREFIX_PATH. When using the DD4hepPlugin
/// both the TGeoPlugin and the DD4hepPlugin components need to be loaded in the user's CMake configuration.
///
/// <B>Prerequisites</B>
///
/// To guarantee a working translation from \a %DD4hep input to Acts geometry the
/// following conditions need to be met:
///
/// * The detector needs to have a barrel-endcap structure:
///   Every hierarchy of subdetectors (e.g. PixelDetector,
///   StripDetector,..)
///   needs to be decomposed of
///   1) {barrel}
///   2) {barrel + 2 endcaps}
///   3) {2 endcaps} - in case there is no barrel at this stage (e.g. forward end caps)
///
/// * If a hierachy is not only a single barrel but is decomposed of a barrel
///   and its corresponding endcaps they need to be grouped together in an assembly
///   using the \c DD4hep_SubdetectorAssembly constructor which is provided by
///   \a %DD4hep.
///   Example of usage in xml file (where \c Barrel0, \c nEndCap0 and \c
///   pEndCap0
///   are sub detectors defined in the file \c PixelTracker.xml):
///   @code
///   <include ref="PixelTracker.xml"/>
///   <detectors>
///     <detector id="1" name="PixelTracker" type="DD4hep_SubdetectorAssembly"
///      vis="BlueVisTrans">
///       <shape name="PixelEnvelope" type="Tube" rmin="Env0_rmin"
///        rmax="Env0_rmax" dz="Env0_dz" material="Air"/>
///         <composite name="Barrel0"/>
///         <composite name="nEndCap0"/>
///         <composite name="pEndCap0"/>
///     </detector>
///   </detectors>
///   @endcode
///
///   If a user wants to create his/her own constructor to group these
///   volumes together the type needs to be set to "compound".
///
///
///	* Since the translation walks trough the \a DetElement tree the following
/// objects
///   need to be declared as a \a %DD4hep \a DetElement:
/// 	- the subvolumes e.g. \b barrel, \b endcap, \b beampipe (they are usually
/// build
/// with
///		  different \a %DD4hep constructors and are therefore \a %DD4hep \a DetElement's
/// per
/// 	  default)
/// 	- \b layers when containing sensitive material and/or the layer should
/// carry
/// 	  material (which will be mapped on the layer if indicated), or the layer is sensitive itself
///       @note the layer does not need to be a direct child of the volume (barrel or endcap), it can be nested in substructures
///
///	    - sensitive detector modules
///		  @note the sensitive detector modules need to be placed in a layer however it can be nested in substructures (can be a component of a modules) i.e. it does not need to be a direct child of the layer
///
///	* The Tracking geometry needs to be built from bottom to top to ensure navigation. Therefore the different hierachies need to be sorted ascending. Per default the sub detectors are sorted by the id of their dd4hep::DetElement. In case another sorting needs to be applied, the users can provide their own function
///
/// * The Acts::ActsExtension's need to be used during the detector construction
///    indicating if a \a DetElement
///		- is a barrel
///		- is an endcap
///		- is the beampipe
///		- is a layer
///
/// There are two modes building the layers around the sensitive detector modules:
/// * The \a DetElements containing the sensitive modules have a geometrical shape
/// 	- the boundaries of the layers in Acts are taken directly from the given shape
/// * The \a DetElements containing the sensitive modules have no specific shape (assembliy)
/// 	- the boundaries of the layers are calculated automatically by adding a tolerance to the geometric extension of the contained surfaces
///		- the tolerances in r and z need to be set for every \a DetElement representing layer using envelopeR and envelopeZ in the Acts::ActsExtension's
///
/// The volumes are automatically build around the layers:
/// 	- the boundaries for the volumes are calculated automatically by adding a tolerance to the geometric extension of the contained layers
///		- the parameters layerEnvelopeR and layerEnvelopeZ (tolerances) need to be set in the Acts::convertDD4hepDetector() function
///
/// Furthermore parameters can be handed over for material mapping or the axes
///   orientation of modules.
///
///
/// Summing up the \a DetElement tree in \a %DD4hep should have the following
/// structure:
/// \image html DD4hepPlugin_DetElementStructure.jpg
///
/// It is also possible to translate a very simple detector geometry, which just
/// consists of cylindrical (for a barrel) or disc (for endcaps) layers which
/// either have material, or, are declared sensitive in dd4hep themselves
/// without containing any detector modules.
///
/// <B>Usage</B>
///
/// To receive the Acts::TrackingGeometry the user should use the global function
/// Acts::convertDD4hepDetector(), where he/she needs to hand over the
/// world \a DetElement of \a %DD4hep.
/// For a valid translation the user needs to make sure, that all prerequisites
/// described above are met and that the right
/// Acts::ActsExtension's are added during the \a %DD4hep construction.


/// @defgroup MaterialPlugins MaterialPlugins
/// @ingroup Plugins
/// @brief Map material onto the Acts geometry.
///
/// The MaterialPlugins allow to map material from a detailed full detector geometry onto the simplfied Acts geometry.
/// The material is mapped onto layers of the tracking geometry which are marked to carry support material. The marking
/// is done during the geometry building process. The material can be mapped onto either, the inner, the outer
/// boundary surface or the central (representing) Acts::Surface of the Acts::Layer. The Acts::Material is described on a two dimensional
/// grid for each layer (Acts::BinnedSurfaceMaterial). The user defines the granularity of the grid during the geometry building process.
/// @note The DD4hepPlugin offers the possiility to mark layers which should carry material and to determine the
/// grid granularity, using the class Acts::ActsExtension.
///
/// Following the Acts philosophy the material mapping is agnostic to any file format and software used to create or store
/// the material maps. The material should be stored in instances of the class Acts::MaterialTrack. This material track record
/// represents a track starting from a certain position, in a certain direction, containing all material along this
/// track. The material along the material track record ist stored as a container of Acts::MaterialStep instances. Each material
/// step contains the material and its thickness at a certain position.
///
/// The material mapping process can be split into two subprocesses:
/// * material assignment
/// * material averaging
///
/// <B>Material Assignment</B>
/// During the material assignment process the decision onto which layer each material step will be assigned is done.
/// To assign a Acts::MaterialTrack the function Acts::MaterialMapping::mapMaterial() should be used.
/// This function extrapolates through the tracking detector, with the start position and direction given by the
/// material track record and collects all layers marked to carry material. Then it loops through all material steps
/// of the material track record and assigns the material of each step to the closest layer :
///
/// \image html MaterialAssignment.jpeg"Example of material assignment onto the inner boundary surface of the layers. The green points are assigned to the current inner layer, the red to the next inner layer."
///
/// <B>Material Averaging</B>
/// During the material mapping the user can decide to average the material whenever he/she prefers by using
/// the function Acts::MaterialMapping::averageLayerMaterial(). In the end when all material track records have been mapped one
/// should use the function Acts::MaterialMapping::finalizeLayerMaterial() in order to finalize the process.
///
///
/// The full material mapping process should be done in the framework of the user.
///
/// Possible workflow:
/// * Create material map(s) of full detector geometry using Acts::MaterialTrack
/// * Read in material map(s) and go through all collected material track records
///		- use Acts::MaterialMapping::mapMaterial() for each material track record
/// * Use Acts::MaterialMapping::averageLayerMaterial() - once per run
/// * In the end of the process use Acts::MaterialMapping::finalizeLayerMaterial() which assigns the final material to the layers

/// @defgroup Contributing Contribution guide

// clang-format on
