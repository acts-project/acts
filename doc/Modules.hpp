// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_MODULES_HPP
#define ACTS_MODULES_HPP 1
// clang-format off

/// @defgroup Design Design and concept descriptions
/// @brief description of general concepts used in ACTS

/// @defgroup Logging Debug output options
/// @ingroup Design
/// @brief description of debug output options
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
/// For your convenience, the macro #ACTS_LOCAL_LOGGER is provided which does
/// the job for you. The ACTS logging facility supports several severity levels
/// which
/// allow you to control the amount of information displayed at run-time. Logger
/// objects can easily be created using the Acts::getDefaultLogger function
/// which should be sufficient to get you started. In case you need more
/// customized debug output, you can make use of the output decorators defined
/// in Acts::Logging or even write your own implementation of
/// Acts::Logging::OutputDecorator.
///
/// @par Code example illustrating the usage
/// @code{.cpp}
/// #include <fstream>
/// #include <memory>
/// #include "ACTS/Utilities/Logger.hpp"
///
/// void myFunction() {
///    // open the logfile
///    std::ofstream logfile("log.txt");
///
///    // setup a logger instance for >= INFO messages, streaming into the log file
///    // make sure you do NOT call the variable 'logger'
///    std::unique_ptr<Acts::Logger> myLogger = Acts::getDefaultLogger("MyLogger", Acts::Logging::INFO, &logfile);
///
///    // make sure the ACTS debug macros can work with your logger
///    ACTS_LOCAL_LOGGER(myLogger);
///
///    ACTS_VERBOSE("This message will not appear in the logfile.");
///    ACTS_INFO("But this one will: Hello World!");
///
///    // do not forget to close the logfile
///    logfile.close();
/// }
/// @endcode

/// @defgroup Core Core classes
/// @brief ACTS core classes

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
/// @brief Description of magnetic field properties

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
/// @brief ACTS Examples

/// @defgroup Plugins Plugins
/// @brief ACTS extensions

/// @defgroup Contributing Contribution guide

// clang-format on
#endif  // ACTS_MODULES_HPP
