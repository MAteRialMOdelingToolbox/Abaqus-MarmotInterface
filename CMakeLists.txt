#/* ---------------------------------------------------------------------
# *                                       _   
# *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_ 
# * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
# * | | | | | | (_| | |  | | | | | | (_) | |_ 
# * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
# * 
# * Unit of Strength of Materials and Structural Analysis
# * University of Innsbruck, 
# * 2020 - today
# * 
# * festigkeitslehre@uibk.ac.at
# * 
# * Matthias Neuner matthias.neuner@uibk.ac.at
# * 
# * This file is part of the MAteRialMOdellingToolbox (marmot).
# *
# * This library is free software; you can redistribute it and/or
# * modify it under the terms of the GNU Lesser General Public
# * License as published by the Free Software Foundation; either
# * version 2.1 of the License, or (at your option) any later version.
# * 
# * The full text of the license can be found in the file LICENSE.md at
# * the top level directory of marmot.
# * ---------------------------------------------------------------------
# */


cmake_minimum_required(VERSION 3.12)

project(Abaqus-MarmotInterfaceTest CXX)

# avoid in-source compilation
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed. Create and/or change to your build folder \
  and run then cmake from there, e.g. mkdir build && cd build && cmake ../ . \\n
  You first may need to remove CMakeCache.txt and CMakeFiles/*.")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_FIND_LIBRARY_PREFIXES lib)

# enable testing
enable_testing()

# set install directories
include(GNUInstallDirs)

## find Eigen library
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})
message("--> found Eigen: ${EIGEN3_INCLUDE_DIR}")

# set C++ PREPROCESSOR DEFINITIONS
add_compile_definitions(NO_ABAQUS)


add_library(${PROJECT_NAME} SHARED user_explicit.cpp)

