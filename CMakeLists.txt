# ============================================================================
# Copyright (c) 2016-2017 Giacomo Resta
#
# This file is part of TightBinding++.
#
# TightBinding++ is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# TightBinding++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ============================================================================

cmake_minimum_required(VERSION 3.2)

#-----------------------------------------------------------------------------
# Require out-of-source build

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR
        "In-source builds are not allowed."
        "\n Remove CMakeCache.txt and run cmake from another directory."
        "\n e.g.:"
        "\n "
        "\n   rm CMakeCache.txt"
        "\n   mkdir build"
        "\n   cd build"
        "\n   cmake ../"
        "\n "
    )
endif()

#-----------------------------------------------------------------------------
# General Items

project(TBPP)
set(TBPP_VERSION_MAJOR 0)
set(TBPP_VERSION_MINOR 1)

# Set build type to release if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Build type empty, setting to 'Release'")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Debug")
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/build_files/cmake_modules")

include(ExternalProject)
set(TBPP_EXTERN_INSTALL_PATH "${CMAKE_BINARY_DIR}/extern_install" CACHE INTERNAL "" FORCE)

#-----------------------------------------------------------------------------
# Output directories

set(TBPP_INSTALL_BIN_DIR bin CACHE INTERNAL "" FORCE)
set(TBPP_INSTALL_LIBRARY_DIR lib CACHE INTERNAL "" FORCE)
set(TBPP_INSTALL_INCLUDE_DIR include/tbpp CACHE INTERNAL "" FORCE)
set(TBPP_INSTALL_DOC_DIR share/doc/tbpp CACHE INTERNAL "" FORCE)

#-----------------------------------------------------------------------------
# Options

option(TBPP_WITH_DOC "Build the Doxygen API Documentation" OFF)
option(TBPP_WITH_PYTHON "Build the python module" OFF)
option(TBPP_WITH_EDITOR "Build the TBPP editor" OFF)

option(TBPP_WITH_OPENMP "Enable OpenMP multithreading (requires compiler support)" OFF)
option(TBPP_WITH_HDF5 "Enable HDF5 file support" ON)

option(TBPP_WITH_SYSTEM_EIGEN "Build using the system Eigen3" ON)
option(TBPP_WITH_SYSTEM_HDF5 "Build using the system HDF5" ON)
option(TBPP_WITH_SYSTEM_VTK "Build using the system VTK" ON)

option(TBPP_NARRAY_RANGE_CHECK "Additional range checking for NArray" OFF)
mark_as_advanced(TBPP_NARRAY_RANGE_CHECK)

#-----------------------------------------------------------------------------
# Require C++11

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#-----------------------------------------------------------------------------
# Build External Libraries

if((NOT TBPP_WITH_SYSTEM_EIGEN)
    OR (TBPP_WITH_HDF5   AND (NOT TBPP_WITH_SYSTEM_HDF5))
    OR (TBPP_WITH_EDITOR AND (NOT TBPP_WITH_SYSTEM_VTK)))
    message(STATUS "================ EXTERNAL PROJECT BUILD ================")
    include(tbpp.cmake)
    return()
endif()

#-----------------------------------------------------------------------------
# EIGEN

# FIXME
#find_package(Eigen3 REQUIRED)
message(STATUS "Using Eigen3 include dir: ${EIGEN3_INCLUDE_DIR}")
include_directories(${EIGEN3_INCLUDE_DIR})

#-----------------------------------------------------------------------------
# HDF5

if(TBPP_WITH_HDF5)
    find_package(HDF5 NAMES hdf5 COMPONENTS C CXX NO_MODULE REQUIRED static)
    message(STATUS "Using HDF5 include dir: ${HDF5_INCLUDE_DIR}")
    message(STATUS "Using HDF5 libraries: ${HDF5_C_STATIC_LIBRARY} ${HDF5_CXX_STATIC_LIBRARY}")
    include_directories(${HDF5_INCLUDE_DIR})
endif()

#-----------------------------------------------------------------------------
# OpenMP

if(TBPP_WITH_OPENMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#-----------------------------------------------------------------------------
# Python and NumPy

if(TBPP_WITH_PYTHON)
    find_package(PythonLibs 3 REQUIRED)
    message(STATUS "Using python include dir: ${PYTHON_INCLUDE_DIRS}")
    message(STATUS "Using python libraries dir: ${PYTHON_LIBRARIES}")

    find_package(PythonInterp 3 REQUIRED)
    find_package(NumPy REQUIRED)
    message(STATUS "Using NumPy include dir: ${PYTHON_NUMPY_INCLUDE_DIR}")
endif()

#-----------------------------------------------------------------------------
# Qt5 and QCustomPlot

if(TBPP_WITH_EDITOR)
    # Do not run moc automatically
    set(CMAKE_AUTOMOC OFF)
    # Do not run uic automatically
    set(CMAKE_AUTOUIC OFF)
    # Find includes in build directories
    set(CMAKE_INCLUDE_CURRENT_DIR on)
    # Find Qt5
    find_package(Qt5 COMPONENTS Core Widgets PrintSupport REQUIRED NO_MODULE)
    # QCustomPlot Library
    add_subdirectory(extern/qcustomplot)
    # Disable RPATH on Mac so that macdeployqt works
    set(CMAKE_MACOSX_RPATH OFF)

    message(STATUS "Using Qt5 include dir: ${Qt5Widgets_INCLUDE_DIRS}")
    message(STATUS "Using Qt5 libraries: ${Qt5Widgets_LIBRARIES}")
endif()

#-----------------------------------------------------------------------------
# VTK

if(TBPP_WITH_EDITOR)
    find_package(VTK 7 REQUIRED NO_MODULE)
    include(${VTK_USE_FILE})
    message(STATUS "Using VTK use file: ${VTK_USE_FILE}")
    message(STATUS "Using VTK include dir: ${VTK_INCLUDE_DIRS}")
    message(STATUS "Using VTK libraries: ${VTK_LIBRARIES}")
endif()

#-----------------------------------------------------------------------------
# Common Include Directories

include_directories("${PROJECT_SOURCE_DIR}/extern")
include_directories("${PROJECT_BINARY_DIR}/extern")
include_directories("${PROJECT_SOURCE_DIR}/source")
include_directories("${PROJECT_BINARY_DIR}/source")

#-----------------------------------------------------------------------------
# Subdirectories

add_subdirectory(source/tbpp)
add_subdirectory(source/tbrun)

if(TBPP_WITH_DOC)
    add_subdirectory(doc/doxygen)
endif()

if(TBPP_WITH_PYTHON)
    add_subdirectory(source/python)
endif()

if(TBPP_WITH_EDITOR)
    add_subdirectory(source/editor)
endif()

#-----------------------------------------------------------------------------
# Packaging

include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${TBPP_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${TBPP_VERSION_MINOR}")
include(CPack)

#-----------------------------------------------------------------------------
# Print all CMake Variables

# get_cmake_property(_var_names VARIABLES)
# foreach (_var_name ${_var_names})
#         message(STATUS "${_var_name}=${${_var_name}}")
# endforeach()

