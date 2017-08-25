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
               
set(TBPP_EXTERN_DEPENDENCIES )

# Whether to download EIGEN
if(NOT TBPP_WITH_SYSTEM_EIGEN)
    list(APPEND TBPP_EXTERN_DEPENDENCIES EIGEN)
    include(extern/eigen.cmake)
endif()

# Whether to build HDF5
if(TBPP_WITH_HDF5 AND (NOT TBPP_WITH_SYSTEM_HDF5))
    list(APPEND TBPP_EXTERN_DEPENDENCIES HDF5)
    include(extern/hdf5.cmake)
endif()

# Whether to build VTK
if(TBPP_WITH_EDITOR AND (NOT TBPP_WITH_SYSTEM_VTK))
    list(APPEND TBPP_EXTERN_DEPENDENCIES VTK)
    include(extern/vtk.cmake)
endif()

ExternalProject_Add(TBPP
    DEPENDS ${TBPP_EXTERN_DEPENDENCIES}
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    BINARY_DIR tbpp-build
    CMAKE_CACHE_ARGS
        -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
        -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
        -DTBPP_WITH_DOC:BOOL=${TBPP_WITH_DOC}
        -DTBPP_WITH_PYTHON:BOOL=${TBPP_WITH_PYTHON}
        -DTBPP_WITH_EDITOR:BOOL=${TBPP_WITH_EDITOR}
        -DTBPP_WITH_OPENMP:BOOL=${TBPP_WITH_OPENMP}
        -DTBPP_WITH_HDF5:BOOL=${TBPP_WITH_HDF5}
        -DTBPP_WITH_SYSTEM_EIGEN:BOOL=ON
        -DTBPP_WITH_SYSTEM_HDF5:BOOL=ON
        -DTBPP_WITH_SYSTEM_VTK:BOOL=ON
        -DTBPP_NARRAY_RANGE_CHECK:BOOL=${TBPP_NARRAY_RANGE_CHECK}
        -DEIGEN3_INCLUDE_DIR:PATH=${EIGEN3_INCLUDE_DIR}
        -DVTK_DIR:PATH=${VTK_DIR}
        -DHDF5_DIR:PATH=${HDF5_DIR}
        -DQt5_DIR:PATH=${Qt5_DIR}
        -DPYTHON_INCLUDE_DIR:PATH=${PYTHON_INCLUDE_DIR}
        -DPYTHON_LIBRARY:PATH=${PYTHON_LIBRARY}
        -DTBPP_INSTALL_BIN_DIR:PATH=${TBPP_INSTALL_BIN_DIR}
        -DTBPP_INSTALL_LIBRARY_DIR:PATH=${TBPP_INSTALL_LIBRARY_DIR}
        -DTBPP_INSTALL_PYTHON_LIBRARY_DIR:PATH=${TBPP_INSTALL_PYTHON_LIBRARY_DIR}
        -DTBPP_INSTALL_INCLUDE_DIR:PATH=${TBPP_INSTALL_INCLUDE_DIR}
        -DTBPP_INSTALL_DOC_DIR:PATH=${TBPP_INSTALL_DOC_DIR}
)

# Force TBPP to recompile if source code is changed
ExternalProject_Add_Step(TBPP forceconfigure
    COMMAND ${CMAKE_COMMAND} -E echo "Force configure TBPP"
    DEPENDEES update
    DEPENDERS configure
    ALWAYS 1
)
