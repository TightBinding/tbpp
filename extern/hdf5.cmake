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

ExternalProject_Add(HDF5
    PREFIX extern/hdf5
    URL http://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.zip
    URL_MD5 0f57dadf13f45edde38d777160e5a7d6
    CMAKE_CACHE_ARGS
        -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
        -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
        -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
        -DBUILD_SHARED_LIBS:BOOL=OFF
        -DHDF5_BUILD_CPP_LIB:BOOL=ON
        -DBUILD_TESTING:BOOL=OFF
        -DHDF5_BUILD_EXAMPLES:BOOL=OFF
        -DHDF5_BUILD_TOOLS:BOOL=OFF
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
)

set(HDF5_DIR ${CMAKE_INSTALL_PREFIX}/share/cmake)
