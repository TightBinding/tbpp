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

ExternalProject_Add(EIGEN
    PREFIX extern/eigen
    URL "http://bitbucket.org/eigen/eigen/get/3.3.3.zip"
    URL_MD5 e99bd8f4363531e23a6cd4f4f5f87d0b 
    CMAKE_CACHE_ARGS
        -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
        -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
        -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
        -DINCLUDE_INSTALL_DIR:PATH=include
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    )

#set(EIGEN3_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include/eigen3 )
set(EIGEN3_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include )
