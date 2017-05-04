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

#-----------------------------------------------------------------------------
# Configure Header

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/common.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/common.h"
)

#-----------------------------------------------------------------------------
# TBPP Core Source Code

set(SRC
    src/context.cpp
    src/node.cpp
    src/lattice.cpp
    src/model.cpp
    src/kpath.cpp
    src/kgrid.cpp
    src/kmodel.cpp
    src/cpa.cpp
    src/cond_dos.cpp
    src/approx_cond_dos.cpp
    src/berryflux.cpp
    src/gen_cond.cpp
)

set(HDR
    context.h
    node.h
    lattice.h
    model.h
    kpath.h
    kgrid.h
    kmodel.h
    cpa.h
    cond_dos.h
    approx_cond_dos.h
    berryflux.h
    gen_cond.h
    matrix.h
    narray.h
)

#-----------------------------------------------------------------------------
# TBPP HDF5 Code

if(TBPP_WITH_HDF5)
    set(SRC ${SRC}
        src/ehfile.cpp
    )
    set(HDR ${HDR}
        ehfile.h
    )
endif()

#-----------------------------------------------------------------------------
# TBPP Library

add_library(tbpp ${SRC})

if(TBPP_WITH_PYTHON)
    set_property(TARGET tbpp PROPERTY POSITION_INDEPENDENT_CODE On)
endif()

#-----------------------------------------------------------------------------
# HDF5

if(TBPP_WITH_HDF5)
    target_link_libraries(tbpp ${HDF5_C_STATIC_LIBRARY} ${HDF5_CXX_STATIC_LIBRARY})
endif()

#-----------------------------------------------------------------------------
# Install

install(TARGETS tbpp DESTINATION ${TBPP_INSTALL_LIBRARY_DIR})
install(FILES ${HDR} DESTINATION ${TBPP_INSTALL_INCLUDE_DIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/common.h DESTINATION ${TBPP_INSTALL_INCLUDE_DIR})
