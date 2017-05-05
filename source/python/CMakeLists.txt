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
# TBPP Python Wrapper Module

set(SRC
    src/tbppmodule.cpp
    src/context.cpp
    src/node.cpp
    src/model.cpp
    src/lattice.cpp
    src/kpath.cpp
    src/kgrid.cpp
    src/kmodel.cpp
    src/cpa.cpp
    src/gen_cond.cpp
    src/cond_dos.cpp
    src/approx_cond_dos.cpp
    src/berryflux.cpp
)

set(SCRIPTS
    scripts/__init__.py
    scripts/tbpp_plot.py
)

add_library(python_tbpp MODULE ${SRC})
set_target_properties(python_tbpp PROPERTIES PREFIX "")
set_target_properties(python_tbpp PROPERTIES OUTPUT_NAME tbpp)

#-----------------------------------------------------------------------------
# Python

include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${PYTHON_NUMPY_INCLUDE_DIR})
target_link_libraries(python_tbpp ${PYTHON_LIBRARIES})
target_link_libraries(python_tbpp tbpp)

# FIXME make compiler independent or fix code to remove warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings")


#-----------------------------------------------------------------------------
# Install

install(TARGETS python_tbpp DESTINATION "${TBPP_INSTALL_LIBRARY_DIR}/pytbpp")
install(TARGETS python_tbpp DESTINATION "${TBPP_INSTALL_LIBRARY_DIR}/pytbpp")
install(FILES ${SCRIPTS} DESTINATION "${TBPP_INSTALL_LIBRARY_DIR}/pytbpp")
