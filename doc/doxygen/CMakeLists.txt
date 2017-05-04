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
# Doxygen documentation

find_package(Doxygen)
if(DOXYGEN_FOUND)
    # C++ API Documentation
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile_api.in
        ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile_api
        @ONLY
    )
    add_custom_target(doc_api ALL
        ${DOXYGEN_EXECUTABLE}
        ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile_api
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT "Generating Doxygen API documentation" VERBATIM
    )
else()
    message(FATAL_ERROR "Doxygen not found")
endif()
