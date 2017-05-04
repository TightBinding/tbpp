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

# This script sets the following two variables
#      PYTHON_NUMPY_FOUND
#      PYTHON_NUMPY_INCLUDE_DIR

if(NOT PYTHON_EXECUTABLE)
    find_package(PythonInterp)
endif()

if(PYTHON_EXECUTABLE)
    #message(STATUS "Searching for numpy with python executable: ${PYTHON_EXECUTABLE}")

    exec_program("${PYTHON_EXECUTABLE}"
        ARGS "-c 'import numpy; print(numpy.get_include())'"
        OUTPUT_VARIABLE NUMPY_PATH)

    find_path(PYTHON_NUMPY_INCLUDE_DIR numpy/arrayobject.h
        HINTS "${NUMPY_PATH}" "${PYTHON_INCLUDE_PATH}")

    if(PYTHON_NUMPY_INCLUDE_DIR)
        set(PYTHON_NUMPY_FOUND 1 CACHE INTERNAL "Python numpy found")
    endif()

else()
    message(STATUS "Python executable not found.")
endif()

