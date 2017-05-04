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

ExternalProject_Add(VTK
    # DEPENDS HDF5
    PREFIX extern/vtk
    URL http://www.vtk.org/files/release/7.1/VTK-7.1.1.zip
    URL_MD5 ee9f921fc5bdc0d2ff8733a28c3e30e6
    CMAKE_CACHE_ARGS
        -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
        -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
        -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX:PATH=${TBPP_EXTERN_INSTALL_PATH}
        -DBUILD_TESTING:BOOL=OFF
        -DBUILD_EXAMPLES:BOOL=OFF
        -DBUILD_SHARED_LIBS:BOOL=OFF
        -DVTK_BUILD_SHARED_LIBS:BOOL=OFF
        # -DVTK_INSTALL_NO_QT_PLUGIN:BOOL=ON
        # -DVTK_INSTALL_NO_PYTHON_EXES:BOOL=ON
        -DVTK_Group_StandAlone:BOOL=OFF 
        -DVTK_Group_Rendering:BOOL=ON 
        -DVTK_RENDERING_BACKEND:STRING=OpenGL2 
        -DVTK_Group_Qt:BOOL=OFF 
        -DVTK_QT_VERSION:STRING=5
        -DVTK_BUILD_QT_DESIGNER_PLUGIN:BOOL=OFF
        -DModule_vtkGUISupportQt:BOOL=ON
        -DModule_vtkRenderingQt:BOOL=ON
        -DModule_vtkViewsQt:BOOL=ON
        # -DHDF5_DIR:PATH=${HDF5_DIR}
        -DQt5_DIR:PATH=${Qt5_DIR}
    INSTALL_DIR ${TBPP_EXTERN_INSTALL_PATH}
)

set(VTK_DIR ${TBPP_EXTERN_INSTALL_PATH}/lib/cmake/vtk-7.1)
