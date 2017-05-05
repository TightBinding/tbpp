/* ===========================================================================
 * Copyright (c) 2016-2017 Giacomo Resta
 *
 * This file is part of TightBinding++.
 *
 * TightBinding++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * TightBinding++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ===========================================================================
 */

#include "python/tbppmodule.h"
#include "tbpp/common.h"
#include <sstream>

#include "python/wrap_numpy.h"
#include "python/context.h"
#include "python/node.h"
#include "python/lattice.h"
#include "python/model.h"
#include "python/kpath.h"
#include "python/kgrid.h"
#include "python/kmodel.h"
#include "python/cpa.h"
#include "python/cond_dos.h"
#include "python/gen_cond.h"
#include "python/approx_cond_dos.h"
#include "python/berryflux.h"

//----------------------------------------------------------------------------

static PyMethodDef TBPPMethods[] = {
    {NULL}
};

PyDoc_STRVAR(tbpp_doc,
"Python Wrapper for TightBinding++");
static struct PyModuleDef tbppmodule = {
    PyModuleDef_HEAD_INIT,
    "tbpp",    /* name of module */
    tbpp_doc,  /* module documentation, may be NULL */
    -1,        /* size of per-interpreter state of the module or -1 if
                  module keeps state in global variables. */
    TBPPMethods,
    NULL, NULL, NULL, NULL
};

//----------------------------------------------------------------------------

PyMODINIT_FUNC
PyInit_tbpp(void)
{
    PyObject *m;

    //------------------------------------------------------------------------

    // Context
    if(PyType_Ready(&ContextPyType) < 0) return NULL;

    // Node
    if(PyType_Ready(&NodePyType) < 0) return NULL;

    // Model
    ModelPyType.tp_base = &NodePyType;
    if(PyType_Ready(&ModelPyType) < 0) return NULL;

    // Lattice
    LatticePyType.tp_base = &NodePyType;
    if(PyType_Ready(&LatticePyType) < 0) return NULL;

    // KPath
    KPathPyType.tp_base = &NodePyType;
    if(PyType_Ready(&KPathPyType) < 0) return NULL;

    // KGrid
    KGridPyType.tp_base = &NodePyType;
    if(PyType_Ready(&KGridPyType) < 0) return NULL;

    // KModel
    KModelPyType.tp_base = &NodePyType;
    if(PyType_Ready(&KModelPyType) < 0) return NULL;

    // CPA
    CPAPyType.tp_base = &NodePyType;
    if(PyType_Ready(&CPAPyType) < 0) return NULL;

    // CondDOS
    CondDOSPyType.tp_base = &NodePyType;
    if(PyType_Ready(&CondDOSPyType) < 0) return NULL;

    // GenCond
    GenCondPyType.tp_base = &NodePyType;
    if(PyType_Ready(&GenCondPyType) < 0) return NULL;

    // ApproxCondDOS
    ApproxCondDOSPyType.tp_base = &NodePyType;
    if(PyType_Ready(&ApproxCondDOSPyType) < 0) return NULL;

    // BerryFlux
    BerryFluxPyType.tp_base = &NodePyType;
    if(PyType_Ready(&BerryFluxPyType) < 0) return NULL;

    //------------------------------------------------------------------------

    // Create Module
    m = PyModule_Create(&tbppmodule);
    if (m == NULL) return NULL;

    init_numpy();

    //------------------------------------------------------------------------

    // Version constants
    PyModule_AddIntConstant(m, "__version_major__", TBPP_VERSION_MAJOR);
    PyModule_AddIntConstant(m, "__version_minor__", TBPP_VERSION_MINOR);

    std::stringstream stream;
    stream << TBPP_VERSION_MAJOR << '.' << TBPP_VERSION_MINOR;
    PyModule_AddStringConstant(m, "__version__", stream.str().c_str());

    //------------------------------------------------------------------------

    // Context
    Py_INCREF(&ContextPyType);
    PyModule_AddObject(m, "Context", (PyObject *)&ContextPyType);

    // Node
    Py_INCREF(&NodePyType);
    PyModule_AddObject(m, "Node", (PyObject *)&NodePyType);

    // Model
    Py_INCREF(&ModelPyType);
    PyModule_AddObject(m, "Model", (PyObject *)&ModelPyType);

    // Lattice
    Py_INCREF(&LatticePyType);
    PyModule_AddObject(m, "Lattice", (PyObject *)&LatticePyType);

    // KPath
    Py_INCREF(&KPathPyType);
    PyModule_AddObject(m, "KPath", (PyObject *)&KPathPyType);

    // KGrid
    Py_INCREF(&KGridPyType);
    PyModule_AddObject(m, "KGrid", (PyObject *)&KGridPyType);

    // KModel
    Py_INCREF(&KModelPyType);
    PyModule_AddObject(m, "KModel", (PyObject *)&KModelPyType);

    // CPA
    Py_INCREF(&CPAPyType);
    PyModule_AddObject(m, "CPA", (PyObject *)&CPAPyType);

    // CondDOS
    Py_INCREF(&CondDOSPyType);
    PyModule_AddObject(m, "CondDOS", (PyObject *)&CondDOSPyType);

    // GenCond
    Py_INCREF(&GenCondPyType);
    PyModule_AddObject(m, "GenCond", (PyObject *)&GenCondPyType);

    // ApproxCondDOS
    Py_INCREF(&ApproxCondDOSPyType);
    PyModule_AddObject(m, "ApproxCondDOS", (PyObject *)&ApproxCondDOSPyType);

    // BerryFlux
    Py_INCREF(&BerryFluxPyType);
    PyModule_AddObject(m, "BerryFlux", (PyObject *)&BerryFluxPyType);

    return m;
}
