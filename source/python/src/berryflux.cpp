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

#include <python/berryflux.h>
#include <python/kgrid.h>
#include <python/wrap.h>
#include <python/wrap_numpy.h>
#include <memory>

using namespace tbpp;

//----------------------------------------------------------------------

static PyObject *
berryflux_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    BerryFluxPy *self  = NULL;
    self = (BerryFluxPy *)type->tp_alloc(type, 0);
    if(self != NULL) {
        self->ptr = nullptr;
    }

    init_numpy();

    return (PyObject *)self;
}

static int
berryflux_init(BerryFluxPy *self, PyObject *args, PyObject *kwds)
{
    self->ptr = std::dynamic_pointer_cast<Node>(std::make_shared<BerryFlux>());
    return 0;
}

//----------------------------------------------------------------------

static PyObject *
berryflux_set_kgrid(BerryFluxPy *self, PyObject *obj) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    BerryFluxPtr ptr = std::dynamic_pointer_cast<BerryFlux>(self->ptr);

    if(!KGridPy_Check(obj)) {
        PyErr_SetString(PyExc_RuntimeError, "KGrid Object Required");
        return NULL;
    }

    NodePy* node = reinterpret_cast<NodePy*>(obj);
    ptr->set_kgrid(std::dynamic_pointer_cast<KGrid>(node->ptr));

    Py_RETURN_NONE;
}

//----------------------------------------------------------------------
//
PyDoc_STRVAR(berryflux_get_flux_k_doc,
"Berry flux for each energy level");
static PyObject *
berryflux_get_flux_k(BerryFluxPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    BerryFluxPtr ptr = std::dynamic_pointer_cast<BerryFlux>(self->ptr);
    return wrap::to_py(ptr->flux_k);
}


PyDoc_STRVAR(berryflux_get_flux_doc,
"Integral of the Berry flux for each energy level");
static PyObject *
berryflux_get_flux(BerryFluxPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    BerryFluxPtr ptr = std::dynamic_pointer_cast<BerryFlux>(self->ptr);
    return wrap::to_py(ptr->flux);
}

static PyObject *
berryflux_get_kgrid(BerryFluxPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    BerryFluxPtr ptr = std::dynamic_pointer_cast<BerryFlux>(self->ptr);
    return kgrid_wrap(ptr->kgrid());
}

static int
berryflux_set_solve_flux_k(BerryFluxPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    BerryFluxPtr ptr = std::dynamic_pointer_cast<BerryFlux>(self->ptr);
    return wrap::to_cpp(value, ptr->solve_flux_k);
}

static PyObject *
berryflux_get_solve_flux_k(BerryFluxPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    BerryFluxPtr ptr = std::dynamic_pointer_cast<BerryFlux>(self->ptr);
    return wrap::to_py(ptr->solve_flux_k);
}


//----------------------------------------------------------------------

static PyMethodDef berryflux_methods[] = {
    {"set_kgrid", (PyCFunction)berryflux_set_kgrid, METH_O, "Set the kgrid"},
    {NULL}
};

static PyGetSetDef berryflux_getseters[] = {
    {"flux_k", (getter)berryflux_get_flux_k, NULL,
        berryflux_get_flux_k_doc, NULL},
    {"flux", (getter)berryflux_get_flux, NULL,
        berryflux_get_flux_doc, NULL},
    {"kgrid", (getter)berryflux_get_kgrid, NULL,
        "Get the KGrid being used", NULL},
    {"solve_flux_k", (getter)berryflux_get_solve_flux_k,
        (setter)berryflux_set_solve_flux_k,
        "Whether to solve for flux_k", NULL},
    {NULL}
};

//----------------------------------------------------------------------

PyTypeObject BerryFluxPyType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "tbpp.BerryFlux",               /* tp_name */
    sizeof(BerryFluxPy),            /* tp_basicsize */
    0,                          /* tp_itemsize */
    0,                          /* tp_dealloc */
    0,                          /* tp_print */
    0,                          /* tp_getattr */
    0,                          /* tp_setattr */
    0,                          /* tp_reserved */
    0,                          /* tp_repr */
    0,                          /* tp_as_number */
    0,                          /* tp_as_sequence */
    0,                          /* tp_as_mapping */
    0,                          /* tp_hash  */
    0,                          /* tp_call */
    0,                          /* tp_str */
    0,                          /* tp_getattro */
    0,                          /* tp_setattro */
    0,                          /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,    /* tp_flags */
    "TBPP BerryFlux object",        /* tp_doc */
    0,                          /* tp_traverse */
    0,                          /* tp_clear */
    0,                          /* tp_richcompare */
    0,                          /* tp_weaklistoffset */
    0,                          /* tp_iter */
    0,                          /* tp_iternext */
    berryflux_methods,              /* tp_methods */
                 0,             /* tp_members */
    berryflux_getseters,            /* tp_getset */
    0,                          /* tp_base */
    0,                          /* tp_dict */
    0,                          /* tp_descr_get */
    0,                          /* tp_descr_set */
    0,                          /* tp_dictoffset */
    (initproc)berryflux_init,       /* tp_init */
    0,                          /* tp_alloc */
    berryflux_new,                  /* tp_new */
};

//----------------------------------------------------------------------

PyObject *
berryflux_wrap(BerryFluxPtr ptr) {
    PyObject *ret = BerryFluxPyType.tp_new((PyTypeObject*)&BerryFluxPyType, NULL, NULL);
    NodePy* obj = reinterpret_cast<NodePy*>(ret);
    obj->ptr = std::dynamic_pointer_cast<Node>(ptr);
    return ret;
}

BerryFluxPtr
berryflux_unwrap(PyObject* obj) {
    if(!BerryFluxPy_Check(obj))
        return nullptr;
    NodePy* m = reinterpret_cast<NodePy*>(obj);
    return std::dynamic_pointer_cast<BerryFlux>(m->ptr);
}

