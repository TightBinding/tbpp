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

#include <python/gen_cond.h>
#include <python/cpa.h>
#include <python/wrap.h>
#include <python/wrap_numpy.h>
#include <memory>

using namespace tbpp;

//----------------------------------------------------------------------

static PyObject *
gen_cond_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    GenCondPy *self  = NULL;
    self = (GenCondPy *)type->tp_alloc(type, 0);
    if(self != NULL) {
        self->ptr = nullptr;
    }

    init_numpy();

    return (PyObject *)self;
}

static int
gen_cond_init(GenCondPy *self, PyObject *args, PyObject *kwds)
{
    self->ptr = std::dynamic_pointer_cast<Node>(std::make_shared<GenCond>());
    return 0;
}

//----------------------------------------------------------------------

static PyObject *
gen_cond_get_zero_j(GenCondPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    GenCondPtr ptr = std::dynamic_pointer_cast<GenCond>(self->ptr);
    return wrap::to_py(ptr->zero_j);
}

static PyObject *
gen_cond_get_solve_sigx(GenCondPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    GenCondPtr ptr = std::dynamic_pointer_cast<GenCond>(self->ptr);
    return wrap::to_py(ptr->solve_sigx);
}

static PyObject *
gen_cond_get_solve_sigy(GenCondPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    GenCondPtr ptr = std::dynamic_pointer_cast<GenCond>(self->ptr);
    return wrap::to_py(ptr->solve_sigy);
}

static PyObject *
gen_cond_get_solve_sigz(GenCondPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    GenCondPtr ptr = std::dynamic_pointer_cast<GenCond>(self->ptr);
    return wrap::to_py(ptr->solve_sigz);
}

static PyObject *
gen_cond_get_sigma(GenCondPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    GenCondPtr ptr = std::dynamic_pointer_cast<GenCond>(self->ptr);
    return cpa_wrap(ptr->sigma());
}

static PyObject *
gen_cond_get_phi(GenCondPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    GenCondPtr ptr = std::dynamic_pointer_cast<GenCond>(self->ptr);
    return wrap::to_py(ptr->phi);
}

static PyObject *
gen_cond_get_cond(GenCondPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    GenCondPtr ptr = std::dynamic_pointer_cast<GenCond>(self->ptr);
    return wrap::to_py(ptr->cond);
}

//----------------------------------------------------------------------

static int
gen_cond_set_zero_j(GenCondPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    GenCondPtr ptr = std::dynamic_pointer_cast<GenCond>(self->ptr);
    return wrap::to_cpp(value, ptr->zero_j);
}

static PyObject *
gen_cond_set_sigma(GenCondPy *self, PyObject *obj) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    GenCondPtr ptr = std::dynamic_pointer_cast<GenCond>(self->ptr);

    if(!CPAPy_Check(obj)) {
        PyErr_SetString(PyExc_RuntimeError, "CPA Object Required");
        return NULL;
    }

    NodePy* node = reinterpret_cast<NodePy*>(obj);
    ptr->set_sigma(std::dynamic_pointer_cast<CPA>(node->ptr));

    Py_RETURN_NONE;
}

static int
gen_cond_set_solve_sigx(GenCondPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    GenCondPtr ptr = std::dynamic_pointer_cast<GenCond>(self->ptr);
    return wrap::to_cpp(value, ptr->solve_sigx);
}

static int
gen_cond_set_solve_sigy(GenCondPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    GenCondPtr ptr = std::dynamic_pointer_cast<GenCond>(self->ptr);
    return wrap::to_cpp(value, ptr->solve_sigy);
}

static int
gen_cond_set_solve_sigz(GenCondPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    GenCondPtr ptr = std::dynamic_pointer_cast<GenCond>(self->ptr);
    return wrap::to_cpp(value, ptr->solve_sigz);
}

//----------------------------------------------------------------------

static PyMethodDef gen_cond_methods[] = {
    {"set_sigma", (PyCFunction)gen_cond_set_sigma, METH_O,
        "Set the self-energy"},
    {NULL}
};

static PyGetSetDef gen_cond_getseters[] = {
    {"sigma", (getter)gen_cond_get_sigma, NULL,
        "Get the self-energy", NULL},
    {"phi", (getter)gen_cond_get_phi, NULL,
        "Get the Phi matrix", NULL},
    {"cond", (getter)gen_cond_get_cond, NULL,
        "Get the computed conductivity", NULL},
    {"solve_sigx", (getter)gen_cond_get_solve_sigx,
        (setter)gen_cond_set_solve_sigx,
        "Solve conductivity along x-axis", NULL},
    {"solve_sigy", (getter)gen_cond_get_solve_sigy,
        (setter)gen_cond_set_solve_sigy,
        "Solve conductivity along y-axis", NULL},
    {"solve_sigz", (getter)gen_cond_get_solve_sigz,
        (setter)gen_cond_set_solve_sigz,
        "Solve conductivity along z-axis", NULL},
    {"zero_j", (getter)gen_cond_get_zero_j,
        (setter)gen_cond_set_zero_j,
        "Minimum value for +0j when computing Greens matrix", NULL},
    {NULL}
};

//----------------------------------------------------------------------

PyTypeObject GenCondPyType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "tbpp.GenCond",             /* tp_name */
    sizeof(GenCondPy),          /* tp_basicsize */
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
    "TBPP GenCond object",      /* tp_doc */
    0,                          /* tp_traverse */
    0,                          /* tp_clear */
    0,                          /* tp_richcompare */
    0,                          /* tp_weaklistoffset */
    0,                          /* tp_iter */
    0,                          /* tp_iternext */
    gen_cond_methods,           /* tp_methods */
                 0,             /* tp_members */
    gen_cond_getseters,         /* tp_getset */
    0,                          /* tp_base */
    0,                          /* tp_dict */
    0,                          /* tp_descr_get */
    0,                          /* tp_descr_set */
    0,                          /* tp_dictoffset */
    (initproc)gen_cond_init,    /* tp_init */
    0,                          /* tp_alloc */
    gen_cond_new,               /* tp_new */
};

//----------------------------------------------------------------------

PyObject *
gen_cond_wrap(GenCondPtr ptr) {
    PyObject *ret = GenCondPyType.tp_new((PyTypeObject*)&GenCondPyType, NULL, NULL);
    NodePy* obj = reinterpret_cast<NodePy*>(ret);
    obj->ptr = std::dynamic_pointer_cast<Node>(ptr);
    return ret;
}

GenCondPtr
gen_cond_unwrap(PyObject* obj) {
    if(!GenCondPy_Check(obj)) return nullptr;
    NodePy* m = reinterpret_cast<NodePy*>(obj);
    return std::dynamic_pointer_cast<GenCond>(m->ptr);
}

