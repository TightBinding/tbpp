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

#include <python/approx_cond_dos.h>
#include <python/kgrid.h>
#include <python/wrap.h>
#include <python/wrap_numpy.h>
#include <memory>

using namespace tbpp;

//----------------------------------------------------------------------

static PyObject *
approx_cond_dos_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    ApproxCondDOSPy *self  = NULL;
    self = (ApproxCondDOSPy *)type->tp_alloc(type, 0);
    if(self != NULL) {
        self->ptr = nullptr;
    }

    init_numpy();

    return (PyObject *)self;
}

static int
approx_cond_dos_init(ApproxCondDOSPy *self, PyObject *args, PyObject *kwds)
{
    self->ptr = std::dynamic_pointer_cast<Node>(std::make_shared<ApproxCondDOS>());
    return 0;
}

//----------------------------------------------------------------------

static PyObject *
approx_cond_dos_set_kgrid(ApproxCondDOSPy *self, PyObject *obj) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ApproxCondDOSPtr ptr = std::dynamic_pointer_cast<ApproxCondDOS>(self->ptr);

    if(!KGridPy_Check(obj)) {
        PyErr_SetString(PyExc_RuntimeError, "KGrid Object Required");
        return NULL;
    }

    NodePy* node = reinterpret_cast<NodePy*>(obj);
    ptr->set_kgrid(std::dynamic_pointer_cast<KGrid>(node->ptr));

    Py_RETURN_NONE;
}

static PyObject *
approx_cond_dos_set_w(ApproxCondDOSPy *self, PyObject *obj) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ApproxCondDOSPtr ptr = std::dynamic_pointer_cast<ApproxCondDOS>(self->ptr);

    if(wrap::to_cpp(obj, ptr->w)) {
        return NULL;
    }

    Py_RETURN_NONE;
}

static PyObject *
approx_cond_dos_set_eta(ApproxCondDOSPy *self, PyObject *obj) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ApproxCondDOSPtr ptr = std::dynamic_pointer_cast<ApproxCondDOS>(self->ptr);

    if(wrap::to_cpp(obj, ptr->eta)) {
        return NULL;
    }

    Py_RETURN_NONE;
}

//----------------------------------------------------------------------

static PyObject *
approx_cond_dos_get_cond(ApproxCondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ApproxCondDOSPtr ptr = std::dynamic_pointer_cast<ApproxCondDOS>(self->ptr);
    return wrap::to_py(ptr->cond);
}

static PyObject *
approx_cond_dos_get_dos_k(ApproxCondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ApproxCondDOSPtr ptr = std::dynamic_pointer_cast<ApproxCondDOS>(self->ptr);
    return wrap::to_py(ptr->dos_k);
}

static PyObject *
approx_cond_dos_get_dos(ApproxCondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ApproxCondDOSPtr ptr = std::dynamic_pointer_cast<ApproxCondDOS>(self->ptr);
    return wrap::to_py(ptr->dos);
}

static PyObject *
approx_cond_dos_get_kgrid(ApproxCondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ApproxCondDOSPtr ptr = std::dynamic_pointer_cast<ApproxCondDOS>(self->ptr);
    return kgrid_wrap(ptr->kgrid());
}

static PyObject *
approx_cond_dos_get_solve_cond(ApproxCondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ApproxCondDOSPtr ptr = std::dynamic_pointer_cast<ApproxCondDOS>(self->ptr);
    return wrap::to_py(ptr->solve_cond);
}

static int
approx_cond_dos_set_solve_cond(ApproxCondDOSPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    ApproxCondDOSPtr ptr = std::dynamic_pointer_cast<ApproxCondDOS>(self->ptr);
    return wrap::to_cpp(value, ptr->solve_cond);
}
static PyObject *
approx_cond_dos_get_solve_dos(ApproxCondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ApproxCondDOSPtr ptr = std::dynamic_pointer_cast<ApproxCondDOS>(self->ptr);
    return wrap::to_py(ptr->solve_dos);
}

static int
approx_cond_dos_set_solve_dos(ApproxCondDOSPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    ApproxCondDOSPtr ptr = std::dynamic_pointer_cast<ApproxCondDOS>(self->ptr);
    return wrap::to_cpp(value, ptr->solve_dos);
}


static PyObject *
approx_cond_dos_get_solve_dos_k(ApproxCondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ApproxCondDOSPtr ptr = std::dynamic_pointer_cast<ApproxCondDOS>(self->ptr);
    return wrap::to_py(ptr->solve_dos_k);
}

static int
approx_cond_dos_set_solve_dos_k(ApproxCondDOSPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    ApproxCondDOSPtr ptr = std::dynamic_pointer_cast<ApproxCondDOS>(self->ptr);
    return wrap::to_cpp(value, ptr->solve_dos_k);
}

static PyObject *
approx_cond_dos_get_w(ApproxCondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ApproxCondDOSPtr ptr = std::dynamic_pointer_cast<ApproxCondDOS>(self->ptr);
    return wrap::to_py(ptr->w);
}

static PyObject *
approx_cond_dos_get_eta(ApproxCondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ApproxCondDOSPtr ptr = std::dynamic_pointer_cast<ApproxCondDOS>(self->ptr);
    return wrap::to_py(ptr->eta);
}


//----------------------------------------------------------------------

static PyMethodDef approx_cond_dos_methods[] = {
    {"set_kgrid", (PyCFunction)approx_cond_dos_set_kgrid, METH_O,
        "Set the self-energy"},
    {"set_w", (PyCFunction)approx_cond_dos_set_w, METH_O,
        "Set the frequencies to solve for"},
    {"set_eta", (PyCFunction)approx_cond_dos_set_eta, METH_O,
        "Broadening to solve for"},
    {NULL}
};

static PyGetSetDef approx_cond_dos_getseters[] = {
    {"kgrid", (getter)approx_cond_dos_get_kgrid, NULL,
        "Get the self-energy", NULL},
    {"w", (getter)approx_cond_dos_get_w, NULL,
        "Frequencies being computed", NULL},
    {"eta", (getter)approx_cond_dos_get_eta, NULL,
        "Broadening to solve for", NULL},
    {"cond", (getter)approx_cond_dos_get_cond, NULL,
        "Get the Conductivity Matrix [ieta,iw,i,j]", NULL},
    {"cond_state", (getter)approx_cond_dos_get_cond, NULL,
        "Get the Conductivity per state [ieta,iw,i,j,s]", NULL},
    {"dos_k", (getter)approx_cond_dos_get_dos_k, NULL,
        "DOS per kpoint dos_k[ieta,iw,ik1,ik2,ik3]", NULL},
    {"dos", (getter)approx_cond_dos_get_dos, NULL,
        "Density of states dos[ieta,iw]", NULL},
    {"solve_cond", (getter)approx_cond_dos_get_solve_cond,
        (setter)approx_cond_dos_set_solve_cond,
        "Whether to solve for conductivity", NULL},
    {"solve_dos", (getter)approx_cond_dos_get_solve_dos,
        (setter)approx_cond_dos_set_solve_dos,
        "Whether to solve for density of states", NULL},
    {"solve_dos_k", (getter)approx_cond_dos_get_solve_dos_k,
        (setter)approx_cond_dos_set_solve_dos_k,
        "Whether to solve for density of states per k-point", NULL},
    {"eta", (getter)approx_cond_dos_get_eta, (setter)approx_cond_dos_set_eta,
        "Minimum value for +0j when computing Green's Matrix", NULL},
    {NULL}
};

//----------------------------------------------------------------------

PyTypeObject ApproxCondDOSPyType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "tbpp.ApproxCondDOS",               /* tp_name */
    sizeof(ApproxCondDOSPy),            /* tp_basicsize */
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
    "TBPP ApproxCondDOS object",        /* tp_doc */
    0,                          /* tp_traverse */
    0,                          /* tp_clear */
    0,                          /* tp_richcompare */
    0,                          /* tp_weaklistoffset */
    0,                          /* tp_iter */
    0,                          /* tp_iternext */
    approx_cond_dos_methods,              /* tp_methods */
                 0,             /* tp_members */
    approx_cond_dos_getseters,            /* tp_getset */
    0,                          /* tp_base */
    0,                          /* tp_dict */
    0,                          /* tp_descr_get */
    0,                          /* tp_descr_set */
    0,                          /* tp_dictoffset */
    (initproc)approx_cond_dos_init,       /* tp_init */
    0,                          /* tp_alloc */
    approx_cond_dos_new,                  /* tp_new */
};

//----------------------------------------------------------------------

PyObject *
approx_cond_dos_wrap(ApproxCondDOSPtr ptr) {
    PyObject *ret = ApproxCondDOSPyType.tp_new((PyTypeObject*)&ApproxCondDOSPyType, NULL, NULL);
    NodePy* obj = reinterpret_cast<NodePy*>(ret);
    obj->ptr = std::dynamic_pointer_cast<Node>(ptr);
    return ret;
}

ApproxCondDOSPtr
approx_cond_dos_unwrap(PyObject* obj) {
    if(!ApproxCondDOSPy_Check(obj))
        return nullptr;
    NodePy* m = reinterpret_cast<NodePy*>(obj);
    return std::dynamic_pointer_cast<ApproxCondDOS>(m->ptr);
}

