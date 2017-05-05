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

#include <python/cond_dos.h>
#include <python/cpa.h>
#include <python/wrap.h>
#include <python/wrap_numpy.h>
#include <memory>

using namespace tbpp;

//----------------------------------------------------------------------

static PyObject *
cond_dos_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    CondDOSPy *self  = NULL;
    self = (CondDOSPy *)type->tp_alloc(type, 0);
    if(self != NULL) {
        self->ptr = nullptr;
    }

    init_numpy();

    return (PyObject *)self;
}

static int
cond_dos_init(CondDOSPy *self, PyObject *args, PyObject *kwds)
{
    self->ptr = std::dynamic_pointer_cast<Node>(std::make_shared<CondDOS>());
    return 0;
}

//----------------------------------------------------------------------

static PyObject *
cond_dos_set_sigma(CondDOSPy *self, PyObject *obj) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);

    if(!CPAPy_Check(obj)) {
        PyErr_SetString(PyExc_RuntimeError, "CPA Object Required");
        return NULL;
    }

    NodePy* node = reinterpret_cast<NodePy*>(obj);
    ptr->set_sigma(std::dynamic_pointer_cast<CPA>(node->ptr));

    Py_RETURN_NONE;
}
//----------------------------------------------------------------------

static PyObject *
cond_dos_get_cond(CondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_py(ptr->cond);
}

static PyObject *
cond_dos_get_cond_state(CondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_py<double, 5>(ptr->cond_state);
}

static PyObject *
cond_dos_get_dos_k(CondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_py(ptr->dos_k);
}

static PyObject *
cond_dos_get_dos(CondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_py(ptr->dos);
}

static PyObject *
cond_dos_get_dos_state(CondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_py(ptr->dos_state);
}

static PyObject *
cond_dos_get_dos_proj(CondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_py(ptr->dos_proj);
}

static PyObject *
cond_dos_get_dos_weights(CondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_py(ptr->dos_weights);
}


static PyObject *
cond_dos_get_sigma(CondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return cpa_wrap(ptr->sigma());
}

static PyObject *
cond_dos_get_solve_cond(CondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_py(ptr->solve_cond);
}

static int
cond_dos_set_solve_cond(CondDOSPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_cpp(value, ptr->solve_cond);
}

static PyObject *
cond_dos_get_solve_sigx(CondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_py(ptr->solve_sigx);
}

static int
cond_dos_set_solve_sigx(CondDOSPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_cpp(value, ptr->solve_sigx);
}

static PyObject *
cond_dos_get_solve_sigy(CondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_py(ptr->solve_sigy);
}

static int
cond_dos_set_solve_sigy(CondDOSPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_cpp(value, ptr->solve_sigy);
}

static PyObject *
cond_dos_get_solve_sigz(CondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_py(ptr->solve_sigz);
}

static int
cond_dos_set_solve_sigz(CondDOSPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_cpp(value, ptr->solve_sigz);
}


static PyObject *
cond_dos_get_solve_dos(CondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_py(ptr->solve_dos);
}

static int
cond_dos_set_solve_dos(CondDOSPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_cpp(value, ptr->solve_dos);
}

static PyObject *
cond_dos_get_solve_dos_state(CondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_py(ptr->solve_dos_state);
}

static int
cond_dos_set_solve_dos_state(CondDOSPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_cpp(value, ptr->solve_dos_state);
}

static PyObject *
cond_dos_get_solve_dos_k(CondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_py(ptr->solve_dos_k);
}

static int
cond_dos_set_solve_dos_k(CondDOSPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_cpp(value, ptr->solve_dos_k);
}

static PyObject *
cond_dos_get_zero_j(CondDOSPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_py(ptr->zero_j);
}

static int
cond_dos_set_zero_j(CondDOSPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);
    return wrap::to_cpp(value, ptr->zero_j);
}

static PyObject *
cond_dos_add_dos_weights(CondDOSPy *self, PyObject *obj) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CondDOSPtr ptr = std::dynamic_pointer_cast<CondDOS>(self->ptr);

    std::vector<double> w;
    if(wrap::to_cpp(obj, w)) {
        return NULL;
    }

    try {
        ptr->add_dos_weights(w);
    } catch(std::exception const& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}


//----------------------------------------------------------------------

static PyMethodDef cond_dos_methods[] = {
    {"set_sigma", (PyCFunction)cond_dos_set_sigma, METH_O,
        "Set the self-energy"},
    {"add_dos_weights", (PyCFunction)cond_dos_add_dos_weights, METH_O,
        "Add a set of weight for DOS projection"},
    {NULL}
};

static PyGetSetDef cond_dos_getseters[] = {
    {"sigma", (getter)cond_dos_get_sigma, NULL,
        "Get the self-energy", NULL},
    {"cond", (getter)cond_dos_get_cond, NULL,
        "Get the Conductivity Matrix [ic,iw,i,j]", NULL},
    {"cond_state", (getter)cond_dos_get_cond_state, NULL,
        "Get the Conductivity per state [ic,iw,i,j,s]", NULL},
    {"dos_k", (getter)cond_dos_get_dos_k, NULL,
        "DOS per kpoint dos_k[ic,iw,ik1,ik2,ik3]", NULL},
    {"dos", (getter)cond_dos_get_dos, NULL,
        "Density of states dos[ic,iw]", NULL},
    {"dos_state", (getter)cond_dos_get_dos_state, NULL,
        "Density of states per state dos_state[ic,iw,is]", NULL},
    {"dos_proj", (getter)cond_dos_get_dos_proj, NULL,
        "Projected density of states dos_proj[ic,iw,ip,ik1,ik2,ik3]", NULL},
    {"dos_weights", (getter)cond_dos_get_dos_weights, NULL,
        "Weights for density of states projection dos_weights[ip,is]", NULL},
    {"solve_cond", (getter)cond_dos_get_solve_cond,
        (setter)cond_dos_set_solve_cond,
        "Whether to solve for conductivity", NULL},
    {"solve_sigx", (getter)cond_dos_get_solve_sigx,
        (setter)cond_dos_set_solve_sigx,
        "Whether to solve for conductivity", NULL},
    {"solve_sigy", (getter)cond_dos_get_solve_sigy,
        (setter)cond_dos_set_solve_sigy,
        "Whether to solve for conductivity", NULL},
    {"solve_sigz", (getter)cond_dos_get_solve_sigz,
        (setter)cond_dos_set_solve_sigz,
        "Whether to solve for conductivity", NULL},
    {"solve_dos", (getter)cond_dos_get_solve_dos,
        (setter)cond_dos_set_solve_dos,
        "Whether to solve for density of states", NULL},
    {"solve_dos_state", (getter)cond_dos_get_solve_dos_state,
        (setter)cond_dos_set_solve_dos_state,
        "Whether to solve for density of states per state", NULL},
    {"solve_dos_k", (getter)cond_dos_get_solve_dos_k,
        (setter)cond_dos_set_solve_dos_k,
        "Whether to solve for density of states per k-point", NULL},
    {"zero_j", (getter)cond_dos_get_zero_j, (setter)cond_dos_set_zero_j,
        "Minimum value for +0j when computing Green's Matrix", NULL},
    {NULL}
};

//----------------------------------------------------------------------

PyTypeObject CondDOSPyType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "tbpp.CondDOS",               /* tp_name */
    sizeof(CondDOSPy),            /* tp_basicsize */
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
    "TBPP CondDOS object",        /* tp_doc */
    0,                          /* tp_traverse */
    0,                          /* tp_clear */
    0,                          /* tp_richcompare */
    0,                          /* tp_weaklistoffset */
    0,                          /* tp_iter */
    0,                          /* tp_iternext */
    cond_dos_methods,              /* tp_methods */
                 0,             /* tp_members */
    cond_dos_getseters,            /* tp_getset */
    0,                          /* tp_base */
    0,                          /* tp_dict */
    0,                          /* tp_descr_get */
    0,                          /* tp_descr_set */
    0,                          /* tp_dictoffset */
    (initproc)cond_dos_init,       /* tp_init */
    0,                          /* tp_alloc */
    cond_dos_new,                  /* tp_new */
};

//----------------------------------------------------------------------

PyObject *
cond_dos_wrap(CondDOSPtr ptr) {
    PyObject *ret = CondDOSPyType.tp_new((PyTypeObject*)&CondDOSPyType, NULL, NULL);
    NodePy* obj = reinterpret_cast<NodePy*>(ret);
    obj->ptr = std::dynamic_pointer_cast<Node>(ptr);
    return ret;
}

CondDOSPtr
cond_dos_unwrap(PyObject* obj) {
    if(!CondDOSPy_Check(obj))
        return nullptr;
    NodePy* m = reinterpret_cast<NodePy*>(obj);
    return std::dynamic_pointer_cast<CondDOS>(m->ptr);
}

