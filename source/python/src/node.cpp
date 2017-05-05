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

#include "python/node.h"
#include "python/wrap.h"
#include <string>

#include "python/model.h"
#include "python/lattice.h"
#include "python/kpath.h"
#include "python/kgrid.h"
#include "python/kmodel.h"
#include "python/cpa.h"
#include "python/gen_cond.h"
#include "python/cond_dos.h"
#include "python/approx_cond_dos.h"
#include "python/berryflux.h"

#include <tbpp/node.h>

using namespace tbpp;

static void
node_dealloc(NodePy* self)
{
    self->ptr = nullptr; //< Required to free shared_ptr reference
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject *
node_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    NodePy *self  = NULL;
    self = (NodePy *)type->tp_alloc(type, 0);
    if(self != NULL) {
        self->ptr = nullptr;
    }
    return (PyObject *)self;
}

static int
node_init(NodePy *self, PyObject *args, PyObject *kwds)
{
    self->ptr = std::make_shared<Node>();
    return 0;
}

//----------------------------------------------------------------------

static PyObject *
node_type(NodePy* self) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    return wrap::to_py(self->ptr->type());
}

static PyObject *
node_name(NodePy* self) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    return wrap::to_py(self->ptr->name());
}

static PyObject *
node_status(NodePy* self) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    std::string status = Node::node_status_to_string(self->ptr->status());
    return wrap::to_py(status);
}

static PyObject *
node_solve(NodePy* self) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    try {
        self->ptr->solve();
    } catch (const std::exception& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject *
node_make_ready(NodePy* self) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    try {
        self->ptr->make_ready();
    } catch (const std::exception& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }
    Py_RETURN_NONE;
}

//----------------------------------------------------------------------

static PyObject *
node_get_parallel(NodePy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    NodePtr ptr = std::dynamic_pointer_cast<Node>(self->ptr);
    return wrap::to_py(ptr->parallel());
}

static int
node_set_parallel(NodePy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    NodePtr ptr = std::dynamic_pointer_cast<Node>(self->ptr);
    bool parallel;
    if(wrap::to_cpp(value, parallel))
        return -1;
    ptr->set_parallel(parallel);
    return 0;
}

//----------------------------------------------------------------------

static PyMethodDef node_methods[] = {
    {"type", (PyCFunction)node_type, METH_NOARGS, "Returns the TBPP type of the node as a string"},
    {"name", (PyCFunction)node_name, METH_NOARGS, "Returns the name of the node"},
    {"status", (PyCFunction)node_status, METH_NOARGS, "Returns the status of the node"},
    {"solve", (PyCFunction)node_solve, METH_NOARGS, "Solve the node regardless of whether it has already been solved"},
    {"make_ready", (PyCFunction)node_make_ready, METH_NOARGS, "Solve the node only if necessary"},
    {NULL}
};

static PyGetSetDef node_getseters[] = {
    {"parallel", (getter)node_get_parallel, (setter)node_set_parallel,
        "Whether to use multithreading", NULL},
    {NULL}
};

//----------------------------------------------------------------------

PyTypeObject NodePyType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "tbpp.Node",          /* tp_name */
    sizeof(NodePy),         /* tp_basicsize */
    0,                          /* tp_itemsize */
    (destructor)node_dealloc, /* tp_dealloc */
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
    Py_TPFLAGS_DEFAULT,         /* tp_flags */
    "Base class for all TBPP Node objects",    /* tp_doc */
    0,                          /* tp_traverse */
    0,                          /* tp_clear */
    0,                          /* tp_richcompare */
    0,                          /* tp_weaklistoffset */
    0,                          /* tp_iter */
    0,                          /* tp_iternext */
    node_methods,               /* tp_methods */
                 0,             /* tp_members */
    node_getseters,             /* tp_getset */
    0,                          /* tp_base */
    0,                          /* tp_dict */
    0,                          /* tp_descr_get */
    0,                          /* tp_descr_set */
    0,                          /* tp_dictoffset */
    (initproc)node_init,        /* tp_init */
    0,                          /* tp_alloc */
    node_new,                   /* tp_new */
};

//----------------------------------------------------------------------

PyObject *
node_wrap(NodePtr ptr) {
    PyObject *ret = node_new((PyTypeObject*)&NodePyType, NULL, NULL);
    NodePy* obj = reinterpret_cast<NodePy*>(ret);
    obj->ptr = ptr;
    return ret;
}

NodePtr
node_unwrap(PyObject* obj) {
    if(!NodePy_Check(obj))
        return nullptr;
    NodePy* m = reinterpret_cast<NodePy*>(obj);
    return m->ptr;
}

PyObject *
node_dynamic_wrap(tbpp::NodePtr ptr) {
    if(ptr == nullptr) return NULL;
    std::string type = ptr->type();

    if(type == "Lattice")
        return lattice_wrap(std::dynamic_pointer_cast<tbpp::Lattice>(ptr));
    if(type == "Model")
        return model_wrap(std::dynamic_pointer_cast<tbpp::Model>(ptr));
    if(type == "KPath")
        return kpath_wrap(std::dynamic_pointer_cast<tbpp::KPath>(ptr));
    if(type == "KGrid")
        return kgrid_wrap(std::dynamic_pointer_cast<tbpp::KGrid>(ptr));
    if(type == "KModel")
        return kmodel_wrap(std::dynamic_pointer_cast<tbpp::KModel>(ptr));
    if(type == "CPA")
        return cpa_wrap(std::dynamic_pointer_cast<tbpp::CPA>(ptr));
    if(type == "CondDOS")
        return cond_dos_wrap(std::dynamic_pointer_cast<tbpp::CondDOS>(ptr));
    if(type == "GenCond")
        return gen_cond_wrap(std::dynamic_pointer_cast<tbpp::GenCond>(ptr));
    if(type == "ApproxCondDOS")
        return approx_cond_dos_wrap(std::dynamic_pointer_cast<tbpp::ApproxCondDOS>(ptr));
    if(type == "BerryFlux")
        return berryflux_wrap(std::dynamic_pointer_cast<tbpp::BerryFlux>(ptr));
    return node_wrap(ptr);
}
