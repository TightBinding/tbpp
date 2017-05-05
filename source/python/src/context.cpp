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

#include "python/context.h"
#include "python/wrap.h"
#include "python/node.h"

#include <tbpp/node.h>
#include <exception>

using namespace tbpp;

static void
context_dealloc(ContextPy* self)
{
    self->ptr = nullptr; //< Required to free shared_ptr reference
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject *
context_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    ContextPy *self  = NULL;
    self = (ContextPy *)type->tp_alloc(type, 0);
    if(self != NULL) {
        self->ptr = nullptr;
    }
    return (PyObject *)self;
}

static int
context_init(ContextPy *self, PyObject *args, PyObject *kwds)
{
    self->ptr = std::make_shared<Context>();
    return 0;
}

//----------------------------------------------------------------------

static PyObject *
context_create(ContextPy* self, PyObject *args, PyObject *kwds) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }

    const char* c_type = NULL;
    const char* c_name = NULL;
    static const char *kwlist[] = {"type","name", NULL};
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|s", (char**)kwlist,
                &c_type, &c_name))
        return NULL;

    std::string type;
    type.assign(c_type);

    std::string name;
    if(c_name) name.assign(c_name);

    NodePtr c = self->ptr->create(type, name);
    return node_dynamic_wrap(c);
}

static PyObject *
context_get(ContextPy* self, PyObject *args) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }

    const char* c_name = NULL;
    if(!PyArg_ParseTuple(args, "s", &c_name))
        return NULL;

    std::string name;
    if(c_name) name.assign(c_name);

    NodePtr c = self->ptr->get(name);
    if (c == nullptr) {
        PyErr_SetString(PyExc_RuntimeError, "No object");
        return NULL;
    }
    return node_dynamic_wrap(c);
}

static PyObject *
context_add(ContextPy* self, PyObject *args) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }

    PyObject* obj;
    const char* c_name = NULL;
    if(!PyArg_ParseTuple(args, "O|s", &obj, &c_name))
        return NULL;

    std::string name;
    if(c_name) name.assign(c_name);

    if(!NodePy_Check(obj)) {
        PyErr_SetString(PyExc_RuntimeError, "Node Object Required");
        return NULL;
    }

    NodePtr ret = self->ptr->add(reinterpret_cast<NodePy*>(obj)->ptr, name);
    if (ret == nullptr) {
        PyErr_SetString(PyExc_RuntimeError, "Could not add object");
        return NULL;
    }
    //return obj;
    Py_RETURN_NONE;
}


static PyObject *
context_rm(ContextPy* self, PyObject *args) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }

    const char* c_name = NULL;
    if(!PyArg_ParseTuple(args, "s", &c_name))
        return NULL;

    std::string name;
    if(c_name) name.assign(c_name);

    NodePtr c = self->ptr->rm(name);
    if (c == nullptr) {
        PyErr_SetString(PyExc_RuntimeError, "No object");
        return NULL;
    }
    return node_dynamic_wrap(c);
}

static PyObject *
context_clear(ContextPy* self) {
    if(self->ptr == nullptr)
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
    else
        self->ptr->clear();
     Py_RETURN_NONE;
}

static PyObject *
context_save(ContextPy* self, PyObject *args) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }

    const char* c_filename = NULL;
    if(!PyArg_ParseTuple(args, "s", &c_filename))
        return NULL;

    std::string filename;
    filename.assign(c_filename);

    self->ptr->save(filename);
    Py_RETURN_NONE;
}

static PyObject *
context_load(ContextPy* self, PyObject *args) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }

    const char* c_filename = NULL;
    if(!PyArg_ParseTuple(args, "s", &c_filename))
        return NULL;

    std::string filename;
    filename.assign(c_filename);

    try {
        self->ptr->load(filename);
    } catch(std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject *
context_make_ready(ContextPy* self) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }

    try {
        self->ptr->make_ready();
    } catch(std::exception& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

static PyObject *
context_solve(ContextPy* self) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }

    try {
        self->ptr->solve();
    } catch(std::exception& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

static PyObject *
context_nodes(ContextPy* self) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }

    PyObject *list = PyList_New(0);
    for(const auto& n : self->ptr->nodes) {
        PyList_Append(list, wrap::to_py(n.first));
    }

    return list;
}



//----------------------------------------------------------------------

static PyMethodDef context_methods[] = {
    {"create", (PyCFunction)context_create, METH_VARARGS | METH_KEYWORDS,
        "Create a new Node"},
    {"get", (PyCFunction)context_get, METH_VARARGS, "Get a Node by name"},
    {"add", (PyCFunction)context_add, METH_VARARGS, "Add an existing Node"},
    {"rm", (PyCFunction)context_rm, METH_VARARGS, "Remove a Node by name"},
    {"solve", (PyCFunction)context_solve, METH_NOARGS, "Solve or resolve all nodes"},
    {"make_ready", (PyCFunction)context_make_ready, METH_NOARGS, "Solve those nodes which are currently not solved"},
    {"clear", (PyCFunction)context_clear, METH_NOARGS, "Delete all nodes"},
    {"save", (PyCFunction)context_save, METH_VARARGS, "Save Context to an HDF5 file"},
    {"load", (PyCFunction)context_load, METH_VARARGS, "Load Context from an HDF5 file"},
    {"nodes", (PyCFunction)context_nodes, METH_NOARGS, "Return a list of names for all nodes in the context"},
    {NULL}
};

//----------------------------------------------------------------------

PyTypeObject ContextPyType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "tbpp.Context",          /* tp_name */
    sizeof(ContextPy),         /* tp_basicsize */
    0,                          /* tp_itemsize */
    (destructor)context_dealloc, /* tp_dealloc */
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
    "Context object",    /* tp_doc */
    0,                          /* tp_traverse */
    0,                          /* tp_clear */
    0,                          /* tp_richcompare */
    0,                          /* tp_weaklistoffset */
    0,                          /* tp_iter */
    0,                          /* tp_iternext */
    context_methods,               /* tp_methods */
                 0,             /* tp_members */
                 0,             /* tp_getset */
    0,                          /* tp_base */
    0,                          /* tp_dict */
    0,                          /* tp_descr_get */
    0,                          /* tp_descr_set */
    0,                          /* tp_dictoffset */
    (initproc)context_init,        /* tp_init */
    0,                          /* tp_alloc */
    context_new,                   /* tp_new */
};

//----------------------------------------------------------------------

PyObject *
context_wrap(ContextPtr ptr) {
    PyObject *ret = context_new((PyTypeObject*)&ContextPyType, NULL, NULL);
    ContextPy* obj = reinterpret_cast<ContextPy*>(ret);
    obj->ptr = ptr;
    return ret;
}

ContextPtr
context_unwrap(PyObject* obj) {
    if(!ContextPy_Check(obj))
        return nullptr;
    ContextPy* m = reinterpret_cast<ContextPy*>(obj);
    return m->ptr;
}

