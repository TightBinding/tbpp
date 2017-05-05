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

#include <python/kpath.h>
#include <python/model.h>
#include <python/wrap.h>
#include <python/wrap_numpy.h>
#include <memory>

using namespace tbpp;

//----------------------------------------------------------------------

static PyObject *
kpath_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    KPathPy *self  = NULL;
    self = (KPathPy *)type->tp_alloc(type, 0);
    if(self != NULL) {
        self->ptr = nullptr;
    }

    init_numpy();

    return (PyObject *)self;
}

static int
kpath_init(KPathPy *self, PyObject *args, PyObject *kwds)
{
    self->ptr = std::dynamic_pointer_cast<Node>(std::make_shared<KPath>());
    return 0;
}

//----------------------------------------------------------------------

static PyObject *
kpath_set_model(KPathPy *self, PyObject *obj) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KPathPtr ptr = std::dynamic_pointer_cast<KPath>(self->ptr);

    if(!ModelPy_Check(obj)) {
        PyErr_SetString(PyExc_RuntimeError, "Model Object Required");
        return NULL;
    }

    NodePy* model = reinterpret_cast<NodePy*>(obj);
    ptr->set_model(std::dynamic_pointer_cast<Model>(model->ptr));

    Py_RETURN_NONE;
}

PyDoc_STRVAR(kpath_add_kpoint_doc,
"Add a kpoint: label, [k1,k2,k3]");
static PyObject *
kpath_add_kpoint(KPathPy *self, PyObject *args) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KPathPtr ptr = std::dynamic_pointer_cast<KPath>(self->ptr);

    char* label;
    double k1,k2,k3;
    if(!PyArg_ParseTuple(args, "s(ddd)", &label, &k1, &k2, &k3))
        return NULL;
    try {
        ptr->add_kpoint(std::string(label), k1, k2, k3);
    } catch (std::exception const& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }
    Py_RETURN_NONE;
}


//----------------------------------------------------------------------

static PyObject *
kpath_get_steps(KPathPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KPathPtr ptr = std::dynamic_pointer_cast<KPath>(self->ptr);
    return wrap::to_py(ptr->steps);
}

static int
kpath_set_steps(KPathPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    KPathPtr ptr = std::dynamic_pointer_cast<KPath>(self->ptr);
    return wrap::to_cpp(value, ptr->steps);
}

PyDoc_STRVAR(kpath_get_k_doc,
"Get the k-point values");
static PyObject *
kpath_get_k(KPathPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KPathPtr ptr = std::dynamic_pointer_cast<KPath>(self->ptr);
    return wrap::to_py(ptr->k);
}

PyDoc_STRVAR(kpath_get_eigval_doc,
"Get the eigenvalues");
static PyObject *
kpath_get_eigval(KPathPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KPathPtr ptr = std::dynamic_pointer_cast<KPath>(self->ptr);
    return wrap::to_py(ptr->eigval);
}

PyDoc_STRVAR(kpath_get_eigvec_doc,
"Get the eigenvectors");
static PyObject *
kpath_get_eigvec(KPathPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KPathPtr ptr = std::dynamic_pointer_cast<KPath>(self->ptr);
    return wrap::to_py(ptr->eigvec);
}

PyDoc_STRVAR(kpath_get_kpoints_doc,
"Get the list of kpoints");
static PyObject *
kpath_get_kpoints(KPathPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KPathPtr ptr = std::dynamic_pointer_cast<KPath>(self->ptr);

    PyObject *array;
    PyObject *dtype_dict;
    PyArray_Descr *dtype;

    dtype_dict = Py_BuildValue(
            "[(s,s),(s,s),(s,s),(s,s)]",
            "label","a64",
               "k1","f8",
               "k2","f8",
               "k3","f8"
            );

    PyArray_DescrConverter(dtype_dict, &dtype);
    Py_DECREF(dtype_dict);

    int nd = 1;
    npy_intp dims[1] = {static_cast<npy_intp>(ptr->kpoints.size())};

    array = PyArray_NewFromDescr((PyTypeObject*)&PyArray_Type,
            dtype, nd, dims, NULL, (void*)ptr->kpoints.data(),
            NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_WRITEABLE, NULL);

    return array;
}

static PyObject *
kpath_get_model(KPathPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KPathPtr ptr = std::dynamic_pointer_cast<KPath>(self->ptr);
    return model_wrap(ptr->model());
}

static PyObject *
kpath_get_solve_eigvec(KPathPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KPathPtr ptr = std::dynamic_pointer_cast<KPath>(self->ptr);
    return wrap::to_py(ptr->solve_eigvec);
}

static int
kpath_set_solve_eigvec(KPathPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    KPathPtr ptr = std::dynamic_pointer_cast<KPath>(self->ptr);
    return wrap::to_cpp(value, ptr->solve_eigvec);
}

//----------------------------------------------------------------------

static PyMethodDef kpath_methods[] = {
    {"set_model", (PyCFunction)kpath_set_model, METH_O, "Set the model"},
    {"add_kpoint", (PyCFunction)kpath_add_kpoint, METH_VARARGS, kpath_add_kpoint_doc},
    {NULL}
};

static PyGetSetDef kpath_getseters[] = {
    {"kpoints", (getter)kpath_get_kpoints, NULL,
        kpath_get_kpoints_doc, NULL},
    {"eigval", (getter)kpath_get_eigval, NULL,
        kpath_get_eigval_doc, NULL},
    {"eigvec", (getter)kpath_get_eigvec, NULL,
        kpath_get_eigvec_doc, NULL},
    {"k", (getter)kpath_get_k, NULL,
        kpath_get_k_doc, NULL},
    {"steps", (getter)kpath_get_steps, (setter)kpath_set_steps,
        "Steps to take between k-points", NULL},
    {"model", (getter)kpath_get_model, NULL,
        "Get the Model being used", NULL},
    {"solve_eigvec", (getter)kpath_get_solve_eigvec, (setter)kpath_set_solve_eigvec,
        "Whether to solve for eigenvectors", NULL},
    {NULL}
};

//----------------------------------------------------------------------

PyTypeObject KPathPyType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "tbpp.KPath",               /* tp_name */
    sizeof(KPathPy),            /* tp_basicsize */
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
    "TBPP KPath object",        /* tp_doc */
    0,                          /* tp_traverse */
    0,                          /* tp_clear */
    0,                          /* tp_richcompare */
    0,                          /* tp_weaklistoffset */
    0,                          /* tp_iter */
    0,                          /* tp_iternext */
    kpath_methods,              /* tp_methods */
                 0,             /* tp_members */
    kpath_getseters,            /* tp_getset */
    0,                          /* tp_base */
    0,                          /* tp_dict */
    0,                          /* tp_descr_get */
    0,                          /* tp_descr_set */
    0,                          /* tp_dictoffset */
    (initproc)kpath_init,       /* tp_init */
    0,                          /* tp_alloc */
    kpath_new,                  /* tp_new */
};

//----------------------------------------------------------------------

PyObject *
kpath_wrap(KPathPtr ptr) {
    PyObject *ret = KPathPyType.tp_new((PyTypeObject*)&KPathPyType, NULL, NULL);
    NodePy* obj = reinterpret_cast<NodePy*>(ret);
    obj->ptr = std::dynamic_pointer_cast<Node>(ptr);
    return ret;
}

KPathPtr
kpath_unwrap(PyObject* obj) {
    if(!KPathPy_Check(obj))
        return nullptr;
    NodePy* m = reinterpret_cast<NodePy*>(obj);
    return std::dynamic_pointer_cast<KPath>(m->ptr);
}

