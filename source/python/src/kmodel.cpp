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

#include <python/kmodel.h>
#include <python/model.h>
#include <python/wrap.h>
#include <python/wrap_numpy.h>
#include <memory>

using namespace tbpp;
using namespace std;

//----------------------------------------------------------------------

static PyObject *
kmodel_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    KModelPy *self  = NULL;
    self = (KModelPy *)type->tp_alloc(type, 0);
    if(self != NULL) {
        self->ptr = nullptr;
    }

    init_numpy();

    return (PyObject *)self;
}

static int
kmodel_init(KModelPy *self, PyObject *args, PyObject *kwds)
{
    self->ptr = std::dynamic_pointer_cast<Node>(std::make_shared<KModel>());
    return 0;
}

//----------------------------------------------------------------------

static PyObject *
kmodel_set_model(KModelPy *self, PyObject *obj) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KModelPtr ptr = std::dynamic_pointer_cast<KModel>(self->ptr);

    if(!ModelPy_Check(obj)) {
        PyErr_SetString(PyExc_RuntimeError, "Model Object Required");
        return NULL;
    }

    NodePy* model = reinterpret_cast<NodePy*>(obj);
    ptr->set_model(std::dynamic_pointer_cast<Model>(model->ptr));

    Py_RETURN_NONE;
}

static PyObject *
kmodel_set_cache(KModelPy *self, PyObject *obj) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KModelPtr ptr = std::dynamic_pointer_cast<KModel>(self->ptr);

    bool value;
    if (wrap::to_cpp(obj, value)) {
        return NULL;
    }
    ptr->set_cache(value);

    Py_RETURN_NONE;
}

static PyObject *
kmodel_set_grid(KModelPy *self, PyObject *args) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KModelPtr ptr = std::dynamic_pointer_cast<KModel>(self->ptr);

    unsigned int nk1=0;
    unsigned int nk2=0;
    unsigned int nk3=0;
    if(!PyArg_ParseTuple(args, "|III", &nk1, &nk2, &nk3))
        return NULL;

    ptr->set_grid(nk1,nk2,nk3);
    Py_RETURN_NONE;
}


//----------------------------------------------------------------------

static PyObject *
kmodel_get_model(KModelPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KModelPtr ptr = std::dynamic_pointer_cast<KModel>(self->ptr);
    return model_wrap(ptr->model());
}

static PyObject *
kmodel_get_Tk(KModelPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KModelPtr ptr = std::dynamic_pointer_cast<KModel>(self->ptr);
    return wrap::to_py(ptr->Tk);
}

static PyObject *
kmodel_get_cache(KModelPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KModelPtr ptr = std::dynamic_pointer_cast<KModel>(self->ptr);
    return wrap::to_py(ptr->cache());
}

static PyObject *
kmodel_get_k1(KModelPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KModelPtr ptr = std::dynamic_pointer_cast<KModel>(self->ptr);
    return wrap::to_py(ptr->k1);
}

static PyObject *
kmodel_get_k2(KModelPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KModelPtr ptr = std::dynamic_pointer_cast<KModel>(self->ptr);
    return wrap::to_py(ptr->k2);
}

static PyObject *
kmodel_get_k3(KModelPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KModelPtr ptr = std::dynamic_pointer_cast<KModel>(self->ptr);
    return wrap::to_py(ptr->k3);
}

PyDoc_STRVAR(kmodel_G_doc,
"G(Ef, zeroj)\n"
"\n Local Green's matrix. Note you must first call set_grid(...) to set the"
"\n k-space grid for integration");
static PyObject *
kmodel_G(ModelPy *self, PyObject *args) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KModelPtr ptr = std::dynamic_pointer_cast<KModel>(self->ptr);

    // Get pointer to model
    ModelPtr mptr = ptr->model();
    if(mptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "No model");
        return NULL;
    }

    double Ef=0, zeroj=0;
    if(!PyArg_ParseTuple(args, "dd", &Ef, &zeroj))
        return NULL;

    npy_intp dims[2];
    dims[0] = mptr->states();
    dims[1] = mptr->states();

    PyObject* array_obj = PyArray_SimpleNew(2, dims, NPY_CDOUBLE);
    PyArrayObject* array = reinterpret_cast<PyArrayObject*>(array_obj);
    std::complex<double>* data = static_cast<std::complex<double>*>(PyArray_DATA(array));
    try {
        ptr->G(data, Ef, zeroj, NULL);
    } catch(std::exception const& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        Py_DECREF(array_obj);
        return NULL;
    }
    return array_obj;
}

//----------------------------------------------------------------------

static PyMethodDef kmodel_methods[] = {
    {"set_model", (PyCFunction)kmodel_set_model, METH_O,
        "Set the model"},
    {"set_cache", (PyCFunction)kmodel_set_cache, METH_O,
        "Set the cache"},
    {"set_grid", (PyCFunction)kmodel_set_grid, METH_VARARGS,
        "Set the k space grid size"},
    {"G", (PyCFunction)kmodel_G, METH_VARARGS, kmodel_G_doc},
    {NULL}
};

static PyGetSetDef kmodel_getseters[] = {
    {"model", (getter)kmodel_get_model, NULL,
        "The model being used.", NULL},
    {"Tk", (getter)kmodel_get_Tk, NULL,
        "Get the Tk[ik1,ik2,ik3,i,j]", NULL},
    {"cache", (getter)kmodel_get_cache, NULL,
        "Whether to cache Tk matrix", NULL},
    {"k1", (getter)kmodel_get_k1, NULL, "Get k1", NULL},
    {"k2", (getter)kmodel_get_k2, NULL, "Get k2", NULL},
    {"k3", (getter)kmodel_get_k3, NULL, "Get k3", NULL},
    {NULL}
};

//----------------------------------------------------------------------

PyTypeObject KModelPyType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "tbpp.KModel",              /* tp_name */
    sizeof(KModelPy),           /* tp_basicsize */
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
    "TBPP KModel object",       /* tp_doc */
    0,                          /* tp_traverse */
    0,                          /* tp_clear */
    0,                          /* tp_richcompare */
    0,                          /* tp_weaklistoffset */
    0,                          /* tp_iter */
    0,                          /* tp_iternext */
    kmodel_methods,             /* tp_methods */
                 0,             /* tp_members */
    kmodel_getseters,           /* tp_getset */
    0,                          /* tp_base */
    0,                          /* tp_dict */
    0,                          /* tp_descr_get */
    0,                          /* tp_descr_set */
    0,                          /* tp_dictoffset */
    (initproc)kmodel_init,      /* tp_init */
    0,                          /* tp_alloc */
    kmodel_new,                 /* tp_new */
};

//----------------------------------------------------------------------

PyObject *
kmodel_wrap(KModelPtr ptr) {
    PyObject *ret = KModelPyType.tp_new((PyTypeObject*)&KModelPyType, NULL, NULL);
    NodePy* obj = reinterpret_cast<NodePy*>(ret);
    obj->ptr = std::dynamic_pointer_cast<Node>(ptr);
    return ret;
}

KModelPtr
kmodel_unwrap(PyObject* obj) {
    if(!KModelPy_Check(obj))
        return nullptr;
    NodePy* m = reinterpret_cast<NodePy*>(obj);
    return std::dynamic_pointer_cast<KModel>(m->ptr);
}

