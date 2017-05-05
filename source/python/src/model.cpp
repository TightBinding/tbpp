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

#include <python/model.h>
#include <python/lattice.h>
#include <python/wrap.h>
#include <python/wrap_numpy.h>
#include <memory>

using namespace tbpp;

//----------------------------------------------------------------------

static PyObject *
model_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    ModelPy *self  = NULL;
    self = (ModelPy *)type->tp_alloc(type, 0);
    if(self != NULL) {
        self->ptr = nullptr;
    }

    init_numpy();

    return (PyObject *)self;
}

static int
model_init(ModelPy *self, PyObject *args, PyObject *kwds)
{
    self->ptr = std::dynamic_pointer_cast<Node>(std::make_shared<Model>());
    return 0;
}

//----------------------------------------------------------------------


static PyObject *
model_uc_size(ModelPy *self) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ModelPtr ptr = std::dynamic_pointer_cast<Model>(self->ptr);

    double ret;
    try {
        ret = ptr->uc_size();
    } catch(std::exception const& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    return wrap::to_py(ret);
}

static PyObject *
model_Tk(ModelPy *self, PyObject *args) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ModelPtr ptr = std::dynamic_pointer_cast<Model>(self->ptr);

    double k1=0, k2=0, k3=0;
    if(!PyArg_ParseTuple(args, "ddd", &k1, &k2, &k3))
        return NULL;

    npy_intp dims[2];
    dims[0] = ptr->states();
    dims[1] = ptr->states();

    PyObject* array = PyArray_SimpleNew(2, dims, NPY_CDOUBLE);
    PyArrayObject* array_obj = reinterpret_cast<PyArrayObject*>(array);
    std::complex<double>* data = static_cast<std::complex<double>*>(PyArray_DATA(array_obj));
    try {
        ptr->Tk(data, k1, k2, k3);
    } catch(std::exception const& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        Py_DECREF(array);
        return NULL;
    }
    return array;
}

static PyObject *
model_Hk(ModelPy *self, PyObject *args) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ModelPtr ptr = std::dynamic_pointer_cast<Model>(self->ptr);

    double k1=0, k2=0, k3=0;
    if(!PyArg_ParseTuple(args, "ddd", &k1, &k2, &k3))
        return NULL;

    npy_intp dims[2];
    dims[0] = ptr->states();
    dims[1] = ptr->states();

    PyObject* array_obj = PyArray_SimpleNew(2, dims, NPY_CDOUBLE);
    PyArrayObject* array = reinterpret_cast<PyArrayObject*>(array_obj);
    std::complex<double>* data = static_cast<std::complex<double>*>(PyArray_DATA(array));
    try {
        ptr->Hk(data, k1, k2, k3);
    } catch(std::exception const& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        Py_DECREF(array_obj);
        return NULL;
    }
    return array_obj;
}

static PyObject *
model_dH_dk(ModelPy *self, PyObject *args) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ModelPtr ptr = std::dynamic_pointer_cast<Model>(self->ptr);

    double k1=0, k2=0, k3=0;
    if(!PyArg_ParseTuple(args, "ddd", &k1, &k2, &k3))
        return NULL;

    npy_intp dims[2];
    dims[0] = ptr->states();
    dims[1] = ptr->states();

    PyObject* dx_obj = PyArray_SimpleNew(2, dims, NPY_CDOUBLE);
    PyObject* dy_obj = PyArray_SimpleNew(2, dims, NPY_CDOUBLE);
    PyObject* dz_obj = PyArray_SimpleNew(2, dims, NPY_CDOUBLE);

    PyArrayObject* dx = reinterpret_cast<PyArrayObject*>(dx_obj);
    PyArrayObject* dy = reinterpret_cast<PyArrayObject*>(dy_obj);
    PyArrayObject* dz = reinterpret_cast<PyArrayObject*>(dz_obj);

    std::complex<double>* data_x = static_cast<std::complex<double>*>(PyArray_DATA(dx));
    std::complex<double>* data_y = static_cast<std::complex<double>*>(PyArray_DATA(dy));
    std::complex<double>* data_z = static_cast<std::complex<double>*>(PyArray_DATA(dz));
    try {
        ptr->dH_dk(data_x, data_y, data_z, k1, k2, k3);
    } catch(std::exception const& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        Py_DECREF(dx_obj);
        Py_DECREF(dy_obj);
        Py_DECREF(dz_obj);
        return NULL;
    }
    return Py_BuildValue("NNN", dx_obj, dy_obj, dz_obj);
}

static PyObject *
model_set_lattice(ModelPy *self, PyObject *obj) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ModelPtr ptr = std::dynamic_pointer_cast<Model>(self->ptr);

    if(!LatticePy_Check(obj)) {
        PyErr_SetString(PyExc_RuntimeError, "Lattice Object Required");
        return NULL;
    }

    NodePy* lattice = reinterpret_cast<NodePy*>(obj);
    ptr->set_lattice(std::dynamic_pointer_cast<Lattice>(lattice->ptr));

    Py_RETURN_NONE;
}

//----------------------------------------------------------------------

static PyObject *
model_get_lattice(ModelPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ModelPtr ptr = std::dynamic_pointer_cast<Model>(self->ptr);
    return lattice_wrap(ptr->lattice());
}

static PyObject *
model_get_a(ModelPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ModelPtr ptr = std::dynamic_pointer_cast<Model>(self->ptr);
    return wrap::to_py<double, 2>(ptr->a);
}

static PyObject *
model_get_b(ModelPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ModelPtr ptr = std::dynamic_pointer_cast<Model>(self->ptr);
    return wrap::to_py<double, 2>(ptr->b);
}

static PyObject *
model_get_V(ModelPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ModelPtr ptr = std::dynamic_pointer_cast<Model>(self->ptr);
    return wrap::to_py<std::complex<double>, 2>(ptr->V);
}

static PyObject *
model_get_T(ModelPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ModelPtr ptr = std::dynamic_pointer_cast<Model>(self->ptr);
    return wrap::to_py<std::complex<double>, 3>(ptr->T);
}

static PyObject *
model_get_R(ModelPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ModelPtr ptr = std::dynamic_pointer_cast<Model>(self->ptr);
    return wrap::to_py<int32_t, 2>(ptr->R);
}

static PyObject *
model_get_kdim(ModelPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ModelPtr ptr = std::dynamic_pointer_cast<Model>(self->ptr);
    return wrap::to_py(ptr->kdim());
}

static PyObject *
model_get_sites(ModelPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ModelPtr ptr = std::dynamic_pointer_cast<Model>(self->ptr);
    return wrap::to_py(ptr->sites());
}

static PyObject *
model_get_states(ModelPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ModelPtr ptr = std::dynamic_pointer_cast<Model>(self->ptr);
    return wrap::to_py(ptr->states());
}

static PyObject *
model_get_hops(ModelPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ModelPtr ptr = std::dynamic_pointer_cast<Model>(self->ptr);
    return wrap::to_py(ptr->hops());
}

static PyObject *
model_get_site_info(ModelPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ModelPtr ptr = std::dynamic_pointer_cast<Model>(self->ptr);

    PyObject *array;
    PyObject *dtype_dict;
    PyArray_Descr *dtype;

    dtype_dict = Py_BuildValue(
            "[(s,s),(s,s),(s,s),(s,s),(s,s),(s,s),(s,s)]",
            "kind","a64",
               "x","f8",
               "y","f8",
               "z","f8",
          "states","u8",
              "si","u8",
              "ei","u8"
            );

    PyArray_DescrConverter(dtype_dict, &dtype);
    Py_DECREF(dtype_dict);

    int nd = 1;
    npy_intp dims[1] = {static_cast<npy_intp>(ptr->site_info.size())};

    array = PyArray_NewFromDescr((PyTypeObject*)&PyArray_Type,
            dtype, nd, dims, NULL, (void*)ptr->site_info.data(),
            NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_WRITEABLE, NULL);

    return array;
}

static PyObject *
model_get_state_info(ModelPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    ModelPtr ptr = std::dynamic_pointer_cast<Model>(self->ptr);

    PyObject *array;
    PyObject *dtype_dict;
    PyArray_Descr *dtype;

    dtype_dict = Py_BuildValue(
            "[(s,s),(s,s),(s,s,i)]",
             "site", "u8",
            "orbit", "u8",
           "spinor","c16",2
            );

    PyArray_DescrConverter(dtype_dict, &dtype);
    Py_DECREF(dtype_dict);

    int nd = 1;
    npy_intp dims[1] = {static_cast<npy_intp>(ptr->state_info.size())};

    array = PyArray_NewFromDescr((PyTypeObject*)&PyArray_Type,
            dtype, nd, dims, NULL, (void*)ptr->state_info.data(),
            NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_WRITEABLE, NULL);

    return array;
}

//----------------------------------------------------------------------

static PyMethodDef model_methods[] = {
    {"set_lattice", (PyCFunction)model_set_lattice, METH_O,
        "Set the lattice"},
    {"uc_size", (PyCFunction)model_uc_size, METH_NOARGS,
        "Size of the unit cell"},
    {"Tk", (PyCFunction)model_Tk, METH_VARARGS,
        "Compute the hopping matrix at a specific k-point"},
    {"Hk", (PyCFunction)model_Hk, METH_VARARGS,
        "Compute the Hamiltonian matrix at a specific k-point"},
    {"dH_dk", (PyCFunction)model_dH_dk, METH_VARARGS,
        "Compute derivatives of the Hamiltonian at a specific k-point"},
    {NULL}
};

static PyGetSetDef model_getseters[] = {
    {"a", (getter)model_get_a, NULL,
        "Matrix with lattice vectors as rows", NULL},
    {"b", (getter)model_get_b, NULL,
        "Matrix with reciprocal lattice vectors as rows", NULL},
    {"V", (getter)model_get_V, NULL,
        "Matrix with on-site energies V[i,j]", NULL},
    {"T", (getter)model_get_T, NULL,
        "Matrix with hopping matrix T[iR,i,j]", NULL},
    {"R", (getter)model_get_R, NULL,
        "Matrix with hopping indices R[iR,i]", NULL},
    {"kdim", (getter)model_get_kdim, NULL,
        "Number of reciprocal lattice vectors", NULL},
    {"sites", (getter)model_get_sites, NULL,
        "Number of sites", NULL},
    {"states", (getter)model_get_states, NULL,
        "Number of states", NULL},
    {"hops", (getter)model_get_hops, NULL,
        "Number of hops", NULL},

    {"site_info", (getter)model_get_site_info, NULL,
        "Structured array with information about sites in the model", NULL},
    {"state_info", (getter)model_get_state_info, NULL,
        "Structured array with information about states in the model", NULL},
    {"lattice", (getter)model_get_lattice, NULL,
        "Get the lattice of the model", NULL},
    {NULL}
};

//----------------------------------------------------------------------

PyTypeObject ModelPyType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "tbpp.Model",               /* tp_name */
    sizeof(ModelPy),            /* tp_basicsize */
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
    "TBPP Model object",        /* tp_doc */
    0,                          /* tp_traverse */
    0,                          /* tp_clear */
    0,                          /* tp_richcompare */
    0,                          /* tp_weaklistoffset */
    0,                          /* tp_iter */
    0,                          /* tp_iternext */
    model_methods,              /* tp_methods */
                 0,             /* tp_members */
    model_getseters,            /* tp_getset */
    0,                          /* tp_base */
    0,                          /* tp_dict */
    0,                          /* tp_descr_get */
    0,                          /* tp_descr_set */
    0,                          /* tp_dictoffset */
    (initproc)model_init,       /* tp_init */
    0,                          /* tp_alloc */
    model_new,                  /* tp_new */
};

//----------------------------------------------------------------------

PyObject *
model_wrap(ModelPtr ptr) {
    PyObject *ret = ModelPyType.tp_new((PyTypeObject*)&ModelPyType, NULL, NULL);
    NodePy* obj = reinterpret_cast<NodePy*>(ret);
    obj->ptr = std::dynamic_pointer_cast<Node>(ptr);
    return ret;
}

ModelPtr
model_unwrap(PyObject* obj) {
    if(!ModelPy_Check(obj))
        return nullptr;
    NodePy* m = reinterpret_cast<NodePy*>(obj);
    return std::dynamic_pointer_cast<Model>(m->ptr);
}

