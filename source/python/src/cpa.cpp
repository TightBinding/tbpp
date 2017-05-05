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

#include <python/cpa.h>
#include <python/kmodel.h>
#include <python/wrap.h>
#include <python/wrap_numpy.h>
#include <memory>

using namespace tbpp;

//----------------------------------------------------------------------

static PyObject *
cpa_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    CPAPy *self  = NULL;
    self = (CPAPy *)type->tp_alloc(type, 0);
    if(self != NULL) {
        self->ptr = nullptr;
    }

    init_numpy();

    return (PyObject *)self;
}

static int
cpa_init(CPAPy *self, PyObject *args, PyObject *kwds)
{
    self->ptr = std::dynamic_pointer_cast<Node>(std::make_shared<CPA>());
    return 0;
}

//----------------------------------------------------------------------

static PyObject *
cpa_set_kmodel(CPAPy *self, PyObject *obj) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CPAPtr ptr = std::dynamic_pointer_cast<CPA>(self->ptr);

    if(!KModelPy_Check(obj)) {
        PyErr_SetString(PyExc_RuntimeError, "KModel Object Required");
        return NULL;
    }

    NodePy* node = reinterpret_cast<NodePy*>(obj);
    ptr->set_kmodel(std::dynamic_pointer_cast<KModel>(node->ptr));

    Py_RETURN_NONE;
}

static PyObject *
cpa_set_w(CPAPy *self, PyObject *obj) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CPAPtr ptr = std::dynamic_pointer_cast<CPA>(self->ptr);

    if(wrap::to_cpp(obj, ptr->w)) {
        return NULL;
    }

    Py_RETURN_NONE;
}

PyDoc_STRVAR(cpa_add_defect_doc,
"Add a site: site_index, concentration, on-site defect potential matrix");
static PyObject *
cpa_add_defect(CPAPy *self, PyObject *args) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CPAPtr ptr = std::dynamic_pointer_cast<CPA>(self->ptr);

    std::vector<uint32_t> sites;
    std::vector<double> c;
    PyObject *onsite_obj, *c_obj, *sites_obj;
    math::DenseMatrix<cxdouble> onsite;
    if(!PyArg_ParseTuple(args, "OOO", &sites_obj, &c_obj, &onsite_obj))
        return NULL;

    if(wrap::to_cpp(onsite_obj, onsite))
        return NULL;
    if(wrap::to_cpp(c_obj, c))
        return NULL;
    if(wrap::to_cpp(sites_obj, sites))
        return NULL;

    try {
        ptr->add_defect(sites, c, onsite);
    } catch (std::exception const& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

//----------------------------------------------------------------------

static PyObject *
cpa_get_kmodel(CPAPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CPAPtr ptr = std::dynamic_pointer_cast<CPA>(self->ptr);
    return kmodel_wrap(ptr->kmodel());
}

static PyObject *
cpa_get_w(CPAPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CPAPtr ptr = std::dynamic_pointer_cast<CPA>(self->ptr);
    return wrap::to_py(ptr->w);
}

static PyObject *
cpa_get_sigma(CPAPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CPAPtr ptr = std::dynamic_pointer_cast<CPA>(self->ptr);
    return wrap::to_py<std::complex<double>, 4>(ptr->sigma);
}

static PyObject *
cpa_get_vca(CPAPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CPAPtr ptr = std::dynamic_pointer_cast<CPA>(self->ptr);
    return wrap::to_py<std::complex<double>, 3>(ptr->vca);
}

static PyObject *
cpa_get_eps(CPAPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CPAPtr ptr = std::dynamic_pointer_cast<CPA>(self->ptr);
    return wrap::to_py(ptr->eps);
}

static int
cpa_set_eps(CPAPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    CPAPtr ptr = std::dynamic_pointer_cast<CPA>(self->ptr);
    return wrap::to_cpp(value, ptr->eps);
}

static PyObject *
cpa_get_zero_j(CPAPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CPAPtr ptr = std::dynamic_pointer_cast<CPA>(self->ptr);
    return wrap::to_py(ptr->zero_j);
}

static int
cpa_set_zero_j(CPAPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    CPAPtr ptr = std::dynamic_pointer_cast<CPA>(self->ptr);
    return wrap::to_cpp(value, ptr->zero_j);
}

static PyObject *
cpa_get_max_iter(CPAPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CPAPtr ptr = std::dynamic_pointer_cast<CPA>(self->ptr);
    return wrap::to_py(ptr->max_iter);
}

static int
cpa_set_max_iter(CPAPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    CPAPtr ptr = std::dynamic_pointer_cast<CPA>(self->ptr);
    return wrap::to_cpp(value, ptr->max_iter);
}


static PyObject *
cpa_get_mix_factor(CPAPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CPAPtr ptr = std::dynamic_pointer_cast<CPA>(self->ptr);
    return wrap::to_py(ptr->mix_factor);
}

static int
cpa_set_mix_factor(CPAPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    CPAPtr ptr = std::dynamic_pointer_cast<CPA>(self->ptr);
    return wrap::to_cpp(value, ptr->mix_factor);
}

static PyObject *
cpa_get_use_vca(CPAPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    CPAPtr ptr = std::dynamic_pointer_cast<CPA>(self->ptr);
    return wrap::to_py(ptr->use_vca);
}

static int
cpa_set_use_vca(CPAPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    CPAPtr ptr = std::dynamic_pointer_cast<CPA>(self->ptr);
    return wrap::to_cpp(value, ptr->use_vca);
}

//----------------------------------------------------------------------

static PyMethodDef cpa_methods[] = {
    {"set_kmodel", (PyCFunction)cpa_set_kmodel, METH_O,
        "Set the KModel"},
    {"set_w", (PyCFunction)cpa_set_w, METH_O,
        "Set the frequencies to calculate"},
    {"add_defect", (PyCFunction)cpa_add_defect, METH_VARARGS,
        "Add a defect"},
    {NULL}
};

static PyGetSetDef cpa_getseters[] = {
    {"kmodel", (getter)cpa_get_kmodel, NULL,
        "Get the kmodel of the model", NULL},
    {"w", (getter)cpa_get_w, NULL,
        "Frequencies being computed", NULL},
    {"sigma", (getter)cpa_get_sigma, NULL,
        "Self-Energy [ic,iw,i,j]", NULL},
    {"vca", (getter)cpa_get_vca, NULL,
        "Virtual Crystal Approximation [ic,i,j]", NULL},
    {"eps", (getter)cpa_get_eps, (setter)cpa_set_eps,
        "Precision to use in CPA self-consistency cycle", NULL},
    {"zero_j", (getter)cpa_get_zero_j, (setter)cpa_set_zero_j,
        "Initial broadening to use for +0j for Green's matrix calculation", NULL},
    {"max_iter", (getter)cpa_get_max_iter, (setter)cpa_set_max_iter,
        "Maximum number of iterations for CPA self-consistent cycle", NULL},
    {"mix_factor", (getter)cpa_get_mix_factor, (setter)cpa_set_mix_factor,
        "Mixing factor to use for CPA (lower values are more stable)", NULL},
    {"use_vca", (getter)cpa_get_use_vca, (setter)cpa_set_use_vca,
        "Whether to use virtual crystal approximation for self-energy", NULL},
    {NULL}
};

//----------------------------------------------------------------------

PyTypeObject CPAPyType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "tbpp.CPA",               /* tp_name */
    sizeof(CPAPy),            /* tp_basicsize */
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
    "TBPP CPA object",        /* tp_doc */
    0,                          /* tp_traverse */
    0,                          /* tp_clear */
    0,                          /* tp_richcompare */
    0,                          /* tp_weaklistoffset */
    0,                          /* tp_iter */
    0,                          /* tp_iternext */
    cpa_methods,              /* tp_methods */
                 0,             /* tp_members */
    cpa_getseters,            /* tp_getset */
    0,                          /* tp_base */
    0,                          /* tp_dict */
    0,                          /* tp_descr_get */
    0,                          /* tp_descr_set */
    0,                          /* tp_dictoffset */
    (initproc)cpa_init,       /* tp_init */
    0,                          /* tp_alloc */
    cpa_new,                  /* tp_new */
};

//----------------------------------------------------------------------

PyObject *
cpa_wrap(CPAPtr ptr) {
    PyObject *ret = CPAPyType.tp_new((PyTypeObject*)&CPAPyType, NULL, NULL);
    NodePy* obj = reinterpret_cast<NodePy*>(ret);
    obj->ptr = std::dynamic_pointer_cast<Node>(ptr);
    return ret;
}

CPAPtr
cpa_unwrap(PyObject* obj) {
    if(!CPAPy_Check(obj))
        return nullptr;
    NodePy* m = reinterpret_cast<NodePy*>(obj);
    return std::dynamic_pointer_cast<CPA>(m->ptr);
}

