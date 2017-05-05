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

#include <python/kgrid.h>
#include <python/kmodel.h>
#include <python/model.h>
#include <python/wrap.h>
#include <python/wrap_numpy.h>
#include <memory>

using namespace tbpp;

//----------------------------------------------------------------------

static PyObject *
kgrid_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    KGridPy *self  = NULL;
    self = (KGridPy *)type->tp_alloc(type, 0);
    if(self != NULL) {
        self->ptr = nullptr;
    }

    init_numpy();

    return (PyObject *)self;
}

static int
kgrid_init(KGridPy *self, PyObject *args, PyObject *kwds)
{
    self->ptr = std::dynamic_pointer_cast<Node>(std::make_shared<KGrid>());
    return 0;
}

//----------------------------------------------------------------------

static PyObject *
kgrid_set_kmodel(KGridPy *self, PyObject *obj) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KGridPtr ptr = std::dynamic_pointer_cast<KGrid>(self->ptr);

    if(!KModelPy_Check(obj)) {
        PyErr_SetString(PyExc_RuntimeError, "KModel Object Required");
        return NULL;
    }

    NodePy* node = reinterpret_cast<NodePy*>(obj);
    ptr->set_kmodel(std::dynamic_pointer_cast<KModel>(node->ptr));

    Py_RETURN_NONE;
}

//----------------------------------------------------------------------

PyDoc_STRVAR(kgrid_get_eigval_doc,
"Get the eigenvalues");
static PyObject *
kgrid_get_eigval(KGridPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KGridPtr ptr = std::dynamic_pointer_cast<KGrid>(self->ptr);
    return wrap::to_py(ptr->eigval);
}

PyDoc_STRVAR(kgrid_get_eigvec_doc,
"Get the eigenvectors");
static PyObject *
kgrid_get_eigvec(KGridPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KGridPtr ptr = std::dynamic_pointer_cast<KGrid>(self->ptr);
    return wrap::to_py(ptr->eigvec);
}

static PyObject *
kgrid_get_kmodel(KGridPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KGridPtr ptr = std::dynamic_pointer_cast<KGrid>(self->ptr);
    return kmodel_wrap(ptr->kmodel());
}

static PyObject *
kgrid_get_solve_eigvec(KGridPy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    KGridPtr ptr = std::dynamic_pointer_cast<KGrid>(self->ptr);
    return wrap::to_py(ptr->solve_eigvec);
}

static int
kgrid_set_solve_eigvec(KGridPy *self, PyObject *value, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return -1;
    }
    KGridPtr ptr = std::dynamic_pointer_cast<KGrid>(self->ptr);
    return wrap::to_cpp(value, ptr->solve_eigvec);
}

//----------------------------------------------------------------------

static PyMethodDef kgrid_methods[] = {
    {"set_kmodel", (PyCFunction)kgrid_set_kmodel, METH_O, "Set the kmodel"},
    {NULL}
};

static PyGetSetDef kgrid_getseters[] = {
    {"eigval", (getter)kgrid_get_eigval, NULL,
        kgrid_get_eigval_doc, NULL},
    {"eigvec", (getter)kgrid_get_eigvec, NULL,
        kgrid_get_eigvec_doc, NULL},
    {"kmodel", (getter)kgrid_get_kmodel, NULL,
        "Get the KModel being used", NULL},
    {"solve_eigvec", (getter)kgrid_get_solve_eigvec, (setter)kgrid_set_solve_eigvec,
        "Whether to solve for eigenvectors", NULL},
    {NULL}
};

//----------------------------------------------------------------------

PyTypeObject KGridPyType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "tbpp.KGrid",               /* tp_name */
    sizeof(KGridPy),            /* tp_basicsize */
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
    "TBPP KGrid object",        /* tp_doc */
    0,                          /* tp_traverse */
    0,                          /* tp_clear */
    0,                          /* tp_richcompare */
    0,                          /* tp_weaklistoffset */
    0,                          /* tp_iter */
    0,                          /* tp_iternext */
    kgrid_methods,              /* tp_methods */
                 0,             /* tp_members */
    kgrid_getseters,            /* tp_getset */
    0,                          /* tp_base */
    0,                          /* tp_dict */
    0,                          /* tp_descr_get */
    0,                          /* tp_descr_set */
    0,                          /* tp_dictoffset */
    (initproc)kgrid_init,       /* tp_init */
    0,                          /* tp_alloc */
    kgrid_new,                  /* tp_new */
};

//----------------------------------------------------------------------

PyObject *
kgrid_wrap(KGridPtr ptr) {
    PyObject *ret = KGridPyType.tp_new((PyTypeObject*)&KGridPyType, NULL, NULL);
    NodePy* obj = reinterpret_cast<NodePy*>(ret);
    obj->ptr = std::dynamic_pointer_cast<Node>(ptr);
    return ret;
}

KGridPtr
kgrid_unwrap(PyObject* obj) {
    if(!KGridPy_Check(obj))
        return nullptr;
    NodePy* m = reinterpret_cast<NodePy*>(obj);
    return std::dynamic_pointer_cast<KGrid>(m->ptr);
}

