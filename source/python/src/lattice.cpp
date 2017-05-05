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

#include <python/lattice.h>
#include <python/wrap_numpy.h>
#include <python/wrap.h>

using namespace tbpp;

//----------------------------------------------------------------------

static PyObject *
lattice_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    LatticePy *self  = NULL;
    self = (LatticePy *)type->tp_alloc(type, 0);
    if(self != NULL) {
        self->ptr = nullptr;
    }

    init_numpy();

    return (PyObject *)self;
}

static int
lattice_init(LatticePy *self, PyObject *args, PyObject *kwds)
{
    self->ptr = std::dynamic_pointer_cast<Node>(std::make_shared<Lattice>());
    return 0;
}

//----------------------------------------------------------------------

PyDoc_STRVAR(lattice_add_trans_doc,
"add_trans(a)\n"
"\n Add a translational symmetry for the lattice. You can only add up to three"
"\n translational symmetries note that the order in which the symmetries are"
"\n added is important. The first added symmetry is the a1 direction, the second"
"\n added symmetry is a2, and the third added is a3."
"\n"
"\n Parameters"
"\n ----------"
"\n a : array_like of length 3"
"\n     Translational symmetry of lattice ([ax,ay,az]).");
static PyObject *
lattice_add_trans(LatticePy *self, PyObject *obj) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    LatticePtr ptr = std::dynamic_pointer_cast<Lattice>(self->ptr);

    std::vector<double> a;
    if(wrap::to_cpp(obj, a)) return NULL;
    if(a.size() != 3) {
        PyErr_SetString(PyExc_ValueError, "Requires vector of length 3");
        return NULL;
    }
    ptr->add_trans(a[0],a[1],a[2]);
    Py_RETURN_NONE;
}

PyDoc_STRVAR(lattice_add_site_doc,
"add_site(kind, r, V, orbit_indices, spinors)\n"
"\n Add a site to the lattice."
"\n"
"\n Parameters"
"\n ----------"
"\n kind : str with a length less that 64"
"\n     A string identifying the kind of site (used to determine hopping matrix)."
"\n r : array_like of length 3"
"\n     Position of the site in real space ([rx,ry,rz])"
"\n V : array_like of dimension 2 and shape (states, states)"
"\n     On-site potential matrix."
"\n orbit_indices : array_like of dimension 1 and length states."
"\n                 For each state must provide an integer identifying which"
"\n                 orbit the state belongs to."
"\n spinors : array_like of dimension 2 and shape (states, 2)"
"\n           For each state must provide a vector containing two complex"
"\n           numbers with the spinor of the state");
static PyObject *
lattice_add_site(LatticePy *self, PyObject *args) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    LatticePtr ptr = std::dynamic_pointer_cast<Lattice>(self->ptr);

    char* kind;
    double x,y,z;
    PyObject *onsite_obj, *orbits_obj, *spinors_obj;
    math::NArray<cxdouble,2> onsite;
    std::vector<uint32_t> orbits;
    math::NArray<std::complex<double>,2> spinors;

    if(!PyArg_ParseTuple(args, "s(ddd)OOO", &kind, &x, &y, &z,
                &onsite_obj, &orbits_obj, &spinors_obj))
        return NULL;

    if(wrap::to_cpp(onsite_obj, onsite))
        return NULL;
    if(wrap::to_cpp(orbits_obj, orbits))
        return NULL;
    if(wrap::to_cpp(spinors_obj, spinors))
        return NULL;

    // TODO Various consistency checks
    if(spinors.size(1) != 2) {
        PyErr_SetString(PyExc_ValueError, "Invalid spinor format");
        return NULL;
    }
    if(spinors.size(0) != orbits.size()) {
        PyErr_SetString(PyExc_ValueError, "Invalid number of spinors");
        return NULL;
    }

    try {
        ptr->add_site(std::string(kind), x, y, z, orbits.size(),
                onsite.data(), orbits.data(), spinors.data());
    } catch (std::exception const& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

PyDoc_STRVAR(lattice_add_hop_doc,
"add_hop(initial_kind, final_kind, dr, T)\n"
"\n Add a real space hopping between two sites separated by a distance dr."
"\n"
"\n Parameters"
"\n ----------"
"\n initial_kind : str with a length less than 64"
"\n     String identifying type of the initial site for the hopping."
"\n final_kind : str with a length less than 64"
"\n     String identifying type of the final site for the hopping."
"\n dr : array like of length 3"
"\n     The (final position) minus (initial position) in real space for the hopping."
"\n T : array_like of shape (initial_states)x(final_states)"
"\n     Hopping matrix.");
static PyObject *
lattice_add_hop(LatticePy *self, PyObject *args) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    LatticePtr ptr = std::dynamic_pointer_cast<Lattice>(self->ptr);

    char *kind1, *kind2;
    double dx,dy,dz;
    PyObject *T_obj;
    math::NArray<std::complex<double>,2> T;

    if(!PyArg_ParseTuple(args, "ss(ddd)O", &kind1, &kind2,
                &dx, &dy, &dz, &T_obj))
        return NULL;

    if(wrap::to_cpp(T_obj, T))
        return NULL;

    try {
        ptr->add_hop(kind1, kind2, dx, dy, dz, T.size(0), T.size(1), T.data());
    } catch (const std::exception& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }
    Py_RETURN_NONE;
}

PyDoc_STRVAR(lattice_copy_along_doc,
"copy_along(v, n)\n"
"\n Create n copies of the lattice along the v direction."
"\n"
"\n Parameters"
"\n ----------"
"\n v : array_like of length 3"
"\n     The direction in which to copy the lattice."
"\n n : int"
"\n     The number of copies to make.");
static PyObject *
lattice_copy_along(LatticePy *self, PyObject *args) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    LatticePtr ptr = std::dynamic_pointer_cast<Lattice>(self->ptr);

    double vx,vy,vz;
    unsigned int n;
    if(!PyArg_ParseTuple(args, "(ddd)I", &vx, &vy, &vz, &n))
        return NULL;

    try {
        ptr->copy_along(vx,vy,vz, n);
    } catch (const std::exception& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }
    Py_RETURN_NONE;
}

PyDoc_STRVAR(lattice_move_along_doc,
"move_along(v)\n"
"\n Move every site in the lattice by the vector v."
"\n"
"\n Parameters"
"\n ----------"
"\n v : array_like of length 3"
"\n     The direction in which to move the lattice.");
static PyObject *
lattice_move_along(LatticePy *self, PyObject *args) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    LatticePtr ptr = std::dynamic_pointer_cast<Lattice>(self->ptr);

    double vx,vy,vz;
    if(!PyArg_ParseTuple(args, "(ddd)", &vx, &vy, &vz))
        return NULL;

    ptr->move_along(vx,vy,vz);
    Py_RETURN_NONE;
}


//----------------------------------------------------------------------

static PyObject *
lattice_get_eps(LatticePy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    LatticePtr ptr = std::dynamic_pointer_cast<Lattice>(self->ptr);
    return wrap::to_py(ptr->eps);
}

static PyObject *
lattice_get_max_hop_range(LatticePy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    LatticePtr ptr = std::dynamic_pointer_cast<Lattice>(self->ptr);
    return wrap::to_py(ptr->max_hop_range);
}

static PyObject *
lattice_get_B(LatticePy *self, void *closure) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    LatticePtr ptr = std::dynamic_pointer_cast<Lattice>(self->ptr);
    return wrap::to_py(ptr->B);
}

PyDoc_STRVAR(lattice_set_B_doc,
"set_B(B)\n"
"\n Set a uniform external magnetic field"
"\n"
"\n Parameters"
"\n ----------"
"\n B : array_like of length 3"
"\n     Magnetic field strength along the x,y,z directions ([Bx,By,Bz])");
static PyObject *
lattice_set_B(LatticePy *self, PyObject *obj) {
    if(self->ptr == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Pointer is nullptr");
        return NULL;
    }
    LatticePtr ptr = std::dynamic_pointer_cast<Lattice>(self->ptr);

    if(wrap::to_cpp(obj, ptr->B)) {
        return NULL;
    }

    Py_RETURN_NONE;
}

//----------------------------------------------------------------------

static PyMethodDef lattice_methods[] = {
    {"add_trans", (PyCFunction)lattice_add_trans, METH_O, lattice_add_trans_doc},
    {"set_B", (PyCFunction)lattice_set_B, METH_O, lattice_set_B_doc},
    {"add_site", (PyCFunction)lattice_add_site, METH_VARARGS, lattice_add_site_doc},
    {"add_hop", (PyCFunction)lattice_add_hop, METH_VARARGS, lattice_add_hop_doc},
    {"copy_along", (PyCFunction)lattice_copy_along, METH_VARARGS, lattice_copy_along_doc},
    {"move_along", (PyCFunction)lattice_move_along, METH_VARARGS, lattice_move_along_doc},
    {NULL}
};

static PyGetSetDef lattice_getseters[] = {
    {"eps", (getter)lattice_get_eps, NULL,
        "Maximum uncertainty in the position of sites. Used as a tolerance when"
        " determining hopping matrix.", NULL},
    {"max_hop_range", (getter)lattice_get_max_hop_range, NULL,
        "Maximum hopping distance (automatically determined)", NULL},
    {"B", (getter)lattice_get_B, (setter)lattice_set_B,
        "External Magnetic Field", NULL},
    {NULL}
};

//----------------------------------------------------------------------

PyTypeObject LatticePyType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "tbpp.Lattice",               /* tp_name */
    sizeof(LatticePy),            /* tp_basicsize */
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
    "TBPP Lattice object",        /* tp_doc */
    0,                          /* tp_traverse */
    0,                          /* tp_clear */
    0,                          /* tp_richcompare */
    0,                          /* tp_weaklistoffset */
    0,                          /* tp_iter */
    0,                          /* tp_iternext */
    lattice_methods,              /* tp_methods */
                 0,             /* tp_members */
    lattice_getseters,            /* tp_getset */
    0,                          /* tp_base */
    0,                          /* tp_dict */
    0,                          /* tp_descr_get */
    0,                          /* tp_descr_set */
    0,                          /* tp_dictoffset */
    (initproc)lattice_init,       /* tp_init */
    0,                          /* tp_alloc */
    lattice_new,                  /* tp_new */
};

//----------------------------------------------------------------------

PyObject *
lattice_wrap(LatticePtr ptr) {
    PyObject *ret = LatticePyType.tp_new((PyTypeObject*)&LatticePyType, NULL, NULL);
    NodePy* obj = reinterpret_cast<NodePy*>(ret);
    obj->ptr = std::dynamic_pointer_cast<Node>(ptr);
    return ret;
}

LatticePtr
lattice_unwrap(PyObject* obj) {
    if(!LatticePy_Check(obj))
        return nullptr;
    NodePy* m = reinterpret_cast<NodePy*>(obj);
    return std::dynamic_pointer_cast<Lattice>(m->ptr);
}

