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

#ifndef __TBPP_PYTHON_WRAP_H__
#define __TBPP_PYTHON_WRAP_H__

/**
 * \file
 * \brief Functions for converting CPP objects to and from NumPy objects
 */

#include <Python.h>
#include <structmember.h>

#include <complex>
#include <string>

namespace wrap {

//----------------------------------------------------------------------
// std::complex<double>

/** \brief Create a new PyComplexObject from an std::complex<double>
 *
 * The data is copied.
 *
 * \warning The returned object must be destroyed with Py_DECREF(...)
 */
inline PyObject* to_py(const std::complex<double> &value) {
    return PyComplex_FromDoubles(real(value), imag(value));
}

/** Copy data from PyComplexObject to an std::complex<double>
 *
 * \return zero if successful, -1 otherwise
 */
inline int to_cpp(PyObject* obj, std::complex<double> &value) {
    if(!PyComplex_Check(obj)) {
        PyErr_SetString(PyExc_TypeError, "A Complex type is required");
        return -1;
    }
    value = std::complex<double>(PyComplex_RealAsDouble(obj),
                                 PyComplex_ImagAsDouble(obj));
    return 0;
}

//----------------------------------------------------------------------
// double

/** \brief Create a new PyFloat from a double
 *
 * The data is copied.
 *
 * \warning The returned object must be destroyed with Py_DECREF(...)
 */
inline PyObject* to_py(double value) {
    return PyFloat_FromDouble(value);
}

/** Copy data from PyFloat to a double
 *
 * \return zero if successful, -1 otherwise
 */
inline int to_cpp(PyObject *obj, double &value) {
    double d = PyFloat_AsDouble(obj);
    if(PyErr_Occurred())
        return -1;
    value = d;
    return 0;
}

//----------------------------------------------------------------------
// int

/** \brief Create a new PyLong from an int
 *
 * The data is copied.
 *
 * \warning The returned object must be destroyed with Py_DECREF(...)
 */
inline PyObject* to_py(int value) {
    return PyLong_FromLong(static_cast<long>(value));
}

/** Copy data from PyLong to an int
 *
 * \return zero if successful, -1 otherwise
 */
inline int to_cpp(PyObject *obj, int &value) {
    long l = PyLong_AsLong(obj);
    if(PyErr_Occurred())
        return -1;
    value = static_cast<int>(l);
    return 0;
}

//----------------------------------------------------------------------
// long

/** \brief Create a new PyLong from a long
 *
 * The data is copied.
 *
 * \warning The returned object must be destroyed with Py_DECREF(...)
 */
inline PyObject* to_py(long value) {
    return PyLong_FromLong(value);
}

/** Copy data from PyLong to a long
 *
 * \return zero if successful, -1 otherwise
 */
inline int to_cpp(PyObject *obj, long &value) {
    long l = PyLong_AsLong(obj);
    if(PyErr_Occurred())
        return -1;
    value = l;
    return 0;
}

//----------------------------------------------------------------------
// unsigned

/** \brief Create a new PyLong from an unsigned
 *
 * The data is copied.
 *
 * \warning The returned object must be destroyed with Py_DECREF(...)
 */
inline PyObject* to_py(unsigned value) {
    return PyLong_FromUnsignedLong(static_cast<unsigned long>(value));
}

/** Copy data from PyLong to an unsigned
 *
 * \return zero if successful, -1 otherwise
 */
inline int to_cpp(PyObject *obj, unsigned &value) {
    unsigned long l = PyLong_AsUnsignedLong(obj);
    if(PyErr_Occurred())
        return -1;
    value = static_cast<unsigned>(l);
    return 0;
}

//----------------------------------------------------------------------
// size_t

/** \brief Create a new PyLong from a size_t
 *
 * The data is copied.
 *
 * \warning The returned object must be destroyed with Py_DECREF(...)
 */
inline PyObject* to_py(size_t value) {
    return PyLong_FromSize_t(value);
}

/** Copy data from PyLong to a size_t
 *
 * \return zero if successful, -1 otherwise
 */
inline int to_cpp(PyObject *obj, size_t &value) {
    size_t l = PyLong_AsSize_t(obj);
    if(PyErr_Occurred())
        return 0;
    value = l;
    return -1;
}

//----------------------------------------------------------------------
// bool

/** \brief Create a new PyBool from a bool
 *
 * The data is copied.
 *
 * \warning The returned object must be destroyed with Py_DECREF(...)
 */
inline PyObject* to_py(bool value) {
    if(value) Py_RETURN_TRUE;
    Py_RETURN_FALSE;
}

/** Copy data from PyBool to a bool
 *
 * \return zero if successful, -1 otherwise
 */
inline int to_cpp(PyObject *obj, bool &value) {
    // FIXME: Use PyObject_IsTrue(obj) instead?
    if(!PyBool_Check(obj)) {
        PyErr_SetString(PyExc_TypeError, "A boolean is required");
        return -1;
    }
    if (obj == Py_True) {
        value = true;
    } else if(obj == Py_False) {
        value = false;
    } else {
        PyErr_SetString(PyExc_TypeError, "Unknown boolean");
        return -1;
    }
    return 0;
}

//----------------------------------------------------------------------
// std::string

/** \brief Create a new Python String from an std::string
 *
 * The data is copied.
 *
 * \warning The returned object must be destroyed with Py_DECREF(...)
 */
inline PyObject* to_py(const std::string &value) {
    return Py_BuildValue("s", value.c_str());
}

/** Copy data from Python String to an std::string
 *
 * \return zero if successful, -1 otherwise
 */
inline int to_cpp(PyObject *obj, std::string &value) {
    PyObject* str = PyUnicode_AsASCIIString(obj);
    if(str == NULL)
        return -1;
    const char* data = PyBytes_AsString(str);
    if(PyErr_Occurred() or data == NULL)
        return -1;
    value.assign(data);
    Py_DECREF(str);
    return 0;
}

//----------------------------------------------------------------------
// char*

/** \brief Create new Python String from const char*
 *
 * The data is copied.
 *
 * \warning The returned object must be destroyed with Py_DECREF(...)
 */
inline PyObject* to_py(const char* value) {
    return Py_BuildValue("s", value);
}

//----------------------------------------------------------------------


} // namespace wrap

#endif /* __TBPP_PYTHON_WRAP_H__ */
