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
#ifndef __TBPP_PYTHON_WRAP_NUMPY_H__
#define __TBPP_PYTHON_WRAP_NUMPY_H__

/**
 * \file
 * \brief Functions for converting CPP objects to and from NumPy objects
 */

#include <Python.h>
#include <structmember.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <complex>
#include <vector>
#include <tbpp/narray.h>
#include <tbpp/matrix.h>

/** Initialize the numpy API
 *
 * This function should be called in every function that uses the NumPy C API
 * before calling any numpy C API functions.
 */
#define init_numpy() {if(PyArray_API == NULL) _import_array();}


namespace wrap {
//----------------------------------------------------------------------

template<typename T> struct npy_type;
template<> struct npy_type<bool> { static const int type = NPY_BOOL; };
template<> struct npy_type<int32_t> { static const int type = NPY_INT32; };
template<> struct npy_type<int64_t> { static const int type = NPY_INT64; };
template<> struct npy_type<uint32_t> { static const int type = NPY_UINT32; };
template<> struct npy_type<uint64_t> { static const int type = NPY_UINT64; };
template<> struct npy_type<float> { static const int type = NPY_FLOAT; };
template<> struct npy_type<double> { static const int type = NPY_DOUBLE; };
template<> struct npy_type<std::complex<float>> { static const int type = NPY_CFLOAT; };
template<> struct npy_type<std::complex<double>> { static const int type = NPY_CDOUBLE; };

//----------------------------------------------------------------------
// NArray

/** \brief Create a NumPy Array which shares the same data as an NArray
 *
 * The NumPy array and NArray share the same data such that changes to either
 * one is reflected in the other. The NArray must exist for the lifetime of
 * the NumPy array.
 *
 * \warning The returned object must be destroyed with Py_DECREF(...)
 */
template<typename T, size_t N>
inline PyObject* to_py(const tbpp::math::NArray<T,N>& value) {
    init_numpy();
    npy_intp dims[N];
    for(size_t i=0; i<N; i++)
        dims[i] = static_cast<npy_intp>(value.size(i));
    PyObject *ret = PyArray_SimpleNewFromData(N, dims, npy_type<T>::type,
            static_cast<void*>(const_cast<T*>(value.data())));
    return ret;
}

/** \brief Copy the data from a NumPy array into an NArray.
 *
 * The NArray must have the same dimension as the NumPy array and is
 * initialized with the same size. The data is copied into the NArray.
 *
 * \return zero if successful, -1 otherwise
 */
template<typename T, size_t N>
inline int to_cpp(PyObject* obj, tbpp::math::NArray<T,N>& value)  {
    init_numpy();
    PyArrayObject *array = NULL;

    array = reinterpret_cast<PyArrayObject*>(
            PyArray_FROM_OTF(obj, npy_type<T>::type, NPY_ARRAY_IN_ARRAY));
    if(array == NULL) return -1;

    int nd = PyArray_NDIM(array);
    npy_intp *dims = PyArray_DIMS(array);
    T *data = static_cast<T*>(PyArray_DATA(array));

    if(nd != N) {
        std::stringstream stream;
        stream << "Must be a " << N << "D array.";
        PyErr_SetString(PyExc_ValueError, stream.str().c_str());
        Py_DECREF(array);
        return -1;
    }

    std::array<size_t, N> array_dims;
    for(size_t i=0; i<N; i++)
        array_dims[i] = static_cast<size_t>(dims[i]);
    value.resize(array_dims);

    T* ptr = value.data();
    for(size_t i=0; i<value.elem(); i++) {
        ptr[i] = data[i];
    }

    Py_XDECREF(array);
    return 0;
}

//----------------------------------------------------------------------
// DenseMatrix

/** \brief Create a NumPy Array which shares the same data as a DenseMatrix
 *
 * The NumPy array and NArray share the same data such that changes to either
 * one is reflected in the other. The DenseMatrix must exist for the lifetime of
 * the NumPy array.
 *
 * \warning The returned object must be destroyed with Py_DECREF(...)
 */
template<typename T>
inline PyObject* to_py(const tbpp::math::DenseMatrix<T>& value) {
    init_numpy();
    npy_intp dims[2];
    dims[0] = static_cast<npy_intp>(value.rows());
    dims[1] = static_cast<npy_intp>(value.cols());

    PyObject *ret = PyArray_SimpleNewFromData(2, dims, npy_type<T>::type,
            static_cast<void*>(const_cast<T*>(value.data())));
    return ret;
}

/** \brief Copy the data from a NumPy array into a DenseMatrix.
 *
 * The NumPy array must have a dimension of 2. The data is copied into the
 * DenseMatrix.
 *
 * \return zero if successful, -1 otherwise
 */
template<typename T>
inline int to_cpp(PyObject* obj, tbpp::math::DenseMatrix<T>& value)  {
    init_numpy();
    PyArrayObject *array = NULL;

    array = reinterpret_cast<PyArrayObject*>(
            PyArray_FROM_OTF(obj, npy_type<T>::type, NPY_ARRAY_IN_ARRAY));
    if(array == NULL) return -1;

    int nd = PyArray_NDIM(array);
    npy_intp *dims = PyArray_DIMS(array);
    T *data = static_cast<T*>(PyArray_DATA(array));

    if(nd != 2) {
        PyErr_SetString(PyExc_ValueError, "Must be a 2D array.");
        Py_DECREF(array);
        return -1;
    }

    value.resize(static_cast<size_t>(dims[0]),
            static_cast<size_t>(dims[1]));

    T* ptr = value.data();
    for(size_t i=0; i<value.size(); i++) {
        ptr[i] = data[i];
    }

    Py_XDECREF(array);
    return 0;
}

//----------------------------------------------------------------------
// std::vector

/** \brief Create a NumPy Array which shares the same data as an std::vector
 *
 * The NumPy array and std::vector share the same data such that changes to either
 * one is reflected in the other. The std::vector must exist for the lifetime of
 * the NumPy array.
 *
 * \warning The returned object must be destroyed with Py_DECREF(...)
 */
template<typename T>
inline PyObject* to_py(const std::vector<T>& value) {
    init_numpy();
    npy_intp dims[1];
    dims[0] = static_cast<npy_intp>(value.size());
    PyObject *ret = PyArray_SimpleNewFromData(1, dims, npy_type<T>::type,
            static_cast<void*>(const_cast<T*>(value.data())));
    return ret;
}

/** \brief Copy the data from a NumPy array into an std::vector
 *
 * The NumPy array must have a dimension of 1.
 *
 * \return zero if successful, -1 otherwise
 */
template<typename T>
inline int to_cpp(PyObject* obj, std::vector<T>& value)  {
    init_numpy();
    PyArrayObject *array = NULL;

    array = reinterpret_cast<PyArrayObject*>(
            PyArray_FROM_OTF(obj, npy_type<T>::type, NPY_ARRAY_IN_ARRAY));
    if(array == NULL) return -1;

    int nd = PyArray_NDIM(array);
    npy_intp *dims = PyArray_DIMS(array);
    T *data = static_cast<T*>(PyArray_DATA(array));

    if(nd != 1) {
        PyErr_SetString(PyExc_ValueError, "Must be 1D array.");
        Py_DECREF(array);
        return -1;
    }

    if (value.size() != static_cast<size_t>(dims[0])) {
        value.resize(static_cast<size_t>(dims[0]));
    }

    for(npy_intp i=0; i<dims[0]; i++) {
        value[i] = data[i];
    }

    Py_XDECREF(array);
    return 0;
}

//----------------------------------------------------------------------
// std::array

/** \brief Create a NumPy Array which shares the same data as an std::array
 *
 * The NumPy array and std::array share the same data such that changes to either
 * one is reflected in the other. The std::array must exist for the lifetime of
 * the NumPy array.
 *
 * \warning The returned object must be destroyed with Py_DECREF(...)
 */
template<typename T, size_t N>
inline PyObject* to_py(const std::array<T,N>& value) {
    init_numpy();
    npy_intp dims[1];
    dims[0] = static_cast<npy_intp>(N);
    PyObject *ret = PyArray_SimpleNewFromData(1, dims, npy_type<T>::type,
            static_cast<void*>(const_cast<T*>(&value[0])));
    return ret;
}

/** \brief Copy the data from a NumPy array into an std::array
 *
 * The NumPy array must have a dimension of 1.
 *
 * \return zero if successful, -1 otherwise
 */
template<typename T, size_t N>
inline int to_cpp(PyObject* obj, std::array<T,N>& value)  {
    init_numpy();
    PyArrayObject *array = NULL;

    array = reinterpret_cast<PyArrayObject*>(
            PyArray_FROM_OTF(obj, npy_type<T>::type, NPY_ARRAY_IN_ARRAY));
    if(array == NULL) return -1;

    int nd = PyArray_NDIM(array);
    npy_intp *dims = PyArray_DIMS(array);
    T *data = static_cast<T*>(PyArray_DATA(array));

    if(nd != 1) {
        PyErr_SetString(PyExc_ValueError, "Must be 1D array.");
        Py_DECREF(array);
        return -1;
    }

    if (N != static_cast<size_t>(dims[0])) {
        PyErr_SetString(PyExc_ValueError, "Incorrect size for array");
        Py_DECREF(array);
        return -1;
    }

    for(npy_intp i=0; i<dims[0]; i++) {
        value[i] = data[i];
    }

    Py_XDECREF(array);
    return 0;
}

//----------------------------------------------------------------------

} // namespace wrap

#endif /* __TBPP_PYTHON_WRAP_NUMPY_H__ */
