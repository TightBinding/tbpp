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

#ifndef __TBPP_MATH_NARRAY_H__
#define __TBPP_MATH_NARRAY_H__

/**
 * \file
 * \brief A multi-dimensional array container.
 */

#include <limits>
#include <fstream>
#include <string>
#include <exception>
#include <vector>
#include <array>
#include <numeric>
#include <cstring>

#include <tbpp/common.h>

namespace tbpp {
namespace math {

/**
 * \brief A multi-dimensional array container.
 *
 * Internally the values are stored using row-major ordering such that the
 * last index is the fastest changing index.
 *
 * Basic usage,
 * \code
 * NArray<double, 3> a; // Create a three dimensional array of doubles
 *
 * // Can also specify the size in constructor
 * //     NArray<double, 3> a(3,4,5);
 *
 * a.resize(3,4,5);     // Change size to 3 x 4 x 5 (all data is erased
 *                        // when changing the size of an array)
 *
 * a.size(2);             // Returns the size of the array along i'th dimension
 *                        // (note that the index is 0-valued). In this case
 *                        // returns 5 since {3,4,5}[2] = 5.
 *
 * a.dim();               // Returns the dimension of the array
 *                        // (returns 3 in this case)
 *
 * a.elem();              // Returns the total number of elements
 *                        // (3x4x5 = 60 in this case)
 *
 * a.fill(1);             // Fill all array elements with ones
 *
 * a(0,0,0) = 10;         // Set the 0,0,0 element to 10
 *
 * double b = a(0,0,0);   // Get the 0,0,0 element
 *
 * a.clear();           // Clear all data and set the dimension and size to 0
 *                      // (data is freed from memory)
 *
 * a.empty();           // Whether array is empty (in this case it returns true
 *                      // since we just ran `clear()`)
 *
 * double* prt;
 * prt = a.data(); // Returns the raw pointer where the data is
 *                   // store (only use this if you know what you are doing).
 *                   // The values are stored in row-major ordering where the
 *                   // last index is the fastest changing index.
 *
 * // When passing NArray as a parameter it is best to use a reference
 * double my_func(const NArray<double,3>& a) {
 *     // code...
 * }
 *
 * // It is possible to get direct access (especially for the last index)
 * double* c = &a(1,2,0); // Get pointer to a(1,2,0)
 * c[0] = 11;             // Set a(1,2,0) to 11
 * c[1] = 12;             // Set a(1,2,1) to 12
 *
 * \endcode
 */
template<typename T, size_t N>
class NArray {
private:
    /// Number of dimensions
    static constexpr size_t _dim = N;

    /// The size for each dimension
    std::array<size_t, N> _sizes;

    /// Number of items to skip to reach next item in a given dimension
    std::array<size_t, N> _strides;

    /// The data
    std::vector<T> _data;

public:

    /** \brief Create a new NArray
     *
     * The resulting array has a size of zero along each dimension and hence
     * has no elements.
     */
    NArray() = default;

    /// Create a new NArray and specify the size of the array.
    template<typename... Size>
    explicit NArray(Size... size) {
        resize(size_t(size)...);
    }

    /// move operator
    NArray(NArray&&) = default;
    /// move operator
    NArray& operator=(NArray&&) = default;
    /// copy operator
    NArray(NArray const&) = default;
    /// copy operator
    NArray& operator=(NArray const&) = default;
    /// destructor
    ~NArray() = default;

    /** Returns the dimension of the array
     *
     * The dimension is the number of indexes that values in the array are
     * accessed by. It is equal to the rank of the resulting tensor. For
     * instance for a vector `dim=1`, for a matrix `dim=2`, etc...
     */
    static constexpr size_t dim() { return _dim; }

    /// Set the size of the array
    template<typename... Size>
    void resize(Size... sizes) {
        static_assert(sizeof...(sizes) == _dim, "Invalid dimension for size");

        // copy new size
        std::array<size_t,_dim> new_sizes{size_t(sizes)...};
        resize(new_sizes);
    }

    /// Set the size of the array using an std::array
    void resize(std::array<size_t,_dim> new_sizes) {
        // copy new size
        _sizes = new_sizes;

        // Strides
        // Example for _dim=5, let n = _sizes
        //     _strides[0] = n[1]*n[2]*n[3]*n[4];
        //     _strides[1] = n[2]*n[3]*n[4];
        //     _strides[2] = n[3]*n[4];
        //     _strides[3] = n[4];
        //     _strides[4] = 1;
        //
        // Loop below generalizes above for arbitrary dimension using:
        //   _strides[i] = _size[data_dim-1]*_size[data_dim-2]..._size[i+1]
        for(size_t i=0; i < _dim; i++) {
            _strides[i] = 1;
            for(size_t j=_dim-1; j > i; j--) {
                _strides[i] *= _sizes[j];
            }
        }

        // Total number of elements
        size_t elem = 1;
        for(size_t i=0; i<_dim; i++)
            elem *= _sizes[i];

        // Allocate new data
        _data.resize(elem);
    }


    /** \brief Returns the size along the i'th dimension
     * \note This is zero-indexed, so first dimension corresponds to \a i=0
     * \throw out_of_range If i >= number of dimensions of array
     */
    size_t size(size_t i) const {
#ifdef TBPP_NARRAY_RANGE_CHECK
        return _sizes.at(i);
#else
        return _sizes[i];
#endif
    }

    /// Returns the total number of elements in the array
    size_t elem() const { return _data.size(); }

    /// Returns true if the array is empty
    bool empty() const { return _data.empty(); }

    /// Clears the data in the array and sets the size to zero
    void clear() {
        _sizes.fill(0);
        _strides.fill(0);
        _data.clear();
    }

    /// Fill all of the entries in the array with value
    void fill(const T& value) {
        for(size_t i=0; i < _data.size(); i++)
            _data[i] = value;
    }

    /** \brief Returns the raw pointer where the data is stored
     *
     * The values are stored in row-major ordering such that the last index is
     * the fastest changing index.
     *
     * \warning It is best to use the operator () for getting and setting
     * values.
     */
    T* data() { return _data.data(); }

    /** \brief Returns the raw pointer where the data is stored
     *
     * The values are stored in row-major ordering such that the last index is
     * the fastest changing index.
     *
     * \warning It is best to use the operator () for getting and setting
     * values.
     */
    const T* data() const { return _data.data(); }

    /// Set the value of an element
    template<typename... Index>
    inline T& operator()(Index... index) {
        static_assert(sizeof...(index) == _dim, "Invalid dimension for index");

        // copy into array
        size_t a_index[_dim]{size_t(index)...};

#ifdef TBPP_NARRAY_RANGE_CHECK
        for(size_t i=0; i<_dim; i++)
            if(a_index[i] >= _sizes[i])
                throw std::out_of_range("Index out of range");
#endif

        size_t raw_index = std::inner_product(a_index, a_index+N,
                _strides.begin(), size_t(0));
        return _data[raw_index];
    }

    /// Get a value in the array
    template<typename... Index>
    inline const T& operator()(Index... index) const {
        static_assert(sizeof...(index) == _dim, "Invalid dimension for index");

        // copy into array
        size_t a_index[_dim]{size_t(index)...};

#ifdef TBPP_NARRAY_RANGE_CHECK
        for(size_t i=0; i<_dim; i++)
            if(a_index[i] >= _sizes[i])
                throw std::out_of_range("Index out of range");
#endif

        size_t raw_index = std::inner_product(a_index, a_index+_dim,
                _strides.begin(), size_t(0));
        return _data[raw_index];
    }

    /// Save array to a file in ASCII dense NArray format
    void save_dense(const std::string& filename) const {
        std::ofstream out(filename, std::ios::out | std::ios::trunc);
        if(!out) throw std::runtime_error("could not open file for writing");

        // Set output precision
        if(std::numeric_limits<T>::is_specialized) {
            out.precision(std::numeric_limits<T>::max_digits10);
        } else {
            out.precision(std::numeric_limits<double>::max_digits10);
        }
        // Set output format
        out.setf(std::ios::scientific);

        // header
        out << "NArray dense\n"
            << "dim " << _dim << '\n';
        out << "sizes";
        for(size_t i=0; i<_dim; i++) {
            out << ' ' << _sizes[i];
        }
        out << '\n';

        // data
        for(size_t i=0; i< elem(); i++) {
            out << _data[i] << '\n';
            if(!out.good()) throw std::runtime_error("error writing to file");
        }
    }

    /// Load array from a file in ASCII dense NArray format
    void load_dense(const std::string& filename) {
        std::ifstream in(filename);
        if(!in) throw std::runtime_error("could not open file for reading");

        std::string buf;
        if(!(in >> buf) or buf != "NArray")
            throw std::runtime_error("[format error]: missing NArray header");
        if(!(in >> buf) or buf != "dense")
            throw std::runtime_error("[format error]: unrecognized data format");
        if(!(in >> buf) or buf != "dim")
            throw std::runtime_error("[format error]: missing dimension information");

        // get dimension
        size_t new_dim;
        if(!(in >> new_dim))
            throw std::runtime_error("[format error]: problem reading array dimension");
        if(new_dim != _dim)
            throw std::runtime_error("[error]: NArray has an incompatible dimension");

        if(!(in >> buf) or buf != "sizes")
            throw std::runtime_error("[format error]: missing sizes information");

        // get size
        std::array<size_t, _dim> new_size;
        for(size_t i=0; i<_dim; i++) {
            if(!(in >> new_size[i]))
                throw std::runtime_error("[format error]: problem reading size");
        }
        resize(new_size);

        // get data
        for(size_t i=0; i<elem(); i++) {
            if(!(in >> _data[i]))
                throw std::runtime_error("[format error]: problem reading data");
        }
    }

    /// Save array to a file in ASCII sparse NArray format
    void save_sparse(const std::string& filename) const {
        std::ofstream out(filename, std::ios::out | std::ios::trunc);
        if(!out) throw std::runtime_error("could not open file for writing");

        // Set output precision
        if(std::numeric_limits<T>::is_specialized) {
            out.precision(std::numeric_limits<T>::max_digits10);
        } else {
            out.precision(std::numeric_limits<double>::max_digits10);
        }
        // Set output format
        out.setf(std::ios::scientific);

        // header
        out << "NArray sparse\n"
            << "dim " << _dim << '\n';
        out << "sizes";
        for(size_t i=0; i<_dim; i++) {
            out << ' ' << _sizes[i];
        }
        out << '\n';

        // data
        for(size_t i=0; i< elem(); i++) {
            if(_data[i] != 0) {
                size_t index = i;
                for(size_t j=0; j<_dim; j++) {
                    size_t id = index/_strides[j];
                    index -= id*_strides[j];
                    out << id << ' ';
                }
                out << _data[i] << '\n';
            }
            if(!out.good()) throw std::runtime_error("error writing to file");
        }
    }

    /// Load array from a file in ASCII sparse NArray format
    void load_sparse(const std::string& filename) {
        std::ifstream in(filename);
        if(!in) throw std::runtime_error("could not open file for reading");

        std::string buf;
        if(!(in >> buf) or buf != "NArray")
            throw std::runtime_error("[format error]: missing NArray header");
        if(!(in >> buf) or buf != "sparse")
            throw std::runtime_error("[format error]: unrecognized data format");
        if(!(in >> buf) or buf != "dim")
            throw std::runtime_error("[format error]: missing dimension information");

        // get dimension
        size_t new_dim;
        if(!(in >> new_dim))
            throw std::runtime_error("[format error]: problem reading array dimension");
        if(new_dim != _dim)
            throw std::runtime_error("[error]: NArray has an incompatible dimension");

        if(!(in >> buf) or buf != "sizes")
            throw std::runtime_error("[format error]: missing sizes information");

        // get size
        std::array<size_t, _dim> new_size;
        for(size_t i=0; i<_dim; i++) {
            if(!(in >> new_size[i]))
                throw std::runtime_error("[format error]: problem reading size");
        }
        resize(new_size);
        fill(0);

        // get data
        while(!in.eof()) {
            size_t index = 0;
            // read index
            for(size_t j=0; j<_dim; j++) {
                size_t id;
                if(!(in >> id))
                    throw std::runtime_error("[format error]: problem reading index");
                index += id*_strides[j];
            }
            T d;
            if(!(in >> d))
                throw std::runtime_error("[format error]: problem reading data");

            _data[index] = d;
        }
    }

    /** \brief Quickly save the contents to the disk in binary format
     *
     * This function is designed to quickly save the array to the disk. The
     * array is saved in native endianness and no attempt is made to convert
     * endianness.
     *
     * \warning Use \ref save_ascii if the file is intended to be shared
     * between computers.
     */
    void save_dump(const std::string& filename) const {
        std::ofstream out(filename, std::ios::trunc | std::ios::binary);
        if(!out) throw std::runtime_error("could not open file for writing");

        // header
        const char header[] = {'N','A','r','r','a','y',' ','d','u','m','p','\0'};
        out.write(header, 12);

        // dimension/sizes/strides
        size_t dim = _dim;
        out.write(reinterpret_cast<const char*>(&dim), sizeof(dim));
        out.write(reinterpret_cast<const char*>(&_sizes[0]), _dim*sizeof(size_t));
        out.write(reinterpret_cast<const char*>(&_strides[0]), _dim*sizeof(size_t));

        // data
        size_t elem = _data.size();
        out.write(reinterpret_cast<const char*>(&elem), sizeof(elem));
        out.write(reinterpret_cast<const char*>(&_data[0]), elem*sizeof(T));

        if(!out.good()) throw std::runtime_error("error writing to file");
    }

    /** \brief Retrieve the contents of a previously dumped array
     *
     * This function is designed to load the contents saved by \ref save_dump.
     *
     * \warning The file is assumed to be in native endianness.
     */
    void load_dump(const std::string& filename) {
        std::ifstream in(filename, std::ios::binary);
        if(!in) throw std::runtime_error("could not open file for reading");

        // header
        const char header[] = {'N','A','r','r','a','y',' ','d','u','m','p','\0'};
        char in_header[12];
        in.read(in_header, 12);
        if(std::strncmp(header, in_header, 12))
            throw std::runtime_error("file does not contain an NArray dump");

        // dimension
        size_t in_dim;
        in.read(reinterpret_cast<char*>(&in_dim), sizeof(in_dim));
        if(in_dim != _dim) throw std::runtime_error("incorrect dimension");

        // sizes/strides
        in.read(reinterpret_cast<char*>(&_sizes[0]), _dim*sizeof(size_t));
        in.read(reinterpret_cast<char*>(&_strides[0]), _dim*sizeof(size_t));

        // elem
        size_t elem;
        in.read(reinterpret_cast<char*>(&elem), sizeof(elem));

        // data
        _data.resize(elem);
        in.read(reinterpret_cast<char*>(&_data[0]), elem*sizeof(T));

        if(!in.good()) throw std::runtime_error("error reading file");
    }

};

}} // namespace tbpp::math

#endif /* __TBPP_MATH_NARRAY_H__ */
