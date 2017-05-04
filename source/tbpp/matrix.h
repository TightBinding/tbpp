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

#ifndef __TBPP_MATH_MATRIX_H__
#define __TBPP_MATH_MATRIX_H__

/**
 * \file
 * \brief Defines matrices with different storage types
 */

#include <vector>
#include <string>
#include <limits>
#include <fstream>
#include <limits>
#include <exception>
#include <cstring>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace tbpp {
namespace math {

//----------------------------------------------------------------------------

template<typename T>
using DenseMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template<typename T>
using DenseMap = Eigen::Map<DenseMatrix<T>>;

//----------------------------------------------------------------------------

/** \brief A Matrix using the band storage format
 *
 * This class provides a more compact method to store matrices which only
 * contain non-zero elements along the diagonal bands.
 */
template<typename T>
class BandMatrix {
private:
    size_t _rows=0;  ///< Number of rows
    size_t _cols=0;  ///< Number of columns
    size_t _band=0;  ///< Number of bands
    size_t _ldab=0;  ///< 2*_band + 1
    std::vector<T> _data; ///< Vector storing data
    T _fill=0; ///< Fill value

public:
    using value_type = T;

    /** \brief Create a new BandMatrix
     *
     * The resulting matrix has a size of zero and no elements.
     */
    BandMatrix() = default;

    /// Create a new BandMatrix and specify the size and band
    explicit BandMatrix(size_t rows, size_t cols, size_t band) {
        resize(rows, cols, band);
    }

    /// move operator
    BandMatrix(BandMatrix&&) = default;
    /// move operator
    BandMatrix& operator=(BandMatrix&&) = default;
    /// copy operator
    BandMatrix(BandMatrix const&) = default;
    /// copy operator
    BandMatrix& operator=(BandMatrix const&) = default;
    // destructor
    ~BandMatrix() = default;

    /// Total number of rows
    size_t rows() const { return _rows; };

    /// Total number of cols
    size_t cols() const { return _cols; };

    /// The band of the matrix
    size_t band() const { return _band; };

    /// Whether the matrix has a size of zero
    bool empty() const { return _data.empty(); };

    /// Clear the matrix and set the size to zero
    void clear() {
        _rows = 0;
        _cols = 0;
        _band = 0;
        _ldab = 0;
        _data.clear();
    }

    /// Fill all entries with value
    void fill(const T& value) {
        for(size_t i=0; i < _data.size(); i++)
            _data[i] = value;
    }

    /// Resize the matrix
    void resize(size_t rows, size_t cols, size_t band) {
        _rows = rows;
        _cols = cols;
        _band = band;
        _ldab = 2*_band+1;
        _data.resize(_ldab*_cols);
    }

    /** \brief Returns the total number of elements
     * \warning This is not equal to rows*cols due to compact storage
     */
    inline size_t elem() const { return _data.size(); }

    /// Returns the raw pointer where the data is stored
    inline T* data() { return _data.data(); }

    /// Returns the raw pointer where the data is stored
    inline const T* data() const { return _data.data(); };

    /// Set the value of the matrix at row \a r and column \a c
    inline T& operator()(size_t r, size_t c) {
        if(std::fabs(r-c) > _band) return _fill;
        return _data[(_band + r - c) + c*_ldab];
    }

    /// Get the value of the matrix at row \a r and column \a c
    inline const T& operator()(size_t r, size_t c) const {
        if(std::fabs(r-c) > _band) return _fill;
        return _data[(_band + r - c) + c*_ldab];
    }

    /// Save the Matrix in ASCII format
    void save_ascii(const std::string& filename) const {
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
        out << "BandMatrix\n"
            << "rows " << _rows << '\n'
            << "cols " << _cols << '\n'
            << "band " << _band << '\n';

        // data
        for(size_t i=0; i< elem(); i++) {
            out << _data[i] << '\n';
            if(!out.good()) throw std::runtime_error("error writing to file");
        }
    }

    /// Load the Matrix in ASCII format
    void load_ascii(const std::string& filename) {
        std::ifstream in(filename);
        if(!in) throw std::runtime_error("could not open file for reading");

        std::string buf;
        if(!(in >> buf) or buf != "BandMatrix")
            throw std::runtime_error("[format error]: missing BandMatrix header");

        // rows
        size_t new_rows;
        if(!(in >> buf) or buf != "rows")
            throw std::runtime_error("[format error]: missing rows information");
        if(!(in >> new_rows))
            throw std::runtime_error("[format error]: problem reading number of rows");

        // cols
        size_t new_cols;
        if(!(in >> buf) or buf != "cols")
            throw std::runtime_error("[format error]: missing cols information");
        if(!(in >> new_cols))
            throw std::runtime_error("[format error]: problem reading number of cols");

        // band
        size_t new_band;
        if(!(in >> buf) or buf != "band")
            throw std::runtime_error("[format error]: missing band information");
        if(!(in >> new_band))
            throw std::runtime_error("[format error]: problem reading number of band");

        resize(new_rows, new_cols, new_band);

        // get data
        for(size_t i=0; i<elem(); i++) {
            if(!(in >> _data[i]))
                throw std::runtime_error("[format error]: problem reading data");
        }
    }

    /** \brief Quickly save the contents to the disk in binary format
     *
     * This function is designed to quickly save the matrix to the disk. The
     * matrix is saved in native endianness and no attempt is made to convert
     * endianness.
     *
     * \warning Use \ref save_ascii if the file is intended to be shared
     * between computers.
     */
    void save_dump(const std::string& filename) const {
        std::ofstream out(filename, std::ios::trunc | std::ios::binary);
        if(!out) throw std::runtime_error("could not open file for writing");

        // header
        constexpr char header[] = "BandMatrix dump";
        constexpr size_t header_size = sizeof(header);
        out.write(header, header_size);

        // rows and cols
        out.write(reinterpret_cast<const char*>(&_rows), sizeof(_rows));
        out.write(reinterpret_cast<const char*>(&_cols), sizeof(_cols));
        out.write(reinterpret_cast<char*>(&_band), sizeof(_band));
        out.write(reinterpret_cast<char*>(&_ldab), sizeof(_ldab));

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
        constexpr char header[] = "BandMatrix dump";
        constexpr size_t header_size = sizeof(header);
        char in_header[header_size];
        in.read(in_header, header_size);
        if(std::strncmp(header, in_header, header_size))
            throw std::runtime_error("file does not contain an BandMatrix dump");

        // rows/cols/band/ldab
        in.read(reinterpret_cast<char*>(&_rows), sizeof(_rows));
        in.read(reinterpret_cast<char*>(&_cols), sizeof(_cols));
        in.read(reinterpret_cast<char*>(&_band), sizeof(_band));
        in.read(reinterpret_cast<char*>(&_ldab), sizeof(_ldab));

        // elem
        size_t elem;
        in.read(reinterpret_cast<char*>(&elem), sizeof(elem));

        // data
        _data.resize(elem);
        in.read(reinterpret_cast<char*>(&_data[0]), elem*sizeof(T));

        if(!in.good()) throw std::runtime_error("error reading file");
    }
};

//----------------------------------------------------------------------------

}} // namespace tbpp::math

#endif /* __TBPP_MATH_MATRIX_H__ */
