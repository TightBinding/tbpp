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

#ifndef __TBPP_EHFILE_H__
#define __TBPP_EHFILE_H__

#include <string>
#include <vector>
#include "H5Cpp.h"
#include <tbpp/narray.h>

/** \file
 * \brief A simple interface to the HDF5 Data Format
 */

// Max length for string attributes
#define EHFILE_MAX_STR_LEN 256

namespace tbpp {

//----------------------------------------------------------------------------

template<typename T> struct h5_type { static H5::DataType type(); };

//template<> struct h5_type<bool> { static H5::DataType type() { return H5::PredType::NATIVE_HBOOL; } };
template<> struct h5_type<int32_t> { static H5::DataType type() { return H5::PredType::NATIVE_INT32; } };
template<> struct h5_type<int64_t> { static H5::DataType type() { return H5::PredType::NATIVE_INT64; } };
template<> struct h5_type<uint32_t> { static H5::DataType type() { return H5::PredType::NATIVE_UINT32; } };
template<> struct h5_type<uint64_t> { static H5::DataType type() { return H5::PredType::NATIVE_UINT64; } };
template<> struct h5_type<float> { static H5::DataType type() { return H5::PredType::NATIVE_FLOAT; } };
template<> struct h5_type<double> { static H5::DataType type() { return H5::PredType::NATIVE_DOUBLE; } };

//----------------------------------------------------------------------

class EHFile {
    H5::Group open_or_create_group(std::vector<std::string> path);

public:
    /** \brief Create a new EHFile
     *
     * Possible values for \a mode,
     *
     *  Mode | Description
     * ------|--------------------------------------
     *    r  | Open in read only. If file does not exist throw exception.
     *   r+  | Open in read-write. If file does not exist throw exception.
     *    w  | Open in read-write. If file already exists overwrite.
     *    x  | Create file. If it already exists, throw exception.
     *
     * Default behavior is `r` open file in read-only.
     *
     * \param filename Path to the file to open or create
     * \param mode Mode to open file in see above for possible modes
     *
     * \throw invalid_argument If mode is invalid
     * \throw H5::FileIException If a file error occurs
     */
    EHFile(const std::string& filename, const std::string& mode="r");

    /// Creates an invalid EHFile()
    EHFile();

    /** \brief Open or create an HDF5 file
     *
     * Possible values for \a mode,
     *
     *  Mode | Description
     * ------|--------------------------------------
     *    r  | Open in read only. If file does not exist throw exception.
     *   r+  | Open in read-write. If file does not exist throw exception.
     *    w  | Open in read-write. If file already exists overwrite.
     *    x  | Create file. If it already exists, throw exception.
     *
     * Default behavior is `r` open file in read-only.
     *
     * \param filename Path to the file to open or create
     * \param mode Mode to open file in see above for possible modes
     *
     * \throw invalid_argument If mode is invalid
     * \throw H5::FileIException If a file error occurs
     */
    void open(const std::string& filename, const std::string& mode="r");

    /** Reopen the file if it was closed
     *
     * \throw H5::FileIException If a file error occurs
     */
    void reopen();

    /// Check whether a dataset exists in the file
    bool has_data(const std::string& path_string);

    /// Check whether a group exists in the file
    bool has_group(const std::string& path);

    // TODO bool has_attr(std::string path, std::string name);
    // TODO void rm_attr(std::string path, std::string name);

    /** \brief Immediately write data to disk and close the file
     *
     * \see reopen()
     * \throw H5::FileIException If a file error occurs
     */
    void close();

    /// Returns the pointer to the H5::H5File
    H5::H5File* get_file_ptr();

    //------------------------------------------------------------------
    // General Types

    void set_data(const std::string& path, const void* data_ptr,
            unsigned dim_n, H5::DataType type,
            const hsize_t* dim);

    void set_attr(const std::string& path, const std::string& name,
            const void* data_ptr, H5::DataType type,
            const hsize_t* dim = NULL, unsigned dim_n=0);

    /** \brief Check the type, dimension and size of a dataset
     *
     * This function allows one to easily insure that a dataset has the
     * correct type, dimension and size for reading. Note that if `dim_n=0`
     * then the dimension checking is disabled and if `size=NULL` or `dim_n=0`
     * then size checking is disabled.
     *
     * \param[in] dataset The dataset to check
     * \param[in] type The type the dataset must have
     * \param[in] dim_n The number of dimensions the data must have
     * \param[in] size Array which contains size dataset must have
     *
     * \throw logic_error if data set does not have required type/dim_n/size
     * \throw logic_error if data set is not simple
     * \throw logic_error if data set is greater than EHFILE_MAX_DIM
     */
    void check_dataset(const H5::DataSet& dataset, H5::DataType type,
            int dim_n=-1, const hsize_t *size=NULL);

    /** \brief Check the type, dimension and size of a dataset
     *
     * This function allows one to easily insure that a dataset has the
     * correct type, dimension and size for reading/writing. Note that if
     * `dim_n<0` then the dimension checking is disabled and if `size=NULL` or
     * `dim_n<=1` then size checking is disabled. If `dim=0` then the dataset
     * is required to be scalar, if `dim > 0` then it is check to have the
     * specified dimension.
     *
     * \param[in] dataset The dataset to check
     * \param[in] type The type the dataset must have
     * \param[in] dim_n The number of dimensions the data must have
     * \param[in] size Array which contains size dataset must have
     *
     * \throw logic_error if data set does not have required type/dim_n/size
     * \throw logic_error if data set is not simple
     * \throw logic_error if data set is greater than EHFILE_MAX_DIM
     */
    void check_attr(const H5::Attribute& attr, H5::DataType type,
            int dim_n=-1, const hsize_t *size=NULL);

    /// Get the number of dimensions a dataset has
    unsigned get_dim_dataset(const H5::DataSet& dataset);

    /// Get the size of the dataset
    void get_size_dataset(const H5::DataSet& dataset, hsize_t *size);

    H5::Attribute open_attr(const std::string& path, const std::string& name);
    void get_size_attr(const H5::Attribute& attr, hsize_t *size);


    //------------------------------------------------------------------
    // Vectors

    template<typename T>
    void set_data(const std::string& path, const std::vector<T>& data) {
        hsize_t dim[1];
        dim[0] = data.size();
        set_data(path, (void*)data.data(), 1, h5_type<T>::type(), dim);
    }

    template<typename T>
    void get_data(const std::string& path, std::vector<T>& data) {
        H5::DataType type = h5_type<T>::type();
        hsize_t size[1];

        H5::DataSet dataset(file.openDataSet(path));
        check_dataset(dataset, type, 1);
        get_size_dataset(dataset, size);

        data.resize(static_cast<size_t>(size[0]));
        dataset.read(data.data(), type);
    }

    void set_data(const std::string& path, const std::vector<std::complex<double>>& data);
    void get_data(const std::string& path, std::vector<std::complex<double>>& data);

    template<typename T>
    void set_attr(const std::string& path, const std::string& name, const std::vector<T>& data) {
        hsize_t dim[1];
        dim[0] = data.size();
        set_attr(path, name, (void*)data.data(), h5_type<T>::type(), dim, 1);
    }

    template<typename T>
    void get_attr(const std::string& path, const std::string& name, std::vector<T>& data) {
        H5::DataType type = h5_type<T>::type();
        hsize_t size[1];
        H5::Attribute attr = open_attr(path, name);
        get_size_attr(attr, size);

        data.resize(static_cast<size_t>(size[0]));
        attr.read(type, data.data());
    }


    //------------------------------------------------------------------
    // NArray

    template<typename T, size_t N>
    void set_data(const std::string& path,
            const tbpp::math::NArray<T,N>& data) {
        hsize_t dim[N];
        for(size_t i=0; i<N; i++)
            dim[i] = data.size(i);
        set_data(path, (void*)data.data(), N, h5_type<T>::type(), dim);
    }

    template<typename T, size_t N>
    void get_data(const std::string& path, tbpp::math::NArray<T,N>& data) {
        H5::DataType type = h5_type<T>::type();
        hsize_t size[N];

        H5::DataSet dataset(file.openDataSet(path));
        check_dataset(dataset, type, N);
        get_size_dataset(dataset, size);

        std::array<size_t,N> new_sizes;
        for(size_t i=0; i<N; i++) {
            new_sizes[i] = size[i];
        }
        data.resize(new_sizes);

        dataset.read(data.data(), type);
    }

    template<size_t N>
    void set_data(const std::string& path,
            const tbpp::math::NArray<std::complex<double>, N>& data) {
        H5::CompType type(sizeof(std::complex<double>));
        type.insertMember("r", 0, H5::PredType::NATIVE_DOUBLE);
        type.insertMember("i", sizeof(double), H5::PredType::NATIVE_DOUBLE);

        hsize_t dim[N];
        for(size_t i=0; i<N; i++)
            dim[i] = data.size(i);
        set_data(path, (void*)data.data(), N, type, dim);
    }

    template<size_t N>
    void get_data(const std::string& path,
            tbpp::math::NArray<std::complex<double>, N>& data) {
        H5::CompType type(sizeof(std::complex<double>));
        type.insertMember("r", 0, H5::PredType::NATIVE_DOUBLE);
        type.insertMember("i", sizeof(double), H5::PredType::NATIVE_DOUBLE);
        hsize_t size[N];

        H5::DataSet dataset(file.openDataSet(path));
        check_dataset(dataset, type, N);
        get_size_dataset(dataset, size);

        std::array<size_t,N> new_sizes;
        for(size_t i=0; i<N; i++) {
            new_sizes[i] = size[i];
        }
        data.resize(new_sizes);

        dataset.read(data.data(), type);
    }

    //------------------------------------------------------------------
    // Scalar Attributes

    void set_attr(const std::string& path, const std::string& name, const char* value);

    void set_attr(const std::string& path, const std::string& name, const std::string& value);
    void get_attr(const std::string& path, const std::string& name, std::string& value);

    void set_attr(const std::string& path, const std::string& name, bool value);
    void get_attr(const std::string& path, const std::string& name, bool& value);

    template<typename T>
    void set_attr(const std::string& path, const std::string& name, T value) {
        set_attr(path, name, (void*)&value, h5_type<T>::type());
    }

    template<typename T>
    void get_attr(const std::string& path, const std::string& name, T &value) {
        H5::DataType type = h5_type<T>::type();
        H5::Attribute attribute = open_attr(path, name);
        check_attr(attribute, type, 0);
        attribute.read(type, &value);
    }

    //------------------------------------------------------------------
    // TODO Low level interface for 1D and N-Dim arrays

    // template<typename T>
    // void set_data(std::string path, T *data, unsigned len);
    // template<typename T>
    // void get_data(std::string path, T *data, unsigned len);
    //
    // template<typename T>
    // void set_data(std::string path, T *data, unsigned* n, unsigned dim);
    // template<typename T>
    // void get_data(std::string path, T *data, unsigned* n, unsigned dim);

private:
    H5::H5File file;
};

//----------------------------------------------------------------------

} // namespace tbpp

#endif /* __TBPP_EHFILE_H__ */

