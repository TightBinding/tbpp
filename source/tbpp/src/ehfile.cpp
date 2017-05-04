#include <tbpp/ehfile.h>
#include <sstream>
#include <stdexcept>
#include <cstring>

using namespace std;
using namespace H5;
using namespace tbpp::math;

namespace tbpp {

//----------------------------------------------------------------------

/** \brief Splits a string into a vector based on delimiter
 *
 * \param[in] s string to split
 * \param[in] m delimiter
 * \return A vector<string> containing the split string
 */
vector<string> split(string s, char m='/') {
    vector<string> ret;
    string frag;
    stringstream stream(s);

    while(std::getline(stream, frag, m)) {
        if(!frag.empty())
            ret.push_back(std::move(frag));
    }
    return ret;
}

EHFile::EHFile() {
    H5::Exception::dontPrint();
}

EHFile::EHFile(const std::string& filename, const std::string& mode) {
    H5::Exception::dontPrint();
    open(filename, mode);
}

void EHFile::open(const std::string& filename, const std::string& mode) {
    unsigned flag;
    if (mode == "r")
        flag = H5F_ACC_RDONLY;
    else if(mode == "r+")
        flag = H5F_ACC_RDWR;
    else if(mode == "w")
        flag = H5F_ACC_TRUNC;
    else if(mode == "x")
        flag = H5F_ACC_EXCL;
    else
        throw std::invalid_argument("Invalid file mode");
    file = H5File(filename, flag);
}


void EHFile::reopen() {
    file.reOpen();
}

void EHFile::close() {
    file.close();
}

H5File* EHFile::get_file_ptr() {
    return &file;
}

H5::Group EHFile::open_or_create_group(std::vector<std::string> path) {
    string absolute_path = "/";
    Group group;
    for(vector<string>::iterator it=path.begin(); it != path.end(); ++it) {
        absolute_path += (*it) + "/";
        try {
            group = file.openGroup(absolute_path);
        } catch( FileIException no_group ) {
            group = file.createGroup(absolute_path);
        }
    }
    return group;
}

bool EHFile::has_group(const std::string& path) {
    try {
        Group group = file.openGroup(path);
    } catch ( FileIException no_group ) {
        return false;
    }
    return true;
}

void EHFile::set_data(const std::string& path_string, const void* data_ptr, unsigned dim_n, H5::DataType type,
        const hsize_t* dim) {

    H5::DataSpace space(dim_n, dim);
    H5::DataSet dataset;

    std::vector<std::string> path = split(path_string);
    std::string data_name = path.back(); // Get name of array to add
    path.pop_back();

    bool check = true;

    if(!path.empty()) {
        // Make sure group exists
        Group group = open_or_create_group(path);
    }

    // Add data to root group
    try {
        dataset = file.openDataSet(path_string);
    } catch( H5::FileIException not_found_error ) {
        dataset = file.createDataSet(path_string, type, space);
        check = false;
    }

    // Check that dataset has correct type and size otherwise raise an error
    if(check) {
        check_dataset(dataset, type, dim_n, dim);
    }
    dataset.write(data_ptr, type);
}


bool EHFile::has_data(const std::string& path_string) {
    try {
        DataSet dataset(file.openDataSet(path_string));
    } catch( H5::FileIException not_found_error ) {
        return false;
    }
    return true;
}

void EHFile::check_dataset(const H5::DataSet& dataset, H5::DataType type,
        int data_dim, const hsize_t *size) {
    // Check type
    if(!(dataset.getDataType() == type))
        throw logic_error("Data set does not have same type");

    // Check dimension
    DataSpace space = dataset.getSpace();
    H5S_class_t class_type = space.getSimpleExtentType();

    if(data_dim == 0 and class_type != H5S_SCALAR)
            throw logic_error("Data set is not scalar");

    if(data_dim > 0) {
        if(!space.isSimple())
            throw logic_error("Data set is not simple");
        int dim_n = space.getSimpleExtentNdims();
        if(data_dim > 0) {
            if(dim_n != data_dim)
                throw logic_error("Data set has different dimension");
            if(size != NULL) {
                // check that size of each dimension is the same
                vector<hsize_t> old_size(dim_n);
                space.getSimpleExtentDims(old_size.data());
                for(int i=0; i < data_dim; i++) {
                    if(old_size[i] != size[i])
                        throw logic_error("Data set has different size");
                }
            }
        }
    }
}

void EHFile::check_attr(const H5::Attribute& attr, H5::DataType type,
        int data_dim, const hsize_t *size) {
    // Check type
    if(!(attr.getDataType() == type))
        throw logic_error("Attribute does not have same type");

    // Check dimension
    DataSpace space = attr.getSpace();
    H5S_class_t class_type = space.getSimpleExtentType();

    if(data_dim == 0 and class_type != H5S_SCALAR)
            throw logic_error("Attribute is not scalar");

    if(data_dim > 0) {
        if(!space.isSimple())
            throw logic_error("Attribute is not simple");
        int dim_n = space.getSimpleExtentNdims();
        if(data_dim > 0) {
            if(dim_n != data_dim)
                throw logic_error("Attribute has different dimension");
            if(size != NULL) {
                // check that size of each dimension is the same
                vector<hsize_t> old_size(dim_n);
                space.getSimpleExtentDims(old_size.data());
                for(int i=0; i < data_dim; i++) {
                    if(old_size[i] != size[i])
                        throw logic_error("Attribute has different size");
                }
            }
        }
    }
}

unsigned EHFile::get_dim_dataset(const H5::DataSet& dataset) {
    DataSpace space = dataset.getSpace();
    return (unsigned)space.getSimpleExtentNdims();
}

void EHFile::get_size_dataset(const H5::DataSet& dataset, hsize_t *size) {
    DataSpace space = dataset.getSpace();
    space.getSimpleExtentDims(size);
}

void EHFile::set_attr(const std::string& path, const std::string& name,
        const void* data_ptr, H5::DataType type, const hsize_t* dim, unsigned dim_n) {
    DataSpace attr_space;
    Attribute attribute;

    if(dim_n == 0)
        attr_space = DataSpace(H5S_SCALAR);
    else
        attr_space = DataSpace(dim_n, dim);

    bool check = true;
    vector<string> path_split = split(path);
    if(path_split.empty()) {
        // Add to root
        if(file.attrExists(name))
            attribute = file.openAttribute(name);
        else {
            attribute = file.createAttribute(name, type, attr_space);
            check = false;
        }
    } else {
        try {
            H5O_type_t class_type = file.childObjType(path);
            if(class_type == H5O_TYPE_GROUP) {
                // Add to group
                Group group = file.openGroup(path);
                if(group.attrExists(name))
                    attribute = group.openAttribute(name);
                else {
                    attribute = group.createAttribute(name, type, attr_space);
                    check = false;
                }
            } else if (class_type== H5O_TYPE_DATASET) {
                // Add to DataSet
                DataSet dataset = file.openDataSet(path);
                if(dataset.attrExists(name))
                    attribute = dataset.openAttribute(name);
                else {
                    attribute = dataset.createAttribute(name, type, attr_space);
                    check = false;
                }
            }
        } catch(FileIException group_does_not_exist) {
            // Create a new group and add attr
            Group group = open_or_create_group(path_split);
            if(group.attrExists(name))
                attribute = group.openAttribute(name);
            else {
                attribute = group.createAttribute(name, type, attr_space);
                check = false;
            }
        }
    }

    if(check) {
        check_attr(attribute, type, dim_n, dim);
    }

    attribute.write(type, data_ptr);
}

H5::Attribute EHFile::open_attr(const std::string& path, const std::string& name) {
    Attribute attribute;

    vector<string> path_split = split(path);
    if(path_split.empty()) {
        // Get from root
        if(!file.attrExists(name))
            throw logic_error("Attribute does not exist");
        attribute = file.openAttribute(name);
    } else {
        H5O_type_t class_type = file.childObjType(path);
        if(class_type == H5O_TYPE_GROUP) {
            // Add to group
            Group group = file.openGroup(path);
            if(!group.attrExists(name))
                throw logic_error("Attribute does not exist");
            attribute = group.openAttribute(name);
        } else if (class_type== H5O_TYPE_DATASET) {
            // Add to DataSet
            DataSet dataset = file.openDataSet(path);
            if(!dataset.attrExists(name))
                throw logic_error("Attribute does not exist");
            attribute = dataset.openAttribute(name);
        }
    }
    return attribute;
}

void EHFile::get_size_attr(const H5::Attribute& attr, hsize_t *size) {
    DataSpace space = attr.getSpace();
    space.getSimpleExtentDims(size);
}

void EHFile::set_attr(const std::string& path, const std::string& name, const std::string& value) {
    unsigned buf_length = EHFILE_MAX_STR_LEN;
    StrType type(PredType::C_S1, buf_length);
    char* buf = new char[buf_length];
    unsigned len = value.copy(buf, buf_length-1);
    buf[len] = '\0';

    set_attr(path, name, (void*)buf, type);

    delete [] buf;
}

void EHFile::set_attr(const std::string& path, const std::string& name, const char* value) {
    unsigned buf_length = EHFILE_MAX_STR_LEN;
    StrType type(PredType::C_S1, buf_length);
    char* buf = new char[buf_length];
    strncpy(buf, value, buf_length);
    buf[buf_length - 1] = '\0';

    set_attr(path, name, (void*)buf, type);

    delete [] buf;
}

void EHFile::get_attr(const std::string& path, const std::string& name, std::string& value) {
    unsigned buf_length = EHFILE_MAX_STR_LEN;
    StrType type(PredType::C_S1, buf_length);
    char* buf = new char[buf_length];

    Attribute attribute = open_attr(path, name);
    check_attr(attribute, type, 0);
    attribute.read(type, buf);
    value = string(buf);

    delete [] buf;
}

void EHFile::set_attr(const std::string& path, const std::string& name, bool value) {
    int ivalue = 0;
    if(value)
        ivalue = 1;
    set_attr(path, name, (void*)&ivalue, PredType::NATIVE_INT);
}

void EHFile::get_attr(const std::string& path, const std::string& name, bool& value) {
    DataType type = PredType::NATIVE_INT;
    Attribute attribute = open_attr(path, name);
    check_attr(attribute, type, 0);
    int ivalue = 0;
    attribute.read(type, &ivalue);
    if(ivalue)
        value = true;
    else
        value = false;
}

void EHFile::set_data(const std::string& path, const std::vector<std::complex<double>>& data) {
    H5::CompType type(sizeof(std::complex<double>));
    type.insertMember("r", 0, H5::PredType::NATIVE_DOUBLE);
    type.insertMember("i", sizeof(double), H5::PredType::NATIVE_DOUBLE);

    hsize_t dim[1];
    dim[0] = data.size();
    set_data(path, (void*)data.data(), 1, type, dim);
}

void EHFile::get_data(const std::string& path, std::vector<std::complex<double>>& data) {
    H5::CompType type(sizeof(std::complex<double>));
    type.insertMember("r", 0, H5::PredType::NATIVE_DOUBLE);
    type.insertMember("i", sizeof(double), H5::PredType::NATIVE_DOUBLE);

    hsize_t size[1];

    H5::DataSet dataset(file.openDataSet(path));
    check_dataset(dataset, type, 1);
    get_size_dataset(dataset, size);

    data.resize(static_cast<size_t>(size[0]));
    dataset.read(data.data(), type);
}

//----------------------------------------------------------------------

} // namespace tbpp
