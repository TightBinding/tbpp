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

#include <tbpp/kpath.h>

using namespace std;

namespace tbpp {

//----------------------------------------------------------------------

const char* KPath::type() const { return "KPath"; };

void KPath::set_model(ModelPtr model) { _model = model; }
ModelPtr KPath::model() const { return _model; }

void KPath::add_kpoint(const string& label, double k1, double k2, double k3) {
    KPoint kp;
    for(size_t i=0; i<TBPP_MAX_SITE_KIND_LEN; i++)
        kp.label[i] = 0;
    label.copy(kp.label, TBPP_MAX_SITE_KIND_LEN-1);
    kp.k[0] = k1;
    kp.k[1] = k2;
    kp.k[2] = k3;
    kpoints.push_back(kp);
}

void KPath::solve() {
    _status = NodeStatus::failed;
    //-----------------------------------------------------------------------

    if (_model == nullptr) throw runtime_error("No Model");
    _model->make_ready();

    header();
    //-----------------------------------------------------------------------

    size_t states = _model->states();
    eigval.resize((kpoints.size()-1)*steps, states);
    k.resize((kpoints.size()-1)*steps, 3);

    if (solve_eigvec)
        eigvec.resize((kpoints.size()-1)*steps, states, states);
    else
        eigvec.clear();

    #pragma omp parallel if(_parallel)
    {
        double *__restrict__  pk;
        math::DenseMatrix<cxdouble> Hk(states, states);
        Eigen::SelfAdjointEigenSolver<math::DenseMatrix<cxdouble>> eig;

        #pragma omp for collapse(2)
        for(size_t i=0; i<kpoints.size()-1; i++) {
            for(size_t j=0; j<steps; j++) {

                // index for current k-point
                size_t index = i*steps + j;

                // interpolate to find k1,k2,k3 for current step
                pk = &k(index,0);
                for(size_t l=0; l<3; l++) {
                    pk[l] = kpoints[i].k[l]
                        + j*(kpoints[i+1].k[l] - kpoints[i].k[l])/(steps-1);
                }

                // get Hamiltonian
                _model->Hk(Hk.data(), pk[0], pk[1], pk[2]);

                // diagonalize
                if(solve_eigvec)
                    eig.compute(Hk, Eigen::ComputeEigenvectors);
                else
                    eig.compute(Hk, Eigen::EigenvaluesOnly);

                // save eigenvalues and eigenvectors
                for(size_t iE=0; iE < states; iE++) {
                    eigval(index, iE) = eig.eigenvalues()[iE];
                    if(solve_eigvec) {
                        for(size_t il=0; il < states; il++) {
                            eigvec(index, iE, il) = eig.eigenvectors().col(iE)[il];
                        }
                    }
                }
            }
        }
    } // #pragma omp parallel

    _status = NodeStatus::done;
}

#ifdef TBPP_WITH_HDF5
void KPath::save(EHFile& file, const string& prefix) const {
    Node::save(file, prefix);
    if(_model != nullptr)
        file.set_attr(prefix, "model", _model->name());
    else
        file.set_attr(prefix, "model", "");

    file.set_attr(prefix, "steps", steps);
    file.set_attr(prefix, "solve_eigvec", solve_eigvec);

    // kpoints
    {
        H5::CompType type(sizeof(KPoint));
        type.insertMember("label", HOFFSET(KPoint, label),
                H5::StrType(H5::PredType::C_S1, TBPP_MAX_SITE_KIND_LEN));
        type.insertMember("k1", HOFFSET(KPoint, k[0]), H5::PredType::NATIVE_DOUBLE);
        type.insertMember("k2", HOFFSET(KPoint, k[1]), H5::PredType::NATIVE_DOUBLE);
        type.insertMember("k3", HOFFSET(KPoint, k[2]), H5::PredType::NATIVE_DOUBLE);

        hsize_t dim = kpoints.size();
        file.set_data(prefix+"/kpoints", (void*)&kpoints[0], 1, type, &dim);
    }

    if(_status == NodeStatus::done) {
        file.set_data(prefix+"/eigval", eigval);
        file.set_data(prefix+"/eigvec", eigvec);
        file.set_data(prefix+"/k", k);
    }
}

void KPath::load(EHFile& file, const string& prefix) {
    Node::load(file, prefix);
    string model_name;
    file.get_attr(prefix, "model", model_name);
    if(model_name.length() > 0 and context() != nullptr)
        _model = dynamic_pointer_cast<Model>(context()->get(model_name));
    else
        _model = nullptr;

    file.get_attr(prefix, "steps", steps);
    file.get_attr(prefix, "solve_eigvec", solve_eigvec);

    // kpoints
    {
        H5::CompType type(sizeof(KPoint));
        type.insertMember("label", HOFFSET(KPoint, label),
                H5::StrType(H5::PredType::C_S1, TBPP_MAX_SITE_KIND_LEN));
        type.insertMember("k1", HOFFSET(KPoint, k[0]), H5::PredType::NATIVE_DOUBLE);
        type.insertMember("k2", HOFFSET(KPoint, k[1]), H5::PredType::NATIVE_DOUBLE);
        type.insertMember("k3", HOFFSET(KPoint, k[2]), H5::PredType::NATIVE_DOUBLE);

        hsize_t size;
        H5::DataSet dataset(file.get_file_ptr()->openDataSet(prefix+"/kpoints"));
        file.check_dataset(dataset, type, 1);
        file.get_size_dataset(dataset, &size);
        kpoints.resize(static_cast<size_t>(size));
        dataset.read(kpoints.data(), type);
    }

    if(_status == NodeStatus::done) {
        file.get_data(prefix+"/eigval", eigval);
        file.get_data(prefix+"/eigvec", eigvec);
        file.get_data(prefix+"/k", k);
    }
}
#endif // TBPP_WITH_HDF5

//----------------------------------------------------------------------

} // namespace tbpp
