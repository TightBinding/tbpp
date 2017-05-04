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

#include <tbpp/kgrid.h>

using namespace std;

namespace tbpp {

//----------------------------------------------------------------------

const char* KGrid::type() const { return "KGrid"; }

void KGrid::set_kmodel(KModelPtr kmodel) { _kmodel = kmodel; }
KModelPtr KGrid::kmodel() const { return _kmodel; }

void KGrid::solve() {
    _status = NodeStatus::failed;
    //-----------------------------------------------------------------------

    if (_kmodel == nullptr) throw runtime_error("No KModel");

    ModelPtr _model = _kmodel->model();
    if (_model == nullptr) throw runtime_error("No Model");

    _model->make_ready();
    _kmodel->make_ready();

    header();
    //-----------------------------------------------------------------------

    const bool cache = _kmodel->cache();

    const vector<double>& k1 = _kmodel->k1;
    const vector<double>& k2 = _kmodel->k2;
    const vector<double>& k3 = _kmodel->k3;

    const size_t kpoints = k1.size()*k2.size()*k3.size();

    size_t states = _model->states();
    eigval.resize(k1.size(), k2.size(), k3.size(), states);

    if(solve_eigvec) {
        eigvec.resize(k1.size(), k2.size(), k3.size(), states, states);
    } else {
        eigvec.clear();
    }

    //-----------------------------------------------------------------------

    #pragma omp parallel if(_parallel)
    {
        math::DenseMatrix<cxdouble> Hk(states, states);
        Eigen::SelfAdjointEigenSolver<math::DenseMatrix<cxdouble>> eig;

        #pragma omp for collapse(3)
        for(size_t ik1=0; ik1<k1.size(); ik1++)
        for(size_t ik2=0; ik2<k2.size(); ik2++)
        for(size_t ik3=0; ik3<k3.size(); ik3++) {

            // get Hamiltonian
            _model->Hk(Hk.data(), k1[ik1], k2[ik2], k3[ik3]);

            // diagonalize
            if(solve_eigvec)
                eig.compute(Hk, Eigen::ComputeEigenvectors);
            else
                eig.compute(Hk, Eigen::EigenvaluesOnly);

            // save eigenvalues and eigenvectors
            for(size_t iE=0; iE < states; iE++) {
                eigval(ik1,ik2,ik3,iE) = eig.eigenvalues()[iE];
            }

            if(solve_eigvec) {
                for(size_t iE=0; iE < states; iE++)
                for(size_t il=0; il < states; il++) {
                    eigvec(ik1,ik2,ik3,iE,il) = eig.eigenvectors().col(iE)[il];
                }
            }
        }
    } // #pragma omp parallel

    _status = NodeStatus::done;
}

#ifdef TBPP_WITH_HDF5
void KGrid::save(EHFile& file, const string& prefix) const {
    Node::save(file, prefix);
    if(_kmodel != nullptr)
        file.set_attr(prefix, "kmodel", _kmodel->name());
    else
        file.set_attr(prefix, "kmodel", "");

    file.set_attr(prefix, "solve_eigvec", solve_eigvec);

    if(_status == NodeStatus::done) {
        file.set_data(prefix+"/eigval", eigval);
        file.set_data(prefix+"/eigvec", eigvec);
    }
}
void KGrid::load(EHFile& file, const string& prefix) {
    Node::load(file, prefix);

    string kmodel_name;
    file.get_attr(prefix, "kmodel", kmodel_name);
    if(kmodel_name.length() > 0 and context() != nullptr)
        _kmodel = dynamic_pointer_cast<KModel>(context()->get(kmodel_name));
    else
        _kmodel = nullptr;

    file.get_attr(prefix, "solve_eigvec", solve_eigvec);

    if(_status == NodeStatus::done) {
        file.get_data(prefix+"/eigval", eigval);
        file.get_data(prefix+"/eigvec", eigvec);
    }
}
#endif // TBPP_WITH_HDF5

//----------------------------------------------------------------------

} // namespace tbpp
