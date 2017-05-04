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

#include <tbpp/berryflux.h>

using namespace std;

namespace tbpp {

//----------------------------------------------------------------------

const char* BerryFlux::type() const { return "BerryFlux"; }

void BerryFlux::set_kgrid(KGridPtr kgrid) { _kgrid = kgrid; }
KGridPtr BerryFlux::kgrid() const { return _kgrid; }

void BerryFlux::solve() {
    _status = NodeStatus::failed;
    //-----------------------------------------------------------------------

    if (_kgrid == nullptr) error("No KGrid");

    KModelPtr _kmodel = _kgrid->kmodel();
    if (_kmodel == nullptr) error("No KModel");

    ModelPtr _model = _kmodel->model();
    if (_model == nullptr) error("No Model");

    _kgrid->make_ready();
    _model->make_ready();
    _kmodel->make_ready();

    if(_kgrid->eigvec.empty() or !_kgrid->solve_eigvec)
        error("KGrid must solve for eigenvectors");

    header();
    //-----------------------------------------------------------------------

    const vector<double>& k1 = _kmodel->k1;
    const vector<double>& k2 = _kmodel->k2;
    const vector<double>& k3 = _kmodel->k3;

    const size_t kpoints = k1.size()*k2.size()*k3.size();
    size_t states = _model->states();

    // Will be set later inside the parallel loop
    size_t threads = 1;
    // reduction for flux [ithread,iE,i]
    math::NArray<double, 3> red_flux;

    flux.resize(states, 3);
    if(solve_flux_k)
        flux_k.resize(k1.size(), k2.size(), k3.size(), states, 3);
    //-----------------------------------------------------------------------

    #pragma omp parallel if(_parallel)
    {
        math::DenseMatrix<cxdouble> dH_dkx(states, states);
        math::DenseMatrix<cxdouble> dH_dky(states, states);
        math::DenseMatrix<cxdouble> dH_dkz(states, states);
        cxdouble nxm, nym, nzm, mxn, myn, mzn;

        // unique thread identifier
        size_t it = omp_get_thread_num();

        #pragma omp master
        {
            threads = omp_get_num_threads();
            red_flux.resize(threads, states, 3);
            red_flux.fill(0);
        }

        #pragma omp barrier

        #pragma omp for collapse(3)
        for(size_t ik1=0; ik1<k1.size(); ik1++)
        for(size_t ik2=0; ik2<k2.size(); ik2++)
        for(size_t ik3=0; ik3<k3.size(); ik3++) {

            // Hamiltonian derivative with respect to reduced coordinates?
            // FIXME should be derivative with respect to k1,k2,k3 not kx,ky,kz
            _model->dH_dk(dH_dkx.data(), dH_dky.data(), dH_dkz.data(),
                    k1[ik1], k2[ik2], k3[ik3]);

            for(size_t n=0; n<states; n++) {
                double vx = 0;
                double vy = 0;
                double vz = 0;
                double En = _kgrid->eigval(ik1,ik2,ik3,n);
                Eigen::Map<Eigen::VectorXcd> un(&_kgrid->eigvec(ik1,ik2,ik3,n,0), states);

                for(size_t m=0; m<states; m++) {
                    if(m == n) continue;
                    Eigen::Map<Eigen::VectorXcd> um(&_kgrid->eigvec(ik1,ik2,ik3,m,0), states);
                    double d = 1.0/pow(abs(En - _kgrid->eigval(ik1,ik2,ik3,m)) + dE, 2);

                    nxm = (un.adjoint()*(dH_dkx*um))(0);
                    nym = (un.adjoint()*(dH_dky*um))(0);
                    nzm = (un.adjoint()*(dH_dkz*um))(0);

                    mxn = (um.adjoint()*(dH_dkx*un))(0);
                    myn = (um.adjoint()*(dH_dky*un))(0);
                    mzn = (um.adjoint()*(dH_dkz*un))(0);

                    vx += d*imag(nym*mzn - nzm*myn);
                    vy += d*imag(nzm*mxn - nxm*mzn);
                    vz += d*imag(nxm*myn - nym*mxn);

                }

                red_flux(it,n,0) += vx;
                red_flux(it,n,1) += vy;
                red_flux(it,n,2) += vz;

                // save berry flux
                if (solve_flux_k) {
                    flux_k(ik1,ik2,ik3,n,0) = vx;
                    flux_k(ik1,ik2,ik3,n,1) = vy;
                    flux_k(ik1,ik2,ik3,n,2) = vz;
                }
            }

        }
    } // #pragma omp parallel

    // Reduction
    for(size_t n=0; n<states; n++) {
        for(size_t it=0; it < threads; it++) {
            flux(n,0) += red_flux(it,n,0)/kpoints;
            flux(n,1) += red_flux(it,n,1)/kpoints;
            flux(n,2) += red_flux(it,n,2)/kpoints;
        }
    }

    _status = NodeStatus::done;
}

#ifdef TBPP_WITH_HDF5
void BerryFlux::save(EHFile& file, const string& prefix) const {
    Node::save(file, prefix);
    file.set_attr(prefix, "solve_flux_k", solve_flux_k);
    if(_kgrid != nullptr)
        file.set_attr(prefix, "kgrid", _kgrid->name());
    else
        file.set_attr(prefix, "kgrid", "");

    if(_status == NodeStatus::done) {
        if (solve_flux_k)
            file.set_data(prefix+"/flux_k", flux_k);
        file.set_data(prefix+"/flux", flux);
    }
}
void BerryFlux::load(EHFile& file, const string& prefix) {
    Node::load(file, prefix);
    file.get_attr(prefix, "solve_flux_k", solve_flux_k);

    string kgrid_name;
    file.get_attr(prefix, "kgrid", kgrid_name);
    if(kgrid_name.length() > 0 and context() != nullptr)
        _kgrid = dynamic_pointer_cast<KGrid>(context()->get(kgrid_name));
    else
        _kgrid = nullptr;

    if(_status == NodeStatus::done) {
        if (solve_flux_k)
            file.get_data(prefix+"/flux_k", flux_k);
        file.get_data(prefix+"/flux", flux);
    }
}
#endif // TBPP_WITH_HDF5

//----------------------------------------------------------------------

} // namespace tbpp
