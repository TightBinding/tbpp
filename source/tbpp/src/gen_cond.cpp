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

#include <tbpp/gen_cond.h>
#include <tbpp/common.h>
#include <cstdio>

// Disable Eigen parallelization
#define EIGEN_DONT_PARALLELIZE

using namespace std;

namespace tbpp {

//----------------------------------------------------------------------

const char* GenCond::type() const { return "GenCond"; }
void GenCond::set_sigma(CPAPtr sigma) { _sigma = sigma; }
CPAPtr GenCond::sigma() const { return _sigma; }

//----------------------------------------------------------------------

void GenCond::solve() {
    _status = NodeStatus::failed;

    //-----------------------------------------------------------------------

    if (_sigma == nullptr) throw runtime_error("No Sigma");

    KModelPtr _kmodel = _sigma->kmodel();
    if (_kmodel == nullptr) throw runtime_error("No KModel");

    ModelPtr _model = _kmodel->model();
    if (_model == nullptr) throw runtime_error("No Model");

    _model->make_ready();
    _kmodel->make_ready();
    _sigma->make_ready();

    header();

    //-----------------------------------------------------------------------

    const size_t states = _model->states();
    const size_t sites = _model->sites();
    const size_t states2 = states*states;
    const size_t kdim = _model->kdim();
    const bool cache = _kmodel->cache();

    const vector<double>& k1 = _kmodel->k1;
    const vector<double>& k2 = _kmodel->k2;
    const vector<double>& k3 = _kmodel->k3;

    const size_t kpoints = k1.size()*k2.size()*k3.size();
    const size_t c_size = _sigma->sigma.size(0);
    const size_t w_size = _sigma->sigma.size(1);

    if((!w_size % 2)) error("Must have odd number of energy samples");

    // Will be set later inside the parallel loop
    size_t threads = 1;
    // reduction for phi [ithread,ic,iw,jw,i,j]
    math::NArray<std::complex<double>,6> red_phi;

    //-----------------------------------------------------------------------
    // Prepare data buffers

    cond.resize(c_size, w_size, 3, 3);
    cond.fill(0);

    phi.resize(c_size, w_size, w_size, 3, 3);
    phi.fill(0);

    #pragma omp parallel if(_parallel)
    {
        //-------------------------------------------------------------------
        // Reduction buffers

        #pragma omp master
        {
            threads = omp_get_num_threads();
            red_phi.resize(threads, c_size, w_size, w_size, 3, 3);
            red_phi.fill(0);
        }

        #pragma omp barrier

        //-------------------------------------------------------------------
        // Temporary data buffers for each thread

        // unique thread identifier
        size_t it = omp_get_thread_num();

        // G(k) = Inv(w*I - T - S)
        math::DenseMatrix<cxdouble> mGk(states, states);
        // P = 1/(2*pi)(conj(Gk.T) - Gk)
        std::deque<math::DenseMatrix<cxdouble>> mP(w_size);

        math::DenseMatrix<cxdouble> mvx(states, states);
        math::DenseMatrix<cxdouble> mvy(states, states);
        math::DenseMatrix<cxdouble> mvz(states, states);

        cxdouble *__restrict__ Gk = mGk.data();

        cxdouble *__restrict__ vx = &mvx(0,0);
        cxdouble *__restrict__ vy = &mvy(0,0);
        cxdouble *__restrict__ vz = &mvz(0,0);

        cxdouble *p1 = nullptr; //< General purpose pointer

        //-------------------------------------------------------------------

        #pragma omp for collapse(3)
        for(size_t ik1=0; ik1<k1.size(); ik1++)
        for(size_t ik2=0; ik2<k2.size(); ik2++)
        for(size_t ik3=0; ik3<k3.size(); ik3++) {
            //-----------------------------------------------------------
            // Get Velocity Matrix

            _model->dH_dk(vx, vy, vz, k1[ik1], k2[ik2], k3[ik3]);

            for(size_t ic=0; ic<c_size; ic++) {
                for(size_t iw=0; iw<w_size; iw++) {
                    //-------------------------------------------------------
                    // Compute Gk

                    cxdouble *SIGMA = &_sigma->sigma(ic,iw,0,0);
                    if(cache) {
                        p1 = &_kmodel->Tk(ik1,ik2,ik3,0,0);
                        for(size_t i=0; i<states2; i++) {
                            Gk[i] = -(p1[i] + SIGMA[i]);
                        }
                    } else {
                        _model->Tk(Gk, k1[ik1],k2[ik2],k3[ik3]);
                        // Now: Gk = T
                        for(size_t i=0; i<states2; i++) {
                            Gk[i] = -(Gk[i] + SIGMA[i]);
                        }
                    }
                    // Now: Gk = -T-S
                    for(size_t i=0; i<states; i++) {
                        cxdouble& v = Gk[i*states + i];
                        v += _sigma->w[iw];
                    }

                    // Now: Gk = w*I - T - S
                    mGk = mGk.inverse();

                    // Now: Gk = Inv(w*I - T - S)
                    mP[iw] = cxdouble(0, 1/(2*data::pi))*(mGk.adjoint() - mGk);
                }

                for(size_t iw1=0; iw1<w_size; iw1++)
                for(size_t iw2=0; iw2<w_size; iw2++) {

                    if(solve_sigx) {
                        // solve_sigxx
                        red_phi(it,ic,iw1,iw2,0,0) += (mvx*mP[iw1]*mvx*mP[iw2]).trace();
                        if(solve_sigy) {
                            // solve_sigxy
                            red_phi(it,ic,iw1,iw2,0,1) += (mvx*mP[iw1]*mvy*mP[iw2]).trace();
                            // solve_sigyx
                            red_phi(it,ic,iw1,iw2,1,0) += (mvy*mP[iw1]*mvx*mP[iw2]).trace();
                        }
                        if(solve_sigz) {
                            // solve_sigxz
                            red_phi(it,ic,iw1,iw2,0,2) += (mvx*mP[iw1]*mvz*mP[iw2]).trace();
                            // solve_sigzx
                            red_phi(it,ic,iw1,iw2,2,0) += (mvz*mP[iw1]*mvx*mP[iw2]).trace();
                        }
                    }

                    if(solve_sigy) {
                        // solve_sigyy
                        red_phi(it,ic,iw1,iw2,1,1) += (mvy*mP[iw1]*mvy*mP[iw2]).trace();
                        if(solve_sigz) {
                            // solve_sigyz
                            red_phi(it,ic,iw1,iw2,1,2) += (mvy*mP[iw1]*mvz*mP[iw2]).trace();
                            // solve_sigzy
                            red_phi(it,ic,iw1,iw2,2,1) += (mvz*mP[iw1]*mvy*mP[iw2]).trace();
                        }
                    }

                    if(solve_sigz) {
                        // solve_sigzz
                        red_phi(it,ic,iw1,iw2,2,2) += (mvz*mP[iw1]*mvz*mP[iw2]).trace();
                    }
                }

                //-----------------------------------------------------------
            }
            #pragma omp critical (stdout)
            {
                if (not (ik1 % 5 or ik2 % 5 or ik3 % 5))
                    cout << "k1=" << ik1 << " k2=" << ik2 << " ik3=" << ik3 << '\n';
            }
        }
    }

    // Reduce phi
    for(unsigned ic=0; ic<c_size; ic++)
    for(unsigned iw1=0; iw1<w_size; iw1++)
    for(unsigned iw2=0; iw2<w_size; iw2++) {
        for(unsigned i=0; i<3; i++)
        for(unsigned j=0; j<3; j++) {
            for(unsigned it=0; it<threads; it++) {
                phi(ic,iw1,iw2,i,j) += red_phi(it,ic,iw1,iw2,i,j);
            }
            phi(ic,iw1,iw2,i,j) /= kpoints;
        }
    }
    red_phi.clear();

    // if (solve_cond) {
    //     double cond_norm = data::pi*pow(data::ec,2)/(
    //             data::hbar*_model->uc_size());
    //
    //     // Normalize cond_state and compute cond
    //     for(size_t ic=0; ic<c_size; ic++)
    //     for(size_t iw=0; iw<w_size; iw++)
    //     for(size_t i=0; i<3; i++)
    //     for(size_t j=0; j<3; j++)
    //     for(size_t s=0; s<states; s++) {
    //         for(size_t it=0; it<threads; it++) {
    //             cond_state(ic,iw,i,j,s) += red_cond_state(it,ic,iw,i,j,s);
    //         }
    //         cond_state(ic,iw,i,j,s) *= cond_norm;
    //         cond(ic,iw,i,j) += cond_state(ic,iw,i,j,s);
    //     }
    // }

    //-----------------------------------------------------------------------

    _status = NodeStatus::done;
}

#ifdef TBPP_WITH_HDF5
void GenCond::save(EHFile& file, const string& prefix) const {
    Node::save(file, prefix);

    if(_sigma != nullptr)
        file.set_attr(prefix, "sigma", _sigma->name());
    else
        file.set_attr(prefix, "sigma", "");

    file.set_attr(prefix, "solve_sigx", solve_sigx);
    file.set_attr(prefix, "solve_sigy", solve_sigy);
    file.set_attr(prefix, "solve_sigz", solve_sigz);
    file.set_attr(prefix, "zero_j", zero_j);

    if(_status == NodeStatus::done) {
        file.set_data(prefix+"/cond", cond);
        file.set_data(prefix+"/phi", phi);
    }
}

void GenCond::load(EHFile& file, const string& prefix) {
    Node::load(file, prefix);

    string sigma_name;
    file.get_attr(prefix, "sigma", sigma_name);
    if(sigma_name.length() > 0 and context() != nullptr)
        _sigma = dynamic_pointer_cast<CPA>(context()->get(sigma_name));
    else
        _sigma = nullptr;

    file.get_attr(prefix, "solve_sigx", solve_sigx);
    file.get_attr(prefix, "solve_sigy", solve_sigy);
    file.get_attr(prefix, "solve_sigz", solve_sigz);
    file.get_attr(prefix, "zero_j", zero_j);

    if(_status == NodeStatus::done) {
        file.get_data(prefix+"/cond", cond);
        file.get_data(prefix+"/phi", phi);
    }
}
#endif // TBPP_WITH_HDF5

//----------------------------------------------------------------------

} // namespace tbpp
