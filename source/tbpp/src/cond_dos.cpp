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

#include <tbpp/cond_dos.h>
#include <tbpp/common.h>

using namespace std;

namespace tbpp {

//----------------------------------------------------------------------

const char* CondDOS::type() const { return "CondDOS"; }

void CondDOS::add_dos_weights(const vector<double>& weights) {
    if(!dos_weights.empty()) {
        if(weights.size() != dos_weights.size(1))
            error("Incorrect length for weights");
    }

    dos_weights.resize(dos_weights.size(0)+1, weights.size());

    for(size_t i=0; i<dos_weights.size(1); i++) {
        dos_weights(dos_weights.size(0)-1, i) = weights[i];
    }
}

void CondDOS::set_sigma(CPAPtr sigma) { _sigma = sigma; }
CPAPtr CondDOS::sigma() const { return _sigma; }

void CondDOS::solve() {
    _status = NodeStatus::failed;

    //-----------------------------------------------------------------------

    if (_sigma == nullptr) throw runtime_error("No Model");

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

    // Will be set later inside the parallel loop
    size_t threads = 1;
    math::NArray<double,6> red_cond_state;
    math::NArray<double,3> red_dos;
    math::NArray<double,4> red_dos_state;

    //-----------------------------------------------------------------------
    // Prepare data buffers

    if(solve_cond) {
        cond.resize(c_size, w_size, 3, 3);
        cond.fill(0);
        cond_state.resize(c_size, w_size, 3, 3, states);
        cond_state.fill(0);
    } else {
        cond.clear();
        cond_state.clear();
    }

    if(solve_dos) {
        dos.resize(c_size, w_size);
        dos.fill(0);
    } else {
        dos.clear();
    }

    if(solve_dos_k) {
        dos_k.resize(c_size, w_size, k1.size(), k2.size(), k3.size());
    } else {
        dos_k.clear();
    }

    if(solve_dos_state) {
        dos_state.resize(c_size, w_size, states);
        dos_state.fill(0);
    } else {
        dos_state.clear();
    }

    bool solve_dos_proj = !dos_weights.empty();
    if(solve_dos_proj) {

        if(dos_weights.size(1) != states) {
            error("Incorrect size for dos_weights");
        }

        dos_proj.resize(c_size, w_size, dos_weights.size(0), k1.size(), k2.size(), k3.size());
        dos_proj.fill(0);

    } else {
        dos_proj.clear();
    }

    if(_parallel) {
        // Disable Eigen Multi-threading
        Eigen::setNbThreads(1);
    } else {
        // Let Eigen determine number of threads to use;
        Eigen::setNbThreads(0);
    }

    #pragma omp parallel if(_parallel)
    {
        //-----------------------------------------------------------------------
        // Reduction buffers

        #pragma omp master
        {
            threads = omp_get_num_threads();

            if(solve_cond) {
                red_cond_state.resize(threads, c_size, w_size, 3, 3, states);
                red_cond_state.fill(0);
            }
            if(solve_dos) {
                red_dos.resize(threads, c_size, w_size);
                red_dos.fill(0);
            }
            if(solve_dos_state) {
                red_dos_state.resize(threads, c_size, w_size, states);
                red_dos_state.fill(0);
            }
        }

        #pragma omp barrier

        //-----------------------------------------------------------------------
        // Temporary data buffers for each thread

        // unique thread identifier
        size_t it = omp_get_thread_num();

        math::DenseMatrix<cxdouble> mGk(states, states); //< G(k) = Inv(w*I - T - S)
        math::DenseMatrix<cxdouble> mP(states, states);  //< P = 1/(2*pi)(conj(Gk.T) - Gk)

        math::DenseMatrix<cxdouble> mvx(states, states);
        math::DenseMatrix<cxdouble> mvy(states, states);
        math::DenseMatrix<cxdouble> mvz(states, states);

        math::DenseMatrix<cxdouble> mvxP(states, states); //< vx*P
        math::DenseMatrix<cxdouble> mvyP(states, states); //< vy*P
        math::DenseMatrix<cxdouble> mvzP(states, states); //< vz*P

        cxdouble *__restrict__ Gk = mGk.data();

        cxdouble *__restrict__ vx = &mvx(0,0);
        cxdouble *__restrict__ vy = &mvy(0,0);
        cxdouble *__restrict__ vz = &mvz(0,0);

        cxdouble *p1 = nullptr; //< General purpose pointer

        //-----------------------------------------------------------------------

        #pragma omp for collapse(3)
        for(size_t ik1=0; ik1<k1.size(); ik1++)
        for(size_t ik2=0; ik2<k2.size(); ik2++)
        for(size_t ik3=0; ik3<k3.size(); ik3++) {
            //-----------------------------------------------------------
            // Get Velocity Matrix

            if(solve_cond) {
                _model->dH_dk(vx, vy, vz, k1[ik1], k2[ik2], k3[ik3]);
            }

            for(size_t ic=0; ic<c_size; ic++)
            for(size_t iw=0; iw<w_size; iw++) {


                //---------------------------------------------------------------
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
                    // if(abs(imag(v)) < zero_j) {
                    //     v = cxdouble(real(v)+_sigma->w[iw], zero_j);
                    // } else {
                        v += _sigma->w[iw];
                    // }
                }

                // Now: Gk = w*I - T - S

                mGk = mGk.inverse();

                // Now: Gk = Inv(w*I - T - S)

                //---------------------------------------------------------------
                // Cond

                if(solve_cond) {

                    //-----------------------------------------------------------
                    // Compute full P

                    mP = cxdouble(0, 1/(2*data::pi))*(mGk.adjoint() - mGk);

                    //-----------------------------------------------------------
                    // Compute mviP

                    if(solve_sigx) mvxP = mvx*mP;
                    if(solve_sigy) mvyP = mvy*mP;
                    if(solve_sigz) mvzP = mvz*mP;

                    //-----------------------------------------------------------
                    // Compute cond_ij

                    if(solve_sigx) {
                        // solve_sigxx
                        for(size_t s=0; s<states; s++) {
                            double v = 0;
                            for(size_t r=0; r<states; r++) {
                                 v += real(mvxP(s,r)*mvxP(r,s));
                            }
                            red_cond_state(it,ic,iw,0,0,s) += v;
                        }
                        if(solve_sigy) {
                            // solve_sigxy
                            for(size_t s=0; s<states; s++) {
                                double v = 0;
                                for(size_t r=0; r<states; r++) {
                                     v += real(mvxP(s,r)*mvyP(r,s));
                                }
                                red_cond_state(it,ic,iw,0,1,s) += v;
                            }
                            // solve_sigyx
                            for(size_t s=0; s<states; s++) {
                                double v = 0;
                                for(size_t r=0; r<states; r++) {
                                     v += real(mvyP(s,r)*mvxP(r,s));
                                }
                                red_cond_state(it,ic,iw,1,0,s) += v;
                            }
                        }
                        if(solve_sigz) {
                            // solve_sigxz
                            for(size_t s=0; s<states; s++) {
                                double v = 0;
                                for(size_t r=0; r<states; r++) {
                                     v += real(mvxP(s,r)*mvzP(r,s));
                                }
                                red_cond_state(it,ic,iw,0,2,s) += v;
                            }
                            // solve_sigzx
                            for(size_t s=0; s<states; s++) {
                                double v = 0;
                                for(size_t r=0; r<states; r++) {
                                     v += real(mvzP(s,r)*mvxP(r,s));
                                }
                                red_cond_state(it,ic,iw,2,0,s) += v;
                            }
                        }
                    }

                    if(solve_sigy) {
                        // solve_sigyy
                        for(size_t s=0; s<states; s++) {
                            double v = 0;
                            for(size_t r=0; r<states; r++) {
                                 v += real(mvyP(s,r)*mvyP(r,s));
                            }
                            red_cond_state(it,ic,iw,1,1,s) += v;
                        }
                        if(solve_sigz) {
                            // solve_sigyz
                            for(size_t s=0; s<states; s++) {
                                double v = 0;
                                for(size_t r=0; r<states; r++) {
                                     v += real(mvyP(s,r)*mvzP(r,s));
                                }
                                red_cond_state(it,ic,iw,1,2,s) += v;
                            }
                            // solve_sigzy
                            for(size_t s=0; s<states; s++) {
                                double v = 0;
                                for(size_t r=0; r<states; r++) {
                                     v += real(mvzP(s,r)*mvyP(r,s));
                                }
                                red_cond_state(it,ic,iw,2,1,s) += v;
                            }
                        }
                    }

                    if(solve_sigz) {
                        // solve_sigzz
                        for(size_t s=0; s<states; s++) {
                            double v = 0;
                            for(size_t r=0; r<states; r++) {
                                 v += real(mvzP(s,r)*mvzP(r,s));
                            }
                            red_cond_state(it,ic,iw,2,2,s) += v;
                        }
                    }

                    //-----------------------------------------------------------

                } else {
                    // If no cond then DOS only requires diagonal terms and no vx,vy,vz
                    for(size_t i=0; i<states; i++) {
                        mP(i,i) = cxdouble(0,1/(2*data::pi))*(conj(mGk(i,i)) - mGk(i,i));
                    }
                }

                //---------------------------------------------------------------
                // DOS

                if(solve_dos or solve_dos_k or solve_dos_state) {
                    double trP = 0;
                    for(size_t i=0; i<states; i++) {
                        double v = -real(mP(i,i));
                        trP += v;
                        if(solve_dos_state) {
                            red_dos_state(it,ic,iw,i) += v;
                        }
                        if(solve_dos_proj) {
                            for(size_t ip=0; ip<dos_weights.size(0); ip++) {
                                dos_proj(ic,iw,ip, ik1,ik2,ik3) += v*dos_weights(ip,i);
                            }
                        }
                    }
                    if(solve_dos)
                        red_dos(it,ic,iw) += trP;
                    if(solve_dos_k)
                        dos_k(ic,iw,ik1,ik2,ik3) = trP;
                }

                //---------------------------------------------------------------
            }
        }
    }


    if (solve_cond) {
        double cond_norm = data::pi*pow(data::ec,2)/(
                data::hbar*kpoints*_model->uc_size());

        // Normalize cond_state and compute cond
        for(size_t ic=0; ic<c_size; ic++)
        for(size_t iw=0; iw<w_size; iw++)
        for(size_t i=0; i<3; i++)
        for(size_t j=0; j<3; j++)
        for(size_t s=0; s<states; s++) {
            for(size_t it=0; it<threads; it++) {
                cond_state(ic,iw,i,j,s) += red_cond_state(it,ic,iw,i,j,s);
            }
            cond_state(ic,iw,i,j,s) *= cond_norm;
            cond(ic,iw,i,j) += cond_state(ic,iw,i,j,s);
        }
    }

    if (solve_dos) {
        const double dos_norm = 1.0f/(kpoints);

        // Normalize dos
        for(size_t ic=0; ic<c_size; ic++)
        for(size_t iw=0; iw<w_size; iw++) {
            for(size_t it=0; it<threads; it++) {
                dos(ic,iw) += red_dos(it,ic,iw);
            }
            dos(ic,iw) *= dos_norm;
        }
    }

    if(solve_dos_state) {
        const double dos_norm = 1.0f/(kpoints);

        // Normalize dos_state
        for(size_t ic=0; ic<c_size; ic++)
        for(size_t iw=0; iw<w_size; iw++)
        for(size_t is=0; is<states; is++) {
            for(size_t it=0; it<threads; it++) {
                dos_state(ic,iw,is) += red_dos_state(it,ic,iw,is);
            }
            dos_state(ic,iw,is) *= dos_norm;
        }
    }

    red_dos.clear();
    red_cond_state.clear();
    red_dos_state.clear();

    //-----------------------------------------------------------------------

    _status = NodeStatus::done;
}

#ifdef TBPP_WITH_HDF5
void CondDOS::save(EHFile& file, const string& prefix) const {
    Node::save(file, prefix);
    if(_sigma != nullptr)
        file.set_attr(prefix, "sigma", _sigma->name());
    else
        file.set_attr(prefix, "sigma", "");

    file.set_attr(prefix, "solve_cond", solve_cond);
    file.set_attr(prefix, "solve_sigx", solve_sigx);
    file.set_attr(prefix, "solve_sigy", solve_sigy);
    file.set_attr(prefix, "solve_sigz", solve_sigz);
    file.set_attr(prefix, "solve_dos", solve_dos);
    file.set_attr(prefix, "solve_dos_state", solve_dos_state);
    file.set_attr(prefix, "solve_dos_k", solve_dos_k);
    file.set_attr(prefix, "zero_j", zero_j);
    file.set_data(prefix+"/dos_weights", dos_weights);

    if(_status == NodeStatus::done) {
        file.set_data(prefix+"/cond", cond);
        file.set_data(prefix+"/cond_state", cond_state);
        file.set_data(prefix+"/dos_k", dos_k);
        file.set_data(prefix+"/dos", dos);
        file.set_data(prefix+"/dos_state", dos_state);
        file.set_data(prefix+"/dos_proj", dos_proj);
    }
}
void CondDOS::load(EHFile& file, const string& prefix) {
    Node::load(file, prefix);

    string sigma_name;
    file.get_attr(prefix, "sigma", sigma_name);
    if(sigma_name.length() > 0 and context() != nullptr)
        _sigma = dynamic_pointer_cast<CPA>(context()->get(sigma_name));
    else
        _sigma = nullptr;

    file.get_attr(prefix, "solve_cond", solve_cond);
    file.get_attr(prefix, "solve_sigx", solve_sigx);
    file.get_attr(prefix, "solve_sigy", solve_sigy);
    file.get_attr(prefix, "solve_sigz", solve_sigz);
    file.get_attr(prefix, "solve_dos_k", solve_dos_k);
    file.get_attr(prefix, "solve_dos", solve_dos);
    file.get_attr(prefix, "zero_j", zero_j);
    file.get_attr(prefix, "solve_dos_state", solve_dos_state);
    file.get_data(prefix+"/dos_weights", dos_weights);

    if(_status == NodeStatus::done) {
        file.get_data(prefix+"/cond", cond);
        file.get_data(prefix+"/cond_state", cond_state);
        file.get_data(prefix+"/dos_k", dos_k);
        file.get_data(prefix+"/dos", dos);
        file.get_data(prefix+"/dos_state", dos_state);
        file.get_data(prefix+"/dos_proj", dos_proj);
    }
}
#endif // TBPP_WITH_HDF5

//----------------------------------------------------------------------

} // namespace tbpp
