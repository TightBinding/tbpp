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

#include <tbpp/approx_cond_dos.h>

using namespace std;

namespace tbpp {

//----------------------------------------------------------------------

const char* ApproxCondDOS::type() const { return "ApproxCondDOS"; }

void ApproxCondDOS::set_w(std::vector<double> _w) { w = _w; }
void ApproxCondDOS::set_eta(std::vector<double> _eta) { eta = _eta; }

void ApproxCondDOS::set_kgrid(KGridPtr kgrid) { _kgrid = kgrid; }
KGridPtr ApproxCondDOS::kgrid() const { return _kgrid; }

void ApproxCondDOS::solve() {
    _status = NodeStatus::failed;

    //-----------------------------------------------------------------------

    if (_kgrid == nullptr) throw runtime_error("No KGrid");

    KModelPtr _kmodel = _kgrid->kmodel();
    if (_kmodel == nullptr) throw runtime_error("No KModel");

    ModelPtr _model = _kmodel->model();
    if (_model == nullptr) throw runtime_error("No Model");

    _model->make_ready();
    _kmodel->make_ready();
    _kgrid->make_ready();

    header();
    //-----------------------------------------------------------------------

    const size_t states = _model->states();
    const size_t sites = _model->sites();
    const size_t kdim = _model->kdim();
    const size_t states2 = states*states;

    const vector<double>& k1 = _kmodel->k1;
    const vector<double>& k2 = _kmodel->k2;
    const vector<double>& k3 = _kmodel->k3;

    const size_t kpoints = k1.size()*k2.size()*k3.size();
    const size_t w_size = w.size();
    const size_t eta_size = eta.size();

    //-----------------------------------------------------------------------
    // Check k-space
    if(kdim >= 1 and k1.size() < 2) error("Incorrect dimension for k1");
    if(kdim >= 2 and k2.size() < 2) error("Incorrect dimension for k2");
    if(kdim >= 3 and k3.size() < 2) error("Incorrect dimension for k3");

    //-----------------------------------------------------------------------
    // Prepare data buffers

    if(solve_cond) {
        cond.resize(eta_size, w_size, 3, 3);
        cond.fill(0);
    } else {
        cond.clear();
    }
    if(solve_dos) {
        dos.resize(eta_size, w_size);
        dos.fill(0);
    } else {
        dos.clear();
    }
    if(solve_dos_k) {
        dos_k.resize(eta_size, w_size, k1.size(), k2.size(), k3.size());
    } else {
        dos_k.clear();
    }

    #pragma omp parallel if(_parallel)
    {
        double dk1 = 1;
        double dk2 = 1;
        double dk3 = 1;
        if(kdim >= 1) dk1 = 2*(k1[1] - k1[0]);
        if(kdim >= 2) dk2 = 2*(k2[1] - k2[0]);
        if(kdim >= 3) dk3 = 2*(k3[1] - k3[0]);

        double v[3] = {0,0,0};  //< velocity with respect to lattice vectors
        double vr[3] = {0,0,0}; //< velocity with respect to x,y,z

        #pragma omp for collapse(2)
        for(size_t ieta=0; ieta<eta_size; ieta++)
        for(size_t iw=0; iw<w_size; iw++) {
            for(size_t ik1=0; ik1<k1.size(); ik1++)
            for(size_t ik2=0; ik2<k2.size(); ik2++)
            for(size_t ik3=0; ik3<k3.size(); ik3++) {
                int ik1p = ik1-1;
                int ik1n = ik1+1;
                int ik2p = ik2-1;
                int ik2n = ik2+1;
                int ik3p = ik3-1;
                int ik3n = ik3+1;
                if(ik1p < 0) ik1p += k1.size();
                if(ik1n >= (int)k1.size()) ik1n -= k1.size();
                if(ik2p < 0) ik2p += k2.size();
                if(ik2n >= (int)k2.size()) ik2n -= k2.size();
                if(ik3p < 0) ik3p += k3.size();
                if(ik3n >= (int)k3.size()) ik3n -= k3.size();

                double val = 0;
                double val_sum = 0;

                for(unsigned iE=0; iE < _kgrid->eigval.size(3); iE++) {
                    val = eta[ieta]/(data::pi*(pow(w[iw] - _kgrid->eigval(ik1,ik2,ik3,iE),2) + pow(eta[ieta],2)));
                    val_sum += val;

                    if(solve_cond) {
                        v[0] = v[1] = v[2] = 0;
                        if (kdim >= 1)
                            v[0] = (_kgrid->eigval(ik1n,ik2,ik3,iE) - _kgrid->eigval(ik1p,ik2,ik3,iE))/dk1;
                        if (kdim >= 2)
                            v[1] = (_kgrid->eigval(ik1,ik2n,ik3,iE) - _kgrid->eigval(ik1,ik2p,ik3,iE))/dk2;
                        if (kdim >= 3)
                            v[2] = (_kgrid->eigval(ik1,ik2,ik3n,iE) - _kgrid->eigval(ik1,ik2,ik3p,iE))/dk3;
                        vr[0] = (v[0]*_model->a(0,0) + v[1]*_model->a(1,0) + v[2]*_model->a(2,0))/(2*data::pi);
                        vr[1] = (v[0]*_model->a(0,1) + v[1]*_model->a(1,1) + v[2]*_model->a(2,1))/(2*data::pi);
                        vr[2] = (v[0]*_model->a(0,2) + v[1]*_model->a(1,2) + v[2]*_model->a(2,2))/(2*data::pi);
                        for(unsigned i=0; i<3; i++)
                        for(unsigned j=0; j<3; j++) {
                            cond(ieta, iw, i, j) += vr[i]*vr[j]*val;
                        }
                    }
                }

                if(solve_dos_k)
                    dos_k(ieta, iw, ik1,ik2,ik3) = val_sum;
                if(solve_dos)
                    dos(ieta, iw) +=  val_sum;
            }
        }
    }


    if (solve_cond) {
        for(size_t ieta=0; ieta<eta_size; ieta++) {
            double cond_norm = pow(data::ec,2)/(2*data::hbar*eta[ieta]*kpoints*_model->uc_size());

            // Normalize cond_state and compute cond
            for(size_t iw=0; iw<w_size; iw++)
            for(size_t i=0; i<3; i++)
            for(size_t j=0; j<3; j++) {
                cond(ieta,iw,i,j) *= cond_norm;
            }
        }
    }

    if (solve_dos) {
        for(size_t ieta=0; ieta<eta_size; ieta++) {
            const double dos_norm = 1.0f/kpoints;

            // Normalize dos
            for(size_t iw=0; iw<w_size; iw++) {
                dos(ieta, iw) *= dos_norm;
            }
        }
    }

    //-----------------------------------------------------------------------

    _status = NodeStatus::done;
}

#ifdef TBPP_WITH_HDF5
void ApproxCondDOS::save(EHFile& file, const string& prefix) const {
    Node::save(file, prefix);
    if(_kgrid != nullptr)
        file.set_attr(prefix, "kgrid", _kgrid->name());
    else
        file.set_attr(prefix, "kgrid", "");

    file.set_attr(prefix, "solve_cond", solve_cond);
    file.set_attr(prefix, "solve_dos_k", solve_dos_k);
    file.set_attr(prefix, "solve_dos", solve_dos);
    file.set_data(prefix+"/w", w);
    file.set_data(prefix+"/eta", eta);

    if(_status == NodeStatus::done) {
        file.set_data(prefix+"/cond", cond);
        file.set_data(prefix+"/dos_k", dos_k);
        file.set_data(prefix+"/dos", dos);
    }
}
void ApproxCondDOS::load(EHFile& file, const string& prefix) {
    Node::load(file, prefix);

    string kgrid_name;
    file.get_attr(prefix, "kgrid", kgrid_name);
    if(kgrid_name.length() > 0 and context() != nullptr)
        _kgrid = dynamic_pointer_cast<KGrid>(context()->get(kgrid_name));
    else
        _kgrid = nullptr;

    file.get_attr(prefix, "solve_cond", solve_cond);
    file.get_attr(prefix, "solve_dos_k", solve_dos_k);
    file.get_attr(prefix, "solve_dos", solve_dos);
    file.get_data(prefix+"/w", w);
    file.get_data(prefix+"/eta", eta);

    if(_status == NodeStatus::done) {
        file.get_data(prefix+"/cond", cond);
        file.get_data(prefix+"/dos_k", dos_k);
        file.get_data(prefix+"/dos", dos);
    }
}
#endif // TBPP_WITH_HDF5


//----------------------------------------------------------------------

} // namespace tbpp
