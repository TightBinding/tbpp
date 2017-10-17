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

#include <iostream>
#include <tbpp/model.h>
#include <cmath>

using namespace std;

namespace tbpp {

//----------------------------------------------------------------------

const char* Model::type() const { return "Model"; }

Model::Model() {
    a.resize(3,3);
    a.fill(0);
    a(0,0) = a(1,1) = a(2,2) = 1;
}

//------------------------------------------------------------------------
// Input type

void Model::set_use_lattice(bool use_lattice) {
    _use_lattice = use_lattice;
    if(not use_lattice) {
        _lattice = nullptr;
    }
}

bool Model::use_lattice() const {
    return _use_lattice;
}

//------------------------------------------------------------------------
// Lattice

bool Model::set_lattice(LatticePtr lattice) {
    if(lattice != nullptr) {
        _use_lattice = true;
        _lattice = lattice;
        return true;
    }
    return false;
}

LatticePtr Model::lattice() const {
    return _lattice;
}

//------------------------------------------------------------------------
// General Settings

unsigned Model::kdim() const {
    return _kdim;
}

bool Model::set_kdim(unsigned kdim) {
    if (kdim > 3)
        return false;
    _kdim = kdim;
    return true;
}

unsigned Model::sites() const {
    return _sites;
}

void Model::set_sites(unsigned sites) {
    _sites = sites;
}

unsigned Model::states() const {
    return _states;
}

void Model::set_states(unsigned states) {
    _states = states;
}

unsigned Model::hops() const {
    return _hops;
}

void Model::set_hops(unsigned hops) {
    _hops = hops;
}

void Model::set_a1(double x, double y, double z) {
    a(0,0) = x;
    a(0,1) = y;
    a(0,2) = z;
}

void Model::set_a2(double x, double y, double z) {
    a(1,0) = x;
    a(1,1) = y;
    a(1,2) = z;
}

void Model::set_a3(double x, double y, double z) {
    a(2,0) = x;
    a(2,1) = y;
    a(2,2) = z;
}

void Model::solve() {
    _status = NodeStatus::failed;
    //-----------------------------------------------------------------------

    if(_use_lattice) {
        if (_lattice == nullptr) error("No Lattice");
        _lattice->make_ready();
    }

    header();

    //-----------------------------------------------------------------------

    if(_use_lattice) {
        _kdim = _lattice->trans.size();
        _sites = _lattice->sites.size();
        _states = 0;
        for(const auto& s : _lattice->sites) {
            _states += s.onsite.rows();
        }
    }

    //-----------------------------------------------------------------------
    // site_info and state_info

    site_info.resize(_sites);
    state_info.resize(_states);

    if(_use_lattice) {
        size_t index=0;
        for(size_t i=0; i<_sites; i++) {

            const Lattice::Site& site = _lattice->sites[i];
            for(size_t j=0; j<TBPP_MAX_SITE_KIND_LEN; j++)
                site_info[i].kind[j] = 0;
            site.kind.copy(site_info[i].kind, TBPP_MAX_SITE_KIND_LEN-1);
            site_info[i].pos[0] = site.pos[0];
            site_info[i].pos[1] = site.pos[1];
            site_info[i].pos[2] = site.pos[2];
            site_info[i].states = site.orbits.size();
            site_info[i].si = index;

            for(size_t j=0; j<site_info[i].states; j++) {
                state_info[index].site = i;
                state_info[index].orbit = site.orbits[j];
                state_info[index].spinor = site.spinors[j];
                index += 1;
            }
            site_info[i].ei = index-1;
        }
    }

    //-----------------------------------------------------------------------
    // a

    if(_use_lattice) {
        for(size_t i=0; i<_kdim; i++) {
            a(i,0) = _lattice->trans[i][0];
            a(i,1) = _lattice->trans[i][1];
            a(i,2) = _lattice->trans[i][2];
        }
    }

    if (_kdim == 0) {
        a.fill(0);
        a(0,0) = a(1,1) = a(2,2) = 1;
    }
    if (_kdim == 1) {
        // gram-schmidt process for a2
        array<double,3> n = {1,0,0};
        if((a(0,1) == 0) and (a(0,2) == 0))
            n = {0,1,0};

        double d = n[0]*a(0,0) + n[1]*a(0,1) + n[2]*a(0,2);
        double da1 = a(0,0)*a(0,0) + a(0,1)*a(0,1) + a(0,2)*a(0,2);
        // a2 = n - dot(n,a1)*a1/dot(a1,a1)
        a(1,0) = n[0] - d*a(0,0)/da1;
        a(1,1) = n[1] - d*a(0,1)/da1;
        a(1,2) = n[2] - d*a(0,2)/da1;

        // normalize a2
        double na2 = sqrt(a(1,0)*a(1,0) + a(1,1)*a(1,1) + a(1,2)*a(1,2));
        a(1,0) /= na2;
        a(1,1) /= na2;
        a(1,2) /= na2;

        // a3 will be set below
    }
    if (_kdim == 2 or _kdim == 1) {
        // a3 = cross(a1, a2)
        a(2,0) =  a(0,1)*a(1,2) - a(1,1)*a(0,2);
        a(2,1) = -a(0,0)*a(1,2) + a(1,0)*a(0,2);
        a(2,2) =  a(0,0)*a(1,1) - a(1,0)*a(0,1);

        // normalize a3
        double na3 = sqrt(a(2,0)*a(2,0) + a(2,1)*a(2,1) + a(2,2)*a(2,2));
        a(2,0) /= na3;
        a(2,1) /= na3;
        a(2,2) /= na3;
    }
    if (_kdim > 3) {
        error("Too many translational symmetries");
    }

    //-----------------------------------------------------------------------
    // b

    b.resize(3,3);
    double twopi_V = 2*data::pi/uc_size();
    double* a1 = &a(0,0);
    double* a2 = &a(1,0);
    double* a3 = &a(2,0);

    // b1 = cross(a2,a3)*twopi_V
    b(0,0) = ( a2[1]*a3[2] - a3[1]*a2[2])*twopi_V;
    b(0,1) = (-a2[0]*a3[2] + a3[0]*a2[2])*twopi_V;
    b(0,2) = ( a2[0]*a3[1] - a3[0]*a2[1])*twopi_V;

    // b2 = cross(a3,a1)*twopi_V
    b(1,0) = ( a3[1]*a1[2] - a1[1]*a3[2])*twopi_V;
    b(1,1) = (-a3[0]*a1[2] + a1[0]*a3[2])*twopi_V;
    b(1,2) = ( a3[0]*a1[1] - a1[0]*a3[1])*twopi_V;

    // b3 = cross(a1,a2)*twopi_V
    b(2,0) = ( a1[1]*a2[2] - a2[1]*a1[2])*twopi_V;
    b(2,1) = (-a1[0]*a2[2] + a2[0]*a1[2])*twopi_V;
    b(2,2) = ( a1[0]*a2[1] - a2[0]*a1[1])*twopi_V;

    //-----------------------------------------------------------------------
    // On site matrix

    V.resize(_states,_states);
    V.fill(0);

    if(_use_lattice) {
        for(size_t site_i=0; site_i<_sites; site_i++)
        for(size_t site_j=0; site_j<_sites; site_j++) {
            size_t si = site_info[site_i].si;

            if(site_i == site_j) {
                for(size_t r=0; r<site_info[site_i].states; r++)
                for(size_t c=0; c<site_info[site_i].states; c++) {
                    V(si+r,si+c) = _lattice->sites[site_i].onsite(r,c);
                }
            } else {
                size_t sj = site_info[site_j].si;

                double dx = site_info[site_j].pos[0] - site_info[site_i].pos[0];
                double dy = site_info[site_j].pos[1] - site_info[site_i].pos[1];
                double dz = site_info[site_j].pos[2] - site_info[site_i].pos[2];

                double px = site_info[site_j].pos[0] + site_info[site_i].pos[0];
                double py = site_info[site_j].pos[1] + site_info[site_i].pos[1];
                double pz = site_info[site_j].pos[2] + site_info[site_i].pos[2];

                // Phase due to magnetic field
                cxdouble exp_phi_b = exp(cxdouble(0, data::ec/(2*data::hbar)*(
                                _lattice->B[1]*pz*dx +
                                _lattice->B[2]*px*dy +
                                _lattice->B[0]*py*dz)));

                math::DenseMatrix<cxdouble> T = _lattice->_T(site_i,site_j,dx,dy,dz);
                if(T.size() > 0) {
                    for(size_t r=0; r<T.rows(); r++)
                    for(size_t c=0; c<T.cols(); c++) {
                        V(sj+r, si+c) = T(r,c)*exp_phi_b;
                        V(si+c, sj+r) = conj(T(r,c))*conj(exp_phi_b);
                    }
                }
            }
        }
    }

    //-----------------------------------------------------------------------
    // Hop matrices

    if(not _use_lattice) {
        T.resize(_hops, _states, _states);
        T.fill(0);

        R.resize(_hops, 3);
        R.fill(0);
    } else {
        _hops = 0;
        T.resize(_hops, _states, _states);
        R.resize(_hops, 3);
        if (_kdim > 0) {
            // hop index range
            array<int,3> n_range = {0,0,0};
            if(_lattice->max_hop_range > 0.0) {
                for(int i=0; i<_kdim; i++) {
                    n_range[i] = ceil(_lattice->max_hop_range/sqrt(a(i,0)*a(i,0)
                                + a(i,1)*a(i,1) + a(i,2)*a(i,2)));
                }
            }

            for(int i=-n_range[0]; i<=n_range[0]; i++)
            for(int j=-n_range[1]; j<=n_range[1]; j++)
            for(int k=-n_range[2]; k<=n_range[2]; k++) {
                if(i == 0 and j == 0 and k == 0) continue;
                bool resize = true;

                double Rx = i*a(0,0) + j*a(1,0) + k*a(2,0);
                double Ry = i*a(0,1) + j*a(1,1) + k*a(2,1);
                double Rz = i*a(0,2) + j*a(1,2) + k*a(2,2);

                for(size_t site_i=0; site_i<_sites; site_i++)
                for(size_t site_j=0; site_j<_sites; site_j++) {
                    double dx = site_info[site_j].pos[0] + Rx - site_info[site_i].pos[0];
                    double dy = site_info[site_j].pos[1] + Ry - site_info[site_i].pos[1];
                    double dz = site_info[site_j].pos[2] + Rz - site_info[site_i].pos[2];

                    math::DenseMatrix<cxdouble> T_sub = _lattice->_T(site_i,site_j,dx,dy,dz);
                    if(T_sub.size() > 0) {
                        if(resize) {
                            // append new {i,j,k} hop matrix and set resize to
                            _hops += 1;
                            T.resize(_hops, _states, _states);
                            R.resize(_hops, 3);
                            resize = false;
                        }

                        size_t si = site_info[site_i].si;
                        size_t sj = site_info[site_j].si;

                        double px = site_info[site_j].pos[0] + Rx + site_info[site_i].pos[0];
                        double py = site_info[site_j].pos[1] + Ry + site_info[site_i].pos[1];
                        double pz = site_info[site_j].pos[2] + Rz + site_info[site_i].pos[2];

                        // Phase due to magnetic field
                        cxdouble exp_phi_b = exp(cxdouble(0, data::ec/(2*data::hbar)*(
                                        _lattice->B[1]*pz*dx +
                                        _lattice->B[2]*px*dy +
                                        _lattice->B[0]*py*dz)));

                        // copy sub matrix
                        for(size_t r=0; r<T_sub.rows(); r++)
                        for(size_t c=0; c<T_sub.rows(); c++) {
                            T(_hops-1, sj+r, si+c) = T_sub(r,c)*exp_phi_b;
                        }
                        R(_hops-1, 0) = i;
                        R(_hops-1, 1) = j;
                        R(_hops-1, 2) = k;
                    }
                }
            }
        }
    }

    //-----------------------------------------------------------------------
    // Print Info

    cout << "kdim = " << _kdim << '\n';
    cout << "sites = " << _sites << '\n';
    cout << "states = " << _states << '\n';
    cout << "hops = " << _hops << "\n\n";

    for(size_t i=0; i<3; i++) {
        cout << "a" << i+1 << " = {" << a(i,0) << ", " << a(i,1) << ", " << a(i,2) << "}\n";
    }
    cout << '\n';

    for(size_t i=0; i<3; i++) {
        cout << "b" << i+1 << " = {" << b(i,0) << ", " << b(i,1) << ", " << b(i,2) << "}\n";
    }
    cout << '\n';

    // for(size_t i=0; i<_sites; i++) {
    //     const auto& site = site_info[i];
    //     cout << "Site " << i << '\n'
    //          << "    " << "kind = " << site.kind << '\n'
    //          << "    " << "states = " << site.states << '\n'
    //          << "    " << "start_index = " << site.si << '\n'
    //          << "    " << "end_index = " << site.ei << '\n'
    //          << "    " << "pos = {"
    //          << site.pos[0] << ", " << site.pos[1] << ", " << site.pos[2] << "}\n"
    //          << "\n\n";
    // }

    cout << "[" << type() << " " << name() << "] Done \n\n";

    //-----------------------------------------------------------------------

    _status = tbpp::NodeStatus::done;
}

double Model::uc_size() const {
    if(a.size(0) != 3 or a.size(1) != 3) {
        error("incorrect size for lattice vectors matrix `a`");
    }

    // Below is: abs(dot(a3, cross(a1,a2))
    double x =  a(0,1)*a(1,2) - a(1,1)*a(0,2);
    double y = -a(0,0)*a(1,2) + a(1,0)*a(0,2);
    double z =  a(0,0)*a(1,1) - a(1,0)*a(0,1);
    return abs(a(2,0)*x + a(2,1)*y + a(2,2)*z);
}

void Model::Tk(complex<double> *d, double k1, double k2, double k3) const {

    switch (_kdim) {
        case 0: k1=k2=k3=0; break;
        case 1:    k2=k3=0; break;
        case 2:       k3=0; break;
        case 3:             break;
        default: error("Invalid value for kdim");
    }

    // k-point in Cartesian coordinates
    double kx = k1*b(0,0) + k2*b(1,0) + k3*b(2,0);
    double ky = k1*b(0,1) + k2*b(1,1) + k3*b(2,1);
    double kz = k1*b(0,2) + k2*b(1,2) + k3*b(2,2);

    // initialize data with zeros
    size_t s2 = _states*_states;
    for(size_t i=0; i<s2; i++) {
        d[i] = 0.0;
    }

    // add hops
    for(size_t iR=0; iR<_hops; iR++) {
        double R_dot_k = 2*M_PI*(R(iR,0)*k1 + R(iR,1)*k2 + R(iR,2)*k3);

        for(size_t si=0; si<_sites; si++)
        for(size_t sj=0; sj<_sites; sj++) {
            // add hop from si to sj
            complex<double> exp_I_phi = exp(complex<double>(0, R_dot_k
                        + (site_info[sj].pos[0] - site_info[si].pos[0])*kx
                        + (site_info[sj].pos[1] - site_info[si].pos[1])*ky
                        + (site_info[sj].pos[2] - site_info[si].pos[2])*kz));

            for(size_t oi=site_info[si].si; oi<=site_info[si].ei; oi++)
            for(size_t oj=site_info[sj].si; oj<=site_info[sj].ei; oj++) {
                // add orbit sub-matrix times phase
                d[oi*_states + oj] += T(iR, oi, oj)*exp_I_phi;
            }
        }
    }
}

void Model::Hk(complex<double> *d, double k1, double k2, double k3) const {

    switch (_kdim) {
        case 0: k1=k2=k3=0; break;
        case 1:    k2=k3=0; break;
        case 2:       k3=0; break;
        case 3:             break;
        default: error("Invalid value for kdim");
    }

    // k-point in Cartesian coordinates
    double kx = k1*b(0,0) + k2*b(1,0) + k3*b(2,0);
    double ky = k1*b(0,1) + k2*b(1,1) + k3*b(2,1);
    double kz = k1*b(0,2) + k2*b(1,2) + k3*b(2,2);

    // initialize data with zeros
    size_t s2 = _states*_states;
    for(size_t i=0; i<s2; i++) {
        d[i] = 0.0;
    }

    // initialize data with on-site
    for(size_t si=0; si<_sites; si++)
    for(size_t sj=0; sj<_sites; sj++) {
        // add hop from si to sj with R=0
        complex<double> exp_I_phi = exp(complex<double>(0,
                      (site_info[sj].pos[0] - site_info[si].pos[0])*kx
                    + (site_info[sj].pos[1] - site_info[si].pos[1])*ky
                    + (site_info[sj].pos[2] - site_info[si].pos[2])*kz));

        for(size_t oi=site_info[si].si; oi<=site_info[si].ei; oi++)
        for(size_t oj=site_info[sj].si; oj<=site_info[sj].ei; oj++) {
            // add orbit sub-matrix times phase
            d[oj*_states + oi] = V(oj, oi)*exp_I_phi;
        }
    }

    // add hops
    for(size_t iR=0; iR<_hops; iR++) {
        double R_dot_k = 2*M_PI*(R(iR,0)*k1 + R(iR,1)*k2 + R(iR,2)*k3);

        for(size_t si=0; si<_sites; si++)
        for(size_t sj=0; sj<_sites; sj++) {
            // add hop from si to sj
            complex<double> exp_I_phi = exp(complex<double>(0, R_dot_k
                        + (site_info[sj].pos[0] - site_info[si].pos[0])*kx
                        + (site_info[sj].pos[1] - site_info[si].pos[1])*ky
                        + (site_info[sj].pos[2] - site_info[si].pos[2])*kz));

            for(size_t oi=site_info[si].si; oi<=site_info[si].ei; oi++)
            for(size_t oj=site_info[sj].si; oj<=site_info[sj].ei; oj++) {
                // add orbit sub-matrix times phase
                d[oj*_states + oi] += T(iR, oj, oi)*exp_I_phi;
            }
        }
    }
}

void Model::dH_dk(complex<double>* dx, complex<double>* dy, complex<double>* dz,
        double k1, double k2, double k3) const {

    switch (_kdim) {
        case 0: k1=k2=k3=0; break;
        case 1:    k2=k3=0; break;
        case 2:       k3=0; break;
        case 3:             break;
        default: error("Invalid value for kdim");
    }

    // k-point in Cartesian coordinates
    double kx = k1*b(0,0) + k2*b(1,0) + k3*b(2,0);
    double ky = k1*b(0,1) + k2*b(1,1) + k3*b(2,1);
    double kz = k1*b(0,2) + k2*b(1,2) + k3*b(2,2);

    // initialize data with zeros
    size_t s2 = _states*_states;
    for(size_t i=0; i<s2; i++) {
        dx[i] = 0.0;
        dy[i] = 0.0;
        dz[i] = 0.0;
    }

    complex<double> iAx,iAy,iAz;

    // initialize data with on-site
    for(size_t si=0; si<_sites; si++)
    for(size_t sj=0; sj<_sites; sj++) {
        // add hop from si to sj with R=0
        complex<double> exp_I_phi = exp(complex<double>(0,
                      (site_info[sj].pos[0] - site_info[si].pos[0])*kx
                    + (site_info[sj].pos[1] - site_info[si].pos[1])*ky
                    + (site_info[sj].pos[2] - site_info[si].pos[2])*kz));

        if(dx) iAx = complex<double>(0, (site_info[sj].pos[0] - site_info[si].pos[0]))*exp_I_phi;
        if(dy) iAy = complex<double>(0, (site_info[sj].pos[1] - site_info[si].pos[1]))*exp_I_phi;
        if(dz) iAz = complex<double>(0, (site_info[sj].pos[2] - site_info[si].pos[2]))*exp_I_phi;

        for(size_t oi=site_info[si].si; oi<=site_info[si].ei; oi++)
        for(size_t oj=site_info[sj].si; oj<=site_info[sj].ei; oj++) {
            // add orbit sub-matrix times phase
            if(dx) dx[oj*_states + oi] = iAx*V(oj, oi);
            if(dy) dy[oj*_states + oi] = iAy*V(oj, oi);
            if(dz) dz[oj*_states + oi] = iAz*V(oj, oi);
        }
    }

    // add hops
    double Rx,Ry,Rz;
    for(size_t iR=0; iR<_hops; iR++) {
        double R_dot_k = 2*M_PI*(R(iR,0)*k1 + R(iR,1)*k2 + R(iR,2)*k3);

        if(dx) Rx = R(iR,0)*a(0,0) + R(iR,1)*a(1,0) + R(iR,2)*a(2,0);
        if(dy) Ry = R(iR,0)*a(0,1) + R(iR,1)*a(1,1) + R(iR,2)*a(2,1);
        if(dz) Rz = R(iR,0)*a(0,2) + R(iR,1)*a(1,2) + R(iR,2)*a(2,2);

        for(size_t si=0; si<_sites; si++)
        for(size_t sj=0; sj<_sites; sj++) {
            // add hop from si to sj
            complex<double> exp_I_phi = exp(complex<double>(0, R_dot_k
                        + (site_info[sj].pos[0] - site_info[si].pos[0])*kx
                        + (site_info[sj].pos[1] - site_info[si].pos[1])*ky
                        + (site_info[sj].pos[2] - site_info[si].pos[2])*kz));

            if(dx) iAx = complex<double>(0, (Rx + site_info[sj].pos[0] - site_info[si].pos[0]))*exp_I_phi;
            if(dy) iAy = complex<double>(0, (Ry + site_info[sj].pos[1] - site_info[si].pos[1]))*exp_I_phi;
            if(dz) iAz = complex<double>(0, (Rz + site_info[sj].pos[2] - site_info[si].pos[2]))*exp_I_phi;

            for(size_t oi=site_info[si].si; oi<=site_info[si].ei; oi++)
            for(size_t oj=site_info[sj].si; oj<=site_info[sj].ei; oj++) {
                // add orbit sub-matrix times phase
                if(dx) dx[oj*_states + oi] += iAx*T(iR, oj, oi);
                if(dy) dy[oj*_states + oi] += iAy*T(iR, oj, oi);
                if(dz) dz[oj*_states + oi] += iAz*T(iR, oj, oi);
            }
        }
    }
}


void Model::Gk(std::complex<double>* d, double k1, double k2, double k3,
        double Ef, double zeroj, std::complex<double>* S) const {
    math::DenseMap<cxdouble> mGk(d, _states, _states);

    Hk(d, k1,k2,k3);
    // Now: Gk = Hk
    mGk *= -1;
    // Now: Gk = -Hk
    complex<double> Ef0j(Ef, fabs(zeroj));
    for(unsigned i=0; i<_states; i++) {
        mGk(i,i) += Ef0j;
    }
    // Now: Gk = (Ef+0j)*I - Hk
    if (S != NULL) {
        math::DenseMap<cxdouble> mS(S, _states, _states);
        mGk -= mS;
    }
    // Now: Gk = (Ef+0j)*I - Hk - S
    mGk = mGk.inverse();
    // Now: Gk = inverse((Ef+0j)*I - Hk - S)
}


#ifdef TBPP_WITH_HDF5
void Model::save(EHFile& file, const string& prefix) const {
    Node::save(file, prefix);

    file.set_data(prefix+"/a", a);
    if(_use_lattice) {
        file.set_attr(prefix, "use_lattice", 1);
    } else {
        file.set_attr(prefix, "use_lattice", 0);
    }

    if(_lattice != nullptr) {
        file.set_attr(prefix, "lattice", _lattice->name());
    } else {
        file.set_attr(prefix, "lattice", "");
    }

    file.set_attr(prefix, "kdim", _kdim);
    file.set_attr(prefix, "sites", _sites);
    file.set_attr(prefix, "states", _states);
    file.set_attr(prefix, "hops", _hops);

    if(_status == NodeStatus::done) {
        file.set_data(prefix+"/b", b);
        file.set_data(prefix+"/V", V);
        file.set_data(prefix+"/T", T);
        file.set_data(prefix+"/R", R);

        // site_info
        {
            H5::CompType type(sizeof(SiteInfo));
            type.insertMember("kind", HOFFSET(SiteInfo, kind),
                    H5::StrType(H5::PredType::C_S1, TBPP_MAX_SITE_KIND_LEN));
            type.insertMember("x", HOFFSET(SiteInfo, pos[0]), H5::PredType::NATIVE_DOUBLE);
            type.insertMember("y", HOFFSET(SiteInfo, pos[1]), H5::PredType::NATIVE_DOUBLE);
            type.insertMember("z", HOFFSET(SiteInfo, pos[2]), H5::PredType::NATIVE_DOUBLE);
            type.insertMember("states", HOFFSET(SiteInfo, states), H5::PredType::NATIVE_UINT64);
            type.insertMember("si", HOFFSET(SiteInfo, si), H5::PredType::NATIVE_UINT64);
            type.insertMember("ei", HOFFSET(SiteInfo, ei), H5::PredType::NATIVE_UINT64);

            hsize_t dim = site_info.size();
            file.set_data(prefix+"/site_info", (void*)&site_info[0], 1, type, &dim);
        }

        // state_info
        {
            H5::CompType ctype(sizeof(std::complex<double>));
            ctype.insertMember("r", 0, H5::PredType::NATIVE_DOUBLE);
            ctype.insertMember("i", sizeof(double), H5::PredType::NATIVE_DOUBLE);

            H5::CompType type(sizeof(StateInfo));
            type.insertMember("site", HOFFSET(StateInfo,site), H5::PredType::NATIVE_UINT64);
            type.insertMember("orbit", HOFFSET(StateInfo,orbit), H5::PredType::NATIVE_UINT64);
            hsize_t spinor_dim = 2;
            type.insertMember("spinor", HOFFSET(StateInfo,spinor), H5::ArrayType(ctype,1,&spinor_dim));

            hsize_t dim = state_info.size();
            file.set_data(prefix+"/state_info", (void*)&state_info[0], 1, type, &dim);
        }
    }
}

void Model::load(EHFile& file, const string& prefix) {
    Node::load(file, prefix);

    a.resize(3,3);
    file.get_data(prefix+"/a", a);

    int i = 0;
    file.get_attr(prefix, "use_lattice", i);
    if(i == 1) {
        _use_lattice = true;
    } else {
        _use_lattice = false;
    }


    string lattice_name;
    file.get_attr(prefix, "lattice", lattice_name);
    if(lattice_name.length() > 0 and context() != nullptr) {
        set_lattice(dynamic_pointer_cast<Lattice>(context()->get(lattice_name)));
    } else {
        set_lattice(nullptr);
    }

    file.get_attr(prefix, "kdim", _kdim);
    file.get_attr(prefix, "sites", _sites);
    file.get_attr(prefix, "states", _states);
    file.get_attr(prefix, "hops", _hops);

    if(_status == NodeStatus::done) {
        file.get_data(prefix+"/b", b);
        file.get_data(prefix+"/V", V);
        file.get_data(prefix+"/T", T);
        file.get_data(prefix+"/R", R);

        // site_info
        {
            H5::CompType type(sizeof(SiteInfo));
            type.insertMember("kind", HOFFSET(SiteInfo, kind),
                    H5::StrType(H5::PredType::C_S1, TBPP_MAX_SITE_KIND_LEN));
            type.insertMember("x", HOFFSET(SiteInfo, pos[0]), H5::PredType::NATIVE_DOUBLE);
            type.insertMember("y", HOFFSET(SiteInfo, pos[1]), H5::PredType::NATIVE_DOUBLE);
            type.insertMember("z", HOFFSET(SiteInfo, pos[2]), H5::PredType::NATIVE_DOUBLE);
            type.insertMember("states", HOFFSET(SiteInfo, states), H5::PredType::NATIVE_UINT64);
            type.insertMember("si", HOFFSET(SiteInfo, si), H5::PredType::NATIVE_UINT64);
            type.insertMember("ei", HOFFSET(SiteInfo, ei), H5::PredType::NATIVE_UINT64);

            hsize_t size;
            H5::DataSet dataset(file.get_file_ptr()->openDataSet(prefix+"/site_info"));
            file.check_dataset(dataset, type, 1);
            file.get_size_dataset(dataset, &size);
            site_info.resize(static_cast<size_t>(size));
            dataset.read(site_info.data(), type);
        }

        // state_info
        {
            H5::CompType ctype(sizeof(std::complex<double>));
            ctype.insertMember("r", 0, H5::PredType::NATIVE_DOUBLE);
            ctype.insertMember("i", sizeof(double), H5::PredType::NATIVE_DOUBLE);

            H5::CompType type(sizeof(StateInfo));
            type.insertMember("site", HOFFSET(StateInfo,site), H5::PredType::NATIVE_UINT64);
            type.insertMember("orbit", HOFFSET(StateInfo,orbit), H5::PredType::NATIVE_UINT64);
            hsize_t spinor_dim = 2;
            type.insertMember("spinor", HOFFSET(StateInfo,spinor), H5::ArrayType(ctype,1,&spinor_dim));

            hsize_t size;
            H5::DataSet dataset(file.get_file_ptr()->openDataSet(prefix+"/state_info"));
            file.check_dataset(dataset, type, 1);
            file.get_size_dataset(dataset, &size);
            state_info.resize(static_cast<size_t>(size));
            dataset.read(state_info.data(), type);
        }
    }
}
#endif // TBPP_WITH_HDF5

//----------------------------------------------------------------------

} // namespace tbpp
