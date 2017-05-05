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

#include <tbpp/lattice.h>

using namespace std;

namespace tbpp {

//----------------------------------------------------------------------

const char* Lattice::type() const { return "Lattice"; }

void Lattice::add_site(const std::string& kind, double x, double y, double z,
            uint32_t states, const std::complex<double> *onsite,
            const uint32_t *orbits, const std::complex<double> *spinors) {
    Site site;
    site.kind = kind;
    site.pos[0] = x;
    site.pos[1] = y;
    site.pos[2] = z;
    site.orbits.resize(states);
    site.spinors.resize(states);
    site.onsite.resize(states, states);

    for(uint32_t i=0; i<states; i++) {
        site.orbits[i] = orbits[i];
    }

    cxdouble* ptr = &site.spinors[0][0];
    for(uint32_t i=0; i<2*states; i++) {
        ptr[i] = spinors[i];
    }

    ptr = site.onsite.data();
    for(uint32_t i=0; i<states*states; i++) {
        ptr[i] = onsite[i];
    }

    sites.push_back(site);
}

void Lattice::add_site(const std::string& kind, double x, double y, double z,
        const math::DenseMatrix<cxdouble>& onsite,
        const std::vector<uint32_t>& orbits,
        const std::vector<Spinor>& spinors) {

    uint32_t states = onsite.rows();
    if(states != onsite.cols()) throw runtime_error("Invalid onsite size");
    if(states != orbits.size()) throw runtime_error("Invalid orbit size");
    if(states != spinors.size()) throw runtime_error("Invalid Spinors size");

    add_site(kind, x, y, z, states, onsite.data(), orbits.data(), &(spinors[0][0]));
}


void Lattice::add_hop(const std::string& initial_kind, const std::string& final_kind,
        double dx, double dy, double dz, uint32_t t_rows, uint32_t t_cols,
        const std::complex<double>* T) {

    double dpos = sqrt(dx*dx + dy*dy + dz*dz);

    if (initial_kind == final_kind and dpos < eps) {
        throw runtime_error("[Error] Hop between onsite states must be provided in the "
                            "site onsite matrix not as a hop term");
    }

    if (initial_kind > final_kind) {
        // conjugate transpose
        math::DenseMatrix<cxdouble> T_conj(t_cols, t_rows);
        for(uint32_t i=0; i<t_rows; i++)
        for(uint32_t j=0; j<t_cols; j++) {
            T_conj(j,i) = conj(T[t_cols*i + j]);
        }

        // add final_kind to initial kind hopping instead
        add_hop(final_kind, initial_kind, -dx, -dy, -dz, T_conj);
        return;
    }

    // check that hop does not already exist
    auto it = hops.find(initial_kind);
    if( it != hops.end() ) {
        // hops[initial_kind] exists
        auto it2 = it->second.find(final_kind);
        if( it2 != it->second.end() ) {
            // hops[initial_kind][final_kind] exists
            for(const auto& h : it2->second) {
                double ddpos = sqrt(pow(dx-h.dpos[0],2) + pow(dy-h.dpos[1], 2) + pow(dz-h.dpos[2], 2));
                if(ddpos < eps) {
                    stringstream stream;
                    stream << "[Error] hopping from " << initial_kind << " to " << final_kind
                           << " for dpos = [ " << dx << ", " << dy << ", " << dz << " ] "
                           << "already exists.";
                    throw runtime_error(stream.str());
                }
            }
        }
    }

    Hop h;
    h.initial_kind = initial_kind;
    h.final_kind = final_kind;
    h.dpos[0] = dx;
    h.dpos[1] = dy;
    h.dpos[2] = dz;
    h.T.resize(t_rows, t_cols);
    cxdouble *ptr = h.T.data();
    for(size_t i=0; i<h.T.size(); i++)
        ptr[i] = T[i];

    // add hop
    max_hop_range = max(dpos, max_hop_range);
    hops[initial_kind][final_kind].push_back(h);

    if(initial_kind == final_kind) {
        // add final_kind to initial kind hopping as well
        h.dpos[0] = -dx;
        h.dpos[1] = -dy;
        h.dpos[2] = -dz;
        h.T.resize(t_cols, t_rows);

        for(uint32_t i=0; i<t_rows; i++)
        for(uint32_t j=0; j<t_cols; j++) {
            h.T(j,i) = conj(T[t_cols*i + j]);
        }

        hops[initial_kind][final_kind].push_back(h);
    }
}

void Lattice::add_hop(const std::string& initial_kind, const std::string& final_kind,
        double dx, double dy, double dz, const math::DenseMatrix<cxdouble>& T) {
    add_hop(initial_kind, final_kind, dx, dy, dz, T.rows(), T.cols(), T.data());
}

void Lattice::add_trans(double x, double y, double z) {
    if(trans.size() >= 3)
        throw runtime_error("Too many translational symmetries");
    trans.push_back({x,y,z});
}

void Lattice::set_B(double Bx, double By, double Bz) {
    B[0] = Bx;
    B[1] = By;
    B[2] = Bz;
}

math::DenseMatrix<cxdouble> Lattice::_T(uint32_t site_i, uint32_t site_j,
        double dx, double dy, double dz) const {
    math::DenseMatrix<cxdouble> ret;

    string kind_i = sites[site_i].kind;
    string kind_j = sites[site_j].kind;
    if (kind_i > kind_j) {
        math::DenseMatrix<cxdouble> t = _T(site_j, site_i, -dx, -dy, -dz);
        ret.resize(t.cols(), t.rows());
        // complex conjugate transpose
        for(uint32_t i=0; i<t.rows(); i++)
        for(uint32_t j=0; j<t.cols(); j++) {
            ret(j,i) = conj(t(i,j));
        }
        return ret;
    }

    auto it = hops.find(kind_i);
    if( it != hops.end() ) {
        // hops[kind_i] exists
        auto it2 = it->second.find(kind_j);
        if( it2 != it->second.end() ) {
            // hops[kind_i][kind_j] exists
            for(const auto& h : it2->second) {
                double ddpos = sqrt(pow(dx-h.dpos[0],2) + pow(dy-h.dpos[1], 2)
                        + pow(dz-h.dpos[2], 2));
                if(ddpos < eps) {
                    ret = h.T;
                    return ret;
                }
            }
        }
    }

    // return empty matrix
    return ret;
}

void Lattice::copy_along(double vx, double vy, double vz, unsigned n){
    size_t sites_n = sites.size();
    for(size_t in=1; in<n; in++) {
        for(size_t i=0; i<sites_n; i++) {
            Site s = sites[i];
            s.pos[0] += in*vx;
            s.pos[1] += in*vy;
            s.pos[2] += in*vz;
            sites.push_back(s);
        }
    }
}

void Lattice::move_along(double vx, double vy, double vz){
    for(auto& s : sites) {
        s.pos[0] += vx;
        s.pos[1] += vy;
        s.pos[2] += vz;
    }
}

//----------------------------------------------------------------------

} // namespace tbpp

