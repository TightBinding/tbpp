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

#include <tbpp/cpa.h>

// Disable Eigen Multi-threading
#define EIGEN_DONT_PARALLELIZE

using namespace std;

namespace tbpp {

//----------------------------------------------------------------------

typedef struct DefectSiteInfo {
    size_t site_index;      ///< index of site
    size_t fsi;             ///< start index for flat layout
    size_t states;          ///< total number of states
    size_t states2;         ///< total number of states
    size_t si;              ///< start index for non-flat layout
    size_t ei;              ///< end index for non-flat layout
    std::vector<double> c_onsite;  ///< probability of original on-site potential
    /// list of pointers to defects
    std::list<const CPA::Defect*> defects_list;
} DefectSiteInfo;

//----------------------------------------------------------------------

const char* CPA::type() const { return "CPA"; }

void CPA::add_defect(const std::vector<size_t>& sites, const std::vector<double>& c,
        const math::DenseMatrix<cxdouble>& V) {
    Defect d;
    d.sites = sites;
    d.c = c;
    d.V = V;

    if(!defects.empty()) {
        if(d.c.size() != defects.front().c.size())
            throw runtime_error("Invalid size of concentration array.");
    }

    defects.push_back(d);
}

void CPA::set_kmodel(KModelPtr kmodel) { _kmodel = kmodel; }
KModelPtr CPA::kmodel() const { return _kmodel; }

void CPA::set_w(const std::vector<double>& _w) { w = _w; }

void CPA::solve() {
    _status = NodeStatus::failed;

    //-----------------------------------------------------------------------
    // General Checks

    if (_kmodel == nullptr ) throw runtime_error("No KModel");
    _kmodel->make_ready();

    ModelPtr _model = _kmodel->model();
    if (_model == nullptr ) throw runtime_error("No Model");

    header();
    //-----------------------------------------------------------------------

    const size_t states = _model->states();
    const size_t sites = _model->sites();
    const size_t states2 = states*states;
    const bool cache = _kmodel->cache();

    const vector<double>& k1 = _kmodel->k1;
    const vector<double>& k2 = _kmodel->k2;
    const vector<double>& k3 = _kmodel->k3;

    const size_t kpoints = k1.size()*k2.size()*k3.size();
    const size_t w_size = w.size();

    size_t c_size = 0;
    if(!defects.empty()) {
        c_size = defects.front().c.size();
        for(const auto& d : defects) {
            if(d.c.size() != c_size) {
                throw runtime_error("Invalid size for list of concentrations");
            }
        }
    } else {
        c_size = 1;
    }

    //-----------------------------------------------------------------------
    // Get Defect Site Info

    /// map from site index to defect site information
    map<size_t, DefectSiteInfo> defect_sites;

    size_t G_size = 0;
    size_t DD_size = 0;

    // Add all defect info to defect_sites
    for(const auto& d : defects) {
        for(const auto& ids : d.sites) {
            if(ids >= _model->site_info.size())
                throw runtime_error("Invalid site index");

            auto it = defect_sites.find(ids);
            DefectSiteInfo& ds = defect_sites[ids];

            if(it == defect_sites.end()) {
                //we have just adding a new site
                const Model::SiteInfo& s = _model->site_info[ids];
                ds.site_index = ids;
                ds.fsi = G_size;
                ds.states = s.states;
                ds.states2 = s.states*s.states;
                ds.si = s.si;
                ds.ei = s.ei;

                G_size += ds.states2;
                DD_size = max(DD_size, ds.states2);
            }

            // add defect info to site
            ds.defects_list.push_back(&d);
            if(d.V.rows() != ds.states or d.V.cols() != ds.states) {
                cout << "row = " << d.V.rows() << " cols = " << d.V.cols() << "\n";
                cout << "states = " << ds.states << "\n";
                throw runtime_error("Invalid size for defect on-site potential");
            }
        }
    }

    // Determine c_onsite for each defect site
    for(auto& it : defect_sites) {
        auto& ds = it.second;
        ds.c_onsite.resize(c_size);
        for(size_t ic=0; ic<c_size; ic++) {
            ds.c_onsite[ic] = 1;
            for(const auto& d : ds.defects_list) {
                ds.c_onsite[ic] -= d->c[ic];
            }
            if(ds.c_onsite[ic] < 0.0)
                throw runtime_error("Invalid defect concentrations");
        }
    }


    //-----------------------------------------------------------------------
    // Virtual Crystal Approximation with initial broadening

    vca.resize(c_size, states, states);

    for(size_t ic=0; ic<c_size; ic++) {

        // Copy Original Potential Matrix
        for(size_t r=0; r<states; r++)
        for(size_t c=0; c<states; c++) {
            vca(ic,r,c) = _model->V(r,c);
        }

        // Adjust defect sites
        for(const auto& it : defect_sites) {
            const auto& ds = it.second;
            // Multiply original potential by c_onsite
            for(size_t r=ds.si; r<=ds.ei; r++)
            for(size_t c=ds.si; c<=ds.ei; c++) {
                vca(ic, r, c) *= ds.c_onsite[ic];
            }
            // Add defects
            for(const auto d : ds.defects_list) {
                for(size_t r=0; r<ds.states; r++)
                for(size_t c=0; c<ds.states; c++) {
                    vca(ic,ds.si+r,ds.si+c) += d->c[ic]*d->V(r,c);
                }
            }
        }

        // Add initial broadening
        for(size_t i=0; i<states; i++)
            vca(ic,i,i) -= cxdouble(0, zero_j);

    }

    //-----------------------------------------------------------------------
    // VCA as self-energy

    sigma.resize(c_size, w_size, states, states);
    if (use_vca) {
        for(size_t ic=0; ic<c_size; ic++)
        for(size_t iw=0; iw<w_size; iw++) {
            for(size_t i=0; i<states; i++)
            for(size_t j=0; j<states; j++) {
                sigma(ic,iw,i,j) = vca(ic,i,j);
            }
        }
        _status = NodeStatus::done;
        return;
    }

    //-----------------------------------------------------------------------
    // Compute Self-Energy using CPA

    #pragma omp parallel if(_parallel)
    {
        math::DenseMatrix<cxdouble> mGk(states, states); //< G(k) = Inv(w*I - T - S)
        math::DenseMatrix<cxdouble> mS1(states, states); //< Self-Energy 1 Buffer
        math::DenseMatrix<cxdouble> mS2(states, states); //< Self-Energy 2 Buffer

        vector<cxdouble> vG(G_size);   //< Local Green matrix (only for defect sites)
        vector<cxdouble> DD(DD_size);  //< Temporary buffer (DD = inv(G) + S)
        vector<cxdouble> DDi(DD_size); //< Temporary buffer (DDi = inv(DD - Vi)
        vector<cxdouble> SD(DD_size);  //< Temporary buffer (SD += c*DDi)

        cxdouble *__restrict__ S = mS1.data();     //< pointer to new self-energy
        cxdouble *__restrict__ Sprev = mS2.data(); //< pointer to previous self-energy
        cxdouble *__restrict__ Gk = mGk.data();    //< pointer to Gk
        cxdouble *__restrict__ G = vG.data();      //< pointer to Local Green
        cxdouble *__restrict__ p1; //< General purpose restricted pointer


        #pragma omp for collapse(2) schedule(dynamic)
        for(size_t ic=0; ic<c_size; ic++)
        for(size_t iw=0; iw<w_size; iw++) {

            // Copy VCA to S and Sprev
            cxdouble *__restrict__ VCA = &vca(ic,0,0);
            for (size_t i=0; i<states2; i++) {
                S[i] = VCA[i];
                Sprev[i] = VCA[i];
            }

            // Self-Consistent Loop
            double dS = 2*eps;
            size_t count = 0;
            while(dS > eps) {
                // Zero Local Greens
                for(size_t i=0; i<G_size; i++) {
                    G[i] = 0;
                }

                // Compute Local Green's Matrix
                for(size_t ik1=0; ik1<k1.size(); ik1++)
                for(size_t ik2=0; ik2<k2.size(); ik2++)
                for(size_t ik3=0; ik3<k3.size(); ik3++) {

                    if(cache) {
                        p1 = &_kmodel->Tk(ik1,ik2,ik3,0,0);
                        for(size_t i=0; i<states2; i++) {
                            Gk[i] = -(p1[i] + S[i]);
                        }
                    } else {
                        _model->Tk(Gk, k1[ik1],k2[ik2],k3[ik3]);
                        // Now: Gk = T
                        for(size_t i=0; i<states2; i++) {
                            Gk[i] = -(Gk[i] + S[i]);
                        }
                    }

                    // Now: Gk = -T-S

                    for(size_t i=0; i<states; i++) {
                        Gk[i*states + i] += w[iw];
                    }

                    // Now: Gk = w*I - T - S

                    mGk = mGk.inverse();

                    // Now: Gk = Inv(w*I - T - S)

                    // Sum only necessary entries of the Local Green's Matrix
                    for(const auto& it: defect_sites) {
                        const auto& ds = it.second;
                        p1 = &G[ds.fsi];
                        for(size_t r=ds.si; r<=ds.ei; r++)
                        for(size_t c=ds.si; c<=ds.ei; c++) {
                            (*p1++) += mGk(r,c);
                        }
                    }
                }

                // Normalize Local Greens with respect to kpoints
                for(size_t i=0; i<G_size; i++) {
                    G[i] /= kpoints;
                }

                // Swap memory pointers (write to S and read from Sprev)
                std::swap(Sprev, S);

                for(const auto& it: defect_sites) {
                    const auto& ds = it.second;
                    math::DenseMap<cxdouble> mG(&G[ds.fsi], ds.states, ds.states);
                    math::DenseMap<cxdouble> mDDi(DDi.data(), ds.states, ds.states);
                    math::DenseMap<cxdouble> mSD(SD.data(), ds.states, ds.states);

                    // Invert local green matrix
                    mG = mG.inverse();

                    // Now: G = inv(G)

                    size_t i=0;
                    for(size_t r=ds.si; r<=ds.ei; r++)
                    for(size_t c=ds.si; c<=ds.ei; c++) {
                        DD[i] = G[ds.fsi + i] + Sprev[r*states + c];
                        SD[i] = 0;
                        i++;
                    }

                    // Now: DD = inv(G) + S
                    //      SD = 0

                    // add original on-site-potential
                    if(ds.c_onsite[ic] > 0) {
                        size_t i=0;
                        for(size_t r=ds.si; r<=ds.ei; r++)
                        for(size_t c=ds.si; c<=ds.ei; c++) {
                            DDi[i] = DD[i] - _model->V(r,c);
                            i++;
                        }

                        // Now: DDi = DD - Vi

                        mDDi = mDDi.inverse();

                        // Now: DDi = inv(DD - Vi)

                        for(size_t i=0; i<ds.states2; i++) {
                            SD[i] += ds.c_onsite[ic]*DDi[i];
                        }
                    }

                    // add defect substitutes
                    for(const auto& d : ds.defects_list) {
                        p1 = const_cast<cxdouble*>(d->V.data());

                        for(size_t i=0; i<ds.states2; i++) {
                            DDi[i] = DD[i] - p1[i];
                        }

                        // Now: DDi = DD - Vi

                        mDDi = mDDi.inverse();

                        // Now: DDi = inv(DD - Vi)

                        for(size_t i=0; i<ds.states2; i++) {
                            SD[i] += d->c[ic]*DDi[i];
                        }
                    }

                    mSD = mSD.inverse();

                    // Now: SD = inv( Sum_i[ c*inv(DD - Vi)] )

                    // Set new Self-Energy for site
                    i=0;
                    for(size_t r=ds.si; r<=ds.ei; r++)
                    for(size_t c=ds.si; c<=ds.ei; c++) {
                        S[r*states + c] = mix_factor*(DD[i] - SD[i])
                            + (1-mix_factor)*Sprev[r*states +c];
                        i++;
                    }
                }

                count++;
                if (count > max_iter) {
                    cerr << "[error] Did not converge for w=" << w[iw] << endl;
                    break;
                }

                dS = 0;
                for(size_t i=0; i<states2; i++) {
                    dS = max(dS, abs(S[i] - Sprev[i]));
                }
            }

            // Store S
            p1 = &sigma(ic,iw,0,0);
            for (size_t i=0; i<states2; i++) {
                p1[i] = S[i];
            }
            #pragma omp critical
            {
                cout << "finished ic = " << ic << " iw = " << iw << endl;
            }
        }
    } // pragma omp parallel


    _status = NodeStatus::done;
}

#ifdef TBPP_WITH_HDF5
void CPA::save(EHFile& file, const std::string& prefix) const {
    Node::save(file, prefix);
    if(_kmodel != nullptr)
        file.set_attr(prefix, "kmodel", _kmodel->name());
    else
        file.set_attr(prefix, "kmodel", "");

    file.set_data(prefix+"/w", w);
    file.set_data(prefix+"/sigma", sigma);
    file.set_data(prefix+"/vca", vca);

    file.set_attr(prefix, "eps", eps);
    file.set_attr(prefix, "zero_j", zero_j);
    file.set_attr(prefix, "max_iter", max_iter);
    file.set_attr(prefix, "mix_factor", mix_factor);
    file.set_attr(prefix, "use_vca", use_vca);

    // save defects data
    size_t max_sites_num = 0;
    size_t max_states_num = 0;
    size_t c_num = 0;
    if(!defects.empty()) {
        c_num = defects.front().c.size();
    }
    for(const auto& d : defects) {
        max_sites_num = max(max_sites_num, d.sites.size());
        max_states_num = max(max_states_num, static_cast<size_t>(d.V.rows()));
        if(c_num != d.c.size())
            error("Invalid size for concentrations");
    }

    // number of sites for each defect
    vector<size_t> sites_num(defects.size());
    // list of sites
    math::NArray<size_t,2> sites(defects.size(), max_sites_num);

    // number of states for each defect
    vector<size_t> states_num(defects.size());
    // V for each defect
    math::NArray<cxdouble,3> V(defects.size(), max_states_num, max_states_num);

    // concentrations for each defect
    math::NArray<double,2> c(defects.size(), c_num);

    size_t i=0;
    for(const auto& d : defects) {

        sites_num[i] = d.sites.size();
        for(size_t j=0; j<d.sites.size(); j++)
            sites(i,j) = d.sites[j];

        states_num[i] = d.V.rows();
        for(size_t r=0; r<states_num[i]; r++)
        for(size_t c=0; c<states_num[i]; c++) {
            V(i,r,c) = d.V(r,c);
        }

        for(size_t j=0; j<d.c.size(); j++)
            c(i,j) = d.c[j];

        i++;
    }

    file.set_data(prefix+"/defects_sites_num", sites_num);
    file.set_data(prefix+"/defects_sites", sites);

    file.set_data(prefix+"/defects_states_num", states_num);
    file.set_data(prefix+"/defects_V", V);

    file.set_data(prefix+"/defects_c", c);
}

void CPA::load(EHFile& file, const std::string& prefix) {
    Node::load(file, prefix);

    string kmodel_name;
    file.get_attr(prefix, "kmodel", kmodel_name);
    if(kmodel_name.length() > 0 and context() != nullptr)
        _kmodel = dynamic_pointer_cast<KModel>(context()->get(kmodel_name));
    else
        _kmodel = nullptr;

    file.get_data(prefix+"/w", w);
    file.get_data(prefix+"/sigma", sigma);
    file.get_data(prefix+"/vca", vca);

    file.get_attr(prefix, "eps", eps);
    file.get_attr(prefix, "zero_j", zero_j);
    file.get_attr(prefix, "max_iter", max_iter);
    file.get_attr(prefix, "mix_factor", mix_factor);
    file.get_attr(prefix, "use_vca", use_vca);

    // load defect info
    vector<size_t> sites_num;
    math::NArray<size_t,2> sites;

    vector<size_t> states_num;
    math::NArray<cxdouble,3> V;

    math::NArray<double,2> c;

    file.get_data(prefix+"/defects_sites_num", sites_num);
    file.get_data(prefix+"/defects_sites", sites);

    file.get_data(prefix+"/defects_states_num", states_num);
    file.get_data(prefix+"/defects_V", V);

    file.get_data(prefix+"/defects_c", c);

    defects.clear();
    for(size_t i=0; i<sites_num.size(); i++) {
        Defect d;
        d.sites.resize(sites_num[i]);
        for(size_t j=0; j<sites_num[i]; j++)
            d.sites[j] = sites(i,j);

        d.V.resize(states_num[i], states_num[i]);
        for(size_t r=0; r<states_num[i]; r++)
        for(size_t c=0; c<states_num[i]; c++)
            d.V(r,c) = V(i,r,c);

        d.c.resize(c.size(1));
        for(size_t j=0; j<c.size(1); j++)
            d.c[j] = c(i,j);

        defects.push_back(d);
    }
}
#endif // TBPP_WITH_HDF5

//----------------------------------------------------------------------

} //namespace tbpp
