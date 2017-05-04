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

#include <tbpp/kmodel.h>

using namespace std;

namespace tbpp {

//----------------------------------------------------------------------

const char* KModel::type() const { return "KModel"; }

KModel::KModel(ModelPtr model, size_t nk1, size_t nk2, size_t nk3,
        bool cache) {
    set_model(model);
    set_grid(nk1, nk2, nk3);
    set_cache(cache);
}

void KModel::set_k1(const std::vector<double>& _k1) { k1 = _k1; }
void KModel::set_k2(const std::vector<double>& _k2) { k2 = _k2; }
void KModel::set_k3(const std::vector<double>& _k3) { k3 = _k3; }

void KModel::set_model(ModelPtr model) { _model = model; }
ModelPtr KModel::model() const { return _model; }

void KModel::set_cache(bool cache) { _cache = cache; }
bool KModel::cache() const { return _cache; }

void KModel::set_grid(size_t nk1, size_t nk2, size_t nk3) {
    k1.resize(nk1);
    k2.resize(nk2);
    k3.resize(nk3);

    for(size_t i=0; i<nk1; i++) k1[i] = static_cast<double>(i)/(nk1-1) - 0.5;
    for(size_t i=0; i<nk2; i++) k2[i] = static_cast<double>(i)/(nk2-1) - 0.5;
    for(size_t i=0; i<nk3; i++) k3[i] = static_cast<double>(i)/(nk3-1) - 0.5;

    if(nk1 == 0 or nk1 == 1) k1 = {0};
    if(nk2 == 0 or nk2 == 1) k2 = {0};
    if(nk3 == 0 or nk3 == 1) k3 = {0};
}

void KModel::solve() {
    _status = NodeStatus::failed;
    //-----------------------------------------------------------------------

    if (_model == nullptr) throw runtime_error("No Model");
    _model->make_ready();

    header();
    //-----------------------------------------------------------------------

    // make k-space consistent with model
    if(_model->kdim() <= 2) k3 = {0};
    if(_model->kdim() <= 1) k2 = {0};
    if(_model->kdim() == 0) k1 = {0};

    if (_cache) {
        size_t nk1 = k1.size();
        size_t nk2 = k2.size();
        size_t nk3 = k3.size();
        size_t states = _model->states();
        size_t states2 = states*states;
        Tk.resize(nk1,nk2,nk3,states,states);

        #pragma omp parallel for collapse(3) if(_parallel)
        for(size_t ik1=0; ik1<nk1; ik1++)
        for(size_t ik2=0; ik2<nk2; ik2++)
        for(size_t ik3=0; ik3<nk3; ik3++) {
            cxdouble *ptr = &Tk(ik1,ik2,ik3,0,0);
            _model->Tk(ptr, k1[ik1], k2[ik2], k3[ik3]);
        }
    } else {
        Tk.clear();
    }

    //-----------------------------------------------------------------------

    _status = NodeStatus::done;
}

#ifdef TBPP_WITH_HDF5
void KModel::save(EHFile& file, const string& prefix) const {
    Node::save(file, prefix);
    if(_model != nullptr)
        file.set_attr(prefix, "model", _model->name());
    else
        file.set_attr(prefix, "model", "");

    file.set_attr(prefix, "cache", _cache);
    file.set_data(prefix+"/k1", k1);
    file.set_data(prefix+"/k2", k2);
    file.set_data(prefix+"/k3", k3);

    if(_status == NodeStatus::done and _cache) {
        file.set_data(prefix+"/Tk", Tk);
    }
}
void KModel::load(EHFile& file, const string& prefix) {
    Node::load(file, prefix);

    string model_name;
    file.get_attr(prefix, "model", model_name);
    if(model_name.length() > 0 and context() != nullptr)
        _model = dynamic_pointer_cast<Model>(context()->get(model_name));
    else
        _model = nullptr;

    file.get_attr(prefix, "cache", _cache);
    file.get_data(prefix+"/k1", k1);
    file.get_data(prefix+"/k2", k2);
    file.get_data(prefix+"/k3", k3);

    if(_status == NodeStatus::done and _cache) {
        file.get_data(prefix+"/Tk", Tk);
    }
}
#endif // TBPP_WITH_HDF5

//----------------------------------------------------------------------

} // namespace tbpp


