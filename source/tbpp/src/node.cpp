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

#include <tbpp/node.h>
#include <iostream>
#include <sstream>

#include <tbpp/lattice.h>
#include <tbpp/model.h>
#include <tbpp/kmodel.h>
#include <tbpp/kpath.h>
#include <tbpp/kgrid.h>
#include <tbpp/cpa.h>
#include <tbpp/cond_dos.h>
#include <tbpp/gen_cond.h>
#include <tbpp/approx_cond_dos.h>
#include <tbpp/berryflux.h>

using namespace std;

namespace tbpp {

//----------------------------------------------------------------------------

shared_ptr<Node> Node::new_node(const string& type) {
    if(type == "Node") return make_shared<Node>();
    if(type == "Lattice") return make_shared<Lattice>();
    if(type == "Model") return make_shared<Model>();
    if(type == "KModel") return make_shared<KModel>();
    if(type == "KPath") return make_shared<KPath>();
    if(type == "KGrid") return make_shared<KGrid>();
    if(type == "CPA") return make_shared<CPA>();
    if(type == "CondDOS") return make_shared<CondDOS>();
    if(type == "GenCond") return make_shared<GenCond>();
    if(type == "ApproxCondDOS") return make_shared<ApproxCondDOS>();
    if(type == "BerryFlux") return make_shared<BerryFlux>();
    throw runtime_error("[error Node::new_node(...)] Unknown node type: " + type);
}

const string& Node::name() const { return _name; }
Context* Node::context() const { return _ct; }
const char* Node::type() const { return "Node"; }

NodeStatus Node::status() const { return _status; }

bool Node::parallel() const { return _parallel; }
void Node::set_parallel(bool parallel) { _parallel = parallel; }

void Node::solve() { header(); _status = NodeStatus::done; }
void Node::clear() { }
void Node::reset() { _status = NodeStatus::not_done; clear(); }

void Node::make_ready() {
    if(_status == NodeStatus::done) return;
    solve();
    if(_status != NodeStatus::done)
        throw runtime_error("Node " + _name + " is not ready");
}

#ifdef TBPP_WITH_HDF5
void Node::save(EHFile& file, const string& prefix) const {
    file.set_attr(prefix, "status", node_status_to_string(_status));
    file.set_attr(prefix, "type", type());
    file.set_attr(prefix, "parallel", _parallel);
}

void Node::load(EHFile& file, const string& prefix) {
    string str;
    file.get_attr(prefix, "status", str);
    _status = node_status_from_string(str);

    string node_type(type());
    file.get_attr(prefix, "type", str);
    if(str != node_type)
        error("Incorrect object type: " + node_type);
    file.get_attr(prefix, "parallel", _parallel);
}
#endif // TBPP_WITH_HDF5

string Node::node_status_to_string(NodeStatus status) {
    switch (status) {
        case NodeStatus::not_done: return "not done";
        case NodeStatus::done: return "done";
        case NodeStatus::failed: return "failed";
    }
    return "unknown";
}

NodeStatus Node::node_status_from_string(const string& status) {
    if(status == "not done") return NodeStatus::not_done;
    if(status == "done") return NodeStatus::done;
    if(status == "failed") return NodeStatus::failed;
    throw runtime_error("unknown node status:" + status);
}

void Node::header() const {
    for(int i=0; i<78; i++) cout << "-";
    cout << '\n';
    cout << type() << " " << name() << '\n';
    for(int i=0; i<78; i++) cout << "-";
    cout << '\n';
}

void Node::debug(const string& msg) const {
    cout << "[debug][" << type() << " " << name() << "] " << msg << '\n';
}

void Node::msg(const string& msg) const {
    cout << "[msg][" << type() << " " << name() << "] " << msg << '\n';
}

void Node::warn(const string& msg) const {
    cout << "[warn]["<< type() << " " << name() << "] " << msg << '\n';
}

void Node::error(const string& msg) const {
    stringstream str;
    str << "[" << type() << " " << name() << "] " << msg << '\n';
    throw runtime_error(str.str());
}

//----------------------------------------------------------------------------

} // namespace tbpp
