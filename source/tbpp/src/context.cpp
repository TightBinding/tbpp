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
#include <iomanip>
#include <sstream>

#include <tbpp/context.h>
#include <tbpp/node.h>

using namespace std;

namespace tbpp {

//----------------------------------------------------------------------------

std::shared_ptr<Node> Context::create(const string& type, const string& name_hint) {
    std::shared_ptr<Node> node = Node::new_node(type);
    return add(node, name_hint);
}

std::shared_ptr<Node> Context::add(std::shared_ptr<Node> node, const string& name_hint) {
    if (!node) return nullptr;
    if (node->_ct == this) {
        cout << "Node already in this context\n";
    } else if (node->_ct != nullptr) {
        cerr << "Node already belongs to another context\n";
        return nullptr;
    }

    // determine name for the node
    string name_prefix = name_hint;
    if(name_prefix.empty()) name_prefix = node->type();

    string name = name_prefix;
    for(int i=1; nodes.count(name); i++) {
        // name already in map so propose a new name
        stringstream stream;
        stream << name_prefix << "." << setfill('0') << setw(3) << i;
        name = stream.str();
    }

    // insert node
    node->_name = name;
    node->_ct = this;
    nodes.insert(make_pair(name, node));
    return node;
}


std::shared_ptr<Node> Context::get(const string& name) {
    auto it = nodes.find(name);
    if(it == nodes.end())
        return nullptr;
    return it->second;
}

std::shared_ptr<Node> Context::rm(const string& name) {
    auto it = nodes.find(name);
    if(it == nodes.end())
        return nullptr;
    std::shared_ptr<Node> node = it->second;
    node->_ct = nullptr;
    nodes.erase(it);
    return node;
}

void Context::solve() {
    for(const auto& n : nodes) {
        n.second->solve();
    }
}

void Context::make_ready() {
    for(const auto& n : nodes) {
        n.second->make_ready();
    }
}

void Context::clear() { nodes.clear(); }

#ifdef TBPP_WITH_HDF5
void Context::save(const string& filename) const {
    EHFile file(filename, "w");
    for(const auto& n : nodes) {
        n.second->save(file, n.first);
    }
}

void Context::load(const string& filename) {
    EHFile file(filename, "r");
    clear();

    H5::H5File *fp = file.get_file_ptr();
    hsize_t obj_count = fp->getNumObjs();

    // first traverse: instantiate all nodes
    for(size_t oi=0; oi<obj_count; oi++) {
        H5G_obj_t type = fp->getObjTypeByIdx(oi);
        if(type == H5G_GROUP) {
            string name = fp->getObjnameByIdx(oi);
            string type;
            try {
                file.get_attr(name, "type", type);
                create(type, name);
            } catch (exception& e) {
                continue;
            }
        }
    }

    // second traverse: let nodes load data and reconnect pointers
    for(auto n : nodes) {
        try {
            n.second->load(file, n.first);
        } catch (...) {
            // FIXME above should not catch all exceptions
            cerr << "Error Loading Node: " << n.first << endl;
        }
    }
}
#else
void Context::save(const string& filename) const {
    throw runtime_error("Must compile with HDF5 to enable saving");
}

void Context::load(const string& filename) {
    throw runtime_error("Must compile with HDF5 to enable loading");
}
#endif // TBPP_WITH_HDF5

} // namespace tbpp
