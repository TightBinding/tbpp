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

#ifndef __TBPP_CONTEXT_H__
#define __TBPP_CONTEXT_H__

/**
 * \file
 * \brief Contains a collection of nodes representing a simulation pipeline
 */


#include <tbpp/common.h>
#include <map>
#include <memory>
#include <string>

namespace tbpp {

//----------------------------------------------------------------------------

class Node;
class Context;
typedef std::shared_ptr<Context> ContextPtr;

//----------------------------------------------------------------------------

/** \brief Contains a collection of nodes representing a simulation pipeline
 *
 * The Context class contains a collection of nodes which are identified by a
 * unique name and facilitates saving a loading the simulation pipeline. The
 * class allows for connections between nodes to be serialized by saving the
 * names of connected nodes and enables connections to be restored when
 * loading by generating all nodes before initializing them.
 */
class Context {
public:
    using NodeMap = typename std::map<std::string, std::shared_ptr<Node>>;

    /// A map containing all nodes in the context
    NodeMap nodes;

    /// Create an empty context
    Context() = default;

    /** Create a new node and add it to the context
     *
     * \param type the type of node to create
     * \param name_hint a hint for the name of the node
     * \return a shared pointer to the newly created node
     */
    std::shared_ptr<Node> create(const std::string& type, const std::string& name_hint="");

    /// Add an existing node to the Context
    std::shared_ptr<Node> add(std::shared_ptr<Node> node, const std::string& name_hint="");

    /// Return a pointer to an existing node in the Context
    std::shared_ptr<Node> get(const std::string& name);

    /// Remove a node from the context
    std::shared_ptr<Node> rm(const std::string& name);

    /// Set whether nodes should save as they finish
    // void set_autosave(bool autosave);

    /// Set filename to use for autosave
    // void set_filename(const std::string& filename,
    //        const std::string& mode="r+");

    /// Solve or resolve all nodes
    void solve();

    /// Solve only those nodes which have not been solved
    void make_ready();

    /// Delete all nodes
    void clear();

    /// Save to file
    void save(const std::string& filename) const;

    /// Load from file
    void load(const std::string& filename);

private:
    std::string _filename;
    std::string _mode;
    // bool _autosave = false;
};

//----------------------------------------------------------------------------

} // namespace tbpp


#endif /* __TBPP_CONTEXT_H__ */
