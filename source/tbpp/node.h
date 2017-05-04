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

#ifndef __TBPP_NODE_H__
#define __TBPP_NODE_H__

/**
 * \file
 * \brief Node class which all TBPP algorithms inherit.
 */

#include <tbpp/common.h>
#include <tbpp/context.h>
#include <string>
#include <memory>

#ifdef TBPP_WITH_HDF5
#include <tbpp/ehfile.h>
#endif // TBPP_WITH_HDF5

namespace tbpp {

//----------------------------------------------------------------------------

class Node;
typedef std::shared_ptr<Node> NodePtr;

//----------------------------------------------------------------------------

enum class NodeStatus {
    not_done,
    done,
    failed
};

//----------------------------------------------------------------------------

/**
 * \brief Node class which all TBPP algorithms inherit
 *
 * The Node class contains basic functionality for running, checking the
 * status of and saving TBPP algorithms. It likewise handles the creation of
 * nodes from a string specifying its type.
 */
class Node {
    friend class Context;
public:

    //------------------------------------------------------------------------

    /// Convert NodeStatus to a string
    static std::string node_status_to_string(NodeStatus status);

    /// Determine NodeStatus from a string
    static NodeStatus node_status_from_string(const std::string& status);

    /// Create a new node from a string describing its type
    static std::shared_ptr<Node> new_node(const std::string& type);

    //------------------------------------------------------------------------

    /// The name of the node (only if the node is in a context)
    const std::string& name() const;

    /// Pointer to the context the node is in (otherwise is nullptr)
    Context* context() const;

    /// The current status of the node
    NodeStatus status() const;

    /// Whether to use multithreading
    bool parallel() const;

    /// Set whether to use multithreading
    void set_parallel(bool);

    /// Solve the node only if it is not solved
    void make_ready();

    /// Clear stored data and set status to not_done
    void reset();

    //------------------------------------------------------------------------
    // To be implemented by each node

    /// The type of the node
    virtual const char* type() const;

    /// Solve the node regardless of whether it has already been solved
    virtual void solve();

    /// Clear stored data
    virtual void clear();

    #ifdef TBPP_WITH_HDF5
    /// Save the node to an EHFile under the provided prefix
    virtual void save(EHFile& file, const std::string& prefix="") const;

    /// Load the node from an EHFile from the provided prefix
    virtual void load(EHFile& file, const std::string& prefix="");
    #endif // TBPP_WITH_HDF5

    //------------------------------------------------------------------------

protected:
    /// Print the header for the node
    void header() const;
    /// Print debug information
    void debug(const std::string&) const;
    /// Print a message
    void msg(const std::string&) const;
    /// Warn of a problem
    void warn(const std::string&) const;
    /// Throw an exception
    void error(const std::string&) const;

    /// Current status of the node
    NodeStatus _status = NodeStatus::not_done;
    /// Whether to use multi-threading
    bool _parallel = true;

private:
    /// Name of node if node is part of context.
    std::string _name;
    /// Pointer to context node it in otherwise is nullptr.
    Context* _ct = nullptr;
};

//----------------------------------------------------------------------------

} // namespace tbpp

#endif /* __TBPP_NODE_H__ */
