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

#ifndef __TBPP_LATTICE_H__
#define __TBPP_LATTICE_H__

/**
 * \file
 * \brief A class representing a list of sites
 */

#include <memory>
#include <tbpp/narray.h>
#include <tbpp/matrix.h>
#include <tbpp/node.h>
#include <cstdint>
#include <iostream>
#include <string>
#include <list>
#include <deque>

namespace tbpp {
//----------------------------------------------------------------------------

class Lattice;
using LatticePtr = std::shared_ptr<Lattice>;
using Spinor = std::array<cxdouble,2>;

//----------------------------------------------------------------------------

/** \brief Contains a list of sites and facilitates their manipulation.
 *
 * The Lattice class represent a list of sites and facilitates the
 * manipulation of the sites prior to determining the Hamiltonian. It is
 * mainly used as an input for the Model class. The Lattice class also encodes
 * the real-space hopping parameters.
 */
class Lattice : public Node {
public:
    virtual const char* type() const;

    //------------------------------------------------------------------------

    /// Information contained for each site
    typedef struct Site {
        std::string kind;                   ///< String identifier
        double pos[3];                      ///< x,y,z position
        math::DenseMatrix<cxdouble> onsite; ///< onsite potential
        std::vector<uint32_t> orbits;       ///< orbit for each state
        std::vector<Spinor> spinors;        ///< spinors for each state
    } Site;

    //------------------------------------------------------------------------

    /// Information contained for each hop
    typedef struct Hop {
        std::string initial_kind;  ///< String identifier for initial site
        std::string final_kind;    ///< String identifier for final site
        double dpos[3];    ///< change in position (pos_{final} - pos_{initial})
        math::DenseMatrix<cxdouble> T; ///< hop matrix
    } Hop;

    //------------------------------------------------------------------------

    using HopMap = std::map<std::string, std::map<std::string, std::list<Hop>>>;
    using TransVec = std::vector<std::array<double,3>>;

    //------------------------------------------------------------------------

    /// List of sites
    std::deque<Site> sites;

    /** \brief Map of hops
     *
     * The hops are stored in a format hops[kind1][kind2][hop_id] where
     * kind1 < kind2 and where hop_id is an integer that iterates through all of the
     * hops from kind1 to kind2.
     */
    HopMap hops;

    /// Translational symmetry directions
    TransVec trans;

    /// Tolerance for computing hops
    double eps = 1e-8;

    /// Maximum range for hop computation (automatically set)
    double max_hop_range = 0.0;

    /// External Uniform Magnetic Field
    std::array<double,3> B;

    //------------------------------------------------------------------------

    /** \brief Low level interface to add a new site
     *
     * \param[in] kind String identifier for the type of the site
     * \param[in] x Position along x axis of the site
     * \param[in] y Position along y axis of the site
     * \param[in] z Position along z axis of the site
     * \param[in] states Total number of states
     * \param[in] onsite Pointer to row-majored onsite matrix of size states by states
     * \param[in] orbits Pointer to vector of site state with orbit number of each state
     * \param[in] spinors Pointer to vector of size 2*state with spinors for each state
     */
    void add_site(const std::string& kind, double x, double y, double z,
            uint32_t states, const std::complex<double> *onsite,
            const uint32_t *orbits, const std::complex<double> *spinors);

    /** \brief Add a new site
     *
     * For example suppose we have a site of type "A" at position
     * [0.0,0.0,0.0] with a total of four states (two orbits and two spins per
     * orbit). Likewise there is a site of type "B" at position [0.5,0.0,0.0]
     * with a total of two states (one orbit and two spins per orbit). We
     * would define this system as follows,
     * \code
     * // Orbit ID for each state: first state belongs to 0 orbit, second
     * // belongs to 0 orbit, third belongs to first orbit, fourth belongs to
     * // first orbit
     * std::vector<uint32_t> A_orbits = {0,0,1,1};
     *
     * // Spinor for each state: first state is [1,0] spin up along z, second is [0,1]
     * // spin down along z, third is [1,1] spin up along x, fourth is [1,-1] spin down along x,
     * std::vector<Spinor> A_spinors = {{1,0}, {0,1}, {1,1}, {1,-1}};
     *
     * // Onsite potential for A
     * cxdouble EA(1,0);
     * DenseMatrix<tbpp::cxdouble> VA(4,4);
     * VA(0,0) = EA;
     * VA(1,1) = EA;
     * VA(2,2) = EA;
     * VA(3,3) = EA;
     *
     * // Onsite potential for B
     * cxdouble EB(2,0);
     * DenseMatrix<tbpp::cxdouble> VB(2,2);
     * VB(0,0) = EB;
     * VB(1,1) = EB;
     *
     * tbpp::Lattice sites;
     *
     * // Add site A
     * sites.add_site("A", 0.0, 0.0, 0.0, VA, A_orbits, A_spinors);
     *
     * // Add site B (can also provide orbits and spinors using initializer list)
     * sites.add_site("B", 0.5, 0.0, 0.0, VB, {0,0}, {{1,0}, {0,1}});
     * \endcode
     *
     * \param[in] kind String identifier for the type of the site
     * \param[in] x Position along x axis of the site
     * \param[in] y Position along y axis of the site
     * \param[in] z Position along z axis of the site
     * \param[in] onsite Onsite potential matrix
     * \param[in] orbits Orbit ID for each state
     * \param[in] spinors Spinor for each state
     */
    void add_site(const std::string& kind, double x, double y, double z,
            const math::DenseMatrix<cxdouble>& onsite,
            const std::vector<uint32_t>& orbits,
            const std::vector<Spinor>& spinors);

    void copy_along(double vx, double vy, double vz, unsigned n);
    void move_along(double vx, double vy, double vz);
    // TODO void rotate(const double* R);
    // TODO void clean_sites();

    // TODO void add_sites(Lattice);
    // TODO void filter(f);

    //------------------------------------------------------------------------

    /** \brief Low level interface to add a new hop
     *
     * \param initial_kind The kind of the initial site
     * \param final_kind The kind of the final site
     * \param dx The change in position x_{final} - x_{initial}
     * \param dy The change in position y_{final} - y_{initial}
     * \param dz The change in position z_{final} - z_{initial}
     * \param t_rows Number of rows of T matrix (number of state in final_kind)
     * \param t_cols Number of columns of T matrix  (number of state in initial_kind)
     * \param T Pointer to row-major ordered hop matrix of size (t_rows, t_cols)
     */
    void add_hop(const std::string& initial_kind, const std::string& final_kind,
            double dx, double dy, double dz, uint32_t t_rows, uint32_t t_cols,
            const std::complex<double>* T);

    /// Add a new hop
    void add_hop(const std::string& initial_kind, const std::string& final_kind,
            double dx, double dy, double dz, const math::DenseMatrix<cxdouble>& T);

    //------------------------------------------------------------------------

    /// Add a translational symmetry direction
    void add_trans(double x, double y, double z);

    /// Set external magnetic field
    void set_B(double Bx, double By, double Bz);

    //------------------------------------------------------------------------

    /** Returns the hop matrix from site_i and site_j where dx,dy,dz is the
     * total distance between the two sites in real space
     *
     * \warning does not include effects due to magnetic field
     * \warning should not be used except by the Model class
     */
    math::DenseMatrix<cxdouble> _T(uint32_t site_i, uint32_t site_j,
            double dx, double dy, double dz) const;
};

//----------------------------------------------------------------------------
} // namespace tbpp

#endif /* __TBPP_LATTICE_H__ */
