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
#ifndef __TBPP_COND_DOS_H__
#define __TBPP_COND_DOS_H__

/**
 * \file
 * \brief Compute DOS and DC electrical conductivity using self-energy
 */

#include <tbpp/cpa.h>

namespace tbpp {

//----------------------------------------------------------------------------

class CondDOS;
using CondDOSPtr = std::shared_ptr<CondDOS>;

//----------------------------------------------------------------------------

/**
 * \brief Compute DOS and DC electrical conductivity using self-energy
 *
 * The CondDOS class allows for computing the Density of States and zero
 * temperature DC electrical conductivity using Green's function methods and
 * the self-energy returned by the CPA module. The electrical conductivity is
 * determined using the Kubo-Greenwood formalism.
 */
class CondDOS : public tbpp::Node {
public:
    virtual const char* type() const;

    //------------------------------------------------------------------------

    /// DC conductivity tensor (cond[ic,iw,i,j]) where sig_ij
    tbpp::math::NArray<double, 4> cond;

    /// DC conductivity contribution per state (cond_state[ic,iw,i,j,is])
    tbpp::math::NArray<double, 5> cond_state;

    /// Density of States per k-point (dos_k[ic,iw,ik1,ik2,ik3])
    tbpp::math::NArray<double, 5> dos_k;

    /// Total Density of States (dos[ic,iw])
    tbpp::math::NArray<double, 2> dos;

    /// Density of states per layer (dos[ic,iw,is])
    tbpp::math::NArray<double, 3> dos_state;

    /// Density of states projected dos_weights (dos_w[ic,iw,ip,ik1,ik2,ik3])
    tbpp::math::NArray<double, 6> dos_proj;

    //------------------------------------------------------------------------

    /// Whether to solve for conductivity
    bool solve_cond = true;
    /// Whether to solve conductivity along x-axis
    bool solve_sigx = true;
    /// Whether to solve conductivity along y-axis
    bool solve_sigy = true;
    /// Whether to solve conductivity along z-axis
    bool solve_sigz = true;

    /// Whether to solve for total density of state
    bool solve_dos = true;
    /// Whether to solve for the density of states per layer
    bool solve_dos_state = true;
    /// Whether to solve for density of state per k-point
    bool solve_dos_k = true;

    /// Minimum value for +0j when computing Greens matrix
    double zero_j = 1e-4;

    /// Weights for DOS projections (dos_weights[ip,is])
    tbpp::math::NArray<double, 2> dos_weights;

    //------------------------------------------------------------------------

    // Add a set of DOS projection weights
    void add_dos_weights(const std::vector<double>& weights);

    /// Set the self-energy
    void set_sigma(CPAPtr sigma);

    /// Get the self-energy
    CPAPtr sigma() const;

    virtual void solve();

#ifdef TBPP_WITH_HDF5
    virtual void save(EHFile& file, const std::string& prefix="") const;
    virtual void load(EHFile& file, const std::string& prefix="");
#endif // TBPP_WITH_HDF5

private:
    /// Self-Energy to use
    CPAPtr _sigma = nullptr;
};

//----------------------------------------------------------------------------

} // namespace tbpp

#endif /* __TBPP_COND_DOS_H__ */
