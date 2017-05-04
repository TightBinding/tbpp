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
#ifndef __TBPP_GEN_COND_H__
#define __TBPP_GEN_COND_H__

/**
 * \file
 * \brief Solve for zero temperature AC electrical conductivity [under development]
 */

#include <tbpp/cpa.h>

namespace tbpp {

//----------------------------------------------------------------------------

class GenCond;
using GenCondPtr = std::shared_ptr<GenCond>;

//----------------------------------------------------------------------------

/** \brief Solve for zero temperature AC electrical conductivity [under development]
 *
 * The GenCond class solves for the AC electrical conductivity at zero
 * temperature using the Kubo-Greenwood formalism.
 *
 * It is currently under development and not ready for production use.
 */
class GenCond : public tbpp::Node {
public:
    virtual const char* type() const;

    //------------------------------------------------------------------------

    /// AC conductivity tensor (cond[ic,iw,i,j]) where sig_ij
    tbpp::math::NArray<double, 4> cond;

    /// phi[ic,iw,jw, i,j]
    tbpp::math::NArray<std::complex<double>, 5> phi;

    //------------------------------------------------------------------------

    /// Whether to solve conductivity along x-axis
    bool solve_sigx = true;
    /// Whether to solve conductivity along y-axis
    bool solve_sigy = true;
    /// Whether to solve conductivity along z-axis
    bool solve_sigz = true;

    /// Minimum value for +0j when computing Greens matrix
    double zero_j = 1e-4;

    //------------------------------------------------------------------------

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

#endif /* __TBPP_GEN_COND_H__ */
