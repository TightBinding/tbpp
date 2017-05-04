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
#ifndef __TBPP_APPROX_COND_DOS_H__
#define __TBPP_APPROX_COND_DOS_H__

/**
 * \file
 * \brief Compute DOS and electrical conductivity using artificial broadening
 */


#include <tbpp/kgrid.h>

namespace tbpp {

//----------------------------------------------------------------------------

class ApproxCondDOS;
using ApproxCondDOSPtr = std::shared_ptr<ApproxCondDOS>;

//----------------------------------------------------------------------------

/** \brief Compute DOS and electrical conductivity using an artificial broadening
 *
 * ApproxCondDOS solves for the Density of States and the zero temperature DC
 * electrical conductivity for a list of Fermi energies and artificial
 * broadenings. The electrical conductivity is computed using the
 * Kubo-Greenwood formalism.
 *
 * \warning Assumes KModel has a uniformly spaced grid in k-space
 */
class ApproxCondDOS : public tbpp::Node {
public:
    virtual const char* type() const;

    //------------------------------------------------------------------------

    /// DC conductivity tensor (cond[ieta,iw,i,j]) where sig_ij
    tbpp::math::NArray<double, 4> cond;

    /// Density of States per k-point (dos_k[ieta,iw,ik1,ik2,ik3])
    tbpp::math::NArray<double, 5> dos_k;

    /// Total Density of States (dos[ieta,iw])
    tbpp::math::NArray<double, 2> dos;

    //------------------------------------------------------------------------

    /// Whether to solve for conductivity
    bool solve_cond = true;

    /// Whether to solve for density of state per k-point
    bool solve_dos_k = true;

    /// Whether to solve for total density of state
    bool solve_dos = true;

    /// Frequencies which to compute
    std::vector<double> w;

    /// Artificial Broadening
    std::vector<double> eta;

    //------------------------------------------------------------------------

    // Set the frequencies to solve for
    void set_w(std::vector<double>);

    // Set the artificial broadenings to solve for
    void set_eta(std::vector<double>);

    /// Set the self-energy
    void set_kgrid(KGridPtr kgrid);

    /// Get the self-energy
    KGridPtr kgrid() const;

    virtual void solve();

#ifdef TBPP_WITH_HDF5
    virtual void save(EHFile& file, const std::string& prefix="") const;
    virtual void load(EHFile& file, const std::string& prefix="");
#endif // TBPP_WITH_HDF5

private:
    /// Self-Energy to use
    KGridPtr _kgrid = nullptr;
};

//----------------------------------------------------------------------------

} // namespace tbpp

#endif /* __TBPP_APPROX_COND_DOS_H__ */
