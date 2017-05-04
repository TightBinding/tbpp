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

#ifndef __TBPP_BERRYFLUX_H__
#define __TBPP_BERRYFLUX_H__

/**
 * \file
 * \brief Compute the Berry Flux in k-space [under development]
 */

#include <tbpp/node.h>
#include <tbpp/kgrid.h>
#include <string>

namespace tbpp {

//----------------------------------------------------------------------------

class BerryFlux;
using BerryFluxPtr = std::shared_ptr<BerryFlux>;

//----------------------------------------------------------------------------

/// Compute the Berry Flux in k-space
class BerryFlux : public Node {
public:
    virtual const char* type() const;

    //------------------------------------------------------------------------

    /// Solved Berry Flux [ik1, ik2, ik3, iE, i]
    tbpp::math::NArray<double, 5> flux_k;
    /// Integral of Berry Flux for each energy level [iE, i]
    tbpp::math::NArray<double, 2> flux;

    //------------------------------------------------------------------------

    bool solve_flux_k = false;
    double dE = 1e-20;

    //------------------------------------------------------------------------

    void set_kgrid(KGridPtr);
    KGridPtr kgrid() const;

    virtual void solve();

#ifdef TBPP_WITH_HDF5
    virtual void save(EHFile& file, const std::string& prefix="") const;
    virtual void load(EHFile& file, const std::string& prefix="");
#endif // TBPP_WITH_HDF5

private:
    KGridPtr _kgrid = nullptr;
};

//----------------------------------------------------------------------------
} // namespace tbpp

#endif /* __TBPP_BERRYFLUX_H__ */
