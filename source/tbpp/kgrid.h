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

#ifndef __TBPP_KGRID_H__
#define __TBPP_KGRID_H__

/**
 * \file
 * \brief Solve for the eigenvalues and eigenvectors on a grid in k-space
 */

#include <tbpp/node.h>
#include <tbpp/kmodel.h>
#include <string>

namespace tbpp {

//----------------------------------------------------------------------------

class KGrid;
using KGridPtr = std::shared_ptr<KGrid>;

//----------------------------------------------------------------------------

/// Solve for the eigenvalues and eigenvectors on a grid in k-space
class KGrid : public Node {
public:
    virtual const char* type() const;

    //------------------------------------------------------------------------

    /// Whether to solve for eigenvectors
    bool solve_eigvec = false;

    //------------------------------------------------------------------------

    /// Solved Eigenvalues [ik1, ik2, ik3, iE]
    tbpp::math::NArray<double, 4> eigval;

    /// Solved Eigenvectors [ik1, ik2, ik3, iE, i]
    tbpp::math::NArray<std::complex<double>, 5> eigvec;

    //------------------------------------------------------------------------

    /// Set the KModel to use
    void set_kmodel(KModelPtr);
    /// Returns the KModel being used (or nullptr)
    KModelPtr kmodel() const;

    virtual void solve();

#ifdef TBPP_WITH_HDF5
    virtual void save(EHFile& file, const std::string& prefix="") const;
    virtual void load(EHFile& file, const std::string& prefix="");
#endif // TBPP_WITH_HDF5

private:
    KModelPtr _kmodel = nullptr;
};

//----------------------------------------------------------------------------
} // namespace tbpp

#endif /* __TBPP_KGRID_H__ */
