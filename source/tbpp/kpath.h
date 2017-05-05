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

#ifndef __TBPP_KPATH_H__
#define __TBPP_KPATH_H__

/**
 * \file
 * \brief Solve for the eigenvalues and eigenvectors along a path in k-space
 */

#include <tbpp/node.h>
#include <tbpp/model.h>
#include <string>

namespace tbpp {

//----------------------------------------------------------------------------

class KPath;
using KPathPtr = std::shared_ptr<KPath>;

//----------------------------------------------------------------------------

/// Solve for the eigenvalues and eigenvectors along a path in k-space
class KPath : public Node {
public:
    virtual const char* type() const;

    //------------------------------------------------------------------------

    typedef struct KPoint {
        char label[TBPP_MAX_SITE_KIND_LEN];
        double k[3];
    } KPoint;

    //------------------------------------------------------------------------

    /// Steps to take between k-points
    unsigned steps = 50;

    /// Whether to solve for eigenvectors
    bool solve_eigvec = false;

    /// List of k-points
    std::vector<KPoint> kpoints;

    //------------------------------------------------------------------------

    /// Solved Eigenvalues [ik, iE]
    tbpp::math::NArray<double, 2> eigval;

    /// Solved Eigenvectors [ik, iE, i]
    tbpp::math::NArray<std::complex<double>, 3> eigvec;

    /// Value of k-point (reduced coordinates)
    tbpp::math::NArray<double, 2> k;

    //------------------------------------------------------------------------

    /// Set the model to use
    void set_model(ModelPtr);

    /// Returns the model being used (or nullptr is model is not set)
    ModelPtr model() const;

    /// Add a k-point (k1,k2,k3  must be in reduced coordinates)
    void add_kpoint(const std::string& label, double k1, double k2, double k3);

    virtual void solve();

#ifdef TBPP_WITH_HDF5
    virtual void save(EHFile& file, const std::string& prefix="") const;
    virtual void load(EHFile& file, const std::string& prefix="");
#endif // TBPP_WITH_HDF5

private:
    ModelPtr _model = nullptr;
};

//----------------------------------------------------------------------------
} // namespace tbpp

#endif /* __TBPP_KPATH_H__ */
