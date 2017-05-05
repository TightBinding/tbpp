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

#ifndef __TBPP_CPA_H__
#define __TBPP_CPA_H__

/**
 * \file
 * \brief Model random substitutions using the Coherent Potential Approximation
 */

#include <tbpp/node.h>
#include <tbpp/kmodel.h>
#include <string>
#include <list>
#include <map>

namespace tbpp {

//----------------------------------------------------------------------------

class CPA;
using CPAPtr = std::shared_ptr<CPA>;

//----------------------------------------------------------------------------

/** \brief Coherent Potential Approximation Solver
 *
 * This class solves for the self-energy of a model using the Coherent
 * Potential Approximation (CPA) for a list of frequencies and
 * concentration distributions.
 */
class CPA : public Node {
public:
    virtual const char* type() const;

    //------------------------------------------------------------------------

    /// Structure containing information for each defect type
    typedef struct Defect {
        /// List of site index where defect occurs
        std::vector<uint32_t> sites;
        /// Concentration of defect (one for each concentration run)
        std::vector<double> c;
        /// On-site potential of the defect
        math::DenseMatrix<cxdouble> V;
    } Defect;

    /// [in] List of defects
    std::list<Defect> defects;

    //------------------------------------------------------------------------
    // Input Parameters

    /// [in] Frequencies to solve the self-energy
    std::vector<double> w;

    /// [in] Precision to use in CPA self-consistency cycle
    double eps = 1e-4;

    /// [in] Value for +0j to use for Green's function calculation
    double zero_j = 1e-3;

    /// [in] Max number of iterations for CPA self-consistent cycle
    unsigned max_iter = 200;

    /// [in] Mixing factor to use for CPA (lower values are more stable)
    double mix_factor = 0.9;

    /// [in] Whether to use Virtual Crystal Approximation for self-energy
    bool use_vca = false;

    //------------------------------------------------------------------------
    // Output Data

    /// [out] Self-Energy Matrix [ic,iw,i,j]
    tbpp::math::NArray<cxdouble, 4> sigma;

    /// [out] Virtual Crystal Approximation [ic,i,j]
    tbpp::math::NArray<cxdouble, 3> vca;

    //------------------------------------------------------------------------
    // Setters

    /// Set the KModel to use
    void set_kmodel(KModelPtr);

    /// Set the frequencies (Fermi energies) to solve for
    void set_w(const std::vector<double>&);

    /** \brief Add a defect
     *
     * For each defect, one must provide a list of sites where the defect
     * occurs, the on-site potential of the defect and the concentration of
     * the defect for each concentration scan.
     *
     * For example, the code below adds a defect with on-site potential Ea to
     * sites 1 and 2 with concentrations [0.1,0.2,0.3] and a defect with
     * on-site potential Eb to sites 1 and 2 with concentration [0.4,0.5,0.6].
     *
     * \code
     *
     * math::DenseMatrix<cxdouble>& mEa(1,1);
     * mEa[0] = Ea;
     *
     * math::DenseMatrix<cxdouble>& mEb(1,1);
     * mEb[0] = Eb;
     *
     * std::vector<double> ca = {0.1, 0.2, 0.3};
     * std::vector<double> cb = {0.4, 0.5, 0.6};
     * std::vector<uint32_t> sites = {1,2};
     *
     * cpa->add_defect(sites, ca, mEa);
     * cpa->add_defect(sites, cb, mEb);
     *
     * \endcode
     *
     * The CPA solver will then run three concentration scans, where on the
     * first concentration scan (ic=0) site 1 and 2 have a 10% probability of
     * having the defect Ea, a 40% probability of having defect Eb and a 50%
     * probability of having the original potential, on the second
     * concentration scan (ic=1) site 1 and 2 have a 20% probability of having
     * defect Ea, a 50% probability of having defect Eb and hence a 30%
     * probability of having the original potential.
     *
     * The length of \a c must be the same for all of the defects and
     * represents the number of concentration scans the CPA solver will
     * perform.
     *
     * The size of \a V must be the same as the on-site potential matrix of
     * the sites were the defect occurs.
     *
     * \param sites List of site indexes where defect occurs
     * \param c List of concentrations (one for each concentration scan)
     * \param V The on-site potential for the defect
     */
    void add_defect(const std::vector<uint32_t>& sites, const std::vector<double>& c,
            const math::DenseMatrix<cxdouble>& V);

    // TODO void add_defect(size_t site_index, double c, cxdouble *V, size_t states);

    //------------------------------------------------------------------------
    // Getters

    /// Returns the KModel being used
    KModelPtr kmodel() const;

    //------------------------------------------------------------------------
    // Methods

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

#endif /* __TBPP_CPA_H__ */
