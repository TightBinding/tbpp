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

#ifndef __TBPP_KMODEL_H__
#define __TBPP_KMODEL_H__

/**
 * \file
 * \brief Represents a Hamiltonian on a discretized k-space
 */

#include <tbpp/model.h>

namespace tbpp {

//----------------------------------------------------------------------------

class KModel;
using KModelPtr = std::shared_ptr<KModel>;

//----------------------------------------------------------------------------

/** \brief Represents a Hamiltonian model on a discretized grid in k-space
 *
 * The KModel class represent a Model on a discretized grid in k-space. The
 * class is mainly used to evaluate a Model Hamiltonian on a grid in k-space
 * such that it can be used with other algorithms which require sampling over
 * the first BZ zone. However the grid does not necessarily need to cover the
 * entire k-space.
 */
class KModel : public tbpp::Node {
public:
    virtual const char* type() const;

    //------------------------------------------------------------------------

    // [out] Hopping matrix Cache Tk[ik1,ik2,ik3,i,j] (only if caching is enabled)
    math::NArray<cxdouble, 5> Tk;

    /// [out] k-space points along b1 (reduced coordinates)
    std::vector<double> k1;

    /// [out] k-space points along b2 (reduced coordinates)
    std::vector<double> k2;

    /// [out] k-space points along b3 (reduced coordinates)
    std::vector<double> k3;

    //------------------------------------------------------------------------

    KModel() = default;
    KModel(ModelPtr model, size_t nk1, size_t nk2, size_t nk3, bool cache=false);

    /// Set the model to use
    void set_model(ModelPtr);

    /// The model being used
    ModelPtr model() const;

    /// Set whether to cache the hopping matrix (Tk)
    void set_cache(bool);

    /// Whether caching of the hopping matrix (Tk) is enabled
    bool cache() const;

    /// Set the grid for the k-space
    void set_grid(size_t nk1=0, size_t nk2=0, size_t nk3=0);

    /// Set points to use along b1 (reduced coordinates)
    void set_k1(const std::vector<double>&);

    /// Set points to use along b2 (reduced coordinates)
    void set_k2(const std::vector<double>&);

    /// Set points to use along b3 (reduced coordinates)
    void set_k3(const std::vector<double>&);

    /// Solve for Tk if caching is enabled
    virtual void solve();

    /** \brief Get the Local Green's Matrix
     *
     * Note that you must first specify the k-space grid for integration by
     * calling set_grid(...) before using this function.
     *
     * \param[out] d Pointer to array of size states*states for G
     * \param[in] Ef Fermi energy to use
     * \param[in] zeroj Value to use for artificial broadening (+0j)
     * \param[in] S Pointer to array of size states*states with self-energy
     */
    void G(std::complex<double>* d, double Ef, double zeroj,
            std::complex<double>* S) const;


#ifdef TBPP_WITH_HDF5
    virtual void save(EHFile& file, const std::string& prefix="") const;
    virtual void load(EHFile& file, const std::string& prefix="");
#endif // TBPP_WITH_HDF5

private:
    ModelPtr _model = nullptr;
    bool _cache = false;
};

//----------------------------------------------------------------------------

} // namespace tbpp

#endif /* __TBPP_KMODEL_H__ */
