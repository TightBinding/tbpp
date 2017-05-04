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

#ifndef __TBPP_MODEL_H__
#define __TBPP_MODEL_H__

/**
 * \file
 * \brief A class representing a tight-binding Hamiltonian
 */

#include <tbpp/lattice.h>

#define TBPP_MAX_SITE_KIND_LEN 64

namespace tbpp {

//----------------------------------------------------------------------------

class Model;
using ModelPtr = std::shared_ptr<Model>;

//----------------------------------------------------------------------------

/**
 * \brief A class representing a tight-binding Hamiltonian
 *
 * The Model class is the central class representing a tight-binding
 * Hamiltonian. The class is usually initialized by using a Lattice class but
 * can also be manually initialized.
 *
 * Once initialized, the class provides methods for determining the Hamiltonian
 * matrix at each k-point as well as derivatives of the Hamiltonian.
 */
class Model : public tbpp::Node {
public:
    virtual const char* type() const;

    //------------------------------------------------------------------------

    /// Structure containing information for each site of the model
    typedef struct SiteInfo {
        char kind[TBPP_MAX_SITE_KIND_LEN]; ///< String identifier
        double pos[3]; ///< Position in real space
        uint64_t states; ///< Total number of states
        uint64_t si;     ///< Start index
        uint64_t ei;     ///< End index
    } SiteInfo;

    /// Structure containing information for each state of the model
    typedef struct StateInfo {
        uint64_t site;
        uint64_t orbit;
        Spinor spinor;
    } StateInfo;

    //------------------------------------------------------------------------

    /// Information for each site
    std::vector<SiteInfo> site_info;

    /// For each state gives the corresponding site, orbit and spinor
    std::vector<StateInfo> state_info;

    /// Matrix with rows containing lattice vectors a[ir,i]
    tbpp::math::NArray<double, 2> a;

    /// Matrix with rows containing reciprocal lattice vectors b[ik,i]
    tbpp::math::NArray<double, 2> b;

    /// On-site potential V[i,j]
    tbpp::math::NArray<std::complex<double>, 2> V;

    /// Hopping matrices T[iR,i,j]
    tbpp::math::NArray<std::complex<double>, 3> T;

    /// Hopping indices in reduced coordinates R[iR,i]
    tbpp::math::NArray<int32_t, 2> R;

    //------------------------------------------------------------------------

    Model();

    //------------------------------------------------------------------------
    // Input type

    /// Whether to use lattice or whether data is manually input
    bool use_lattice() const;

    /// Set whether to use lattice or whether data is manually input
    void set_use_lattice(bool);

    //------------------------------------------------------------------------
    // Lattice

    /// Get the lattice
    LatticePtr lattice() const;

    /// Set the lattice (also runs set_use_lattice)
    bool set_lattice(LatticePtr lattice);

    //------------------------------------------------------------------------
    // General Settings

    /// Dimension of k-space
    size_t kdim() const;

    /// Set the dimensionality of k-space
    bool set_kdim(unsigned);

    /// Number of sites
    size_t sites() const;

    /// Set the total number of sites
    void set_sites(unsigned);

    /// Total Number of states
    size_t states() const;

    /// Set the total number of states in the system
    void set_states(unsigned);

    /// Number of hopping matrices
    size_t hops() const;

    /// Set the total number of hopping matrices
    void set_hops(unsigned);

    /// Set a1 lattice vector
    void set_a1(double x, double y, double z);

    /// Set a2 lattice vector
    void set_a2(double x, double y, double z);

    /// Set a3 lattice vector
    void set_a3(double x, double y, double z);

    //------------------------------------------------------------------------

    /// Solve the model
    virtual void solve();

    //------------------------------------------------------------------------
    // Methods after being solved

    /// Returns the volume of the unit cell in real-space
    double uc_size() const;

    /** Returns the hopping matrix for a specific k-point
     *
     * Note that k1,k2,k3 are in reduced coordinates. This means that the
     * first BZ is from -0.5 < k1,k2,k3 <0.5 and likewise that
     *
     *          k = k1*b1 + k2*b2 + k3*b3
     *
     * \param[out] d Pointer to array of size states*states
     * \param[in] k1 k-point in reduced coordinates along b1
     * \param[in] k2 k-point in reduced coordinates along b2
     * \param[in] k3 k-point in reduced coordinates along b3
     */
    void Tk(std::complex<double> *d, double k1, double k2, double k3) const;

    /** \brief Get the Hamiltonian at a specific k-point
     *
     * Note that k1,k2,k3 are in reduced coordinates. This means that the
     * first BZ is from -0.5 < k1,k2,k3 <0.5.
     *
     * \param[out] d Pointer to array of size states*states
     * \param[in] k1 k-point in reduced coordinates along b1
     * \param[in] k2 k-point in reduced coordinates along b2
     * \param[in] k3 k-point in reduced coordinates along b3
     */
    void Hk(std::complex<double>* d, double k1, double k2, double k3) const;

    /** \brief Get the derivative of Hamiltonian along x,y,z
     *
     * Note that k1,k2,k3 are in reduced coordinates. This means that the
     * first BZ is from -0.5 < k1,k2,k3 <0.5.
     *
     * If dx,dy or dz is NULL, then the derivative is not computed along that
     * axis.
     *
     * \param[out] dx Pointer to array of size states*states for dH_dkx
     * \param[out] dy Pointer to array of size states*states for dH_dky
     * \param[out] dz Pointer to array of size states*states for dH_dkz
     * \param[in] k1 k-point in reduced coordinates along b1
     * \param[in] k2 k-point in reduced coordinates along b2
     * \param[in] k3 k-point in reduced coordinates along b3
     */
    void dH_dk(std::complex<double>* dx, std::complex<double>* dy, std::complex<double>* dz,
            double k1, double k2, double k3) const;

#ifdef TBPP_WITH_HDF5
    virtual void save(EHFile& file, const std::string& prefix="") const;
    virtual void load(EHFile& file, const std::string& prefix="");
#endif // TBPP_WITH_HDF5

private:
    /// Whether to use Lattice
    bool _use_lattice = true;

    /// Lattice
    LatticePtr _lattice = nullptr;

    /// Number of k-space axis
    size_t _kdim = 0;

    /// Number of sites
    size_t _sites = 0;

    /// Number of states
    size_t _states = 0;

    /// Number of hopping matrices
    size_t _hops = 0;
};

//----------------------------------------------------------------------------

} // namespace tbpp


#endif /* __TBPP_MODEL_H__ */
