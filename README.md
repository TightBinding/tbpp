# TightBinding++

TightBinding++ is a framework for simulating quantum tight-binding models.

*The software is currently under heavy development.*

The main abilities are:
- Automated generation of the Hamiltonian matrix.
- Calculate band structure and density of states.
- Simulating random on-site disorder using Coherent Potential Approximation (CPA).
- Compute linear response electrical conductivity using Kubo-Greenwood Formalism.

The main features include:
- Python 3 wrapper allowing for all algorithms to be used from python.
- Algorithms are entirely written in C++11 and can be used in new C++ code by
  linking against the TBPP library.
- Saving and loading parameters and results using the HDF5 file format.
- OpenMP multi-threading (requires compiler support).
- Graphical interface simplifying the creation, inspection and post-processing
  of simulations.

## License
TightBinding++ is licensed under the GNU General Public License version 3.
See [license](LICENSE) for the full license.
