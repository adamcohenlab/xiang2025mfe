# xiang2025mfe

### *Supporting Code*

#### Mechanism of giant magnetic field effect in a red fluorescent protein

Katherine M. Xiang<sup>1</sup>, Hana Lampson<sup>1</sup>, Rebecca Frank Hayward<sup>2</sup>, Andrew G. York<sup>4</sup>, Maria Ingaramo<sup>4</sup>, Adam E. Cohen<sup>1,3</sup> * 

<sup>1</sup>Department of Physics, Harvard University, Cambridge, MA 02138, USA

<sup>2</sup>School of Engineering and Applied Sciences, Harvard University, Cambridge, MA 02138, USA

<sup>3</sup>Department of Chemistry and Chemical Biology, Harvard University, Cambridge, MA 02138, USA

<sup>4</sup>Calico Life Sciences LLC, South San Francisco, CA 94080, USA

---

#### Description of files:

`MFE_kinetic_model.m` : MATLAB script for simulating the kinetic photophysical model of the mScarlet3/FMN/FMNH~2~ magnetic field effect. 

`SCRP_spin_sim.m` : MATLAB script for simulating spin dynamics. Includes code for simulating a single spin precessing in an arbitrary static magnetic field, and for an uninteracting two-spin system in separate arbitrary static magnetic fields. 

`SCRP_spin_and_chem_sim.m` : MATLAB script for incorporating chemical reactions to the spin simulations to simulate SCRP separation and recombination. The function `spin_chem_sim` numerically computes the time evolution, while the function `spin_chem_sim_analytical` uses the analytical expressions derived in `spin_dynamics.nb` which is faster.

`spin_dynamics.nb` : Mathematica notebook with derivations for the time evolution of  $T_0$ and $T_+$-born SCRPs in arbitrary static magnetic fields, and the corresponding singlet and triplet probabilities as a function of time. 
