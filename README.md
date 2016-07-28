# <a name="top"></a>axionCAMB
*[axionCAMB and CAMB](#intro)
*[Basic Usage](#basics)
*[Physics](#physics)
*[Known Issues and Warnings](#warnings)

----------------------------------------------------------------------
#### <a name="intro"></a>axionCAMB and CAMB

CODE Authors:

* Daniel Grin
* David JE Marsh
* Renee Hlozek

axionCAMB is a modified version of the publicly available code "CAMB," which is available at http://camb.info/ and on GitHub at https://github.com/cmbant/CAMB
axionCAMB computes cosmological observables for comparison with data. This is normally the CMB power spectra (T,E,B,\phi in auto and cross power), but also includes the matter power spectrum. 
We try and keep axionCAMB up to date with developments in CAMB, but the codes are independent, so they may not always be equivalent. The "base" version of CAMB that axionCAMB is built off is Nov13. 
The physics of axionCAMB is described in detail in http://arxiv.org/abs/1410.2896 (HGMF). When using this code, you should cite HGMF, and the original version of CAMB.
For parameter estimation, axionCAMB is compatible with cosmomc, with minimal modifications. However, you should be careful with sampling and degeneracies, as described in HGMF.
axionCAMB is also compatible as a module in cosmosis (https://bitbucket.org/joezuntz/cosmosis/wiki/Home) when made with the MakefileCosmosis Makefile.

----------------------------------------------------------------------
#### <a name="basics"></a>Basic Usage

Refer to CAMB documentation for all non-axion related things.
axionCAMB adds one new fluid component: the axion fluid. This is included in equations_ppf.f90.
The axion fluid is consistently incorporated into the evolution of all background and perturbation variables, including, for example, recombination in recfast_axion.f90.
There are two new parameters: the axion mass, ma (measured in eV), and the axion energy density, omaxh2. These are entered in params.ini as usual.
Sensible values for the axion mass range between 1e-33 and 1e-18, though the code works for values outside of this range.
The axion fluid can function as either a dark matter component or a dark energy component depending on the mass.

axionCAMB can be used in conjunction with the additional neutrino parameters of neutrino mass and number of species.
axionCAMB cannot (currently) be used in conjunction with additional dark energy parameters, such as the equation of state.

----------------------------------------------------------------------
#### <a name="physics"></a>Physics

At early times, the fluid equations solved are equivalent to the first order perturbed Klein-Gordon equation in synchronous gauge.
At late times, the WKB approximation matches the axion field to an effective fluid with equation of state zero, and a scale-dependent sound speed.
The energy density is found using a shooting method from an initial guess for the axion field value. 
The background evolution is computed in axion_background.f90, which creates an array for the axion equation of state. This array is carried around everywhere the energy densities are needed.

----------------------------------------------------------------------
#### <a name="warnings"></a>Known Issues and Warnings

You should be careful if you use any non-linear options (setting do_nonlinear \= 0). axionCAMB uses halofit_ppf.f90 with various modifications. The treatment of axions is not expected to be quantitatively correct, and various approximations are made (these will documented at a later date). The default version of halofit is set to the original "Smith" version, as this is expected to be most stable.
You should also be careful of galaxy bias if you use galaxy survey data, as discussed in HGMF.
Also please not that modern recombination codes (HyREC and CosmoREC) have not been appropriately modified  to include axions, but RecFAST has. This would be a useful fork for 
someone to pursue.
