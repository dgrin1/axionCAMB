# <a name="top"></a>AxionCamb 2.0
*[AxionCamb and CAMB](#intro)
*[Basic Usage](#basics)
*[Physics](#physics)
*[Known Issues and Warnings](#warnings)

----------------------------------------------------------------------
#### <a name="intro"></a>AxionCamb and CAMB

AxionCamb is a modified version of the publicly available code "CAMB," which is available at http://camb.info/ and on GitHub at https://github.com/cmbant/CAMB
AxionCamb computes cosmological observables for comparison with data. This is normally the CMB power spectra (T,E,B,\phi in auto and cross power), but also includes the matter power spectrum. AxionCamb was ued to obtain the axion constraints in http://arxiv.org/abs/1410.2896, the forecasts of, http://arxiv.org/abs/1708.05681, and the axion + isocurvature constraints of http://arxiv.org/abs/1607.08208.
The "base" version of CAMB that AxionCamb is built off is Nov13. 
The physics of AxionCamb is described in detail in http://arxiv.org/abs/1410.2896 (HGMF). When using this code, you should cite HGMF, and the original version of CAMB.
For parameter estimation, AxionCamb is compatible with cosmomc, with minimal modifications. However, you should be careful with sampling and degeneracies, as described in HGMF.
AxionCamb is also compatible as a module in cosmosis (https://bitbucket.org/joezuntz/cosmosis/wiki/Home)

Posted here is axionCAMB 2.0, incorporating changes made by Rayne Liu, Wayne Hu, and Daniel Grin. 



----------------------------------------------------------------------
#### <a name="basics"></a>Compiling AxionCAMB

To compile for use in cosmosis (as in arxiv/1706.05681 or arXiv/1607.08208), copy the file MakefileCosmosis to Makefile and run Make. To plot and explore power spectra, use the Makefile MakefileCAMB. Do not erase the file Makefile_main, but modify it as needed if you include new physics. Then run make.

#### <a name="basics"></a>Basic Usage

Refer to CAMB documentation for all non-axion related things. 

AxionCamb adds one new fluid component: the axion fluid. This is included in equations_ppf.f90.
The axion fluid is consistently incoroporated into the evolution of all background and perturbation variables, including, for example, recombination in recfast_axion.f90.
There are two new parameters: the axion mass, ma (measured in eV), and the axion energy density, omaxh2. These are entered in params.ini as usual.
Sensible values for the axion mass range between 1e-33 and 1e-18, though the code works for values outside of this range.
The axion fluid can function as either a dark matter component or a dark energy component depending on the mass.

Isocurvature perturbations in the axion fluid can be turned on using initial_condition=6. 
For isocurvature as described in http://arxiv.org/abs/1708.05681, set axion_isocurvature=T. Hinf is log10(H_inf), where H_inf is the inflationary Hubble parameter in GeV. In arxiv/1708.05681, this is the parameter varied and the tensor/scalar ratio and isocurvature amplitude alpha_ax are derived parameters.  When axion_isocurvature=T, the code outputs the sum of the adiabatic and isocurvature power spectrum in the CMB output files.
AxionCamb can be used in conjunction with the additional neutrino parameters of neutrino mass and number of species.
AxionCamb cannot (currently) be used in conjunction with additional dark energy parameters, such as the equation of state.

----------------------------------------------------------------------
#### <a name="physics"></a>Physics

At early times, the fluid equations solved are equivalent to the first order perturbed Klein-Gordon equation in synchronous gauge.
At late times, the WKB approximation matches the axion field to an effectieve fluid with equation of state zero, and a scale-dependent sound-speed.
The energy density is found using a shooting method from an initial guess for the axion field value. 
The background evolution is computed in axion_background.f90, which creates an array for the axion equation of state. This array is carried around everywhere the energy densities are needed.

----------------------------------------------------------------------
#### <a name="warnings"></a>Known Issues and Warnings

You should be careful if you use any non-linear options (setting do_nonlinear \= 0). AxionCamb uses halofit_ppf.f90 with various modifications. The treatment of axions is not expected to be quantitatively correct, and various approximations are made (these will documented at a later date). The default version of halofit is set to the original "Smith" version, as this is expected to be most stable.
You should also be careful of galaxy bias if you use galaxy survey data, as discussed in HGMF.

#### <a name="warnings"></a>Bug Fix

4 bugs were identified in axionCAMB and systematically investigated by Rayne Liu and Wayne Hu (with help from axionCAMB authors):

1) In the background module (axion_background.f90), a number of factors of little h were identified in the first-order casting of the KG equation, the Friedmann equation, and the adiabatic sound speed.
2) Neutrino variables were normalized in a way that led to incorrect (at the several perecent level) neutrino contributions to the energy density when mnu!=0 -- tests of the code for non-zero neutrino mass helped uncover this issue.
3) One coefficient in the 8th order RK solved used in the axion_background.f90 module was off in the 3rd decimal place.
4) Scattered use of the COBE CMB temperature in the code (2.7255K) instead of the input value.

These bugs have now been fixed. 





