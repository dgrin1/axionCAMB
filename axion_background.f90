!!!!!!!!!!!!!!!! axion_background !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module computes the evolution of the background in the presence of an axion field
! We use exact equations up until m=3H, and then treat as a w=0 fluid
! For more details, see Hlozek et al, 2014, arXiv:1410.2896
! 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Input cosmological parameters, and receive evolution of equation of state, density, 
!adiabatic sound speed
! Also receive initial field value, oscillation scale factor, equality scale factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This code The Klein-Gordon equation for a scalar field in an expanding universe (including
!photons, massive and massless neutrinos, cold dark matter, and baryons)
! \ddot{\phi}+2H \phi+m^2 a \phi=0 is transformed using the variable definitions
!where dots denote derivatives with respect to conformal time and H
!is the conformal Hubble parameter and a is the cosmological scale factor
!phi=\sqrt(3/(4\pi G)) v_{1}^{twiddle}
!phi_dot=u_2, u_2=H_{0} u_{2}^twiddle, dimensionless conformal time t^twiddle=H_{0} t
!H_{0} is the usual Hubble parameter, and finally u_{2}^{twiddle}=\sqrt{3/(4\pi G)} v_{2}^{twiddle}
!A dimensionless Hubble parameter littleh is defined as H/(100 km /s/ Mpc)
!as is a dimensionless axion mass maxion_twiddle=m/(H_0), where h
!is the dimensionless Hubble parameter
!This transforms the second order ODE into a pair of first order ODES
!In terms of these variables, there is a nice expression for the adiabatic
! sound speed c_ad^{2}=Pdot/\rho_dot=1+(2/3) maxion_twiddle^{2} v_{1}^{twiddle}a^{2}/(v_{2}^{twiddle}*littleh)





module axion_background

contains

subroutine w_evolve(Params, badflag)
!Params is the overall parameter structure to be passed back to camb, including the estimated a_osc 
!value, table for interpolation of the adiabatic sound speed, axion EOS and axion energy density
!flag for failed histories where H becomes negative (axion drives collapse). Flag is raised in derivs routine or
!in w_evolve if a history has a turnaround within 10 log(a) bin steps of the last time step
!Give subroutine access to model parameters from camb, constants from constants.f90, massive neutrino routines

use ModelParams
use constants
use Precision
use MassiveNu
implicit none
type(CAMBparams) :: Params 
integer badflag ! flag for failed histories where H becomes negative (axion drives collapse). Flag is raised in derivs routine.
integer i 
integer k,j !general indices



!Internal: Neutrino-related constants from CAMB input
!contribution to H^2/(100 km /s/ Mpc)^2 of massive and massless neutrinos, respectively
real(dl) lhsqcont_massive(Params%Nu_mass_eigenstates),lhsqcont_massless
!zeta(3), conversion factor to get neutrino masses in eV for check comparison with rest of CAMB
real(dl) zeta3, conv
!Correction factor from massless to massive neutrino energy densities, constant to go from massive neutrino
!mass fractions to massive neutrino massive values
real(dl)  rhonu,nu_constant
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! non-trivial aspects to sharing neutrino data structures in new subroutines of CAMB
! so we recompute somethings twice, but these are single numbers
! should not be a significant slow down



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal:: Cosmological Parameters and useful algebraic combinations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!omega h^2 of various species
!Matter, Cosmological constant, massive neutrino density
real(dl) omegah2_m, omegah2_lambda,omnuh2
!baryons + cold dark matter, Hubble parameter in units of eV, omega_curvature
real(dl) omegah2_regm,H_ev,rhocrit,omk
!baryons, dark matter separately, axions, axion mass in units of Hubble, dimensionless Hubble Parameter
real(dl) omegah2_b,omegah2_dm,omegah2_ax,maxion_twiddle,hnot
!Dimensionless Hubble^2
real (dl) hsq
!Desired axion fraction omega_a/(omega_a+omega_b+omegac)
real(dl) fax
! scale factor at equality estimated under assumption axions redshift purely as CDM
real(dl) regzeq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal: Evolution control parameters (other than those defined in modules.f90)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!dfac sets m=dfac*H condition required for code to switch from evolving KG+Einstein equations recast in terms of fluid variables
!to evolving WKB equations
!Details in Hlozek et al 2014. arXiv:1410.2896
real(dl) dfac
!how much past aosc to go in the final a integration at the end
real(dl),parameter::eps=0.001
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Internal: Evolution variable arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(dl) a_arr(ntable) ! scale factor
!scalar field and its derivative (in appropriate units), conformal Hubble parameter /(100 km /s/ Mpc) as a function of scale factor
real(dl) v_vec(2,ntable),littlehfunc(ntable) 
!m/(dfac*littlehfunc*a) diagnostic parameter that is used to identify approximate aosc
real(dl) diagnostic(ntable)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Internal: Fluid quantities of interest, computed after evolution generated
!Scalar field energy density
real(dl) grhoax_table_internal(ntable)
!log_e(m/3H) used for array of this quantity to find a_osc, the time of scalar field oscillation
!!arrays used to find aosc once final scalar field initial condition is determined
real(dl) f_arr(1:ntable)
!Array used to find estimate of scale factor at matter-radiation equality (axions counter as matter)
real(dl) eq_arr(ntable)
!!!!!!!!!!!!!!!!!!!!!!!!

! NB some of these quantities a_arr, and grhoax_table_internal are here so that
! we don't have to keep taking logarithms and exponentials of lists
! ultimately we just want the logs of these for cubic spline interpolation in the rest of CAMB
! but to compute Hubble, field derivatives, we need the quantities themselves and not their logs
! By keeping these quantities we can compute them once inside the scalar field routine (per initial condition)
! and keep the information to use for when it is needed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Internal storage arrays for spline interpolation of scalar field solution
!for fixed scalar initial condition
real(dl) v_buff(1:ntable) ! spline buffer scalar field and its derivative with respect to a
real (dl) abuff(1:ntable) !spline buffer for scale factor interpolation (ie finding 'when
! various quantities like rho_rad/rho_m, m/3H, etc, reach values of interest)
!
real (dl) eq_arr_buff(ntable)
!nstop is a control index that excludes really bad scalar field evolution histories and
!prevents them from screwing up the spline
integer  nstop
!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Quantities of interest that are generated for each scalar field evolution 
! and sometimes read into arrays
! energy fraction output (omega_ax/(omega_ax+omega_c+omega_b))
real(dl) fout
!Estimated scalar factor when m=dfac*H, scalar field v(1,i) when this value is reached,
!Estimated scalar field derivative at same point i time
real(dl) aosc,phiosc,phidosc
!log(aosc)
real(dl) laosc
!output value of scalar field value needed to get desired axion density
!obtained via cubic splines, final fout/fax for the corresponding scalar field history
real(dl) vtw,final_fout_check
!critical densities, joining value of rho when axion starts to act like CDM, better log initial value so derivatives can be taken t offset, v2 value at aosc joining point
!used to set full density evolution of final chosen scalar field history
real(dl) drefp
!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Internal storage arrays to sweep through quantities of interest for different 
!scalar field evolutions corresponding to different initial conditions
!
! omega_ax/(omega_ax+omega_c+omega_b) for different scalar field initial conditions
real (dl) fout_arr(nphi)
!array of aosc values for different initial conditions
! for the scalar field
real(dl) aosc_arr(nphi)
!log of aosc values for different scalar field histories
real(dl) laosc_arr(nphi)
!array of fout_arr(nphi)/fax (actual axion mass fraction vs desired axion mass fraction)
!logarithm later taken for root finding
real(dl) fout_check_arr(nphi)
! stepsize in log of initial scalar field value
real(dl) dlogvtwiddle_init
!array of initial v1 values
real(dl) vtwiddle_init_arr(nphi)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!WORK NEEDED ON COMMENTS
!!!!!! Simple analytic estimates for required scalar field initial condition
!!! used to set a broad range of scalar initial conditions 
!Scalar field evolution is generated for each initial condition and
!cubic spline interpolation is then used to find the one which generates
!desired axion energy fraction today
!initial value max and min for phi_init, and final output 
real(dl) vtwiddle_initmax,vtwiddle_initmin,vtwiddle_init
!list of initial values for phi_init
real(dl) vtwiddle_initlist(3)
!!!!!!!!
!WORK NEEDED ON COMMENTS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Internal storage arrays for spline interpolation of quantities of interest
! for different scalar initial conditions
!array of axion energy fractions, buffer array for spline fitting in space of initial phi values
real (dl) vtwiddle_init_arr_buff(nphi)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Spline derivatives at endpoints (used for all splines)
real(dl) d1,d2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!
!Timing variables
!real :: clock_start, clock_stop   
!!!!!!!!!!

!WORK NEEDED ON OMMENTS
! numbers for initial and final value of a, and stepsize
! initial needed to ensure axion sub-dominance to everything else at a_init
real(dl) a_init,a_m, a_lambda, a_rel,as_scalar,as_rad,as_matt,a_final,dloga,log_a_final,log_a_init
real(dl) cmat(1:16,1:16),kvec(1:16,1:2),kfinal(1:2),svec(1:16),avec(1:16)
!END WORK NEEDED ON COMMENTS


!Time Code
!clock_start = 0.0
!call cpu_time(clock_start)
!Set value of neutrino constants as elsewhre in camb
zeta3=1.2020569031595942853997d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DM cosmological parameters etc. in units for integrator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
hnot = Params%H0/100.d0
hsq=hnot**2.0d0
omegah2_dm=Params%omegac*(hnot**2.0d0)
omegah2_b=Params%omegab*(hnot**2.0d0) 
omnuh2=Params%omegan*(hnot**2.0d0)   

omegah2_lambda=Params%omegav*(hnot**2.0d0)  
omegah2_ax=Params%omegaax*(hnot**2.0d0)   
maxion_twiddle= Params%ma
!Hubble parameter in eV
H_eV=1.d14*6.5821d-25*hnot/(MPC_in_sec*c)
!convert axion mass units from eV to H
maxion_twiddle = maxion_twiddle/H_eV 
!Total densities to juggle and use later
omegah2_regm=omegah2_dm+omegah2_b
omegah2_m=omegah2_regm+omegah2_ax
!Axion Mass Fraction
fax=omegah2_ax/omegah2_m

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Initialize massive and massless neutrino variables
call Nu_init
grhom = 3.0d0*(hsq*1.d10)/(c**2.0d0) 


!!!!!4/8 RL TCMB changes again (COBE->input TCMB)
! grhog = ((kappa/(c**2.0d0)*4.0d0*sigma_boltz)/(c**3.0d0))*(COBE_CMBTemp**4.0d0)*(Mpc**2.0d0) 

grhog = ((kappa/(c**2.0d0)*4.0d0*sigma_boltz)/(c**3.0d0))*(Params%TCMB**4.0d0)*(Mpc**2.0d0)



grhor = (7.0d0/8.0d0)*((4.0d0/11.0d0)**(4.0d0/3.0d0))*grhog 
!calculate critical density
rhocrit=(8.0d0*const_pi*G*1.d3/(3.0d0*((1.d7/(MPC_in_sec*c*1.d2))**(2.0d0))))**(-1.0d0)

!4/8 RL correct COBE-> regular CMB temperature

! Params%omegah2_rad=((COBE_CMBTemp**4.0d0)/(rhocrit))/(c**2.0d0)

Params%omegah2_rad=((Params%TCMB**4.0d0)/(rhocrit))/(c**2.0d0)


Params%omegah2_rad=Params%omegah2_rad*a_rad*1.d1/(1.d4)
!calculate omega rad using starndard formula
!Contribution of photons and massless neutrinos to H/(100 km /s/Mpc)

!!!DG 4/8/2023 -- again, deal with neutrino error
!lhsqcont_massless=(Params%Num_Nu_massless*grhor*(c**2.0d0)/((1.d5**2.0d0)))/3.0d0
lhsqcont_massless=(Params%nu_massless_degeneracy*grhor*(c**2.0d0)/((1.d5**2.0d0)))/3.0d0
!!!!!



Params%omegah2_rad=Params%omegah2_rad+lhsqcont_massless
!print*, 'hi renee', Params%omegah2_rad
!const from modules.f90
nu_constant=(7.0d0/120.0d0)*(const_pi**4.0d0)/(zeta3*1.5d0)
nu_constant=nu_constant*omnuh2*(grhom/grhor)/hsq

do k=1,Params%Nu_mass_eigenstates,1
!Compute neutrino masses corresponding to input mass fractions and degeneracies
	Nu_masses(k)=nu_constant*Params%Nu_mass_fractions(k)/(Params%Nu_mass_degeneracies(k))
!Compute contribution of massive neutrinos to Hubble parameter/100 km/s/Mpc
 	lhsqcont_massive(k)=Params%Nu_mass_degeneracies(k)*(grhor*(c**2.0d0)/((1.d5**2.0d0)))/3.0d0
!Print some useful checks for neutrinos
 !	call Nu_rho(Nu_Masses(k),rhonu)
!	print*,Nu_masses(k),omegah2_rad+lhsqcont_massive(1)*rhonu
 !	conv = k_B*(8.0d0*grhor/grhog/7.0d0)**0.25d0*COBE_CMBTemp/elecV * &
 !   &(Params%Nu_mass_degeneracies(k)/dble(Params%Nu_mass_numbers(k)))**0.25d0!
!	call Nu_rho(1.0d0,rhonu)
!	print*,Nu_masses(k)*conv
enddo
!print*,lhsqcont_massive(1)/Params%Nu_mass_degeneracies(1)

!! Params%omegah2_rad=Params%omegah2_rad+

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Switch between  evolution regimes
!definition of m=nH pmatch point, use 3 for convention
dfac=3.0d0
!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Constants for integrator from
!This is an implementation of the 8th order formulas on page 75 of 
!Classical Fifth-, Sixth-, Seventh-, and Eighth-Order Runge-Kutta Formulas with Stepsize Control
!by E Fehlberg
!http://hdl.handle.net/2060/19680027281
!(NASA Huntsville, 1968)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Constant matrix
avec=0.0d0
kvec=0.0d0
avec(1)=0.4436894037649818d0
avec(2)=0.6655341056474727d0
avec(3)=0.9983011584712091d0
avec(4)=0.31550d0
avec(5)=0.5054410094816906d0
avec(6)=0.1714285714285714d0
avec(7)=0.8285714285714285d0
avec(8)=0.6654396612101156d0
avec(9)=0.2487831796806265d0
avec(10)=0.1090d0
avec(11)=0.8910d0
avec(12)=0.3995d0
avec(13)=0.6005d0
avec(14)=1.0d0
avec(15)=0.0d0
avec(16)=1.0d0
cmat=0.0d0
cmat(1,1)=avec(1)
cmat(2,1)=0.1663835264118681d0
cmat(2,2)=0.49915057d0
cmat(3,1)=0.24957528d0
cmat(3,2)=0.0d0
cmat(3,3)=0.74872586d0
cmat(4,1)=0.20661891d0
cmat(4,2)=0.0d0
cmat(4,3)=0.17707880d0
cmat(4,4)=-0.68197715d-1
cmat(5,1)=0.10927823d0
cmat(5,2)=0.0d0
cmat(5,3)=0.0d0
cmat(5,4)=0.40215962d-2
cmat(5,5)=0.39214118d0
cmat(6,1)=0.98899281d-1
cmat(6,2)=0.0d0
cmat(6,3)=0.0d0
cmat(6,4)=0.35138370d-1
cmat(6,5)=0.12476099d0
cmat(6,6)=-0.55745546d-1
cmat(7,1)=-0.36806865d0
cmat(7,2)=0.0d0
cmat(7,3)=0.0d0
cmat(7,4)=0.0d0
cmat(7,5)=-0.22273897d1
cmat(7,6)=0.13742908d1
cmat(7,7)=0.20497390d1
cmat(8,1)=0.45467962d-1
cmat(8,2)=0.0d0
cmat(8,3)=0.0d0
cmat(8,4)=0.0d0
cmat(8,5)=0.0d0
cmat(8,6)=0.32542131d0
cmat(8,7)=0.28476660d0
cmat(8,8)=0.97837801d-2
cmat(9,1)=0.60842071d-1
cmat(9,2)=0.0d0
cmat(9,3)=0.0d0
cmat(9,4)=0.0d0
cmat(9,5)=0.0d0
cmat(9,6)=-0.21184565d-1
cmat(9,7)=0.19596557d0
cmat(9,8)=-0.42742640d-2
cmat(9,9)=0.17434365d-1
cmat(10,1)=0.54059783d-1
cmat(10,2)=0.0d0
cmat(10,3)=0.0d0
cmat(10,4)=0.0d0
cmat(10,5)=0.0d0
cmat(10,6)=0.0d0
cmat(10,7)=.11029325d0
cmat(10,8)=-.12565008d-2
cmat(10,9)=0.36790043d-2
cmat(10,10)=-.57780542d-1
cmat(11,1)=.12732477d0
cmat(11,2)=0.0d0
cmat(11,3)=0.0d0
cmat(11,4)=0.0d0
cmat(11,5)=0.0d0
cmat(11,6)=0.0d0
cmat(11,7)=0.0d0
cmat(11,8)=0.11448805
cmat(11,9)=0.28773020
cmat(11,10)=0.50945379d0
cmat(11,11)=-0.14799682d0
cmat(12,1)=-0.36526793d-2
cmat(12,2)=0.0d0
cmat(12,3)=0.0d0
cmat(12,4)=0.0d0
cmat(12,5)=0.0d0
cmat(12,6)=0.81629896d-1
cmat(12,7)=-0.38607735d0
cmat(12,8)=0.30862242d-1
cmat(12,9)=-0.58077254d-1
cmat(12,10)=0.33598659d0
cmat(12,11)=0.41066880d0
cmat(12,12)=-0.11840245d-1
cmat(13,1)=-0.12375357d1
cmat(13,2)=0.0d0
cmat(13,3)=0.0d0
cmat(13,4)=0.0d0
cmat(13,5)=0.0d0
cmat(13,6)=-0.24430768d2
cmat(13,7)=0.54779568d0
cmat(13,8)=-0.44413863d1
cmat(13,9)=0.10013104d2
cmat(13,10)=-0.14995773d2
cmat(13,11)=0.58946948d1
cmat(13,12)=0.17380377d1
cmat(13,13)=0.27512330d2
cmat(14,1)=-0.35260859d0
cmat(14,2)=0.0d0
cmat(14,3)=0.0d0
cmat(14,4)=0.0d0
cmat(14,5)=0.0d0
cmat(14,6)=-0.18396103d0
cmat(14,7)=-0.65570189d0
cmat(14,8)=-.39086144d0
cmat(14,9)=0.26794646d0
cmat(14,10)=-0.10383022d1
cmat(14,11)=0.16672327d1
cmat(14,12)=0.49551925d0
cmat(14,13)=.11394001d1
cmat(14,14)=0.51336696d-1
cmat(15,1)=0.10464847d-2
cmat(15,2)=0.0d0
cmat(15,3)=0.0d0
cmat(15,4)=0.0d0
cmat(15,5)=0.0d0
cmat(15,6)=0.0d0
cmat(15,7)=0.0d0
cmat(15,8)=0.0d0
cmat(15,9)=-0.67163886d-2
cmat(15,10)=0.81828762d-2
cmat(15,11)=-0.42640342d-2
cmat(15,12)=0.280090294741d-3
cmat(15,13)=-0.87835333d-2
cmat(15,14)=0.10254505d-1
cmat(15,15)=0.0d0
cmat(16,1)=-0.13536550d1
cmat(16,2)=0.0d0
cmat(16,3)=0.0d0
cmat(16,4)=0.0d0
cmat(16,5)=0.0d0
cmat(16,6)=-0.18396103d0
cmat(16,7)=-0.65570189d0
cmat(16,8)=-0.39086144d0
cmat(16,9)=0.27466285d0
cmat(16,10)=-0.10464851d1
cmat(16,11)=0.16714967d1
cmat(16,12)=0.49523916d0
cmat(16,13)=0.11481836d1
cmat(16,14)=0.41082191d-1
cmat(16,15)=0.0d0
cmat(16,16)=1.0d0

svec(1)=0.32256083d-1
svec(2)=0.0d0
svec(3)=0.0d0
svec(4)=0.0d0
svec(5)=0.0d0
svec(6)=0.0d0
svec(7)=0.0d0
svec(8)=0.0d0
svec(9)=0.25983725d0
svec(10)=0.92847805d-1
svec(11)=.16452330d0
svec(12)=0.176659510d0
svec(13)=0.23920102d0
svec(14)=0.39484274d-2

!4/8 DG typo in RK78 coefficient (see jupyter notebook for reference)
!svec(15)=0.3082649547580d-1
svec(15)=0.3072649547580d-1




!END NEEDS WORK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! begin initialization procedure, shoot for best phi_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!NEEDS WORK IN COMMNENTS
!Initial aosc guess to tell computer to try to find aosc, also a value that flags when axion doesnt start oscillating by today
aosc=15.0d0
!compute curvature parameter for this set of omegas
omk=1.0d0-(omegah2_m+Params%omegah2_rad+omegah2_lambda+omnuh2)/hsq

!!!!!!!!!
! Determine a range of epochs to be well before in starting axion field evolution
!scale factor at time of axions becoming irrelevant in competition with different species
!and other important transitionary epochs
as_matt=(omegah2_regm/(maxion_twiddle**2.0d0))**(1.0d0/3.0d0)
as_rad=(Params%omegah2_rad/(maxion_twiddle**2.0d0))**(1.0d0/4.0d0)
!scale factor at equality of other matter with radiation
a_m=(Params%omegah2_rad/(omegah2_regm))
!This tims 10^-7 will be another maximum initial scale factor (see below)
a_rel=10.0d0
!Subdominance of dark energy to ordinary radiation
a_lambda=(Params%omegah2_rad/omegah2_lambda)**(0.25d0)
!aosc (start well before coherent oscillation of the axion field)
as_scalar=(omegah2_ax/(maxion_twiddle**2.0d0))**(1.0d0/3.0d0)
!find safe initial a (take analytic estimates for when axions are negligible to everything else
! and look at scalar factors 7 orders of magnitude smaller as initial scale factors!!!
 a_init=min(a_rel,a_lambda,a_m,as_matt,as_rad,as_scalar)*1.d-7
!maximum a
a_final=1.0d0
   !log of various a values, dlog a, step back one
   log_a_init=dlog(a_init)
   log_a_final=dlog(a_final)   
  !set log a bin size
   dloga=(log_a_final-log_a_init)/(dble(ntable-1)) ! RH: DG's original
	a_arr(1)=dexp(log_a_init)
!Feed to output parameter table
   !initial loga and a value
   Params%loga_table(1)=log_a_init
   !Generate arrays of log scale factor for output and internal scale factor
   forall(i=2:ntable)
      Params%loga_table(i)=log_a_init+dloga*dble(i-1)
      a_arr(i)=dexp(Params%loga_table(i))
   end forall
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


!!!!!!!
!assuming axion acts like cosmological constant then matter, 
!analytically obtain
!initial v_1^twiddle value to get right density today,
! under different assumptions
! Expressions will agree with arXiv:1410.2896 in appropriate asymptotic limits
!During Axion domination
vtwiddle_initlist(1)=fax
!Oscillation starts during matter domination
vtwiddle_initlist(2)=fax*omegah2_m/(dsqrt(maxion_twiddle)*((Params%omegah2_rad**0.750d0)))
!Oscillation starts during radiation domination
vtwiddle_initlist(3)=fax*omegah2_m/(dsqrt(maxion_twiddle)*(Params%omegah2_rad))
!Concatenate these initial v_1^twiddle
vtwiddle_initlist=dsqrt(vtwiddle_initlist)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Set range of initial scalar field initial conditions to try
!bracket range of initial vtwiddle 1 values (basically phi 1)
!Spanning many orders of magnitude beyond the simples guesses above to bracket a wide range of initial conditions
!and make sure to choose the correct one
vtwiddle_initmin=min(vtwiddle_initlist(1),vtwiddle_initlist(2),vtwiddle_initlist(3))/1.0d2
vtwiddle_initmax=max(vtwiddle_initlist(1),vtwiddle_initlist(2),vtwiddle_initlist(3))*100.0d1
!set size of log step in v1 initial
dlogvtwiddle_init=(dlog(vtwiddle_initmax)-dlog(vtwiddle_initmin))/(dble(nphi)-1.0d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Generate array of initial phi values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
forall(j=1:nphi)
   vtwiddle_init_arr(j)=dexp(dble(j-1)*dlogvtwiddle_init+dlog(vtwiddle_initmin))
end forall
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!LOOP OVER SCALAR FIELD INITIAL CONDITIONS
do j=1,nphi,1
!Control bad value of aosc -- if a history never gets to m=dfac*H the value of this parameter alerts the code
! that this is a history for which coherent oscillation never begins

   aosc=15.0d0
!!!!initial phi value
   vtwiddle_init=vtwiddle_init_arr(j)
  
   !start axion at rest at top of hill, with phi guessed at current guess
   v_vec(1,1)=vtwiddle_init
   v_vec(2,1)=0.0d0
   !calculate dimensionless hubble including standardmatter and axion potential and kinetic energy
   call lh(omegah2_regm,Params%omegah2_rad,omegah2_lambda,omk,hsq,maxion_twiddle,a_arr(1),v_vec(1:2,1),littlehfunc(1),badflag,&
        lhsqcont_massless,lhsqcont_massive,Params%Nu_mass_eigenstates,Nu_masses)
   
   !calculate initial value of dphi/da and dphidot/da
    !Compute first step parameters for the pair of first-order ODEs being solved
   kvec=0.0d0
   kfinal=0.0d0
   call next_step(a_arr(1),v_vec(1:2,1),kvec(1:16,1:2),kfinal(1:2),avec(1:16),&
   	&omegah2_regm,Params%omegah2_rad,&
    &omegah2_lambda,omk,hsq,&
	&maxion_twiddle,badflag,dloga,16,cmat(1:16,1:16),lhsqcont_massless,lhsqcont_massive,&
	&Params%Nu_mass_eigenstates,Nu_masses)
   
   diagnostic(1)=dfac*littlehfunc(1)/a_arr(1)
   
!!!!

!!!!!!!!!!!!!!!!!!!!!!
!Integration performed using eight order Runge-Kutta method
!Derivatives at 8 points in advance of the one in question (with trial function values)
!Added to current value using optimal quadrature coefficients under certain assumptions
!This is an implementation of the 8th order formulas on page 75 of 
!Classical Fifth-, Sixth-, Seventh-, and Eighth-Order Runge-Kutta Formulas with Stepsize Control
!by E Fehlberg
!http://hdl.handle.net/2060/19680027281
!(NASA Huntsville, 1968)
!It turns out a reasonably accurate integrator is required to accurately obtain the adiabatic sound speed at earlier times
!Using super-horizon analytic solutions, we found that this integrator + 5000 grid points
!was necessary to avoid exciting a non-physically large low-l ISW effect
! DM please add comments on which condition led us to this realiztion???


   do i=2,ntable,1
!!!!integrate ODE using 16 pt (8th order Runge-Kutta) rule
      !increment fluid (homogeneous values) using precomputed steps
      v_vec(:,i)=v_vec(:,i-1)+(svec(1)*kvec(1,:)+svec(9)*kvec(9,:)+svec(10)*kvec(10,:)&
           &+svec(11)*kvec(11,:)+svec(12)*kvec(12,:)+svec(13)*kvec(13,:)&
           &+svec(14)*kvec(14,:)+svec(15)*kvec(15,:))
!!!!
      !calculate hubble for next step
      call lh(omegah2_regm,Params%omegah2_rad,omegah2_lambda,omk,hsq,maxion_twiddle,a_arr(i),v_vec(1:2,i),littlehfunc(i),badflag,&
           &lhsqcont_massless,lhsqcont_massive,Params%Nu_mass_eigenstates,Nu_masses)
      
      kvec=0.0d0
      kfinal=0.0d0
      
      !Compute next steps in scalar field and its derivative
      call next_step(a_arr(i),v_vec(1:2,i),kvec(1:16,1:2),kfinal(1:2),&
      &avec(1:16),omegah2_regm,Params%omegah2_rad,omegah2_lambda,omk,hsq,&
           &maxion_twiddle,badflag,dloga,16,cmat(1:16,1:16),&
           &lhsqcont_massless,lhsqcont_massive,Params%Nu_mass_eigenstates,Nu_masses)
      
      
      !compare m to nH and identify a reasonable guess for when a given history crosses the m=3H 
      ! condition (and thus when the code will switch from pure scalar evolution to coherent oscillation)
      diagnostic(i)=dfac*littlehfunc(i)/a_arr(i)
      ! find first guess for aosc
      if (i .gt. 1) then
         if (aosc .eq. 15.0d0) then
            if (maxion_twiddle .lt.diagnostic(i-1)) then
            if (maxion_twiddle .ge.diagnostic(i)) then
               aosc=a_arr(i)
            endif
         endif
      endif
   endif
enddo
!!!!now there's an array of dlog(diagnostic), use a simple spline to find aosc

f_arr(1:ntable)=maxion_twiddle/(diagnostic(1:ntable))
f_arr=dlog(f_arr)



if (aosc .ne. 15.0d0) then

d1=(Params%loga_table(2)-Params%loga_table(1))/(f_arr(2)-f_arr(1))
d2=(Params%loga_table(ntable)-Params%loga_table(ntable-1))/(f_arr(ntable)-f_arr(ntable-1))

call spline(f_arr(1:ntable),(Params%loga_table(1:ntable)),ntable,d1,d2,abuff(1:ntable))
call spline_out(f_arr(1:ntable),Params%loga_table(1:ntable),&
	&abuff(1:ntable)&
	&,ntable,0.0d0,laosc)

!use a simple spline to calculate phi at time oscillation begins (important to get relic density later, which will be off by a factor of order unity depending on how aosc is defined)
!Finds root of equation log(m/3H)=0 uisng a cubic spline
d1=(v_vec(1,2)-v_vec(1,1))/(Params%loga_table(2)-Params%loga_table(1))
d2=(v_vec(1,ntable)-v_vec(1,ntable-1))/(Params%loga_table(ntable)-Params%loga_table(ntable-1))
do i=1,ntable-10,1
!!       print*,i,isnan(f_arr(i))
!
do k=0,10,1
call lh(omegah2_regm,Params%omegah2_rad,omegah2_lambda,omk,hsq,maxion_twiddle,a_arr(i+k),v_vec(1:2,i+k),littlehfunc(i+k),badflag,&
	&lhsqcont_massless,lhsqcont_massive,Params%Nu_mass_eigenstates,Nu_masses)
        if (((isnan(f_arr(i+k))).or.(isnan(v_vec(1,i+k))) ).or.(isnan(v_vec(2,i+k))))then
                print*,i,k,littlehfunc(i+k),v_vec(1,i+k),v_vec(2,i+k),a_arr(i+k),dexp(laosc)
        endif
enddo
enddo

!Find scalar field value at moment of m=nH (fiducial onset of coherent oscillation)
call spline(Params%loga_table(1:ntable),v_vec(1,1:ntable),ntable,d1,d2,&
      v_buff(1:ntable))
call spline_out(Params%loga_table(1:ntable),v_vec(1,1:ntable),&
	&v_buff(1:ntable)&
	&,ntable,laosc,phiosc)  

!Find scalar field derivative with scale factor at moment of m=nH (fiducal onset of coherent
!oscillation)  
!use a cubic spline to calculate phidot at time rolling begins, we can see perfect equipartition between kinetic and potential scalar field energy is not yet achieved at m=3H, but this is good enough for government work
d1=(v_vec(2,2)-v_vec(2,1))/(Params%loga_table(2)-Params%loga_table(1))
d2=(v_vec(2,ntable)-v_vec(2,ntable-1))/(Params%loga_table(ntable)-Params%loga_table(ntable-1))      
call spline(Params%loga_table(1:ntable),v_vec(2,1:ntable),ntable,d1,d2,v_buff(1:ntable))   
call spline_out(Params%loga_table(1:ntable),v_vec(2,1:ntable),&
     &v_buff(1:ntable)&
     &,ntable,laosc,phidosc) 


aosc=dexp(laosc)  
!calculate axion energy density fraction (compared to total matter+axion)
!parcel these numbers of interest into an array which sweeps through different v1 values  
!Beyond a>aosc, a simple a^-3 energy density scaling is applied                
fout=((maxion_twiddle*phiosc)**2.0d0+(phidosc/aosc)**2.0d0)*(aosc**3.0d0)
fout=fout/(omegah2_regm+((maxion_twiddle*phiosc)**2.0d0+(phidosc/aosc)**2.0d0)*(aosc**3.0d0))
!print*,vtwiddle_init,aosc,fout,fax
fout_arr(j)=fout
fout_check_arr(j)=fout/fax
aosc_arr(j)=(aosc)
f_arr=dexp(f_arr)
else
!(case of no oscillation, just use result from full scalar field evolution)
f_arr=dexp(f_arr)
fout=(v_vec(2,ntable)/a_arr(ntable))**2.0d0+(maxion_twiddle*v_vec(1,ntable))**2.0d0
fout=fout/(fout+omegah2_regm)
fout_arr(j)=fout
fout_check_arr(j)=fout/fax
aosc_arr(j)=aosc
endif



enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End of sweep through different evolutions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!DM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Spline results of shooting method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

vtwiddle_init_arr=dlog(vtwiddle_init_arr)
fout_check_arr=dlog(fout_check_arr)
fout_arr=dlog(fout_arr)


!Exclude histories that crashed
nstop=0
do j=1,nphi,1
   if (j .lt. nphi) then
      if ((dexp(fout_check_arr(j+1)).eq.dexp(fout_check_arr(j)))) then
         if (nstop.eq.0) then
            nstop=j
         endif
      endif
   endif
!   print*,j,fout_check_arr(j)
enddo
if (nstop .eq. 0) then
   nstop=nphi
endif
!
!!!cubic spline to find best initial phi value to get axion energy fraction that is desire
d1=(vtwiddle_init_arr(2)-vtwiddle_init_arr(1))/(fout_check_arr(2)-fout_check_arr(1))
d2=(vtwiddle_init_arr(nstop)-vtwiddle_init_arr(nstop-1))
d2=d2/(fout_check_arr(nstop)-fout_check_arr(nstop-1))
call spline(fout_check_arr(1:nstop),vtwiddle_init_arr(1:nstop),nstop,d1,d2,vtwiddle_init_arr_buff(1:nstop))
call spline_out(fout_check_arr(1:nstop),vtwiddle_init_arr(1:nstop),&
	&vtwiddle_init_arr_buff(1:nstop)&
	&,nstop,0.0d0,vtw) 



!!!use a cubic spline to check that correct axion energy fraction is indeed achieved here
d1=(vtwiddle_init_arr(2)-vtwiddle_init_arr(1))/(fout_arr(2)-fout_arr(1))
d1=1.d0/d1
d2=(vtwiddle_init_arr(nstop)-vtwiddle_init_arr(nstop-1))
d2=d2/(fout_arr(nstop)-fout_arr(nstop-1))       
d2=1.0d0/d2
!!!             
call spline(vtwiddle_init_arr(1:nstop),fout_arr(1:nstop),nstop,d1,d2,vtwiddle_init_arr_buff(1:nstop))   
call spline_out(vtwiddle_init_arr(1:nstop),fout_arr(1:nstop),&
     &vtwiddle_init_arr_buff(1:nstop)&
     &,nstop,vtw,final_fout_check) 
!print*,dexp(final_fout_check)*omegah2_regm/(1.0d0-dexp(final_fout_check)),fax*omegah2_regm/(1.0d0-(fax))

!Use cubic spline to get aosc for this best (and chosen history)
laosc_arr=dlog(aosc_arr)
d1=(vtwiddle_init_arr(2)-vtwiddle_init_arr(1))/(aosc_arr(2)-aosc_arr(1))
d1=1.d0/d1
d2=(vtwiddle_init_arr(nstop)-vtwiddle_init_arr(nstop-1))
d2=d2/(laosc_arr(nstop)-laosc_arr(nstop-1))       
d2=1.0d0/d2
call spline(vtwiddle_init_arr(1:nstop),laosc_arr(1:nstop),nstop,d1,d2,vtwiddle_init_arr_buff(1:nstop))   
call spline_out(vtwiddle_init_arr(1:nstop),laosc_arr(1:nstop),&
     &vtwiddle_init_arr_buff(1:nstop)&
     &,nstop,vtw,Params%a_osc)
vtw=dexp(vtw)
Params%a_osc=dexp(Params%a_osc)
if (Params%a_osc .ge. 1.0d0) then
   Params%a_osc=1.0d0
endif



!!!!!!!!!!!!!!!!!!!!!!!!!
! DM log a integration at best fit initial field value for output
!!!!!!!!!!!!!!!!!!!!!!!!


!Do same log scale factor integration at bestfit phi value to get history of fluid variables
!!Set scale factor integration range
a_init=dexp(log_a_init)
if (Params%a_osc .le. 1.0d0) then
a_final=Params%a_osc*(1.0d0+eps)
else
a_final=1.0d0
endif
log_a_final=dlog(a_final)
dloga=(log_a_final-log_a_init)/dble(ntable-1)

!Set intial value
vtwiddle_init=vtw


!initial loga and a value
Params%loga_table(i)=log_a_init
a_arr(1)=dexp(Params%loga_table(i))
!start axion at rest at top of hill, with phi guessed at current guess
v_vec(1,1)=vtwiddle_init
v_vec(2,1)=0.0d0

!Compute initial Hubble
call lh(omegah2_regm,Params%omegah2_rad,omegah2_lambda,omk,hsq,maxion_twiddle,a_arr(1),v_vec(1:2,1),littlehfunc(1),badflag,&
	&lhsqcont_massless,lhsqcont_massive,Params%Nu_mass_eigenstates,Nu_masses)

!Compute first step parameters
kvec=0.0d0
kfinal=0.0d0
call next_step(a_arr(1),v_vec(1:2,1),kvec(1:16,1:2),kfinal(1:2),avec(1:16),omegah2_regm,Params%omegah2_rad,omegah2_lambda,omk,hsq,&
	&maxion_twiddle,badflag,dloga,16,cmat(1:16,1:16),lhsqcont_massless,lhsqcont_massive,Params%Nu_mass_eigenstates,Nu_masses)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DM: do loop for pulling out values at correct initial phi
!time integration, pull out same values at each time step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

badflag=0

forall(i=2:ntable)
   Params%loga_table(i)=log_a_init+dloga*dble(i-1)
   a_arr(i)=dexp(Params%loga_table(i))
end forall

do i=2,ntable,1
   !Take step using previously computed step parameters   
   v_vec(1:2,i)=v_vec(1:2,i-1)+(svec(1)*kvec(1,1:2)+svec(9)*kvec(9,1:2)+svec(10)*kvec(10,1:2)&
	&+svec(11)*kvec(11,1:2)+svec(12)*kvec(12,1:2)+svec(13)*kvec(13,1:2)&
	&+svec(14)*kvec(14,1:2)+svec(15)*kvec(15,1:2))
   
	!Compute Hubble/(km/s/Mpc)
   call lh(omegah2_regm,Params%omegah2_rad,omegah2_lambda,omk,hsq,maxion_twiddle,a_arr(i),v_vec(1:2,i),littlehfunc(i),badflag,&
	&lhsqcont_massless,lhsqcont_massive,Params%Nu_mass_eigenstates,Nu_masses)

kvec=0.0d0
kfinal=0.0d0
!compute next step parameters

call next_step(a_arr(i),v_vec(1:2,i),kvec(1:16,1:2),kfinal(1:2),avec(1:16),omegah2_regm,&
	&Params%omegah2_rad,omegah2_lambda,omk,hsq,&
	&maxion_twiddle,badflag,dloga,16,cmat(1:16,1:16),&
	&lhsqcont_massless,lhsqcont_massive,Params%Nu_mass_eigenstates,Nu_masses)
enddo

!Compute m/3H as a function of scale factor, as well as the scalar field equation of state
!and appropriately normalized energy density
forall(i=1:ntable)
diagnostic(i)=dfac*littlehfunc(i)/a_arr(i)
Params%wax_table(i)=(((v_vec(2,i)/(a_arr(i)))**2.0d0)-((v_vec(1,i)*maxion_twiddle)**2.0d0))
Params%wax_table(i)=Params%wax_table(i)/(((v_vec(2,i)/(a_arr(i)))**2.0d0)+((v_vec(1,i)*maxion_twiddle)**2.0d0))
!!tabulate axion energy density
grhoax_table_internal(i)=((v_vec(2,i)/a_arr(i))**2.0d0+(maxion_twiddle*v_vec(1,i))**2.0d0)
end forall
!This is the axion energy density *h^2/(3H_0^2/8\pi G), where h is the dimensionless Hubble Parameter
!so in camb's definitions grhoa2_axion=grhom*grhox_table_internal(i)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Spline final arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!with complete history in hand, spline to find axion energy density at a=aosc
Params%a_osc=dlog(Params%a_osc)
Params%grhoax_table=dlog(grhoax_table_internal)
   d1=(Params%grhoax_table(2)-Params%grhoax_table(1))/(Params%loga_table(2)-Params%loga_table(1))
d2=(Params%grhoax_table(ntable)-Params%grhoax_table(ntable-1))/(Params%loga_table(ntable)-Params%loga_table(ntable-1))      
call spline(Params%loga_table(1:ntable),Params%grhoax_table(1:ntable),ntable,d1,d2,Params%grhoax_table_buff)
call spline_out(Params%loga_table(1:ntable),Params%grhoax_table(1:ntable),&
	&Params%grhoax_table_buff(1:ntable)&
	&,ntable,Params%a_osc,drefp)


Params%a_osc=dexp(Params%a_osc)
drefp=dexp(drefp)

!Replace putative scalar field solution with a^-3 solution at match point
!for CAMB+CosmoMC (or Cosmosis) adequate time griddings, scalar field EOMs cannot be solved exactly
!Thus we switch to this scaling
do i=1,ntable,1
if (a_arr(i) .gt. Params%a_osc) then
grhoax_table_internal(i)=drefp*((Params%a_osc/a_arr(i))**3.0d0)
endif
enddo
Params%grhoax_table=dlog(grhoax_table_internal)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!calculate axion adiabatic sound speed Pdot/rhodot, asymptoting to constant (physical) value at low a, when small machine numbers start to behave badly
!!Apply analytic expression for Pdot/rdhot so that
!!noisy numerical derivatives need no be used
forall(i=1:ntable)

!little h error in adiabatic sound speed
!4/28 RL + DG
! Params%cs2_table(i)=1.0d0+2.0d0*((maxion_twiddle*a_arr(i))**2.0d0)*v_vec(1,i)/(3.0d0*v_vec(2,i)*littlehfunc(i))
!Corrected expression
Params%cs2_table(i)=1.0d0+2.0d0*((maxion_twiddle*a_arr(i))**2.0d0)*v_vec(1,i)*dsqrt(hsq)/(3.0d0*v_vec(2,i)*littlehfunc(i))
end forall
Params%cs2_table(1)=Params%cs2_table(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!Use WKB averages of cs2_adiabatic and w_axion after a>a_osc
do i=1,ntable,1
if (a_arr(i) .gt. Params%a_osc) then
Params%cs2_table(i)=0.0d0
Params%wax_table(i)=0.0d0
endif
enddo
!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!
!use same spline method used in rest of subroutine to find new aeq value in this cosmological history (will be needed in initial condition code to normalize, time, scale factor, etc, in particular for isocurvature mode
forall(i=1:ntable)
eq_arr(i)=((omegah2_regm/(a_arr(i)**3.0d0)+grhoax_table_internal(i))/(Params%omegah2_rad/(a_arr(i)**4.0d0)))
end forall
eq_arr=dlog(eq_arr)
d1=(Params%loga_table(2)-Params%loga_table(1))/(eq_arr(2)-eq_arr(1))
d2=(Params%loga_table(ntable)-Params%loga_table(ntable-1))/(eq_arr(ntable)-eq_arr(ntable-1))      
call spline(eq_arr(1:ntable),(Params%loga_table(1:ntable)),ntable,d1,d2,eq_arr_buff(1:ntable))
call spline_out(eq_arr(1:ntable),Params%loga_table(1:ntable),&
     &eq_arr_buff(1:ntable)&
     &,ntable,0.0d0,Params%aeq)
     Params%aeq=dexp(Params%aeq)
!Sometimes this spline breaks if a_osc<a_eq, in that case simpler expressions can be used
	regzeq=(Params%omegah2_rad+sum(lhsqcont_massive))/(omegah2_b+omegah2_dm+omegah2_ax)
	if (Params%a_osc.lt.regzeq) then
		Params%aeq=regzeq
	endif
!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!put axion energy density in units used in camb time integrations
forall(i=1:ntable)
	Params%grhoax_table(i)=Params%grhoax_table(i)-dlog(hsq)
end forall
Params%drefp_hsq=drefp/hsq
Params%grhoax_table=dlog10(dexp(Params%grhoax_table))

! Now moved to equations_ppf.f90
!!! For the fisher code for later:
!open(unit=983, file="/Users/reneehlozek/Code/OxFishDec15_axion/results/cambOutput/grhoax.dat", action="write", status="replace")
!do i=1,ntable
!   write(983,*) dexp(Params%loga_table(i)), Params%grhoax_table(i)
!end do

!close(983)
!!! RH
Params%loga_table=dlog10(dexp(Params%loga_table))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!create spline buffer arrays for all quatntities of interest as a function of time in code, this will allow camb to calculate via interpolation any of the quantities of interest at any time
!EOS of axions
d1=(Params%wax_table(2)-Params%wax_table(1))/(Params%loga_table(2)-Params%loga_table(1))
d2=(Params%wax_table(ntable)-Params%wax_table(ntable-1))/(Params%loga_table(ntable)-Params%loga_table(ntable-1))      
call spline(Params%loga_table(1:ntable),Params%wax_table(1:ntable),ntable,d1,d2,Params%wax_table_buff(1:ntable))
!adiabatic sound speed of axions=Pdot/rhodot
d1=(Params%cs2_table(3)-Params%cs2_table(2))/(Params%loga_table(3)-Params%loga_table(2))
d2=(Params%cs2_table(ntable)-Params%cs2_table(ntable-1))/(Params%loga_table(ntable)-Params%loga_table(ntable-1))
call spline(Params%loga_table(2:ntable),(Params%cs2_table(2:ntable)),ntable-1,d1,d2,Params%cs2_table_buff(2:ntable))
!SPINE SCALAR FIELD ENERGY DENSITY FOR LATER USE IN CAMB
d1=(Params%grhoax_table(2)-Params%grhoax_table(1))/(Params%loga_table(2)-Params%loga_table(1))
d2=(Params%grhoax_table(ntable)-Params%grhoax_table(ntable-1))/(Params%loga_table(ntable)-Params%loga_table(ntable-1))
call spline(Params%loga_table(1:ntable),Params%grhoax_table(1:ntable),ntable,d1,d2,Params%grhoax_table_buff(1:ntable))
!!!!!!!!!!!!!!!!!!!!!!!

!!feed out real (dl) valued version of all this stuff for output
!put scalar field initial condition in planck units
Params%phiinit=vtw*sqrt(6.0d0)
if (Params%use_axfrac) then
   Params%axfrac = Params%axfrac
else
   Params%axfrac = fax
end if

if (Params%axion_isocurvature) then

Params%amp_i = Params%Hinf**2/(pi**2*Params%phiinit**2)
Params%r_val  = 2*(Params%Hinf**2/(pi**2.*Params%InitPower%ScalarPowerAmp(1)))
Params%alpha_ax = Params%amp_i/Params%InitPower%ScalarPowerAmp(1)
!print*, 'computing isocurvature', Params%amp_i, Params%r_val, Params%axfrac**2*Params%amp_i/Params%InitPower%ScalarPowerAmp(1), Params%Hinf
!print*, 'computing isocurvature, Params%amp_i, Params%r_val, Params%axfrac**2*Params%amp_i/Params%InitPower%ScalarPowerAmp(1), Params%Hinf'
end if

!output omega_r
Params%omegar=Params%omegah2_rad/hsq

!!!!Timing Stuff
!clock_stop = 0.0
!call cpu_time(clock_stop)
!print*,'axion bg subroutine timing:', clock_stop - clock_start
end subroutine w_evolve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!
! Begin derivative routine
!!!!!!!!!!!!!!!!!!!!!!!!

subroutine derivs(a,v,dvt_dloga,omegah2_regm,omegah2_rad,omegah2_lambda,omk,hsq,maxion_twiddle,badflag,&
	&lhsqcont_massless,lhsqcont_massive,Nu_mass_eigenstates,Nu_masses)
use constants
use Precision
implicit none
integer badflag,Nu_mass_eigenstates
real(dl) a
real(dl) dvt_dloga(1:2),dvt_da(1:2),lhr
real(dl) v(1:2)
real(dl) omegah2_regm,omegah2_rad,omegah2_lambda
real(dl) maxion_twiddle,Nu_masses(Nu_mass_eigenstates)
real(dl) omk,hsq,lhsqcont_massless,lhsqcont_massive(Nu_mass_eigenstates)


!Compute hubble dimensionless 
call lh(omegah2_regm,omegah2_rad,omegah2_lambda,omk,hsq,maxion_twiddle,a,v,lhr,badflag,&
	&lhsqcont_massless,lhsqcont_massive,Nu_mass_eigenstates,Nu_masses)

!calculate fluid derivatives (d/dloga )for next step
!Solving Equation described and defined in top part of this fortran file

! First reported by astralsight5 in https://github.com/dgrin1/axionCAMB/issues/6.

!Confirmed by RL and DG

!Correction 4/8 -- dimensionless H normalization error -- see notebook and or latex note

! dvt_da(1)=v(2)/(a*(lhr))
! dvt_da(2)=-2.0d0*v(2)/(a)-(maxion_twiddle**2.0d0)*a&
!       *v(1)/(lhr)
      

!correct code 4/8
dvt_da(1)=v(2)*dsqrt(hsq)/(a*(lhr))
dvt_da(2)=-2.0d0*v(2)/(a)-(maxion_twiddle**2.0d0)*a&
      *v(1)*dsqrt(hsq)/(lhr)      

dvt_dloga(1:2)=a*dvt_da(1:2)
end subroutine derivs





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Compute H/(100 km/s/Mpc) at any moment in our cosmological history
subroutine lh(omegah2_regm,omegah2_rad,omegah2_lambda,omk,hsq,maxion_twiddle,a,v,littlehfunc,&
	&badflag,&
	&lhsqcont_massless,lhsqcont_massive,Nu_mass_eigenstates,nu_masses)
use constants
use Precision
use MassiveNu
implicit none
integer badflag,Nu_mass_eigenstates,i
real(dl) omegah2_regm,omegah2_rad,omegah2_lambda,omk,hsq,maxion_twiddle,a
real(dl) v(1:2),littlehfunc,lhsqcont_massless,lhsqcont_massive(Nu_mass_eigenstates)
real(dl) mass_correctors(Nu_mass_eigenstates),rhonu
real(dl) nu_masses(Nu_mass_eigenstates)

!Compute corrections to neutrino energy density from massless to massive case
do i=1,Nu_mass_eigenstates,1
call Nu_rho(a*nu_masses(i),rhonu)
mass_correctors(i)=rhonu
enddo
!!!!!!!


!Compute H/(100 km /s/Mpc), contributions from normal stuff
littlehfunc=(omegah2_regm/(a**3.0d0)+omegah2_rad/(a**4.0d0))+&
	&sum(lhsqcont_massive*mass_correctors,Nu_mass_eigenstates)/(a**4.0d0)

!Contributions from cosmological constant and scalar field
littlehfunc=littlehfunc+omegah2_lambda+((v(2)/a)**2.0d0)
littlehfunc=littlehfunc+(maxion_twiddle*v(1))**2.0d0
littlehfunc=littlehfunc*(a**2.0d0)+omk*hsq
! DM flag histories where h goes negative, i.e. do not want collapsing universe
if (littlehfunc .le. 0.0d0) then
        badflag=1
endif
if (isnan(littlehfunc)) then 
badflag = 1
endif
!

littlehfunc=dsqrt(littlehfunc)

end subroutine


!!!!! Take next step in integration using method described above

subroutine next_step(a,v,kvec,kfinal,avec,omegah2_regm,omegah2_rad,omegah2_lambda,omk,&
	&hsq,maxion_twiddle,badflag,dloga,nstep,cmat,&
	&lhsqcont_massless,lhsqcont_massive,Nu_mass_eigenstates,Nu_masses)
use constants
use Precision
implicit none

real(dl) hsq,a,v(1:2),kvec(1:nstep,1:2),omegah2_regm,omegah2_rad,omegah2_lambda,maxion_twiddle,dloga,omk
real(dl) vfeed(1:2),cmat(1:nstep,1:nstep),kfinal(1:2),avec(1:nstep)
integer Nu_mass_eigenstates
integer nstep,cp,m,badflag
real(dl) lhsqcont_massless,lhsqcont_massive(Nu_mass_eigenstates)
real(dl) Nu_masses(Nu_mass_eigenstates)

	kvec=0.0d0
kfinal=0.0d0
			call derivs(a,v(1:2),kvec(1,1:2),omegah2_regm,omegah2_rad,omegah2_lambda,omk,hsq,maxion_twiddle,badflag,&
	&lhsqcont_massless,lhsqcont_massive,Nu_mass_eigenstates,Nu_masses)
			kvec(1,1:2)=kvec(1,1:2)*dloga
			do m=1,nstep,1
			do cp=1,2,1
	vfeed(cp)=dot_product(cmat(m,1:m),kvec(1:m,cp))

	enddo
     if (m .le. (nstep-1)) then
	call derivs(a*dexp(dloga*avec(m)),v(1:2)+vfeed(1:2),&
	&kvec(m+1,1:2),omegah2_regm,omegah2_rad,omegah2_lambda,omk,hsq&
     &,maxion_twiddle,badflag,lhsqcont_massless,lhsqcont_massive,Nu_mass_eigenstates,Nu_masses)
kvec(m+1,1:2)=kvec(m+1,1:2)*dloga
else
call derivs(a*dexp(dloga*avec(m)),v(1:2)+vfeed(1:2),&
	&kfinal(1:2),omegah2_regm,omegah2_rad,omegah2_lambda,omk,hsq&
	&,maxion_twiddle,badflag,lhsqcont_massless,lhsqcont_massive,Nu_mass_eigenstates,Nu_masses)
	kfinal(1:2)=kfinal(1:2)*dloga
endif
enddo
end subroutine

!!!!!!!!








end module axion_background
