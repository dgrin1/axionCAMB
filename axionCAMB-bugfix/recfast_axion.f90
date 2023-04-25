!Recombination module for CAMB, using RECFAST

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!C Integrator for Cosmic Recombination of Hydrogen and Helium,
!C developed by Douglas Scott (dscott@astro.ubc.ca)
!C based on calculations in the paper Seager, Sasselov & Scott
!C (ApJ, 523, L1, 1999).
!and "fudge" updates in Wong, Moss & Scott (2008).
!C
!C Permission to use, copy, modify and distribute without fee or royalty at
!C any tier, this software and its documentation, for any purpose and without
!C fee or royalty is hereby granted, provided that you agree to comply with
!C the following copyright notice and statements, including the disclaimer,
!C and that the same appear on ALL copies of the software and documentation,
!C including modifications that you make for internal use or for distribution:
!C
!C Copyright 1999-2010 by University of British Columbia.  All rights reserved.
!C
!C THIS SOFTWARE IS PROVIDED "AS IS", AND U.B.C. MAKES NO 
!C REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.  
!C BY WAY OF EXAMPLE, BUT NOT LIMITATION,
!c U.B.C. MAKES NO REPRESENTATIONS OR WARRANTIES OF 
!C MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT 
!C THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE 
!C ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.   
!C
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!CN     Name:        RECFAST
!CV     Version:     1.5.2
!C 
!CP     Purpose:  Calculate ionised fraction as a function of redshift.
!CP            Solves for H and He simultaneously, and includes
!CP           H "fudge factor" for low z effect, as well as
!CP           HeI fudge factor.
!C
!CD     Description: Solves for ionisation history since recombination
!CD     using the equations in Seager, Sasselov & Scott (ApJ, 1999).
!CD     The Cosmological model can be flat or open.
!CD	 The matter temperature is also followed, with an update from
!CD	 Scott & Scott (2009).
!CD	 The values for \alpha_B for H are from Hummer (1994).
!CD	 The singlet HeI coefficient is a fit from the full code.
!CD	 Additional He "fudge factors" are as described in Wong, Moss
!CD	 and Scott (2008).
!CD	 Extra fitting function included (in optical depth) to account
!CD	 for extra H physics described in Rubino-Martin et al. (2010).
!CD	 Care is taken to use the most accurate constants.
!C            
!CA     Arguments:
!CA     Name, Description
!CA     real(dl) throughout
!CA
!CA     z is redshift - W is sqrt(1+z), like conformal time
!CA     x is total ionised fraction, relative to H
!CA     x_H is ionized fraction of H - y(1) in R-K routine
!CA     x_He is ionized fraction of He - y(2) in R-K routine
!CA       (note that x_He=n_He+/n_He here and not n_He+/n_H)
!CA     Tmat is matter temperature - y(3) in R-K routine
!CA     f's are the derivatives of the Y's
!CA     alphaB is case B recombination rate
!CA     alpHe is the singlet only HeII recombination rate
!CA     a_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
!CA     b_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
!CA     c_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
!CA     d_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
!CA     a_VF is Verner and Ferland type fitting parameter for Helium
!CA     b_VF is Verner and Ferland type fitting parameter for Helium
!CA     T_0 is Verner and Ferland type fitting parameter for Helium
!CA     T_1 is Verner and Ferland type fitting parameter for Helium
!CA     Tnow is the observed CMB temperature today
!CA     Yp is the primordial helium abundace
!CA     fHe is He/H number ratio = Yp/4(1-Yp)
!CA     Trad and Tmat are radiation and matter temperatures
!CA	    epsilon is the approximate difference (=Trad-Tmat) at high z
!CA     OmegaB is Omega in baryons today
!CA     H is Hubble constant in units of 100 km/s/Mpc
!CA     HO is Hubble constant in SI units
!CA     bigH is 100 km/s/Mpc in SI units
!CA	    Hz is the value of H at the specific z (in ION)
!CA     G is grvitational constant
!CA     n is number density of hydrogen
!CA     Nnow is number density today
!CA     x0 is initial ionized fraction
!CA     x_H0 is initial ionized fraction of Hydrogen
!CA     x_He0 is initial ionized fraction of Helium
!CA     rhs is dummy for calculating x0
!CA     zinitial and zfinal are starting and ending redshifts
!CA     zeq is the redshift of matter-radiation equality
!CA     zstart and zend are for each pass to the integrator
!CA     C,k_B,h_P: speed of light, Boltzmann's and Planck's constants
!CA     m_e,m_H: electron mass and mass of H atom in SI
!CA     not4: ratio of 4He atomic mass to 1H atomic mass
!CA     sigma: Thomson cross-section
!CA     a_rad: radiation constant for u=aT^4
!CA     Lambda: 2s-1s two photon rate for Hydrogen
!CA     Lambda_He: 2s-1s two photon rate for Helium
!CA     DeltaB: energy of first excited state from continuum = 3.4eV
!CA     DeltaB_He: energy of first excited state from cont. for He = 3.4eV
!CA     L_H_ion: level for H ionization in m^-1
!CA     L_H_alpha: level for H Ly alpha in m^-1
!CA     L_He1_ion: level for HeI ionization
!CA     L_He2_ion: level for HeII ionization
!CA     L_He_2s: level for HeI 2s
!CA     L_He_2p: level for HeI 2p (21P1-11S0) in m^-1
!CA     Lalpha: Ly alpha wavelength in SI
!CA     Lalpha_He: Helium I 2p-1s wavelength in SI
!CA     mu_H,mu_T: mass per H atom and mass per particle
!CA     H_frac: follow Tmat when t_Compton / t_Hubble > H_frac
!CA     CDB=DeltaB/k_B                     Constants derived from B1,B2,R
!CA     CDB_He=DeltaB_He/k_B  n=2-infinity for He in Kelvin
!CA     CB1=CDB*4.         Lalpha and sigma_Th, calculated
!CA     CB1_He1: CB1 for HeI ionization potential
!CA     CB1_He2: CB1 for HeII ionization potential
!CA     CR=2*Pi*(m_e/h_P)*(k_B/h_P)  once and passed in a common block
!CA     CK=Lalpha**3/(8.*Pi)
!CA     CK_He=Lalpha_He**3/(8.*Pi)
!CA     CL=C*h_P/(k_B*Lalpha)
!CA     CL_He=C*h_P/(k_B*Lalpha_He)
!CA     CT=(8./3.)*(sigma/(m_e*C))*a
!CA     Bfact=exp((E_2p-E_2s)/kT)    Extra Boltzmann factor
!CA b_He= "fudge factor" for HeI, to approximate higher z behaviour
!CA Heswitch=integer for modifying HeI recombination
!CA Parameters and quantities to describe the extra triplet states
!CA  and also the continuum opacity of H, with a fitting function
!CA  suggested by KIV, astro-ph/0703438
!CA a_trip: used to fit HeI triplet recombination rate
!CA b_trip: used to fit HeI triplet recombination rate
!CA L_He_2Pt: level for 23P012-11S0 in m^-1
!CA L_He_2St: level for 23S1-11S0 in m^-1
!CA L_He2St_ion: level for 23S1-continuum in m^-1
!CA A2P_s: Einstein A coefficient for He 21P1-11S0
!CA A2P_t: Einstein A coefficient for He 23P1-11S0    
!CA sigma_He_2Ps: H ionization x-section at HeI 21P1-11S0 freq. in m^2
!CA sigma_He_2Pt: H ionization x-section at HeI 23P1-11S0 freq. in m^2
!CA CL_PSt = h_P*C*(L_He_2Pt - L_He_2st)/k_B
!CA CfHe_t: triplet statistical correction
!CA	Hswitch is an boolean for modifying the H recombination
!CA	AGauss1 is the amplitude of the 1st Gaussian for the H fudging
!CA	AGauss2 is the amplitude of the 2nd Gaussian for the H fudging
!CA	zGauss1 is the ln(1+z) central value of the 1st Gaussian
!CA	zGauss2 is the ln(1+z) central value of the 2nd Gaussian
!CA	wGauss1 is the width of the 1st Gaussian
!CA	wGauss2 is the width of the 2nd Gaussian


!CA     tol: tolerance for the integrator
!CA     cw(24),w(3,9): work space for DVERK
!CA     Ndim: number of d.e.'s to solve (integer)
!CA     Nz: number of output redshitf (integer)
!CA     I: loop index (integer)
!CA     ind,nw: work-space for DVERK (integer)
!C
!CF     File & device access:
!CF     Unit /I,IO,O  /Name (if known)
!C
!CM     Modules called:
!CM     DVERK (numerical integrator)
!CM     GET_INIT (initial values for ionization fractions)
!CM     ION (ionization and Temp derivatives)
!C
!CC     Comments:
!CC     none
!C
!CH     History:
!CH     CREATED            (simplest version) 19th March 1989
!CH     RECREATED    11th January 1995
!CH               includes variable Cosmology
!CH               uses DVERK integrator
!CH               initial conditions are Saha
!CH     TESTED              a bunch, well, OK, not really
!CH     MODIFIED     January 1995 (include Hummer's 1994 alpha table)
!CH               January 1995 (include new value for 2s-1s rate)
!CH               January 1995 (expand comments)
!CH               March 1995 (add Saha for Helium)
!CH               August 1997 (add HeII alpha table)
!CH               July 1998 (include OmegaT correction and H fudge factor)
!CH               Nov 1998 (change Trad to Tmat in Rup)
!CH               Jan 1999 (tidied up for public consumption)
!CH               Sept 1999 (switch to formula for alpha's, fix glitch)
!CH                  Sept 1999 modified to CMBFAST by US & MZ          
!CH                     Nov 1999 modified for F90 and CAMB (AML)
!CH                     Aug 2000 modified to prevent overflow erorr in He_Boltz (AML)
!CH                     Feb 2001 corrected fix of Aug 2000 (AML)
!CH                     Oct 2001 fixed error in hubble parameter, now uses global function (AML)
!                       March 2003 fixed bugs reported by savita gahlaut
!                       March 2005 added option for corrections from astro-ph/0501672.
!                                  thanks to V.K.Dubrovich, S.I.Grachev
!                       June 2006 defined RECFAST_fudge as free parameter (AML)
!                       October 2006 (included new value for G)
!                       October 2006 (improved m_He/m_H to be "not4")
!                       October 2006 (fixed error, x for x_H in part of f(1))
!CH              January 2008 (improved HeI recombination effects,
!CH                       including HeI rec. fudge factor)
!                Feb 2008   Recfast 1.4 changes above added (AML)
!                           removed Dubrovich option (wrong anyway)
!CH   			 Sept 2008 (added extra term to make transition, smoother for Tmat evolution)
!                Sept 2008 Recfast 1.4.2 changes above added (AML) 
!                          General recombination module structure, fix to make He x_e smooth also in recfast (AML)
!CH		      	 Jan 2010 (added fitting function to modify K
!CH			 	 to match x_e(z) for new H physics)
!AL             June 2012 updated fudge parameters to match HyRec and CosmoRec (AML)                    
!AL             Sept 2012 changes now in public recfast, version number changed to match Recfast 1.5.2. 
!!      ===============================================================

       module RECDATA
        use constants
    !    use ModelParams
        implicit none
         

        real(dl) Lambda,DeltaB,DeltaB_He,Lalpha,mu_H,mu_T,H_frac
        real(dl) Lambda_He,Lalpha_He,Bfact,CK_He,CL_He
        real(dl) L_H_ion,L_H_alpha,L_He1_ion,L_He2_ion,L_He_2s,L_He_2p
        real(dl) CB1,CDB,CR,CK,CL,CT,CB1_He1,CB1_He2,CDB_He,fu
        real(dl) A2P_s,A2P_t,sigma_He_2Ps,sigma_He_2Pt
        real(dl)  L_He_2Pt,L_He_2St,L_He2St_ion


        real(dl), parameter :: bigH=100.0D3/Mpc !Ho in s-1
        real(dl), parameter :: sigma = sigma_thomson
        real(dl), parameter :: not4  = mass_ratio_He_H    !mass He/H atom

        real(dl) Tnow,HO
        integer :: n_eq = 3

!The following only used for approximations where small effect
       real(dl) OmegaK, OmegaT, z_eq


!Fundamental constants in SI units
!      ("not4" pointed out by Gary Steigman)

        data    Lambda      /8.2245809d0/
        data    Lambda_He   /51.3d0/    !new value from Dalgarno
        data    L_H_ion     /1.096787737D7/ !level for H ion. (in m^-1)
        data    L_H_alpha   /8.225916453D6/ !averaged over 2 levels
        data    L_He1_ion   /1.98310772D7/  !from Drake (1993)
        data    L_He2_ion   /4.389088863D7/ !from JPhysChemRefData (1987)
        data    L_He_2s     /1.66277434D7/  !from Drake (1993)
        data    L_He_2p     /1.71134891D7/  !from Drake (1993)
!   2 photon rates and atomic levels in SI units

        data    A2P_s       /1.798287D9/    !Morton, Wu & Drake (2006)
        data    A2P_t       /177.58D0/      !Lach & Pachuski (2001)
        data    L_He_2Pt    /1.690871466D7/ !Drake & Morton (2007)
        data    L_He_2St    /1.5985597526D7/ !Drake & Morton (2007)
        data    L_He2St_ion /3.8454693845D6/ !Drake & Morton (2007)
        data    sigma_He_2Ps    /1.436289D-22/  !Hummer & Storey (1998)
        data    sigma_He_2Pt    /1.484872D-22/  !Hummer & Storey (1998)
!    Atomic data for HeI 


       end module RECDATA


        module Recombination
        use constants
        use AMLUtils
        implicit none
        private

        real(dl), parameter ::  zinitial = 1e4_dl !highest redshift
        real(dl), parameter ::  zfinal=0._dl
        integer,  parameter :: Nz=10000  
        real(dl), parameter :: delta_z = (zinitial-zfinal)/Nz

        integer, parameter ::  RECFAST_Heswitch_default = 6
        real(dl), parameter :: RECFAST_fudge_He_default = 0.86_dl !Helium fudge parameter
        logical, parameter  :: RECFAST_Hswitch_default = .true. !include H corrections (v1.5, 2010)
        real(dl), parameter :: RECFAST_fudge_default = 1.14_dl
        real(dl), parameter :: RECFAST_fudge_default2 = 1.105d0 + 0.02d0
                      !fudge parameter if RECFAST_Hswitch
        real(dl) :: AGauss1 =      -0.14D0  !Amplitude of 1st Gaussian
        real(dl) :: AGauss2 =       0.079D0 ! 0.05D0  !Amplitude of 2nd Gaussian
        real(dl) :: zGauss1 =       7.28D0  !ln(1+z) of 1st Gaussian
        real(dl) :: zGauss2=        6.73D0  !ln(1+z) of 2nd Gaussian
        real(dl) :: wGauss1=        0.18D0  !Width of 1st Gaussian
        real(dl) :: wGauss2=        0.33D0  !Width of 2nd Gaussian
        !Gaussian fits for extra H physics (fit by Adam Moss , modified by Antony Lewis)

        type RecombinationParams

          real(dl) :: RECFAST_fudge 
          real(dl) :: RECFAST_fudge_He 
          integer  :: RECFAST_Heswitch
          logical  :: RECFAST_Hswitch  
         !0) no change from old Recfast'
         !1) full expression for escape probability for singlet'
         !'   1P-1S transition'
         !2) also including effect of contiuum opacity of H on HeI'
         !'   singlet (based in fitting formula suggested by'
         !'   Kholupenko, Ivanchik & Varshalovich, 2007)'
         !3) only including recombination through the triplets'
         !4) including 3 and the effect of the contiuum '
         !'   (although this is probably negligible)' 
         !5) including only 1, 2 and 3'
         !6) including all of 1 to 4'
   
        end  type RecombinationParams

        character(LEN=*), parameter :: Recombination_Name = 'Recfast_1.5'
      
        real(dl) zrec(Nz),xrec(Nz),dxrec(Nz), Tsrec(Nz) ,dTsrec(Nz), tmrec(Nz),dtmrec(Nz)

        real(dl), parameter :: Do21cm_mina = 1/(1+900.) !at which to start evolving Delta_TM
        logical, parameter :: evolve_Ts = .false. !local equilibrium is very accurate
        real(dl), parameter :: Do21cm_minev = 1/(1+400.) !at which to evolve T_s
       
 
        real(dl), parameter :: B01 = 3*B10
        real(dl) :: NNow, fHe 
               

        logical :: Do21cm = .false.
        logical :: doTmatTspin = .false.

        real(dl) :: recombination_saha_z !Redshift at which saha OK
        real(dl) :: recombination_saha_tau !set externally


        public RecombinationParams, Recombination_xe, Recombination_tm,Recombination_ts ,Recombination_init,   &
               Recombination_ReadParams, Recombination_SetDefParams, Recombination_Validate, Recombination_Name, &
               kappa_HH_21cm,kappa_eH_21cm,kappa_pH_21cm, &
               Do21cm, doTmatTspin, Do21cm_mina, dDeltaxe_dtau, &
               recombination_saha_tau, recombination_saha_z
 
       contains



         subroutine Recombination_ReadParams(R, Ini)
          use IniFile
          Type(RecombinationParams) :: R
          Type(TIniFile) :: Ini


             R%RECFAST_fudge_He = Ini_Read_Double_File(Ini,'RECFAST_fudge_He',RECFAST_fudge_He_default)
             R%RECFAST_Heswitch = Ini_Read_Int_File(Ini, 'RECFAST_Heswitch',RECFAST_Heswitch_default)
             R%RECFAST_Hswitch = Ini_Read_Logical_File(Ini, 'RECFAST_Hswitch',RECFAST_Hswitch_default)
             R%RECFAST_fudge = Ini_Read_Double_File(Ini,'RECFAST_fudge',RECFAST_fudge_default)
             AGauss1 = Ini_REad_Double_File(Ini,'AGauss1',AGauss1)
             AGauss2 = Ini_REad_Double_File(Ini,'AGauss2',AGauss2)
             zGauss1 = Ini_REad_Double_File(Ini,'zGauss1',zGauss1)
             zGauss2 = Ini_REad_Double_File(Ini,'zGauss2',zGauss2)
             wGauss1 = Ini_REad_Double_File(Ini,'wGauss1',wGauss1)
             wGauss2 = Ini_REad_Double_File(Ini,'wGauss2',wGauss2)
  
             if (R%RECFAST_Hswitch) then
                R%RECFAST_fudge = R%RECFAST_fudge - (RECFAST_fudge_default - RECFAST_fudge_default2)
             end if   
         end subroutine Recombination_ReadParams 

        subroutine Recombination_SetDefParams(R)
         type (RecombinationParams) ::R

      
          R%RECFAST_fudge = RECFAST_fudge_default
          R%RECFAST_fudge_He = RECFAST_fudge_He_default !Helium fudge parameter
          R%RECFAST_Heswitch = RECFAST_Heswitch_default
          R%RECFAST_Hswitch =  RECFAST_Hswitch_default
          if (R%RECFAST_Hswitch) then
                R%RECFAST_fudge = RECFAST_fudge_default2
          end if   
    
        end subroutine Recombination_SetDefParams


        subroutine Recombination_Validate(R, OK)
          Type(RecombinationParams),intent(in) :: R
          logical, intent(inout) :: OK
 
              if (R%RECFAST_Heswitch<0 .or. R%RECFAST_Heswitch > 6) then
                     OK = .false.
                     write(*,*) 'RECFAST_Heswitch unknown'
               end if
             
         end subroutine Recombination_Validate


        function Recombination_tm(a)
        use RECDATA, only : Tnow
        real(dl) zst,a,z,az,bz,Recombination_tm
        integer ilo,ihi
        
        if (.not. doTmatTspin) stop 'RECFAST: Recombination_tm not stored'
        z=1/a-1
        if (z >= zrec(1)) then
            Recombination_tm=Tnow/a
        else
         if (z <=zrec(nz)) then
          Recombination_tm=Tmrec(nz)
         else
          zst=(zinitial-z)/delta_z
          ihi= int(zst)
          ilo = ihi+1
          az=zst - int(zst)
          bz=1-az     
          Recombination_tm=az*Tmrec(ilo)+bz*Tmrec(ihi)+ &
           ((az**3-az)*dTmrec(ilo)+(bz**3-bz)*dTmrec(ihi))/6._dl
         endif
        endif

        end function Recombination_tm


        function Recombination_ts(a)
        !zrec(1) is zinitial-delta_z
        real(dl), intent(in) :: a
        real(dl) zst,z,az,bz,Recombination_ts
        integer ilo,ihi
        
        z=1/a-1
        if (z.ge.zrec(1)) then
          Recombination_ts=tsrec(1)
        else
         if (z.le.zrec(nz)) then
          Recombination_ts=tsrec(nz)
         else
          zst=(zinitial-z)/delta_z
          ihi= int(zst)
          ilo = ihi+1
          az=zst - int(zst)
          bz=1-az     

          Recombination_ts=az*tsrec(ilo)+bz*tsrec(ihi)+ &
           ((az**3-az)*dtsrec(ilo)+(bz**3-bz)*dtsrec(ihi))/6._dl
         endif
        endif

        end function Recombination_ts


        function Recombination_xe(a)
        real(dl), intent(in) :: a
        real(dl) zst,z,az,bz,Recombination_xe
        integer ilo,ihi
        
        z=1/a-1
        if (z.ge.zrec(1)) then
          Recombination_xe=xrec(1)
        else
         if (z.le.zrec(nz)) then
          Recombination_xe=xrec(nz)
         else
          zst=(zinitial-z)/delta_z
          ihi= int(zst)
          ilo = ihi+1
          az=zst - int(zst)
          bz=1-az     
          Recombination_xe=az*xrec(ilo)+bz*xrec(ihi)+ &
           ((az**3-az)*dxrec(ilo)+(bz**3-bz)*dxrec(ihi))/6._dl
         endif
        endif

        end function Recombination_xe



!!!!!!!!!!!
!Modified master subroutine for recast
!includes now size of spline tables, aeq, aosc (time when m=dfac*H from axion_background
!.f90), omega_ah^2*aosc^3 for a^-3 scaling of density when a>aosc
!logarithm of scale factor, energy density of axions and spilne buffer all from axion_background.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine Recombination_init(Recomb, OmegaC, OmegaB, Omegan, Omegav, &
      h0inp,tcmb,yp,OmegaAx,omegar,ntable,aeq,aosc,drefp_hsq,loga_table,grhoax_table,grhoax_table_buff)
        !Would love to pass structure as arguments, but F90 would give circular reference...
        !hence mess passing parameters explcitly and non-generally
        !Note recfast only uses OmegaB, h0inp, tcmb and yp - others used only for Tmat approximation where effect small
        use RECDATA
        use AMLUtils
        use constants
        implicit none
        Type (RecombinationParams) :: Recomb
 
        real(dl), save :: last_OmB =0, Last_YHe=0, Last_H0=0, Last_dtauda=0, last_fudge, last_fudgeHe 

        real(dl) Trad,Tmat,Tspin,d0hi,d0lo
        integer I

        real(dl) OmegaB,OmegaC, Omegan, Omegav, H
         
       
        real(dl) z,x,x0,rhs,x_H,x_He,x_H0,x_He0,h0inp
!   Unecessary variable
        real(dl) n
        real(dl) zstart,zend,tcmb
        real(dl) cw(24)
        real(dl), dimension(:,:), allocatable :: w
        real(dl) y(4)
        real(dl) yp
        real(dl) C10, tau_21Ts
        !un-necessary fi find aeq using axion code
  !      real(dl) fnu
       integer ind,nw
 
 		!!!!!!!Some more Axion definitions
       	!Omega_axion, m=dfac*H transition point (scale factor aosc), energy density at that time
		real(dl) OmegaAx,aosc,drefp_hsq
       	!aeq from axion_background.f90
        real(dl) aeq
  	     integer ntable
      	!log scale factor, log energy density and spline bufffer for axion inclusion
      	 real(dl) loga_table(ntable),grhoax_table(ntable),grhoax_table_buff(ntable)
		!!!!!!!!!!!!



!       --- Parameter statements
        real(dl), parameter :: tol=1.D-3                !Tolerance for R-K

        real(dl) dtauda,omegar
        external dtauda, dverk

!       ===============================================================

        if (Last_OmB==OmegaB .and. Last_H0 == h0inp .and. yp == Last_YHe .and. & 
             dtauda(0.2352375823_dl) == Last_dtauda .and. last_fudge == Recomb%RECFAST_fudge &
              .and. last_fudgeHe==Recomb%RECFAST_fudge_He) return
           !This takes up most of the single thread time, so cache if at all possible
           !For example if called with different reionization, or tensor rather than scalar
        
        Last_dtauda =  dtauda(0.2352375823_dl) !Just get it at a random scale factor
        Last_OmB = OmegaB
        Last_H0 = h0inp
        Last_YHe=yp
        last_fudge = Recomb%RECFAST_FUDGE
        last_fudgeHe = Recomb%RECFAST_FUDGE_He

        if (Do21cm) doTmatTspin = .true.


!       write(*,*)'recfast version 1.0'
!       write(*,*)'Using Hummer''s case B recombination rates for H'
!       write(*,*)' with fudge factor = 1.14'
!       write(*,*)'and tabulated HeII singlet recombination rates'
!       write(*,*)

        n_eq = 3
        if (Evolve_Ts) n_eq=4
        allocate(w(n_eq,9))

        recombination_saha_z=0.d0

        Tnow=tcmb
!       These are easy to inquire as input, but let's use simple values
        z = zinitial
!       will output every 1 in z, but this is easily changed also

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !Not general, but only for approx 
 !DM11: notes that the neutrinos are not included here. Gives wrong curvature.
 !DG15: Have added massive and massless neutrino+photon contribution'
 ! Does not make a huge difference
 ! See comment in halofit also.
        OmegaT=OmegaC+OmegaB        !total dark matter + baryons + axions
        OmegaK=1.d0-OmegaT-OmegaAx-OmegaV-omegar-Omegan      !curvature
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  
!       convert the Hubble constant units
        H = H0inp/100._dl
        HO = H*bigH


!       sort out the helium abundance parameters
        mu_H = 1.d0/(1.d0-Yp)           !Mass per H atom
        mu_T = not4/(not4-(not4-1.d0)*Yp)   !Mass per atom
        fHe = Yp/(not4*(1.d0-Yp))       !n_He_tot / n_H_tot


        Nnow = 3._dl*HO*HO*OmegaB/(8._dl*Pi*G*mu_H*m_H)

        n = Nnow * (1._dl+z)**3
 !       fnu = (21.d0/8.d0)*(4.d0/11.d0)**(4.d0/3.d0)
!	(this is explictly for 3 massless neutrinos - change if N_nu.ne.3; but only used for approximation so not critical)
        !z_eq = !(3.d0*(HO*C)**2/(8.d0*Pi*G*a_rad*(1.d0+fnu)*Tnow**4))*(OmegaB+OmegaC)
        z_eq=aeq**(-1.0d0)
        z_eq = z_eq - 1.d0
!print*,z_eq
      
!       Set up some constants so they don't have to be calculated later
        Lalpha = 1.d0/L_H_alpha
        Lalpha_He = 1.d0/L_He_2p
        DeltaB = h_P*C*(L_H_ion-L_H_alpha)
        CDB = DeltaB/k_B
        DeltaB_He = h_P*C*(L_He1_ion-L_He_2s)   !2s, not 2p
        CDB_He = DeltaB_He/k_B
        CB1 = h_P*C*L_H_ion/k_B
        CB1_He1 = h_P*C*L_He1_ion/k_B   !ionization for HeI
        CB1_He2 = h_P*C*L_He2_ion/k_B   !ionization for HeII
        CR = 2.d0*Pi*(m_e/h_P)*(k_B/h_P)
        CK = Lalpha**3/(8.d0*Pi)
        CK_He = Lalpha_He**3/(8.d0*Pi)
        CL = C*h_P/(k_B*Lalpha)
        CL_He = C*h_P/(k_B/L_He_2s) !comes from det.bal. of 2s-1s
        CT = Compton_CT / MPC_in_sec

        Bfact = h_P*C*(L_He_2p-L_He_2s)/k_B

        
!       Matter departs from radiation when t(Th) > H_frac * t(H)
!       choose some safely small number
        H_frac = 1D-3

!       Fudge factor to approximate for low z out of equilibrium effect
        fu=Recomb%RECFAST_fudge

!       Set initial matter temperature
        y(3) = Tnow*(1._dl+z)            !Initial rad. & mat. temperature
        Tmat = y(3)
        y(4) = Tmat
        Tspin = Tmat

        call get_init(z,x_H0,x_He0,x0)
    
        y(1) = x_H0
        y(2) = x_He0

!       OK that's the initial conditions, now start writing output file


!       Set up work-space stuff for DVERK
        ind  = 1
        nw   = n_eq
        do i = 1,24
          cw(i) = 0._dl
        end do

        do i = 1,Nz
!       calculate the start and end redshift for the interval at each z
!       or just at each z
          zstart = zinitial  - real(i-1,dl)*delta_z
          zend   = zinitial  - real(i,dl)*delta_z

! Use Saha to get x_e, using the equation for x_e for ionized helium
! and for neutral helium.
! Everything ionized above z=8000.  First ionization over by z=5000.
! Assume He all singly ionized down to z=3500, then use He Saha until
! He is 99% singly ionized, and *then* switch to joint H/He recombination.

          z = zend
        
          if (zend > 8000._dl) then

            x_H0 = 1._dl
            x_He0 = 1._dl
            x0 = 1._dl+2._dl*fHe
            y(1) = x_H0
            y(2) = x_He0
            y(3) = Tnow*(1._dl+z)
            y(4) = y(3)

          else if(z > 5000._dl)then

            x_H0 = 1._dl
            x_He0 = 1._dl
            rhs = exp( 1.5d0 * log(CR*Tnow/(1._dl+z)) &
                - CB1_He2/(Tnow*(1._dl+z)) ) / Nnow
            rhs = rhs*1._dl            !ratio of g's is 1 for He++ <-> He+
            x0 = 0.5d0 * ( sqrt( (rhs-1._dl-fHe)**2 &
                + 4._dl*(1._dl+2._dl*fHe)*rhs) - (rhs-1._dl-fHe) )
            y(1) = x_H0
            y(2) = x_He0
            y(3) = Tnow*(1._dl+z)
            y(4) = y(3)

          else if(z > 3500._dl)then

            x_H0 = 1._dl
            x_He0 = 1._dl
            x0 = x_H0 + fHe*x_He0
            y(1) = x_H0
            y(2) = x_He0
            y(3) = Tnow*(1._dl+z)
            y(4) = y(3)

          else if(y(2) > 0.99)then

            x_H0 = 1._dl
            rhs = exp( 1.5d0 * log(CR*Tnow/(1._dl+z)) &
                - CB1_He1/(Tnow*(1._dl+z)) ) / Nnow
            rhs = rhs*4._dl            !ratio of g's is 4 for He+ <-> He0
            x_He0 = 0.5d0 * ( sqrt( (rhs-1._dl)**2 &
                + 4._dl*(1._dl+fHe)*rhs )- (rhs-1._dl))
            x0 = x_He0
            x_He0 = (x0 - 1._dl)/fHe
            y(1) = x_H0
            y(2) = x_He0
            y(3) = Tnow*(1._dl+z)
            y(4) = y(3)

          else if (y(1) > 0.99d0) then

            rhs = exp( 1.5d0 * log(CR*Tnow/(1._dl+z)) &
                - CB1/(Tnow*(1._dl+z)) ) / Nnow
            x_H0 = 0.5d0 * (sqrt( rhs**2+4._dl*rhs ) - rhs )
!recdverk is same runge kutte integrator used ins ubroutines, just added pasisng of splines
            call recdverk(Recomb,OmegaAx,omegar,ntable,aeq,aosc,drefp_hsq,loga_table,grhoax_table,&      
      grhoax_table_buff,3,ION,zstart,y,zend,tol,ind,cw,nw,w)
            y(1) = x_H0
            x0 = y(1) + fHe*y(2)
            y(4)=y(3)
          else
            
            call recdverk(Recomb,OmegaAx,omegar,ntable,aeq,aosc,drefp_hsq,loga_table,grhoax_table,&
      grhoax_table_buff,nw,ION,zstart,y,zend,tol,ind,cw,nw,w)
          
            x0 = y(1) + fHe*y(2)
          
          end if
          
          Trad = Tnow * (1._dl+zend)
          Tmat = y(3)
          x_H = y(1)
          x_He = y(2)
          x = x0

          zrec(i)=zend
          xrec(i)=x

    
          if (doTmatTspin) then
              if (Evolve_Ts .and. zend< 1/Do21cm_minev-1 ) then
               Tspin = y(4)
              else 
               C10 = Nnow * (1._dl+zend)**3*(kappa_HH_21cm(Tmat,.false.)*(1-x_H) + kappa_eH_21cm(Tmat,.false.)*x)
               tau_21Ts = line21_const*NNow*(1+zend)*dtauda(1/(1+zend))/1000
        
               Tspin = Trad*( C10/Trad + A10/T_21cm)/(C10/Tmat + A10/T_21cm) + &
                     tau_21Ts/2*A10*( 1/(C10*T_21cm/Tmat+A10) -  1/(C10*T_21cm/Trad+A10) )
          
               y(4) = Tspin
              end if

              tsrec(i) = Tspin
              tmrec(i) = Tmat
           
          end if   

!          write (*,'(5E15.5)') zend, Trad, Tmat, Tspin, x
     
        end do
        
        d0hi=1.0d40
        d0lo=1.0d40
        call spline(zrec,xrec,nz,d0lo,d0hi,dxrec)
        if (doTmatTspin) then
         call spline(zrec,tsrec,nz,d0lo,d0hi,dtsrec)
         call spline(zrec,tmrec,nz,d0lo,d0hi,dtmrec)
        end if
        deallocate(w)
        
        end subroutine Recombination_init

!       ===============================================================
        subroutine GET_INIT(z,x_H0,x_He0,x0)

!       Set up the initial conditions so it will work for general,
!       but not pathological choices of zstart
!       Initial ionization fraction using Saha for relevant species
        use RECDATA
        implicit none
  
        
        real(dl) z,x0,rhs,x_H0,x_He0
  

        if(z > 8000._dl)then

            x_H0 = 1._dl
            x_He0 = 1._dl
            x0 = 1._dl+2._dl*fHe

        else if(z > 3500._dl)then

            x_H0 = 1._dl
            x_He0 = 1._dl
            rhs = exp( 1.5d0 * log(CR*Tnow/(1._dl+z)) &
                - CB1_He2/(Tnow*(1._dl+z)) ) / Nnow
        rhs = rhs*1._dl    !ratio of g's is 1 for He++ <-> He+
        x0 = 0.5d0 * ( sqrt( (rhs-1._dl-fHe)**2 &
                + 4._dl*(1._dl+2._dl*fHe)*rhs) - (rhs-1._dl-fHe) )

        else if(z > 2000._dl)then

        x_H0 = 1._dl
            rhs = exp( 1.5d0 * log(CR*Tnow/(1._dl+z)) &
                - CB1_He1/(Tnow*(1._dl+z)) ) / Nnow
        rhs = rhs*4._dl    !ratio of g's is 4 for He+ <-> He0
            x_He0 = 0.5d0  * ( sqrt( (rhs-1._dl)**2 + 4._dl*(1._dl+fHe)*rhs )- (rhs-1._dl))
            x0 = x_He0
            x_He0 = (x0 - 1._dl)/fHe

        else

            rhs = exp( 1.5d0 * log(CR*Tnow/(1._dl+z)) &
                - CB1/(Tnow*(1._dl+z)) ) / Nnow
            x_H0 = 0.5d0 * (sqrt( rhs**2+4._dl*rhs ) - rhs )
            x_He0 = 0._dl
            x0 = x_H0

        end if

        
        end subroutine GET_INIT



        subroutine ION(Recomb,OmegaAx,omegar,ntable,aeq,aosc,drefp_hsq,loga_table,grhoax_table,&
      grhoax_table_buff,Ndim,z,Y,f)
        use RECDATA
     !   use ModelParams
        implicit none

        integer Ndim
        Type (RecombinationParams) :: Recomb

        real(dl) z,x,n,n_He,Trad,Tmat,Tspin,x_H,x_He, Hz
        real(dl) y(Ndim),f(Ndim)
        real(dl) Rup,Rdown,K,K_He,Rup_He,Rdown_He,He_Boltz
        real(dl) timeTh,timeH
        real(dl) a_VF,b_VF,T_0,T_1,sq_0,sq_1,a_PPB,b_PPB,c_PPB,d_PPB
        real(dl) tauHe_s,pHe_s
        real(dl) a_trip,b_trip,Rdown_trip,Rup_trip
        real(dl) Doppler,gamma_2Ps,pb,qb,AHcon
        real(dl) tauHe_t,pHe_t,CL_PSt,CfHe_t,gamma_2Pt
        real(dl) epsilon, omegar
        integer Heflag,ntable
        real(dl) dtauda
        real (dl) gr,sfac,aosc,drefp_hsq
        real(dl) C10, dHdz
        real (dl) OmegaAx,aeq,loga_table(ntable),grhoax_table(ntable),grhoax_table_buff(ntable),dorp,dorpa,deriv_eps
        external dtauda
        integer i

!       the Pequignot, Petitjean & Boisson fitting parameters for Hydrogen    
        a_PPB = 4.309d0
        b_PPB = -0.6166d0
        c_PPB = 0.6703d0
        d_PPB = 0.5300d0
!       the Verner and Ferland type fitting parameters for Helium
!       fixed to match those in the SSS papers, and now correct
        a_VF = 10.d0**(-16.744d0)
        b_VF = 0.711d0
        T_0 = 10.d0**(0.477121d0)   !3K
        T_1 = 10.d0**(5.114d0)
!      fitting parameters for HeI triplets
!      (matches Hummer's table with <1% error for 10^2.8 < T/K < 10^4)

        a_trip = 10.d0**(-16.306d0)
        b_trip = 0.761D0

       
        x_H = y(1)
        x_He = y(2)
        x = x_H + fHe * x_He
        Tmat = y(3)
!        Tspin = y(4)

        n = Nnow * (1._dl+z)**3
        n_He = fHe * Nnow * (1._dl+z)**3
        Trad = Tnow * (1._dl+z)

        Hz = 1/dtauda(1/(1._dl+z))*(1._dl+z)**2/MPC_in_sec       
     

!       Get the radiative rates using PPQ fit, identical to Hummer's table
        
        Rdown=1.d-19*a_PPB*(Tmat/1.d4)**b_PPB &
                /(1._dl+c_PPB*(Tmat/1.d4)**d_PPB)
        Rup = Rdown * (CR*Tmat)**(1.5d0)*exp(-CDB/Tmat)
      
!       calculate He using a fit to a Verner & Ferland type formula
        sq_0 = sqrt(Tmat/T_0)
        sq_1 = sqrt(Tmat/T_1)
!       typo here corrected by Wayne Hu and Savita Gahlaut
        Rdown_He = a_VF/(sq_0*(1.d0+sq_0)**(1.d0-b_VF))
        Rdown_He = Rdown_He/(1.d0+sq_1)**(1.d0+b_VF)
        Rup_He = Rdown_He*(CR*Tmat)**(1.5d0)*exp(-CDB_He/Tmat)
        Rup_He = 4.d0*Rup_He    !statistical weights factor for HeI
!       Avoid overflow (pointed out by Jacques Roland)
        if((Bfact/Tmat) > 680.d0)then
          He_Boltz = exp(680.d0)
        else
          He_Boltz = exp(Bfact/Tmat)
        end if
!	now deal with H and its fudges
    if (.not. Recomb%RECFAST_Hswitch) then
      K = CK/Hz		!Peebles coefficient K=lambda_a^3/8piH
    else
!c	fit a double Gaussian correction function
      K = CK/Hz*(1.0d0 &
        +AGauss1*exp(-((log(1.0d0+z)-zGauss1)/wGauss1)**2.d0) &
        +AGauss2*exp(-((log(1.0d0+z)-zGauss2)/wGauss2)**2.d0))
    end if        
        
        
 !  add the HeI part, using same T_0 and T_1 values
    Rdown_trip = a_trip/(sq_0*(1.d0+sq_0)**(1.0-b_trip))
    Rdown_trip = Rdown_trip/((1.d0+sq_1)**(1.d0+b_trip))
    Rup_trip = Rdown_trip*dexp(-h_P*C*L_He2St_ion/(k_B*Tmat))
    Rup_trip = Rup_trip*((CR*Tmat)**(1.5d0))*(4.d0/3.d0)
!   last factor here is the statistical weight

!       try to avoid "NaN" when x_He gets too small
    if ((x_He.lt.5.d-9) .or. (x_He.gt.0.98d0)) then 
      Heflag = 0
    else
      Heflag = Recomb%RECFAST_Heswitch
    end if
    if (Heflag.eq.0)then        !use Peebles coeff. for He
      K_He = CK_He/Hz
    else    !for Heflag>0       !use Sobolev escape probability
      tauHe_s = A2P_s*CK_He*3.d0*n_He*(1.d0-x_He)/Hz
      pHe_s = (1.d0 - dexp(-tauHe_s))/tauHe_s
      K_He = 1.d0/(A2P_s*pHe_s*3.d0*n_He*(1.d0-x_He))
!      if (((Heflag.eq.2) .or. (Heflag.ge.5)) .and. x_H < 0.99999d0) then 
      if (((Heflag.eq.2) .or. (Heflag.ge.5)) .and. x_H < 0.9999999d0) then
       !AL changed July 08 to get smoother Helium

!   use fitting formula for continuum opacity of H
!   first get the Doppler width parameter
        Doppler = 2.D0*k_B*Tmat/(m_H*not4*C*C)
        Doppler = C*L_He_2p*dsqrt(Doppler)
        gamma_2Ps = 3.d0*A2P_s*fHe*(1.d0-x_He)*C*C &
            /(dsqrt(Pi)*sigma_He_2Ps*8.d0*Pi*Doppler*(1.d0-x_H)) &
            /((C*L_He_2p)**2.d0)
        pb = 0.36d0  !value from KIV (2007)
        qb = Recomb%RECFAST_fudge_He
!   calculate AHcon, the value of A*p_(con,H) for H continuum opacity
        AHcon = A2P_s/(1.d0+pb*(gamma_2Ps**qb))
        K_He=1.d0/((A2P_s*pHe_s+AHcon)*3.d0*n_He*(1.d0-x_He))
      end if
      if (Heflag.ge.3) then     !include triplet effects
        tauHe_t = A2P_t*n_He*(1.d0-x_He)*3.d0
        tauHe_t = tauHe_t /(8.d0*Pi*Hz*L_He_2Pt**(3.d0))
        pHe_t = (1.d0 - dexp(-tauHe_t))/tauHe_t
        CL_PSt = h_P*C*(L_He_2Pt - L_He_2st)/k_B
        if ((Heflag.eq.3) .or. (Heflag.eq.5).or.(x_H.gt.0.99999d0)) then !Recfast 1.4.2 (?)
!        if ((Heflag.eq.3) .or. (Heflag.eq.5) .or. x_H >= 0.9999999d0) then    !no H cont. effect
            CfHe_t = A2P_t*pHe_t*dexp(-CL_PSt/Tmat)
            CfHe_t = CfHe_t/(Rup_trip+CfHe_t)   !"C" factor for triplets
        else                  !include H cont. effect
            Doppler = 2.d0*k_B*Tmat/(m_H*not4*C*C)
            Doppler = C*L_He_2Pt*dsqrt(Doppler)
            gamma_2Pt = 3.d0*A2P_t*fHe*(1.d0-x_He)*C*C &
            /(dsqrt(Pi)*sigma_He_2Pt*8.d0*Pi*Doppler*(1.d0-x_H)) &
            /((C*L_He_2Pt)**2.d0)
    !   use the fitting parameters from KIV (2007) in this case
            pb = 0.66d0
            qb = 0.9d0
            AHcon = A2P_t/(1.d0+pb*gamma_2Pt**qb)/3.d0
            CfHe_t = (A2P_t*pHe_t+AHcon)*dexp(-CL_PSt/Tmat)
            CfHe_t = CfHe_t/(Rup_trip+CfHe_t)   !"C" factor for triplets
        end if
      end if
    end if
        
        
!       Estimates of Thomson scattering time and Hubble time
        timeTh=(1._dl/(CT*Trad**4))*(1._dl+x+fHe)/x       !Thomson time
        timeH=2./(3.*HO*(1._dl+z)**1.5)      !Hubble time

!       calculate the derivatives
!       turn on H only for x_H<0.99, and use Saha derivative for 0.98<x_H<0.99
!       (clunky, but seems to work)
        if (x_H > 0.99) then   !don't change at all
                f(1) = 0._dl
!!        else if (x_H > 0.98_dl) then
      else if (x_H.gt.0.985d0) then     !use Saha rate for Hydrogen
            f(1) = (x*x_H*n*Rdown - Rup*(1.d0-x_H)*dexp(-CL/Tmat)) /(Hz*(1.d0+z))
            recombination_saha_z = z  
!AL: following commented as not used
!   for interest, calculate the correction factor compared to Saha
!   (without the fudge)
!       factor=(1.d0 + K*Lambda*n*(1.d0-x_H))
!       /(Hz*(1.d0+z)*(1.d0+K*Lambda*n*(1.d0-x)
!       +K*Rup*n*(1.d0-x)))       
      else !use full rate for H

        f(1) = ((x*x_H*n*Rdown - Rup*(1.d0-x_H)*exp(-CL/Tmat)) &
                *(1.d0 + K*Lambda*n*(1.d0-x_H))) &
                /(Hz*(1.d0+z)*(1.d0/fu+K*Lambda*n*(1.d0-x_H)/fu &
                +K*Rup*n*(1.d0-x_H)))

      end if
   
!       turn off the He once it is small
      if (x_He < 1.e-15) then 
                f(2)=0.d0
      else

        f(2) = ((x*x_He*n*Rdown_He &
            - Rup_He*(1-x_He)*exp(-CL_He/Tmat)) &
                *(1 + K_He*Lambda_He*n_He*(1.d0-x_He)*He_Boltz)) &
                /(Hz*(1+z) &
                * (1 + K_He*(Lambda_He+Rup_He)*n_He*(1.d0-x_He)*He_Boltz))
                
!   Modification to HeI recombination including channel via triplets
          if (Heflag.ge.3) then
            f(2) = f(2)+ (x*x_He*n*Rdown_trip & 
             - (1.d0-x_He)*3.d0*Rup_trip*dexp(-h_P*C*L_He_2st/(k_B*Tmat))) &
             *CfHe_t/(Hz*(1.d0+z))
          end if

        end if

        if (timeTh < H_frac*timeH) then
!                f(3)=Tmat/(1._dl+z)      !Tmat follows Trad
!	additional term to smooth transition to Tmat evolution,
!	(suggested by Adam Moss)
!model params used here added 9.19/2012, omega r put i in explicitly , omega m separated from omega axion, and spline table introduced for omegaaxion



!OLD RECFAST COMPUTATION OF DHDZ
  !       dHdz = (HO**2/2.d0/Hz)*(4.d0*((1.d0+z)**3)*omegar &
   !      + 3.d0*OmegaT*(1.d0+z)**2 + 2.d0*OmegaK*(1.d0+z) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Ours here includes omega_k computed self consistently including radiation neutrinos etc
!!!! as well as of course ULAS
sfac=1.0d0/(1.0d0+z)
	if (sfac .lt. aosc) then
        call spline_out(loga_table,grhoax_table,grhoax_table_buff,ntable,dlog10(sfac),gr)
        if (gr .eq. gr) then
           !delog it and multiply by physical units
           dorp=(1.0d1**(dble(gr)))
        else
           dorp=0.0d0
        endif
        else
        	dorp=drefp_hsq*((aosc/sfac)**3.0d0)
        endif
        !this is grhoax_table_internal=rho_ax/rho_crit^today so already 
        !in the form needed for recfast's normalization conventions
        !(unlike equations_ppf.f90, no multiplication by grhom needed)!        
        !Obtain H(z)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Obtain H(z) at slightly offset time period (lower redshift)
!

deriv_eps=1.d-3*real(sfac)
call spline_out(loga_table,grhoax_table,grhoax_table_buff,ntable,dlog10(sfac+deriv_eps),gr)
!
	if ((sfac+deriv_eps) .lt. aosc) then
        call spline_out(loga_table,grhoax_table,grhoax_table_buff,ntable,dlog10(sfac+deriv_eps),gr)
!        !write(*,*) a, gr
!        !compute log10 of density use of interpolation
!        !print*,a,grhoc,grhom*(1.0d1**(dble(gr)))/((grhoc+grhob)/(a**3.0d0))
        if (gr .eq. gr) then
           !delog it and multiply by physical units
           dorpa=(1.0d1**(dble(gr)))
        else
           dorpa=0.0d0
        endif
        else
        	dorpa=drefp_hsq*((aosc/(sfac+deriv_eps))**3.0d0)
        endif
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!simple derivative of axion energy density
 dorpa=(dorpa-dorp)/(deriv_eps)
!above calculate derivative of dimensionless axion density
 dHdz = (HO**2/2.d0/Hz)*(4.d0*((1.d0+z)**3)*omegar &
         + 3.d0*OmegaT*(1.d0+z)**2 + 2.d0*OmegaK*(1.d0+z)&
      &-dorpa/((1.0d0+z)**2.0d0))
!!!!! Include self consistently in dH/dz
      

        epsilon = Hz*(1.d0+x+fHe)/(CT*Trad**3*x)
        f(3) = Tnow &
        + epsilon*((1.d0+fHe)/(1.d0+fHe+x))*((f(1)+fHe*f(2))/x) &
        - epsilon* dHdz/Hz + 3.0d0*epsilon/(1.d0+z) 
                
        else
                f(3)= CT * (Trad**4) * x / (1._dl+x+fHe) &
                        * (Tmat-Trad) / (Hz*(1._dl+z)) + 2._dl*Tmat/(1._dl+z)    
        end if

         ! print *, z, f(3)*(1+z)/Tmat
         
        if (Do21cm .and. evolve_Ts) then

    !       follow the matter temperature once it has a chance of diverging
            if (timeTh < H_frac*timeH) then
                f(4) = Tnow !spin follows Trad and Tmat
            else
                if (z< 1/Do21cm_minev-1) then
       
                 Tspin = y(4) 
                 C10 = n*(kappa_HH_21cm(Tmat,.false.)*(1-x_H) + kappa_eH_21cm(Tmat,.false.)*x)
       
                 f(4) = 4*Tspin/Hz/(1+z)*( (Tspin/Tmat-1._dl)*C10 + Trad/T_21cm*(Tspin/Trad-1._dl)*A10) - f(1)*Tspin/(1-x_H)
                else
                 f(4)=f(3)
                end if
            end if
      
        end if

        end subroutine ION



        function dDeltaxe_dtau(a, Delta_xe,Delta_nH, Delta_Tm, hdot, kvb)
        !d x_e/d tau assuming Helium all neutral and temperature perturbations negligible
        !it is not accurate for x_e of order 1
        use RECDATA
        implicit none
        real(dl) dDeltaxe_dtau
        real(dl), intent(in):: a, Delta_xe,Delta_nH, Delta_Tm, hdot, kvb
        real(dl) Delta_Tg
        real(dl) xedot,z,x,n,n_He,Trad,Tmat,x_H,Hz, C_r, dlnC_r
        real(dl) Rup,Rdown,K
        real(dl) a_PPB,b_PPB,c_PPB,d_PPB
        real(dl) delta_alpha, delta_beta, delta_K, clh
        real(dl) dtauda
        external dtauda

       
        Delta_tg =Delta_Tm
        x_H = min(1._dl,Recombination_xe(a))

!       the Pequignot, Petitjean & Boisson fitting parameters for Hydrogen    
        a_PPB = 4.309d0
        b_PPB = -0.6166d0
        c_PPB = 0.6703d0
        d_PPB = 0.5300d0
 
        z=1/a-1

        x = x_H 

        n = Nnow /a**3
        n_He = fHe * n
        Trad = Tnow /a
        clh = 1/dtauda(a)/a !conformal time
        Hz = clh/a/MPC_in_sec !normal time in seconds

        Tmat = Recombination_tm(a)

!       Get the radiative rates using PPQ fit, identical to Hummer's table
        
        Rdown=1.d-19*a_PPB*(Tmat/1.d4)**b_PPB &
                /(1._dl+c_PPB*(Tmat/1.d4)**d_PPB)   !alpha
        Rup = Rdown * (CR*Tmat)**(1.5d0)*exp(-CDB/Tmat)
      
        K = CK/Hz              !Peebles coefficient K=lambda_a^3/8piH


        Rdown = Rdown*fu
        Rup = Rup*fu
        C_r =  a*(1.d0 + K*Lambda*n*(1.d0-x_H)) /( 1.d0+K*(Lambda+Rup)*n*(1.d0-x_H) )*MPC_in_sec
 
        xedot = -(x*x_H*n*Rdown - Rup*(1.d0-x_H)*exp(-CL/Tmat))*C_r
               
        delta_alpha = (b_PPB + c_PPB*(Tmat/1d4)**d_PPB*(b_PPB-d_PPB))/(1+c_PPB*(Tmat/1d4)**d_PPB)*Delta_Tg
        delta_beta = delta_alpha + (3./2 + CDB/Tmat)*delta_Tg !(Rup = beta)
        delta_K = - hdot/clh - kvb/clh/3


        dlnC_r = -Rup*K*n*( (Delta_nH+Delta_K + Delta_beta*(1+K*Lambda*n*(1-x_H)))*(1-x_H) - x_H*Delta_xe) &
          / ( 1.d0+K*(Lambda+Rup)*n*(1.d0-x_H) ) /(1.d0 + K*Lambda*n*(1.d0-x_H)) 
             
        dDeltaxe_dtau= xedot/x_H*(dlnC_r +Delta_alpha - Delta_xe) &
         - C_r*( (2*Delta_xe + Delta_nH)*x_H*n*Rdown + (Delta_xe - (3./2+ CB1/Tmat)*(1/x_H-1)*Delta_Tg)*Rup*exp(-CL/Tmat))
        
 
!Approximate form valid at late times
!        dDeltaxe_dtau= xedot/x_H*(Delta_alpha + Delta_xe + Delta_nH)


        end function dDeltaxe_dtau

!  ===============================================================


          function polevl(x,coef,N)
          implicit none
          integer N
          real(dl) polevl
          real(dl) x,ans
          real(dl) coef(N+1)

          integer i

          ans=coef(1)
          do i=2,N+1
             ans=ans*x+coef(i)
          end do
          polevl=ans
      
          end function polevl


          function derivpolevl(x,coef,N)
          implicit none
          integer N
          real(dl) derivpolevl
          real(dl) x,ans
          real(dl) coef(N+1)
          integer i

          ans=coef(1)*N
          do i=2,N
             ans=ans*x+coef(i)*(N-i+1)
          end do
          derivpolevl=ans
      
          end function derivpolevl


        function kappa_HH_21cm(T, deriv)
        !Polynomail fit to Hydrogen-Hydrogen collision rate as function of Tmatter, from astro-ph/0608032
        !if deriv return d log kappa / d log T
         real(dl), intent(in) :: T
         logical, intent(in) :: deriv
 !        real(dl), dimension(8), parameter :: fit = &
 !         (/ 0.00120402_dl, -0.0322247_dl,0.339581_dl, -1.75094_dl,4.3528_dl,-4.03562_dl, 1.26899_dl, -29.6113_dl /)
         integer, parameter :: n_table = 27
         integer, dimension(n_table), parameter :: Temps = &
          (/ 1, 2, 4, 6,8,10,15,20,25,30,40,50,60,70,80,90,100,200,300,500,700,1000,2000,3000,5000,7000,10000/) 
         real, dimension(n_table), parameter :: rates = &
          (/ 1.38e-13, 1.43e-13,2.71e-13, 6.60e-13,1.47e-12,2.88e-12,9.10e-12,1.78e-11,2.73e-11,&
           3.67e-11,5.38e-11,6.86e-11,8.14e-11,9.25e-11, &
           1.02e-10,1.11e-10,1.19e-10,1.75e-10,2.09e-10,2.56e-10,2.91e-10,3.31e-10,4.27e-10,&
           4.97e-10,6.03e-10,6.87e-10,7.87e-10/) 
          
         real(dl) kappa_HH_21cm, logT, logRate
         real(dl), save, dimension(:), allocatable :: logRates, logTemps, ddlogRates
         integer xlo, xhi
         real(dl) :: a0, b0, ho

         if (.not. allocated(logRates)) then

             allocate(logRates(n_table),logTemps(n_table),ddlogRates(n_table))
             logRates = log(real(rates,dl)*0.01**3)
             logTemps = log(real(Temps,dl))
             call spline(logTemps,logRates,n_table,1d30,1d30,ddlogRates)
         end if

         if (T<=Temps(1)) then
             if (deriv) then
              kappa_HH_21cm = 0
             else
              kappa_HH_21cm = rates(1)*0.01**3
             end if
            return
         elseif (T >=Temps(n_table)) then
             if (deriv) then
              kappa_HH_21cm = 0
             else
              kappa_HH_21cm = rates(n_table)*0.01**3
             end if
             return
         end if
         
         logT = log(T)
         xlo=0
         do xhi=2, n_table
          if (logT < logTemps(xhi)) then
            xlo = xhi-1
            exit
          end  if
         end do
         xhi = xlo+1
 
         ho=logTemps(xhi)-logTemps(xlo) 
         a0=(logTemps(xhi)-logT)/ho
         b0=1-a0
 
         if (deriv) then
          kappa_HH_21cm  = (logRates(xhi) - logRates(xlo))/ho + &
              ( ddlogRates(xhi)*(3*b0**2-1) - ddlogRates(xlo)*(3*a0**2-1))*ho/6
!          kappa_HH_21cm = derivpolevl(logT,fit,7)     
         else
          logRate = a0*logRates(xlo)+ b0*logRates(xhi)+ ((a0**3-a0)* ddlogRates(xlo) +(b0**3-b0)*ddlogRates(xhi))*ho**2/6
          kappa_HH_21cm = exp(logRate)
!          kappa_HH_21cm = exp(polevl(logT,fit,7))*0.01**3          
  
         end if

        end function kappa_HH_21cm


        function kappa_eH_21cm(T, deriv)
        !Polynomail fit to electron-Hydrogen collision rate as function of Tmatter; from astro-ph/0608032
        !if deriv return d log kappa / d log T
        ! from astro-ph/0608032
        !    1 2.39e-10
        !    2 3.37e-10
        !    5 5.3e-10
        !    10 7.46e-10
        !    20 1.05e-9
        !    50 1.63e-9
        !    100 2.26e-9
        !    200 3.11e-9
        !    500 4.59e-9
        !    1000 5.92e-9
        !    2000 7.15e-9
        !    5000 8.17e-9
        !    10000 8.37e-9
        !    15000 8.29e-9
        !    20000 8.11e-9
         real(dl), intent(in) :: T
         logical, intent(in) :: deriv
         real(dl), dimension(6), parameter :: fit = &
          (/5.86236d-005,  -0.00171375_dl, 0.0137303_dl, -0.0435277_dl, 0.540905_dl,-22.1596_dl /)
       
         real(dl) kappa_eH_21cm, logT

         logT = log(T)
         if (deriv) then
          kappa_eH_21cm = derivpolevl(logT,fit,5)      
         else
          kappa_eH_21cm = exp(polevl(logT,fit,5))*0.01**3          
         end if

        end function kappa_eH_21cm




        function kappa_pH_21cm(T, deriv) ! from astro-ph/0702487
        !Not actually used
        !Polynomail fit to proton-Hydrogen collision rate as function of Tmatter
        !if deriv return d log kappa / d log T
         real(dl), intent(in) :: T
         logical, intent(in) :: deriv
         integer, parameter :: n_table = 17
         integer, dimension(n_table), parameter :: Temps = &
          (/ 1, 2, 5, 10,20,50,100,200,500,1000,2000,3000,5000,7000,10000,15000,20000/) 
         real, dimension(n_table), parameter :: rates = &
          (/ 0.4028, 0.4517,0.4301,0.3699,0.3172,0.3047, 0.3379, 0.4043, 0.5471, 0.7051, 0.9167, 1.070, &
              1.301, 1.48,1.695,1.975,2.201/) 
          
         real(dl) kappa_pH_21cm, logT, logRate
         real(dl), save, dimension(:), allocatable :: logRates, logTemps, ddlogRates
         integer xlo, xhi
         real(dl) :: a0, b0, ho
         real(dl):: factor = 0.01**3*1e-9

         if (.not. allocated(logRates)) then

             allocate(logRates(n_table),logTemps(n_table),ddlogRates(n_table))
             logRates = log(real(rates,dl)*factor)
             logTemps = log(real(Temps,dl))
             call spline(logTemps,logRates,n_table,1d30,1d30,ddlogRates)
         end if

         if (T<=Temps(1)) then
             if (deriv) then
              kappa_pH_21cm = 0
             else
              kappa_pH_21cm = rates(1)*factor
             end if
            return
         elseif (T >=Temps(n_table)) then
             if (deriv) then
              kappa_pH_21cm = 0
             else
              kappa_pH_21cm = rates(n_table)*factor
             end if
             return
         end if
         
         logT = log(T)
         xlo=0
         do xhi=2, n_table
          if (logT < logTemps(xhi)) then
            xlo = xhi-1
            exit
          end  if
         end do
         xhi = xlo+1
 
         ho=logTemps(xhi)-logTemps(xlo) 
         a0=(logTemps(xhi)-logT)/ho
         b0=1-a0
 
         if (deriv) then
          kappa_pH_21cm  = (logRates(xhi) - logRates(xlo))/ho + &
           ( ddlogRates(xhi)*(3*b0**2-1) - ddlogRates(xlo)*(3*a0**2-1))*ho/6
         else
          logRate = a0*logRates(xlo)+ b0*logRates(xhi)+ ((a0**3-a0)* ddlogRates(xlo) +(b0**3-b0)*ddlogRates(xhi))*ho**2/6
          kappa_pH_21cm = exp(logRate)
         end if

        end function kappa_pH_21cm
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !recombination tailored verk that properly passes around splines
        !put in to avoid a a make catastrophe, variables as defined earlier
        !use recfast's integrator instead of my own for minimal catastrophe
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
      subroutine recdverk (EV,OmegaAx,omegar,ntable,aeq,aosc,drefp_hsq,loga_table,grhoax_table,&
      grhoax_table_buff,n, fcn, x, y, xend, tol, ind, c, nw, w)
      use Precision
      use AMLUtils
      integer n, ind, nw, k,ntable
      real(dl) x, y(n), xend, tol, c(*), w(nw,9), temp
  !    real EV !It isn't, but as long as it maintains it as a pointer we are OK
                Type(RecombinationParams) ::EV
 real(dl) loga_table(ntable),grhoax_table(ntable),grhoax_table_buff(ntable)
 real(dl) OmegaAx,omegar,aeq,aosc,drefp_hsq
      
!
!***********************************************************************
!                                                                      *
! note added 11/14/85.                                                 *
!                                                                      *
! if you discover any errors in this subroutine, please contact        *
!                                                                      *
!        kenneth r. jackson                                            *
!        department of computer science                                *
!        university of toronto                                         *
!        toronto, ontario,                                             *
!        canada   m5s 1a4                                              *
!                                                                      *
!        phone: 416-978-7075                                           *
!                                                                      *
!        electronic mail:                                              *
!        uucp:   {cornell,decvax,ihnp4,linus,uw-beaver}!utcsri!krj     *
!        csnet:  krj@toronto                                           *
!        arpa:   krj.toronto@csnet-relay                               *
!        bitnet: krj%toronto@csnet-relay.arpa                          *
!                                                                      *
! dverk is written in fortran 66.                                      *
!                                                                      *
! the constants dwarf and rreb -- c(10) and c(11), respectively -- are *
! set for a  vax  in  double  precision.  they  should  be  reset,  as *
! described below, if this program is run on another machine.          *
!                                                                      *
! the c array is declared in this subroutine to have one element only, *
! although  more  elements  are  referenced  in this subroutine.  this *
! causes some compilers to issue warning messages.  there is,  though, *
! no  error  provided  c is declared sufficiently large in the calling *
! program, as described below.                                         *
!                                                                      *
! the following external statement  for  fcn  was  added  to  avoid  a *
! warning  message  from  the  unix  f77 compiler.  the original dverk *
! comments and code follow it.                                         *
!                                                                      *
!***********************************************************************
!
      external fcn
!
!***********************************************************************
!                                                                      *
!     purpose - this is a runge-kutta  subroutine  based  on  verner's *
! fifth and sixth order pair of formulas for finding approximations to *
! the solution of  a  system  of  first  order  ordinary  differential *
! equations  with  initial  conditions. it attempts to keep the global *
! error proportional to  a  tolerance  specified  by  the  user.  (the *
! proportionality  depends  on the kind of error control that is used, *
! as well as the differential equation and the range of integration.)  *
!                                                                      *
!     various options are available to the user,  including  different *
! kinds  of  error control, restrictions on step sizes, and interrupts *
! which permit the user to examine the state of the  calculation  (and *
! perhaps make modifications) during intermediate stages.              *
!                                                                      *
!     the program is efficient for non-stiff systems.  however, a good *
! variable-order-adams  method  will probably be more efficient if the *
! function evaluations are very costly.  such a method would  also  be *
! more suitable if one wanted to obtain a large number of intermediate *
! solution values by interpolation, as might be the case  for  example *
! with graphical output.                                               *
!                                                                      *
!                                    hull-enright-jackson   1/10/76    *
!                                                                      *
!***********************************************************************
!                                                                      *
!     use - the user must specify each of the following                *
!                                                                      *
!     n  number of equations                                           *
!                                                                      *
!   fcn  name of subroutine for evaluating functions - the  subroutine *
!           itself must also be provided by the user - it should be of *
!           the following form                                         *
!              subroutine fcn(n, x, y, yprime)                         *
!              integer n                                               *
!              real(dl) x, y(n), yprime(n)                     *
!                      *** etc ***                                     *
!           and it should evaluate yprime, given n, x and y            *
!                                                                      *
!     x  independent variable - initial value supplied by user         *
!                                                                      *
!     y  dependent variable - initial values of components y(1), y(2), *
!           ..., y(n) supplied by user                                 *
!                                                                      *
!  xend  value of x to which integration is to be carried out - it may *
!           be less than the initial value of x                        *
!                                                                      *
!   tol  tolerance - the subroutine attempts to control a norm of  the *
!           local  error  in  such  a  way  that  the  global error is *
!           proportional to tol. in some problems there will be enough *
!           damping  of  errors, as well as some cancellation, so that *
!           the global error will be less than tol. alternatively, the *
!           control   can   be  viewed  as  attempting  to  provide  a *
!           calculated value of y at xend which is the exact  solution *
!           to  the  problem y' = f(x,y) + e(x) where the norm of e(x) *
!           is proportional to tol.  (the norm  is  a  max  norm  with *
!           weights  that  depend on the error control strategy chosen *
!           by the user.  the default weight for the k-th component is *
!           1/max(1,abs(y(k))),  which therefore provides a mixture of *
!           absolute and relative error control.)                      *
!                                                                      *
!   ind  indicator - on initial entry ind must be set equal to  either *
!           1  or  2. if the user does not wish to use any options, he *
!           should set ind to 1 - all that remains for the user to  do *
!           then  is  to  declare c and w, and to specify nw. the user *
!           may also  select  various  options  on  initial  entry  by *
!           setting ind = 2 and initializing the first 9 components of *
!           c as described in the next section.  he may also  re-enter *
!           the  subroutine  with ind = 3 as mentioned again below. in *
!           any event, the subroutine returns with ind equal to        *
!              3 after a normal return                                 *
!              4, 5, or 6 after an interrupt (see options c(8), c(9))  *
!              -1, -2, or -3 after an error condition (see below)      *
!                                                                      *
!     c  communications vector - the dimension must be greater than or *
!           equal to 24, unless option c(1) = 4 or 5 is used, in which *
!           case the dimension must be greater than or equal to n+30   *
!                                                                      *
!    nw  first dimension of workspace w -  must  be  greater  than  or *
!           equal to n                                                 *
!                                                                      *
!     w  workspace matrix - first dimension must be nw and second must *
!           be greater than or equal to 9                              *
!                                                                      *
!     the subroutine  will  normally  return  with  ind  =  3,  having *
! replaced the initial values of x and y with, respectively, the value *
! of xend and an approximation to y at xend.  the  subroutine  can  be *
! called  repeatedly  with new values of xend without having to change *
! any other argument.  however, changes in tol, or any of the  options *
! described below, may also be made on such a re-entry if desired.     *
!                                                                      *
!     three error returns are also possible, in which  case  x  and  y *
! will be the most recently accepted values -                          *
!     with ind = -3 the subroutine was unable  to  satisfy  the  error *
!        requirement  with a particular step-size that is less than or *
!        equal to hmin, which may mean that tol is too small           *
!     with ind = -2 the value of hmin  is  greater  than  hmax,  which *
!        probably  means  that the requested tol (which is used in the *
!        calculation of hmin) is too small                             *
!     with ind = -1 the allowed maximum number of fcn evaluations  has *
!        been  exceeded,  but  this  can only occur if option c(7), as *
!        described in the next section, has been used                  *
!                                                                      *
!     there are several circumstances that will cause the calculations *
! to  be  terminated,  along with output of information that will help *
! the user determine the cause of  the  trouble.  these  circumstances *
! involve  entry with illegal or inconsistent values of the arguments, *
! such as attempting a normal  re-entry  without  first  changing  the *
! value of xend, or attempting to re-enter with ind less than zero.    *
!                                                                      *
!***********************************************************************
!                                                                      *
!     options - if the subroutine is entered with ind = 1, the first 9 *
! components of the communications vector are initialized to zero, and *
! the subroutine uses only default values  for  each  option.  if  the *
! subroutine  is  entered  with ind = 2, the user must specify each of *
! these 9 components - normally he would first set them all  to  zero, *
! and  then  make  non-zero  those  that  correspond to the particular *
! options he wishes to select. in any event, options may be changed on *
! re-entry  to  the  subroutine  -  but if the user changes any of the *
! options, or tol, in the course of a calculation he should be careful *
! about  how  such changes affect the subroutine - it may be better to *
! restart with ind = 1 or 2. (components 10 to 24 of c are used by the *
! program  -  the information is available to the user, but should not *
! normally be changed by him.)                                         *
!                                                                      *
!  c(1)  error control indicator - the norm of the local error is  the *
!           max  norm  of  the  weighted  error  estimate  vector, the *
!           weights being determined according to the value of c(1) -  *
!              if c(1)=1 the weights are 1 (absolute error control)    *
!              if c(1)=2 the weights are 1/abs(y(k))  (relative  error *
!                 control)                                             *
!              if c(1)=3 the  weights  are  1/max(abs(c(2)),abs(y(k))) *
!                 (relative  error  control,  unless abs(y(k)) is less *
!                 than the floor value, abs(c(2)) )                    *
!              if c(1)=4 the weights are 1/max(abs(c(k+30)),abs(y(k))) *
!                 (here individual floor values are used)              *
!              if c(1)=5 the weights are 1/abs(c(k+30))                *
!              for all other values of c(1), including  c(1) = 0,  the *
!                 default  values  of  the  weights  are  taken  to be *
!                 1/max(1,abs(y(k))), as mentioned earlier             *
!           (in the two cases c(1) = 4 or 5 the user must declare  the *
!           dimension of c to be at least n+30 and must initialize the *
!           components c(31), c(32), ..., c(n+30).)                    *
!                                                                      *
!  c(2)  floor value - used when the indicator c(1) has the value 3    *
!                                                                      *
!  c(3)  hmin specification - if not zero, the subroutine chooses hmin *
!           to be abs(c(3)) - otherwise it uses the default value      *
!              10*max(dwarf,rreb*max(weighted norm y/tol,abs(x))),     *
!           where dwarf is a very small positive  machine  number  and *
!           rreb is the relative roundoff error bound                  *
!                                                                      *
!  c(4)  hstart specification - if not zero, the subroutine  will  use *
!           an  initial  hmag equal to abs(c(4)), except of course for *
!           the restrictions imposed by hmin and hmax  -  otherwise it *
!           uses the default value of hmax*(tol)**(1/6)                *
!                                                                      *
!  c(5)  scale specification - this is intended to be a measure of the *
!           scale of the problem - larger values of scale tend to make *
!           the method more reliable, first  by  possibly  restricting *
!           hmax  (as  described  below) and second, by tightening the *
!           acceptance requirement - if c(5) is zero, a default  value *
!           of  1  is  used.  for  linear  homogeneous  problems  with *
!           constant coefficients, an appropriate value for scale is a *
!           norm  of  the  associated  matrix.  for other problems, an *
!           approximation to  an  average  value  of  a  norm  of  the *
!           jacobian along the trajectory may be appropriate           *
!                                                                      *
!  c(6)  hmax specification - four cases are possible                  *
!           if c(6).ne.0 and c(5).ne.0, hmax is taken to be            *
!              min(abs(c(6)),2/abs(c(5)))                              *
!           if c(6).ne.0 and c(5).eq.0, hmax is taken to be  abs(c(6)) *
!           if c(6).eq.0 and c(5).ne.0, hmax is taken to be            *
!              2/abs(c(5))                                             *
!           if c(6).eq.0 and c(5).eq.0, hmax is given a default  value *
!              of 2                                                    *
!                                                                      *
!  c(7)  maximum number of function evaluations  -  if  not  zero,  an *
!           error  return with ind = -1 will be caused when the number *
!           of function evaluations exceeds abs(c(7))                  *
!                                                                      *
!  c(8)  interrupt number  1  -  if  not  zero,  the  subroutine  will *
!           interrupt   the  calculations  after  it  has  chosen  its *
!           preliminary value of hmag, and just before choosing htrial *
!           and  xtrial  in  preparation for taking a step (htrial may *
!           differ from hmag in sign, and may  require  adjustment  if *
!           xend  is  near) - the subroutine returns with ind = 4, and *
!           will resume calculation at the point  of  interruption  if *
!           re-entered with ind = 4                                    *
!                                                                      *
!  c(9)  interrupt number  2  -  if  not  zero,  the  subroutine  will *
!           interrupt   the  calculations  immediately  after  it  has *
!           decided whether or not to accept the result  of  the  most *
!           recent  trial step, with ind = 5 if it plans to accept, or *
!           ind = 6 if it plans to reject -  y(*)  is  the  previously *
!           accepted  result, while w(*,9) is the newly computed trial *
!           value, and w(*,2) is the unweighted error estimate vector. *
!           the  subroutine  will  resume calculations at the point of *
!           interruption on re-entry with ind = 5 or 6. (the user  may *
!           change ind in this case if he wishes, for example to force *
!           acceptance of a step that would otherwise be rejected,  or *
!           vice versa. he can also restart with ind = 1 or 2.)        *
!                                                                      *
!***********************************************************************
!                                                                      *
!  summary of the components of the communications vector              *
!                                                                      *
!     prescribed at the option       determined by the program         *
!           of the user                                                *
!                                                                      *
!                                    c(10) rreb(rel roundoff err bnd)  *
!     c(1) error control indicator   c(11) dwarf (very small mach no)  *
!     c(2) floor value               c(12) weighted norm y             *
!     c(3) hmin specification        c(13) hmin                        *
!     c(4) hstart specification      c(14) hmag                        *
!     c(5) scale specification       c(15) scale                       *
!     c(6) hmax specification        c(16) hmax                        *
!     c(7) max no of fcn evals       c(17) xtrial                      *
!     c(8) interrupt no 1            c(18) htrial                      *
!     c(9) interrupt no 2            c(19) est                         *
!                                    c(20) previous xend               *
!                                    c(21) flag for xend               *
!                                    c(22) no of successful steps      *
!                                    c(23) no of successive failures   *
!                                    c(24) no of fcn evals             *
!                                                                      *
!  if c(1) = 4 or 5, c(31), c(32), ... c(n+30) are floor values        *
!                                                                      *
!***********************************************************************
!                                                                      *
!  an overview of the program                                          *
!                                                                      *
!     begin initialization, parameter checking, interrupt re-entries   *
!  ......abort if ind out of range 1 to 6                              *
!  .     cases - initial entry, normal re-entry, interrupt re-entries  *
!  .     case 1 - initial entry (ind .eq. 1 or 2)                      *
!  v........abort if n.gt.nw or tol.le.0                               *
!  .        if initial entry without options (ind .eq. 1)              *
!  .           set c(1) to c(9) equal to zero                          *
!  .        else initial entry with options (ind .eq. 2)               *
!  .           make c(1) to c(9) non-negative                          *
!  .           make floor values non-negative if they are to be used   *
!  .        end if                                                     *
!  .        initialize rreb, dwarf, prev xend, flag, counts            *
!  .     case 2 - normal re-entry (ind .eq. 3)                         *
!  .........abort if xend reached, and either x changed or xend not    *
!  .        re-initialize flag                                         *
!  .     case 3 - re-entry following an interrupt (ind .eq. 4 to 6)    *
!  v        transfer control to the appropriate re-entry point.......  *
!  .     end cases                                                  .  *
!  .  end initialization, etc.                                      .  *
!  .                                                                v  *
!  .  loop through the following 4 stages, once for each trial step .  *
!  .     stage 1 - prepare                                          .  *
!***********error return (with ind=-1) if no of fcn evals too great .  *
!  .        calc slope (adding 1 to no of fcn evals) if ind .ne. 6  .  *
!  .        calc hmin, scale, hmax                                  .  *
!***********error return (with ind=-2) if hmin .gt. hmax            .  *
!  .        calc preliminary hmag                                   .  *
!***********interrupt no 1 (with ind=4) if requested.......re-entry.v  *
!  .        calc hmag, xtrial and htrial                            .  *
!  .     end stage 1                                                .  *
!  v     stage 2 - calc ytrial (adding 7 to no of fcn evals)        .  *
!  .     stage 3 - calc the error estimate                          .  *
!  .     stage 4 - make decisions                                   .  *
!  .        set ind=5 if step acceptable, else set ind=6            .  *
!***********interrupt no 2 if requested....................re-entry.v  *
!  .        if step accepted (ind .eq. 5)                              *
!  .           update x, y from xtrial, ytrial                         *
!  .           add 1 to no of successful steps                         *
!  .           set no of successive failures to zero                   *
!**************return(with ind=3, xend saved, flag set) if x .eq. xend *
!  .        else step not accepted (ind .eq. 6)                        *
!  .           add 1 to no of successive failures                      *
!**************error return (with ind=-3) if hmag .le. hmin            *
!  .        end if                                                     *
!  .     end stage 4                                                   *
!  .  end loop                                                         *
!  .                                                                   *
!  begin abort action                                                  *
!     output appropriate  message  about  stopping  the  calculations, *
!        along with values of ind, n, nw, tol, hmin,  hmax,  x,  xend, *
!        previous xend,  no of  successful  steps,  no  of  successive *
!        failures, no of fcn evals, and the components of y            *
!     stop                                                             *
!  end abort action                                                    *
!                                                                      *
!***********************************************************************
!
!     ******************************************************************
!     * begin initialization, parameter checking, interrupt re-entries *
!     ******************************************************************
!
!  ......abort if ind out of range 1 to 6
         if (ind.lt.1 .or. ind.gt.6) go to 500
!
!        cases - initial entry, normal re-entry, interrupt re-entries
!         go to (5, 5, 45, 1111, 2222, 2222), ind
         if (ind==3) goto 45
         if (ind==4) goto 1111
         if (ind==5 .or. ind==6) goto 2222

!        case 1 - initial entry (ind .eq. 1 or 2)
!  .........abort if n.gt.nw or tol.le.0
            if (n.gt.nw .or. tol.le.0._dl) go to 500
            if (ind.eq. 2) go to 15
!              initial entry without options (ind .eq. 1)
!              set c(1) to c(9) equal to 0
               do k = 1, 9
                  c(k) = 0._dl
               end do
               go to 35
   15       continue
!              initial entry with options (ind .eq. 2)
!              make c(1) to c(9) non-negative
               do k = 1, 9
                  c(k) = dabs(c(k))
               end do
!              make floor values non-negative if they are to be used
               if (c(1).ne.4._dl .and. c(1).ne.5._dl) go to 30
                  do k = 1, n
                     c(k+30) = dabs(c(k+30))
                  end do
   30          continue
   35       continue
!           initialize rreb, dwarf, prev xend, flag, counts
            c(10) = 2._dl**(-56)
            c(11) = 1.d-35
!           set previous xend initially to initial value of x
            c(20) = x
            do k = 21, 24
               c(k) = 0._dl
            end do
            go to 50
!        case 2 - normal re-entry (ind .eq. 3)
!  .........abort if xend reached, and either x changed or xend not
   45       if (c(21).ne.0._dl .and. &
                              (x.ne.c(20) .or. xend.eq.c(20))) go to 500
!           re-initialize flag
            c(21) = 0._dl
            go to 50
!        case 3 - re-entry following an interrupt (ind .eq. 4 to 6)
!           transfer control to the appropriate re-entry point..........
!           this has already been handled by the computed go to        .
!        end cases                                                     v
   50    continue
!
!     end initialization, etc.
!
!     ******************************************************************
!     * loop through the following 4 stages, once for each trial  step *
!     * until the occurrence of one of the following                   *
!     *    (a) the normal return (with ind .eq. 3) on reaching xend in *
!     *        stage 4                                                 *
!     *    (b) an error return (with ind .lt. 0) in stage 1 or stage 4 *
!     *    (c) an interrupt return (with ind  .eq.  4,  5  or  6),  if *
!     *        requested, in stage 1 or stage 4                        *
!     ******************************************************************
!
99999 continue
!
!        ***************************************************************
!        * stage 1 - prepare - do calculations of  hmin,  hmax,  etc., *
!        * and some parameter  checking,  and  end  up  with  suitable *
!        * values of hmag, xtrial and htrial in preparation for taking *
!        * an integration step.                                        *
!        ***************************************************************
!
!***********error return (with ind=-1) if no of fcn evals too great
            if (c(7).eq.0._dl .or. c(24).lt.c(7)) go to 100
               ind = -1
               return
  100       continue
!
!           calculate slope (adding 1 to no of fcn evals) if ind .ne. 6
            if (ind .eq. 6) go to 105
               call fcn(EV,OmegaAx,omegar,ntable,aeq,aosc,drefp_hsq,loga_table,grhoax_table,&
      grhoax_table_buff,n, x, y, w(1,1))
               c(24) = c(24) + 1._dl
  105       continue
!
!           calculate hmin - use default unless value prescribed
            c(13) = c(3)
            if (c(3) .ne. 0._dl) go to 165
!              calculate default value of hmin
!              first calculate weighted norm y - c(12) - as specified
!              by the error control indicator c(1)
               temp = 0._dl
               if (c(1) .ne. 1._dl) go to 115
!                 absolute error control - weights are 1
                  do 110 k = 1, n
                     temp = dmax1(temp, dabs(y(k)))
  110             continue
                  c(12) = temp
                  go to 160
  115          if (c(1) .ne. 2._dl) go to 120
!                 relative error control - weights are 1/dabs(y(k)) so
!                 weighted norm y is 1
                  c(12) = 1._dl
                  go to 160
  120          if (c(1) .ne. 3._dl) go to 130
!                 weights are 1/max(c(2),abs(y(k)))
                  do 125 k = 1, n
                     temp = dmax1(temp, dabs(y(k))/c(2))
  125             continue
                  c(12) = dmin1(temp, 1._dl)
                  go to 160
  130          if (c(1) .ne. 4._dl) go to 140
!                 weights are 1/max(c(k+30),abs(y(k)))
                  do 135 k = 1, n
                     temp = dmax1(temp, dabs(y(k))/c(k+30))
  135             continue
                  c(12) = dmin1(temp, 1._dl)
                  go to 160
  140          if (c(1) .ne. 5._dl) go to 150
!                 weights are 1/c(k+30)
                  do 145 k = 1, n
                     temp = dmax1(temp, dabs(y(k))/c(k+30))
  145             continue
                  c(12) = temp
                  go to 160
  150          continue
!                 default case - weights are 1/max(1,abs(y(k)))
                  do 155 k = 1, n
                     temp = dmax1(temp, dabs(y(k)))
  155             continue
                  c(12) = dmin1(temp, 1._dl)
  160          continue
               c(13) = 10._dl*dmax1(c(11),c(10)*dmax1(c(12)/tol,dabs(x)))
  165       continue
!
!           calculate scale - use default unless value prescribed
            c(15) = c(5)
            if (c(5) .eq. 0._dl) c(15) = 1._dl
!
!           calculate hmax - consider 4 cases
!           case 1 both hmax and scale prescribed
               if (c(6).ne.0._dl .and. c(5).ne.0._dl) &
                                          c(16) = dmin1(c(6), 2._dl/c(5))
!           case 2 - hmax prescribed, but scale not
               if (c(6).ne.0._dl .and. c(5).eq.0._dl) c(16) = c(6)
!           case 3 - hmax not prescribed, but scale is
               if (c(6).eq.0._dl .and. c(5).ne.0._dl) c(16) = 2._dl/c(5)
!           case 4 - neither hmax nor scale is provided
               if (c(6).eq.0._dl .and. c(5).eq.0._dl) c(16) = 2._dl
!
!***********error return (with ind=-2) if hmin .gt. hmax
            if (c(13) .le. c(16)) go to 170
               ind = -2
               return
  170       continue
!
!           calculate preliminary hmag - consider 3 cases
            if (ind .gt. 2) go to 175
!           case 1 - initial entry - use prescribed value of hstart, if
!              any, else default
               c(14) = c(4)
               if (c(4) .eq. 0._dl) c(14) = c(16)*tol**(1._dl/6._dl)
               go to 185
  175       if (c(23) .gt. 1._dl) go to 180
!           case 2 - after a successful step, or at most  one  failure,
!              use min(2, .9*(tol/est)**(1/6))*hmag, but avoid possible
!              overflow. then avoid reduction by more than half.
               temp = 2._dl*c(14)
               if (tol .lt. (2._dl/.9d0)**6*c(19)) &
                       temp = .9d0*(tol/c(19))**(1._dl/6._dl)*c(14)
               c(14) = dmax1(temp, .5d0*c(14))
               go to 185
  180       continue
!           case 3 - after two or more successive failures
               c(14) = .5d0*c(14)
  185       continue
!
!           check against hmax
            c(14) = dmin1(c(14), c(16))
!
!           check against hmin
            c(14) = dmax1(c(14), c(13))
!
!***********interrupt no 1 (with ind=4) if requested
            if (c(8) .eq. 0._dl) go to 1111
               ind = 4
               return
!           resume here on re-entry with ind .eq. 4   ........re-entry..
 1111       continue
!
!           calculate hmag, xtrial - depending on preliminary hmag, xend
            if (c(14) .ge. dabs(xend - x)) go to 190
!              do not step more than half way to xend
               c(14) = dmin1(c(14), .5d0*dabs(xend - x))
               c(17) = x + dsign(c(14), xend - x)
               go to 195
  190       continue
!              hit xend exactly
               c(14) = dabs(xend - x)
               c(17) = xend
  195       continue
!
!           calculate htrial
            c(18) = c(17) - x
!
!        end stage 1
!
!        ***************************************************************
!        * stage 2 - calculate ytrial (adding 7 to no of  fcn  evals). *
!        * w(*,2), ... w(*,8)  hold  intermediate  results  needed  in *
!        * stage 3. w(*,9) is temporary storage until finally it holds *
!        * ytrial.                                                     *
!        ***************************************************************
!
            temp = c(18)/1398169080000._dl
!
            do 200 k = 1, n
               w(k,9) = y(k) + temp*w(k,1)*233028180000._dl
  200       continue
            call fcn(EV,real(OmegaAx),omegar,ntable,aeq,aosc,drefp_hsq,loga_table,grhoax_table,&
      grhoax_table_buff,n, x + c(18)/6._dl, w(1,9), w(1,2))
!
            do 205 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*74569017600._dl &
                                      + w(k,2)*298276070400._dl  )
  205       continue
            call fcn(EV,real(OmegaAx),omegar,ntable,aeq,aosc,drefp_hsq,loga_table,grhoax_table,&
      grhoax_table_buff,n, x + c(18)*(4._dl/15._dl), w(1,9), w(1,3))
!
            do 210 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*1165140900000._dl &
                                      - w(k,2)*3728450880000._dl &
                                      + w(k,3)*3495422700000._dl )
  210       continue
            call fcn(EV,real(OmegaAx),omegar,ntable,aeq,aosc,drefp_hsq,loga_table,grhoax_table,&
      grhoax_table_buff,n, x + c(18)*(2._dl/3._dl), w(1,9), w(1,4))
!
            do 215 k = 1, n
               w(k,9) = y(k) + temp*( - w(k,1)*3604654659375._dl &
                                      + w(k,2)*12816549900000._dl &
                                      - w(k,3)*9284716546875._dl &
                                      + w(k,4)*1237962206250._dl )
  215       continue
            call fcn(EV,real(OmegaAx),omegar,ntable,aeq,aosc,drefp_hsq,loga_table,grhoax_table,&
      grhoax_table_buff,n, x + c(18)*(5._dl/6._dl), w(1,9), w(1,5))
!
            do 220 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*3355605792000._dl &
                                      - w(k,2)*11185352640000._dl &
                                      + w(k,3)*9172628850000._dl &
                                      - w(k,4)*427218330000._dl &
                                      + w(k,5)*482505408000._dl  )
  220       continue
            call fcn(EV,real(OmegaAx),omegar,ntable,aeq,aosc,drefp_hsq,loga_table,grhoax_table,&
      grhoax_table_buff,n, x + c(18), w(1,9), w(1,6))
!
            do 225 k = 1, n
               w(k,9) = y(k) + temp*( - w(k,1)*770204740536._dl &
                                      + w(k,2)*2311639545600._dl &
                                      - w(k,3)*1322092233000._dl &
                                      - w(k,4)*453006781920._dl &
                                      + w(k,5)*326875481856._dl  )
  225       continue
            call fcn(EV,real(OmegaAx),omegar,ntable,aeq,aosc,drefp_hsq,loga_table,grhoax_table,&
      grhoax_table_buff,n, x + c(18)/15._dl, w(1,9), w(1,7))
!
            do 230 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*2845924389000._dl &
                                      - w(k,2)*9754668000000._dl &
                                      + w(k,3)*7897110375000._dl &
                                      - w(k,4)*192082660000._dl &
                                      + w(k,5)*400298976000._dl &
                                      + w(k,7)*201586000000._dl  )
  230       continue
            call fcn(EV,real(OmegaAx),omegar,ntable,aeq,aosc,drefp_hsq,loga_table,grhoax_table,&
      grhoax_table_buff,n, x + c(18), w(1,9), w(1,8))
!
!           calculate ytrial, the extrapolated approximation and store
!              in w(*,9)
            do 235 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*104862681000._dl &
                                      + w(k,3)*545186250000._dl &
                                      + w(k,4)*446637345000._dl &
                                      + w(k,5)*188806464000._dl &
                                      + w(k,7)*15076875000._dl &
                                      + w(k,8)*97599465000._dl   )
  235       continue
!
!           add 7 to the no of fcn evals
            c(24) = c(24) + 7._dl
!
!        end stage 2
!
!        ***************************************************************
!        * stage 3 - calculate the error estimate est. first calculate *
!        * the  unweighted  absolute  error  estimate vector (per unit *
!        * step) for the unextrapolated approximation and store it  in *
!        * w(*,2).  then  calculate the weighted max norm of w(*,2) as *
!        * specified by the error  control  indicator  c(1).  finally, *
!        * modify  this result to produce est, the error estimate (per *
!        * unit step) for the extrapolated approximation ytrial.       *
!        ***************************************************************
!
!           calculate the unweighted absolute error estimate vector
            do 300 k = 1, n
               w(k,2) = (   w(k,1)*8738556750._dl &
                          + w(k,3)*9735468750._dl &
                          - w(k,4)*9709507500._dl &
                          + w(k,5)*8582112000._dl &
                          + w(k,6)*95329710000._dl &
                          - w(k,7)*15076875000._dl &
                          - w(k,8)*97599465000._dl)/1398169080000._dl
  300       continue
!
!           calculate the weighted max norm of w(*,2) as specified by
!           the error control indicator c(1)
            temp = 0._dl
            if (c(1) .ne. 1._dl) go to 310
!              absolute error control
               do 305 k = 1, n
                  temp = dmax1(temp,dabs(w(k,2)))
  305          continue
               go to 360
  310       if (c(1) .ne. 2._dl) go to 320
!              relative error control
               do 315 k = 1, n
                  temp = dmax1(temp, dabs(w(k,2)/y(k)))
  315          continue
               go to 360
  320       if (c(1) .ne. 3._dl) go to 330
!              weights are 1/max(c(2),abs(y(k)))
               do 325 k = 1, n
                  temp = dmax1(temp, dabs(w(k,2)) &
                                   / dmax1(c(2), dabs(y(k))) )
  325          continue
               go to 360
  330       if (c(1) .ne. 4._dl) go to 340
!              weights are 1/max(c(k+30),abs(y(k)))
               do 335 k = 1, n
                  temp = dmax1(temp, dabs(w(k,2)) &
                                   / dmax1(c(k+30), dabs(y(k))) )
  335          continue
               go to 360
  340       if (c(1) .ne. 5._dl) go to 350
!              weights are 1/c(k+30)
               do 345 k = 1, n
                  temp = dmax1(temp, dabs(w(k,2)/c(k+30)))
  345          continue
               go to 360
  350       continue
!              default case - weights are 1/max(1,abs(y(k)))
               do 355 k = 1, n
                  temp = dmax1(temp, dabs(w(k,2)) &
                                   / dmax1(1._dl, dabs(y(k))) )
  355          continue
  360       continue
!
!           calculate est - (the weighted max norm of w(*,2))*hmag*scale
!              - est is intended to be a measure of the error  per  unit
!              step in ytrial
            c(19) = temp*c(14)*c(15)
!
!        end stage 3
!
!        ***************************************************************
!        * stage 4 - make decisions.                                   *
!        ***************************************************************
!
!           set ind=5 if step acceptable, else set ind=6
            ind = 5
            if (c(19) .gt. tol) ind = 6
!
!***********interrupt no 2 if requested
            if (c(9) .eq. 0._dl) go to 2222
               return
!           resume here on re-entry with ind .eq. 5 or 6   ...re-entry..
 2222       continue
!
            if (ind .eq. 6) go to 410
!              step accepted (ind .eq. 5), so update x, y from xtrial,
!                 ytrial, add 1 to the no of successful steps, and set
!                 the no of successive failures to zero
               x = c(17)
               do 400 k = 1, n
                  y(k) = w(k,9)
  400          continue
               c(22) = c(22) + 1._dl
               c(23) = 0._dl
!**************return(with ind=3, xend saved, flag set) if x .eq. xend
               if (x .ne. xend) go to 405
                  ind = 3
                  c(20) = xend
                  c(21) = 1._dl
                  return
  405          continue
               go to 420
  410       continue
!              step not accepted (ind .eq. 6), so add 1 to the no of
!                 successive failures
               c(23) = c(23) + 1._dl
!**************error return (with ind=-3) if hmag .le. hmin
               if (c(14) .gt. c(13)) go to 415
                  ind = -3
                  return
  415          continue
  420       continue
!
!        end stage 4
!
      go to 99999
!     end loop
!
!  begin abort action
  500 continue
!

      write (*,*) 'Error in dverk, x =',x, 'xend=', xend !, 'not aborting'
      call MpiStop()
!
!  end abort action
!
      end subroutine recdverk
        end module Recombination

