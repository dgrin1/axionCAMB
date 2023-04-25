    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! The `halofit' code models the nonlinear evolution of cold matter
    ! cosmological power spectra. The full details of the way in which
    ! this is done are presented in Smith et al. (2002), MNRAS, ?, ?.
    !
    ! The code `halofit' was written by R. E. Smith & J. A. Peacock.
    ! See http://www.astro.upenn.edu/~res,
    !
    ! Subsequent updates as below
    ! Only tested for basic models with power law initial power spectra

    ! Adapted for F90 and CAMB, AL March 2005
    !!BR09 Oct 09: generalized expressions for om(z) and ol(z) to include w

    ! RT12 Oct: update some fitting parameters in the code to enhance
    !           the power spectrum at small scales (arXiv:1208.2701)

    !!JD 08/13: generalized expressions for om(z) and ol(z) to include
    !           w_0 and w_a
    ! SPB14 Feb: update the fitting parameters for neutrinos to work with RT12
    !           modifications
    ! AL Sept 14: added halofit_version parameter to change approximation used;
    !   separate halofit.f90 is no longer needed as equations.f90 defined fixed wa_ppf
    ! Jan 15: Suggested change from Simeon Bird to avoid issues with very large Omm and neutrinos

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module NonLinear
    use ModelParams
    use transfer
    use LambdaGeneral
    implicit none
    private

    real, parameter :: Min_kh_nonlinear = 0.005
    real(dl):: om_m,om_v,fnu,omm0, acur

    integer, parameter :: halofit_original = 1, halofit_bird=2, halofit_peacock=3, halofit_takahashi=4
    integer, parameter :: halofit_default = halofit_original ! DM15: Takahashi is not stable for axion models. Other versions agree well and are sensible to percent level for lensing \ell<4000.
    integer :: halofit_version = halofit_default
    public Min_kh_nonlinear,NonLinear_GetNonLinRatios, NonLinear_ReadParams
    public halofit_version,halofit_default, halofit_original, halofit_bird, halofit_peacock, halofit_takahashi
    contains

    subroutine NonLinear_ReadParams(Ini)
    use IniFile
    Type(TIniFile) :: Ini

    halofit_version = Ini_Read_Int_File(Ini, 'halofit_version', halofit_default)

    end subroutine NonLinear_ReadParams

    subroutine NonLinear_GetNonLinRatios(CAMB_Pk)
    !Fill the CAMB_Pk%nonlin_scaling array with sqrt(non-linear power/linear power)
    !for each redshift and wavenumber
    !This implementation uses Halofit
    type(MatterPowerData) :: CAMB_Pk
    integer itf
    real(dl) a,plin,pq,ph,pnl,rk
    real(dl) sig,rknl,rneff,rncur,d1,d2
    real(dl) diff,xlogr1,xlogr2,rmid
    integer i

    !!BR09 putting neutrinos into the matter as well, not sure if this is correct, but at least one will get a consisent omk.
    
    !! DM16: modification to treat axions in non-linear lensing S4 fiducial model.
    ! Axion mass kluge: include in computation on non-linear ratio only
    ! for masses that are non-linear at z>2.
    ! Otherwise they are treated as quintessence, i.e. ignored in this.
    ! Boundary found by hand for fiducial "test field" approximation,
    ! valid for low axion density.
    ! Make sure to set the same boundary mass, 1.e-25 eV, in subroutine outtransf.
    if (CP%ma.ge.1.e-25) then
       omm0 = CP%omegac+CP%omegab+CP%omegan+CP%omegaax
    else
       omm0 = CP%omegac+CP%omegab+CP%omegan
    end if
    fnu = CP%omegan/omm0
    
    
    CAMB_Pk%nonlin_ratio = 1
    
    do itf = 1, CAMB_Pk%num_z
       
       ! calculate nonlinear wavenumber (rknl), effective spectral index (rneff) and
       ! curvature (rncur) of the power spectrum at the desired redshift, using method
       ! described in Smith et al (2002).
       a = 1/real(1+CAMB_Pk%Redshifts(itf),dl)
       om_m = omega_m(a, omm0, CP%omegav, w_lam,wa_ppf)
       om_v = omega_v(a, omm0, CP%omegav, w_lam,wa_ppf)
       acur = a
       xlogr1=-2.0
       xlogr2=3.5
       do
          rmid=(xlogr2+xlogr1)/2.0
          rmid=10**rmid
          call wint(CAMB_Pk,itf,rmid,sig,d1,d2)
          diff=sig-1.0
          if (abs(diff).le.0.001) then
             rknl=1./rmid
             rneff=-3-d1
             rncur=-d2
             exit
          elseif (diff.gt.0.001) then
             xlogr1=log10(rmid)
          elseif (diff.lt.-0.001) then
             xlogr2=log10(rmid)
          endif
          if (xlogr2 < -1.9999) then
             !is still linear, exit
             goto 101
          else if (xlogr2>3.4999) then
             ! Totally crazy non-linear
             global_error_flag=349
                write(*,*) 'Error in halofit'
                goto 101
             end if
          end do
          
          
          ! now calculate power spectra for a logarithmic range of wavenumbers (rk)
          
          do i=1, CAMB_PK%num_k
             rk = exp(CAMB_Pk%log_kh(i))
             
             if (rk > Min_kh_nonlinear) then
                
                ! linear power spectrum !! Remeber => plin = k^3 * P(k) * constant
                ! constant = 4*pi*V/(2*pi)^3
                
                plin= MatterPowerData_k(CAMB_PK, rk, itf)*(rk**3/(2*pi**2))
                
                ! calculate nonlinear power according to halofit: pnl = pq + ph,
                ! where pq represents the quasi-linear (halo-halo) power and
                ! where ph is represents the self-correlation halo term.
                call halofit(rk,rneff,rncur,rknl,plin,pnl,pq,ph)   ! halo fitting formula
                CAMB_Pk%nonlin_ratio(i,itf) = sqrt(pnl/plin)
                !write(9,*) rk,plin/(rk**3/(2*pi**2)),pnl/(rk**3/(2*pi**2)) ! DM: write out pnl inside halofit
             end if
             
          enddo
          
101       continue
       end do
       
     end subroutine NonLinear_GetNonLinRatios
     
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    subroutine halofit(rk,rn,rncur,rknl,plin,pnl,pq,ph)
    implicit none

    real(dl) gam,a,b,c,xmu,xnu,alpha,beta,f1,f2,f3
    real(dl) rk,rn,plin,pnl,pq,ph,plinaa
    real(dl) rknl,y,rncur
    real(dl) f1a,f2a,f3a,f1b,f2b,f3b,frac
    real(dl) extragam, peacock_fudge

    if (halofit_version ==halofit_original .or. halofit_version ==halofit_bird &
        .or. halofit_version == halofit_peacock) then
    ! halo model nonlinear fitting formula as described in
    ! Appendix C of Smith et al. (2002)
    !SPB11: Standard halofit underestimates the power on the smallest scales by a
    !factor of two. Add an extra correction from the simulations in Bird, Viel,
    !Haehnelt 2011 which partially accounts for this.
    if (halofit_version ==halofit_bird) then
        extragam = 0.3159 -0.0765*rn -0.8350*rncur
        gam=extragam+0.86485+0.2989*rn+0.1631*rncur
    else
        gam=0.86485+0.2989*rn+0.1631*rncur
    end if
    a=1.4861+1.83693*rn+1.67618*rn*rn+0.7940*rn*rn*rn+ &
        0.1670756*rn*rn*rn*rn-0.620695*rncur
    a=10**a
    b=10**(0.9463+0.9466*rn+0.3084*rn*rn-0.940*rncur)
    c=10**(-0.2807+0.6669*rn+0.3214*rn*rn-0.0793*rncur)
    xmu=10**(-3.54419+0.19086*rn)
    xnu=10**(0.95897+1.2857*rn)
    alpha=1.38848+0.3701*rn-0.1452*rn*rn
    beta=0.8291+0.9854*rn+0.3400*rn**2+fnu*(-6.4868+1.4373*rn**2)
    elseif (halofit_version == halofit_takahashi) then
        !RT12 Oct: the halofit in Smith+ 2003 predicts a smaller power
        !than latest N-body simulations at small scales.
        !Update the following fitting parameters of gam,a,b,c,xmu,xnu,
        !alpha & beta from the simulations in Takahashi+ 2012.
        !The improved halofit accurately provide the power spectra for WMAP
        !cosmological models with constant w.
        gam=0.1971-0.0843*rn+0.8460*rncur
        a=1.5222+2.8553*rn+2.3706*rn*rn+0.9903*rn*rn*rn+ &
            0.2250*rn*rn*rn*rn-0.6038*rncur+0.1749*om_v*(1.+w_lam+wa_ppf*(1-acur))
        a=10**a
        b=10**(-0.5642+0.5864*rn+0.5716*rn*rn-1.5474*rncur+ &
            0.2279*om_v*(1.+w_lam+wa_ppf*(1-acur)))
        c=10**(0.3698+2.0404*rn+0.8161*rn*rn+0.5869*rncur)
        xmu=0.
        xnu=10**(5.2105+3.6902*rn)
        alpha=abs(6.0835+1.3373*rn-0.1959*rn*rn-5.5274*rncur)
        beta=2.0379-0.7354*rn+0.3157*rn**2+1.2490*rn**3+ &
            0.3980*rn**4-0.1682*rncur + fnu*(1.081 + 0.395*rn**2)
    else
        stop 'Unknown halofit_version'
    end if

    if(abs(1-om_m).gt.0.01) then ! omega evolution
        f1a=om_m**(-0.0732)
        f2a=om_m**(-0.1423)
        f3a=om_m**(0.0725)
        f1b=om_m**(-0.0307)
        f2b=om_m**(-0.0585)
        f3b=om_m**(0.0743)
        frac=om_v/(1.-om_m)
        f1=frac*f1b + (1-frac)*f1a
        f2=frac*f2b + (1-frac)*f2a
        f3=frac*f3b + (1-frac)*f3a
    else
        f1=1.0
        f2=1.
        f3=1.
    endif

    y=(rk/rknl)


    ph=a*y**(f1*3)/(1+b*y**(f2)+(f3*c*y)**(3-gam))
    ph=ph/(1+xmu*y**(-1)+xnu*y**(-2))*(1+fnu*0.977)
    plinaa=plin*(1+fnu*47.48*rk**2/(1+1.5*rk**2))
    pq=plin*(1+plinaa)**beta/(1+plinaa*alpha)*exp(-y/4.0-y**2/8.0)

    pnl=pq+ph

    if (halofit_version == halofit_peacock) then
        !From http://www.roe.ac.uk/~jap/haloes/
        !(P-P_linear) -> (P-P_linear) * (1+2y^2)/(1+y^2), where y = k/10 h Mpc^(-1).
        peacock_fudge = rk/10
        pnl = plin + (pnl-plin)*(1+2*peacock_fudge**2)/(1+peacock_fudge**2)
    end if

    end subroutine halofit


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! The subroutine wint, finds the effective spectral quantities
    ! rknl, rneff & rncur. This it does by calculating the radius of
    ! the Gaussian filter at which the variance is unity = rknl.
    ! rneff is defined as the first derivative of the variance, calculated
    ! at the nonlinear wavenumber and similarly the rncur is the second
    ! derivative at the nonlinear wavenumber.

    subroutine wint(CAMB_Pk,itf,r,sig,d1,d2)
    implicit none
    integer, intent(in) :: itf
    type(MatterPowerData) :: CAMB_Pk
    real(dl) sum1,sum2,sum3,t,y,x,w1,w2,w3
    real(dl) x2,rk, fac,r, sig, d1,d2, anorm
    integer i,nint

    nint=3000
    sum1=0.d0
    sum2=0.d0
    sum3=0.d0
    anorm = 1/(2*pi**2)
    do i=1,nint
        t=(i-0.5_dl)/nint
        y=-1.d0+1.d0/t
        rk=y
        d2=MatterPowerData_k(CAMB_PK, rk, itf)*(rk**3*anorm)
        x=y*r
        x2=x*x
        w1=exp(-x2)
        w2=2*x2*w1
        w3=4*x2*(1-x2)*w1
        fac=d2/y/t/t
        sum1=sum1+w1*fac
        sum2=sum2+w2*fac
        sum3=sum3+w3*fac
    enddo
    sum1=sum1/nint
    sum2=sum2/nint
    sum3=sum3/nint
    sig=sqrt(sum1)
    d1=-sum2/sum1
    d2=-sum2*sum2/sum1/sum1 - sum3/sum1

    end subroutine wint

    !!JD 08/13 generalize to variable w

    function omega_m(aa,om_m0,om_v0,wval,waval)
    implicit none
    real(dl) omega_m,omega_t,om_m0,om_v0,aa,wval,waval,Qa2
    Qa2= aa**(-1.0-3.0*(wval+waval))*dexp(-3.0*(1-aa)*waval)
    omega_t=1.0+(om_m0+om_v0-1.0)/(1-om_m0-om_v0+om_v0*Qa2+om_m0/aa)
    omega_m=omega_t*om_m0/(om_m0+om_v0*aa*Qa2)
    end function omega_m

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! evolution of omega lambda with expansion factor

    function omega_v(aa,om_m0,om_v0,wval,waval)
    implicit none
    real(dl) aa,omega_v,om_m0,om_v0,omega_t,wval,waval,Qa2
    Qa2= aa**(-1.0-3.0*(wval+waval))*dexp(-3.0*(1-aa)*waval)
    omega_t=1.0+(om_m0+om_v0-1.0)/(1-om_m0-om_v0+om_v0*Qa2+om_m0/aa)
    omega_v=omega_t*om_v0*Qa2/(om_v0*Qa2+om_m0/aa)
    end function omega_v

    !!JD end generalize to variable w

    end module NonLinear


    !workaround for f90 circular-module reference
    subroutine NonLinear_GetRatios(CAMB_Pk)
    use Transfer
    use NonLinear
    type(MatterPowerData) :: CAMB_Pk

    call NonLinear_GetNonLinRatios(CAMB_Pk)

    end subroutine NonLinear_GetRatios



    subroutine NonLinear_GetRatios_all(CAMB_Pk)
    use Transfer
    use NonLinear
    type(MatterPowerData) :: CAMB_Pk

    stop 'Halofit module doesn''t support non-linear velocities'

    end subroutine NonLinear_GetRatios_All

