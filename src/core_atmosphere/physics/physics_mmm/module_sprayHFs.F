MODULE module_sprayHFs

IMPLICIT NONE
PUBLIC :: sprayHFs
PRIVATE

!===============================================================================
!
!        SEASTATE-DEPENDENT AIR-SEA HEAT FLUXES WITH SPRAY IN HIGH WINDS 
!
!     This module contains subroutines to incorporate seastate-dependent sea 
! spray heat flux physics into an existing bulk surface layer scheme in coupled 
! regional or global Earth system models.  This code is maintained at 
! https://github.com/bwbarr/sprayHFs, where documentation and the latest version
! are available.  Please send any questions about model physics or this module's
! implementation to Ben Barr at benjamin.barr@whoi.edu.
!
!===============================================================================
!
! Implemented in MPAS/MONAN models by S.R.Freitas (INPE @ Jan 2025)
!===============================================================================


CONTAINS


!===============================================================================

SUBROUTINE sprayHFs(z_1,t_1,q_1,U_1,gf,p_0,t_0,eps,dcp,swh,mss,L,z0,z0t,z0q,&
               z_ref,whichSSGF,paramWaves,dHS1_spr,dHL1_spr,dthstar_spr,&
               dqstar_spr,dthvstar_spr,dthref_spr,dtref_spr,dqref_spr,dsref_spr,&
               tau,H_S0pr,H_L0pr,M_spr,H_Tspr,H_Rspr,H_SNspr,H_Lspr,alpha_S,&
               beta_S,beta_L,gamma_S,gamma_L)

!-------------------------------------------------------------------------------
!
! Subroutine to calculate spray modifications to heat fluxes and surface layer 
! properties.  Implemented per Barr et al. (2023) (hereafter BCF23).  Heat 
! fluxes from the ocean to the atmosphere are defined as positive, and momentum
! fluxes (stress) from the atmosphere to the ocean are defined as positive.  
! Positive heat fluxes correspond to negative values of the turbulent flux 
! scales thstar, qstar, and thvstar.
!
! Inputs:
!     z_1 - height of lowest atmospheric model mass level [m]
!     t_1 - air temperature at z_1 [K]
!     q_1 - air specific humidity at z_1 [kg kg-1]
!     U_1 - current-relative windspeed magnitude at z_1 [m s-1].  Does not
!         include gustiness.
!     gf - gust factor accounting for gustiness at low winds [-]. gf = Sg_1/U_1,
!         where Sg_1 is the windspeed magnitude at z_1 including gustiness.  If 
!         not accounting for gustiness, pass gf = 1.
!     p_0 - surface pressure [Pa]
!     t_0 - sea surface temperature (interfacial temp, not bulk layer temp) [K]
!     eps - wave energy dissipation flux [W m-2]
!     dcp - dominant wave phase speed [m s-1]
!     swh - significant wave height [m]
!     mss - mean squared waveslope [-]
!     L - Obukhov stability length [m].  This is considered "known" for this
!         subroutine's computations but may be determined iteratively by the 
!         existing bulk algorithm code.
!     z0 - momentum roughness length [m].  If surface stress comes from a bulk
!         algorithm (e.g., COARE), pass the parameterized z0.  If surface stress
!         comes from a wave model, invert and pass an equivalent z0 based on MO 
!         theory using stress, U_1, gf, and L.
!     z0t - thermal roughness length [m]
!     z0q - moisture roughness length [m].  Does not have to be equal to z0t.
!     z_ref - reference height [m] for calculating spray changes to subgrid surface 
!         layer potential temperature, temperature, specific humidity, and 
!         saturation ratio.  This may be used to update the existing code's 
!         calculations of thermodynamic variables at a reference height (e.g., 
!         2 m).  The physics of this parameterization do not produce direct 
!         changes to the subgrid wind profile (e.g., 10-m windspeed components).  
!         z_ref should not be larger than z_1.  Pass -1 for z_ref to calculate 
!         these changes at the mid-spray-layer height.
!     whichSSGF - string naming which SSGF to use.  Options are:
!         'BCF23_Seastate': seastate-dependent, per BCF23 Eq. 1.
!         'F94_MOM80': widely used F94 wind-based SSGF, using original whitecap 
!             fraction per Monahan and O'Muircheartaigh (1980).
!         'F94_BCF23': updated F94 wind-based SSGF given as BCF23 Eq. 3.
!     paramWaves - True to parameterize wave processes, False to use input values
!         for wave properties (eps, dcp, swh, mss).  If True, dummy values
!         (e.g., zeroes) should be passed in the subroutine call for all four
!         wave properties.  If False, externally-defined values (e.g., from a 
!         wave model) should be passed in the subroutine call.
! Outputs:
!     dHS1_spr - change to total surface sensible heat flux due to spray [W m-2]
!     dHL1_spr - change to total surface latent heat flux due to spray [W m-2]
!     dthstar_spr - change to total theta flux scale due to spray [K]
!     dqstar_spr - change to total q flux scale due to spray [kg kg-1]
!     dthvstar_spr - change to total thetav (i.e., buoyancy) flux scale due to 
!         spray [K]
!     dthref_spr - change to potential temperature at user-specified reference
!         height due to spray heat flux feedback [K], (+) = warming
!     dtref_spr - change to temperature at user-specified reference height due 
!         to spray heat flux feedback [K], (+) = warming
!     dqref_spr - change to q at user-specified reference height due to spray 
!         heat flux feedback [kg kg-1], (+) = moistening
!     dsref_spr - change to saturation ratio at user-specified reference height 
!         due to spray heat flux feedback [-], (+) = s increases
!     tau - bulk stress [Pa].  Provided for comparison to bulk stress from 
!         existing bulk code.  Should not be used to override existing bulk stress.
!     H_S0pr - bulk SHF without spray [W m-2].  Provided for comparison to bulk
!         SHF from existing code.  Should not be used to override existing bulk
!         SHF.
!     H_L0pr - bulk LHF without spray [W m-2].  Provided for comparison to bulk
!         LHF from existing code.  Should not be used to override existing bulk
!         LHF.
!     M_spr - spray generated mass flux [kg m-2 s-1]
!     H_Tspr - spray heat flux due to temperature change [W m-2].  This is
!         exactly equal to the spray enthalpy flux, H_Kspr.  H_Kspr is not a 
!         variable used in this code, but it is discussed in the BCF23 paper.
!     H_Rspr - spray heat flux due to size change [W m-2]
!     H_SNspr - net spray sensible heat flux [W m-2], including both the direct
!         heating flux (H_Sspr) and the evaporative cooling flux (-H_Rspr), 
!         i.e., H_SNspr = H_Sspr - H_Rspr.
!     H_Lspr - spray latent heat flux [W m-2]
!     alpha_S - feedback coefficient between H_Sspr and surface layer 
!         temp/humidity profiles [-]
!     beta_S - feedback coefficient between H_Rspr and surface layer
!         temp/humidity profiles [-]
!     beta_L - feedback coefficient between H_Lspr and surface layer
!         temp/humidity profiles [-]
!     gamma_S - feedback coefficient between spray heat fluxes and interfacial
!         sensible heat flux [-]
!     gamma_L - feedback coefficient between spray heat fluxes and interfacial
!         latent heat flux [-]
!
!-------------------------------------------------------------------------------

REAL,INTENT(IN) :: z_1,t_1,q_1,U_1,gf,p_0,t_0,L,z0,z0t,z0q,z_ref
CHARACTER(LEN=*),INTENT(IN) :: whichSSGF
LOGICAL,INTENT(IN) :: paramWaves
REAL,INTENT(INOUT) :: eps,dcp,swh,mss
REAL,INTENT(OUT) :: dHS1_spr
REAL,INTENT(OUT) :: dHL1_spr
REAL,INTENT(OUT) :: dthstar_spr
REAL,INTENT(OUT) :: dqstar_spr
REAL,INTENT(OUT) :: dthvstar_spr
REAL,INTENT(OUT) :: dthref_spr
REAL,INTENT(OUT) :: dtref_spr
REAL,INTENT(OUT) :: dqref_spr
REAL,INTENT(OUT) :: dsref_spr
REAL,INTENT(OUT) :: tau
REAL,INTENT(OUT) :: H_S0pr
REAL,INTENT(OUT) :: H_L0pr
REAL,INTENT(OUT) :: M_spr
REAL,INTENT(OUT) :: H_Tspr
REAL,INTENT(OUT) :: H_Rspr
REAL,INTENT(OUT) :: H_SNspr
REAL,INTENT(OUT) :: H_Lspr
REAL,INTENT(OUT) :: alpha_S
REAL,INTENT(OUT) :: beta_S
REAL,INTENT(OUT) :: beta_L
REAL,INTENT(OUT) :: gamma_S
REAL,INTENT(OUT) :: gamma_L

REAL,PARAMETER :: sprayLB = 10.0    ! Lower bound on U_10 for spray production. 10 m s-1 is typical for whitecapping
REAL,PARAMETER :: kappa = 0.4    ! von Karman constant [-]
REAL,PARAMETER :: g = 9.81    ! Acceleration due to gravity [m s-2]
REAL,PARAMETER :: Rdry = 287.1    ! Dry air gas constant [J kg-1 K-1]
REAL,PARAMETER :: rho_sw = 1030.    ! Density of seawater [kg m-3]
REAL,PARAMETER :: rho_dry = 2160.    ! Density of chrystalline salt [kg m-3]
REAL,PARAMETER :: cp_sw = 4200.    ! Specific heat capacity of seawater [J kg-1 K-1]
REAL,PARAMETER :: cp_a = 1004.67    ! Specific heat capacity of air [J kg-1 K-1]
REAL,PARAMETER :: nu = 2.    ! Number of ions into which NaCl dissociates [-]
REAL,PARAMETER :: Phi_s = 0.924    ! Practical osmotic coefficient at molality of 0.6 [-]
REAL,PARAMETER :: Mw = 18.02    ! Molecular weight of water [g mol-1]
REAL,PARAMETER :: Ms = 58.44    ! Molecular weight of salt [g mol-1]
REAL,PARAMETER :: xs = 0.035    ! Mass fraction of salt in seawater [-]

REAL :: Sg_1,ustar,U_10
REAL :: Lv,delspr,zref,rdryBYr0,y0,q_0,tv_1,rho_a,p_1,p_delsprD2,p_zref,th_0,&
        th_1,tC_1,k_a,nu_a,Dv_a,gammaWB
REAL :: psiH_1,psiH_delspr,psiH_delsprD2,psiH_zref,phisprH_delspr,phisprH_zref,thstar_pr,&
        qstar_pr,G_S,G_L,t_delsprD2pr,th_zref_pr,t_zref_pr,q_delsprD2pr,q_zref_pr,&
        s_delsprD2pr,s_zref_pr
REAL :: zR,H_Sspr,H_Tsprpr,H_Ssprpr,H_Rsprpr,H_Lsprpr
REAL :: etaT,Cs_pr,C_HIG,H_IG,Psi,Chi,Lambda,A,B,C,B2M4AC,s_hatPOS,H_Rspr_IG
REAL :: H_S1,H_L1,H_S0,H_L0,th_zref,t_zref,q_zref,s_zref,thstar,qstar,thvstar_pr,thvstar
CHARACTER(LEN=100) :: Wform

! Droplet radius vector [m]
REAL,DIMENSION(25) :: r0       = (/ 10.,   20.,   30.,   40.,   50.,   60.,&
        70.,   80.,   90., 102.5, 122.5, 157.5,  215.,  300.,  400.,  500.,&
       600.,  700.,  800.,  900.,1037.5, 1250., 1500., 1750., 2000. /)*1e-6
! Droplet bin width vector [m]
REAL,DIMENSION(25) :: delta_r0 = (/ 10.,   10.,   10.,   10.,   10.,   10.,&
        10.,   10.,   10.,   15.,   25.,   45.,   70.,  100.,  100.,  100.,&
       100.,  100.,  100.,  100.,  175.,  250.,  250.,  250.,  250. /)*1e-6
REAL,DIMENSION(25) :: v_g,dmdr0,tauf,Fp,tauT,zT


! 1. Check if U_10 >= sprayLB.  If so, perform spray calculations --------------
Sg_1 = U_1*gf    ! Windspeed magnitude at z_1 including gustiness [m s-1]
ustar = Sg_1*kappa/(LOG(z_1/z0) - stabIntM(z_1/L))    ! Local turbulence intensity (not bulk fric velocity) [m s-1]
U_10 = ustar/kappa/gf*(LOG(10./z0) - stabIntM(10./L))    ! 10m windspeed [m s-1], gustiness removed
IF (U_10 < sprayLB) THEN    ! Set returned fields to zero or dummy (999) values

    dHS1_spr = 0.; dHL1_spr = 0.; dthstar_spr = 0.; dqstar_spr = 0.; dthvstar_spr = 0.;
    dthref_spr = 0.; dtref_spr = 0.; dqref_spr = 0.; dsref_spr = 0.
    tau = 999.; H_S0pr = 999.; H_L0pr = 999.; M_spr = 999.
    H_Tspr = 999.; H_Rspr = 999.; H_SNspr = 999.; H_Lspr = 999.
    alpha_S = 999.; beta_S = 999.; beta_L = 999.; gamma_S = 999.; gamma_L = 999.

ELSE IF (U_10 >= sprayLB) THEN    ! Perform spray calculations

    ! 2. Get properties and environmental variables ----------------------------
    Lv = (2.501-0.00237*(t_0-273.15))*1e6    ! Latent heat of vap for water [J kg-1]
    rdryBYr0 = (xs*rho_sw/rho_dry)**0.33333    ! Ratio of dry salt radius to r0 [-]
    y0 = -nu*Phi_s*Mw/Ms*rho_dry/rho_sw*rdryBYr0**3/(1. - rho_dry/rho_sw*rdryBYr0**3)    ! y for surface seawater [-]
    q_0 = qsat0(t_0,p_0)*(1. + y0)    ! Specific humidity at surface (accounting for salt) [kg kg-1]
    tv_1 = t_1*(1.+0.608*q_1)    ! Virtual temperature at z_1 [K]
    rho_a = (p_0 - 1.25*g*z_1)/(Rdry*tv_1)    ! Air density at z_1 [kg m-3], adjusting pressure using rho_a~1.25 kg m-3
    IF (paramWaves) THEN    ! Parameterize eps, dcp, swh, and mss
        CALL waveProps(U_10,ustar/gf**0.5,rho_a,eps,dcp,swh,mss)    ! Uses bulk friction velocity without gustiness
    END IF
    delspr = MIN(swh,z_1)    ! Spray layer thickness [m], nominally one swh per M&V2014a.  Limited to z_1.
    IF (z_ref < 0.) THEN    ! Set reference height to user-defined value or mid-spray layer height
        zref = delspr/2.
    ELSE
        zref = z_ref
    END IF
    p_1        = p_0 - rho_a*g*z_1    ! Pressure at z_1 [Pa].  All pressure adjustments assume hydrostatic gradient.
    p_delsprD2 = p_0 - rho_a*g*delspr/2.    ! Pressure at delspr/2 [Pa]
    p_zref     = p_0 - rho_a*g*zref    ! Pressure at zref [Pa]
    th_0 = t_0*(1.e5/p_0)**0.286    ! Potential temperature at surface [K]
    th_1 = t_1*(1.e5/p_1)**0.286    ! Potential temperature at z_1 [K]
    tC_1 = t_1 - 273.15    ! Convert t_1 to [C] for property calculations
    k_a = 2.411e-2*(1.+3.309e-3*tC_1-1.441e-6*tC_1**2)    ! Thermal conductivity of air [W m-1 K-1]
    nu_a = 1.326e-5*(1.+6.542e-3*tC_1+8.301e-6*tC_1**2-4.84e-9*tC_1**3)    ! Kin visc of air [m2 s-1]
    Dv_a = 2.11e-5*((tC_1+273.)/273.)**1.94    ! Water vapor diffusivity in air [m2 s-1]
    gammaWB = 240.97*17.502/(tC_1+240.97)**2    ! gamma = (dqsat/dT)/qsat [K-1], per Buck (1981)

    ! 3. Calculate fluxes w/o spray and some stability/feedback items ----------
    psiH_1        = stabIntH(z_1/L)    ! Stability integral for heat at z_1 [-]
    psiH_delspr   = stabIntH(delspr/L)    ! Stability integral for heat at delspr [-]
    psiH_delsprD2 = stabIntH(delspr/2./L)    ! Stability integral for heat at delspr/2 [-]
    psiH_zref     = stabIntH(zref/L)    ! Stability integral for heat at zref [-]
    phisprH_delspr = stabIntSprayH(delspr/L)    ! Stability integral for heat with spray at delspr [-]
    phisprH_zref   = stabIntSprayH(zref/L)    ! Stability integral for heat with spray at zref [-]
    thstar_pr = -(th_0 - th_1)*kappa/(LOG(z_1/z0t) - psiH_1)    ! theta flux scale without spray [K]
    qstar_pr  = -( q_0 -  q_1)*kappa/(LOG(z_1/z0q) - psiH_1)    ! q flux scale without spray [kg kg-1]
    G_S = rho_a*cp_a*kappa*ustar    ! Dimensional group for sensible heat [W m-2 K-1]
    G_L = rho_a*Lv*kappa*ustar    ! Dimensional group for latent heat [W m-2]
    tau = rho_a*ustar**2/gf    ! Bulk stress [Pa]
    H_S0pr = -G_S/kappa*thstar_pr    ! Bulk SHF without spray [W m-2]
    H_L0pr = -G_L/kappa*qstar_pr    ! Bulk LHF without spray [W m-2]
    t_delsprD2pr = (th_0 - H_S0pr/G_S*(LOG(delspr/2./z0t) - psiH_delsprD2))*(p_delsprD2/1.e5)**0.286   ! t at mid-layer w/o fdbk [K]
    th_zref_pr   =  th_0 - H_S0pr/G_S*(LOG(zref/z0t)      - psiH_zref    )   ! Potential temp at zref w/o fdbk [K]
    q_delsprD2pr =   q_0 - H_L0pr/G_L*(LOG(delspr/2./z0q) - psiH_delsprD2)    ! q at mid-layer w/o fdbk [kg kg-1]
    q_zref_pr    =   q_0 - H_L0pr/G_L*(LOG(zref/z0q)      - psiH_zref    )    ! q at zref w/o fdbk [kg kg-1]
    t_zref_pr    = th_zref_pr*(p_zref/1.e5)**0.286    ! t at zref w/o fdbk [K]
    s_delsprD2pr = satratio(t_delsprD2pr,p_delsprD2,q_delsprD2pr,0.99999)    ! s at mid-layer w/o fdbk [-]
    s_zref_pr    = satratio(t_zref_pr   ,p_zref    ,q_zref_pr   ,0.99999)    ! s at zref w/o fdbk [-]
    gamma_S = (LOG(delspr/z0t) - psiH_delspr - 1. + phisprH_delspr)/(LOG(z_1/z0t) - psiH_1)    ! Interfacial SHF fdbk coeff [-]
    gamma_L = (LOG(delspr/z0q) - psiH_delspr - 1. + phisprH_delspr)/(LOG(z_1/z0q) - psiH_1)    ! Interfacial LHF fdbk coeff [-]

    ! 4. Calculate spray generation and some droplet properties ----------------
    v_g = fall_velocity_PK97(r0)    ! Droplet settling velocity [m s-1]
    IF (whichSSGF == 'BCF23_Seastate') THEN
        CALL ssgf_BCF23(eps,swh,dcp,mss,ustar,z0,L,gf,r0,delta_r0,v_g,M_spr,dmdr0)
    ELSE IF ((whichSSGF == 'F94_MOM80') .OR. (whichSSGF == 'F94_BCF23')) THEN
        IF (whichSSGF == 'F94_MOM80') THEN
            Wform = 'MOM80'
        ELSE IF (whichSSGF == 'F94_BCF23') THEN
            Wform = 'BCF23'
        END IF
        CALL ssgf_F94(U_10,Wform,r0,delta_r0,M_spr,dmdr0)
    END IF
    tauf = delspr/v_g    ! Characteristic droplet settling time [s]
    Fp = 1. + 0.25*(2.*v_g*r0/nu_a)**0.5    ! Slip factor (Pr&Kl) [-]
    tauT = rho_sw*cp_sw*r0**2/3./k_a/Fp    ! Characteristic droplet cooling timescale [s]
    zT = MIN(0.5*delspr,0.5*v_g*tauT)    ! H_Tspr ambient property height [m]
    zR = 0.5*delspr    ! H_Rspr ambient property height [m]

    ! 5. Calculate spray heat fluxes without feedback --------------------------
    ! Initialize spray heat fluxes to zero
    H_Sspr = 0.
    H_Rspr = 0.
    H_Lspr = 0.
    ! Calculate spray heat fluxes without feedback
    CALL update_HFs(H_Sspr,H_Rspr,H_Lspr,r0,delta_r0,zT,zR,tauf,tauT,Fp,&
            dmdr0,p_0,gammaWB,y0,t_0,rho_a,Dv_a,xs,delspr,z0t,z0q,L,th_0,q_0,&
            G_S,G_L,H_S0pr,H_L0pr,gamma_S,gamma_L,Lv,cp_a,rho_sw,nu,Phi_s,Mw,&
            Ms,cp_sw,g,H_Tspr)
    ! Store spray heat fluxes without feedback
    H_Tsprpr = H_Tspr    ! Spray HF due to temp change, without feedback [W m-2]
    H_Ssprpr = H_Sspr    ! Spray SHF, without feedback [W m-2]
    H_Rsprpr = H_Rspr    ! Spray HF due to size change, without feedback [W m-2]
    H_Lsprpr = H_Lspr    ! Spray LHF, without feedback [W m-2]

    ! 6. Calculate spray heat fluxes with feedback -----------------------------
    ! Use simple model to calculate initial guess for H_Rspr
    etaT = 17.502*240.97/(t_delsprD2pr-273.15+240.97)**2    ! [K-1]
    Cs_pr = (1. + y0 - s_delsprD2pr)**2/(1. - s_delsprD2pr)    ! [-]
    IF (delspr < 4.) THEN
        C_HIG = 1.0    ! Tuneable constant for equivalent height [-]
    ELSE IF (delspr > 10.) THEN
        C_HIG = 0.7
    ELSE
        C_HIG = -0.05*delspr + 1.2
    END IF
    H_IG = MIN(C_HIG*delspr,z_1)    ! Equivalent height for heating in simple model [m]
    Psi = etaT*LOG(z_1/H_IG)/G_S*H_Tsprpr    ! [-]
    Chi = etaT*LOG(z_1/H_IG)/G_S + LOG(z_1/H_IG)/G_L/q_delsprD2pr    ! [m2 W-1]
    Lambda = (Psi - 1. + (1. + y0)/s_delsprD2pr)/Chi    ! [W m-2]
    A = H_Rsprpr/Cs_pr + 1./s_delsprD2pr/Chi    ! [W m-2]
    B = -Lambda - y0/s_delsprD2pr/Chi    ! [W m-2]
    C = y0*Lambda    ! [W m-2]
    B2M4AC = B**2 - 4.*A*C    ! Argument of radical in quadratic formula [W2 m-4]
    IF (B2M4AC >= 0.) THEN    ! Real solutions
        s_hatPOS = (-B + SQRT(B2M4AC))/2./A    ! Larger root [-], seems physical
        H_Rspr_IG = s_hatPOS**2/(s_hatPOS - y0)*H_Rsprpr/Cs_pr    ! IG for H_Rspr [W m-2]
    ELSE    ! For imaginary solutions set H_Rspr_IG to zero
        H_Rspr_IG = 0.
    END IF
    ! Adjust heat fluxes using IG
    H_Rspr = H_Rspr_IG
    H_Lspr = H_Rspr_IG + H_Tspr - H_Sspr
    ! Calculate spray heat fluxes with feedback using H_Rspr_IG
    CALL update_HFs(H_Sspr,H_Rspr,H_Lspr,r0,delta_r0,zT,zR,tauf,tauT,Fp,&
            dmdr0,p_0,gammaWB,y0,t_0,rho_a,Dv_a,xs,delspr,z0t,z0q,L,th_0,q_0,&
            G_S,G_L,H_S0pr,H_L0pr,gamma_S,gamma_L,Lv,cp_a,rho_sw,nu,Phi_s,Mw,&
            Ms,cp_sw,g,H_Tspr)

    ! 7. Calculate final diagnosed quantities ----------------------------------
    H_SNspr = H_Sspr - H_Rspr    ! Net spray sensible heat flux [W m-2]
    dHS1_spr = gamma_S*H_SNspr    ! Change to total SHF due to spray [W m-2]
    dHL1_spr = gamma_L*H_Lspr    ! Change to total LHF due to spray [W m-2]
    dthstar_spr = -kappa/G_S*dHS1_spr    ! Change to thstar due to spray [K]
    dqstar_spr  = -kappa/G_L*dHL1_spr    ! Change to qstar due to spray [kg kg-1]
    
    !--srf: got floating-point exception - erroneous arithmetic operation.
    !alpha_S = H_Sspr/H_Ssprpr    ! Feedback coefficient for H_Sspr [-]
    !beta_S = H_Rspr/H_Rsprpr    ! Feedback coefficient for H_Rspr [-]
    !beta_L = H_Lspr/H_Lsprpr    ! Feedback coefficient for H_Lspr [-]
    !----

    H_S1 = H_S0pr + dHS1_spr    ! Total SHF [W m-2]
    H_L1 = H_L0pr + dHL1_spr    ! Total LHF [W m-2]
    H_S0 = H_S1 - H_SNspr    ! Surface SHF (interfacial with feedback) [W m-2]
    H_L0 = H_L1 - H_Lspr    ! Surface LHF (interfacial with feedback) [W m-2]
    IF (zref < delspr) THEN    ! zref is within spray layer
        th_zref = th_0 - 1./G_S*(H_S0*(LOG(zref/z0t) - psiH_zref) &
                + zref/delspr*(1. - phisprH_zref)*H_SNspr)    ! th at zref w/ fdbk [K]
        q_zref  =  q_0 - 1./G_L*(H_L0*(LOG(zref/z0q) - psiH_zref) &
                + zref/delspr*(1. - phisprH_zref)*H_Lspr)    ! q at zref w/ fdbk [kg kg-1]
    ELSE    ! zref is above spray layer
        th_zref = th_1 + H_S1/G_S*(LOG(z_1/zref) - psiH_1 + psiH_zref)    ! [K]
        q_zref  =  q_1 + H_L1/G_L*(LOG(z_1/zref) - psiH_1 + psiH_zref)    ! [kg kg-1]
    END IF
    t_zref = th_zref*(p_zref/1.e5)**0.286    ! t at zref w/ fdbk [K]
    s_zref = satratio(t_zref,p_zref,q_zref,0.99999)    ! s at zref w/ fdbk [-]
    dthref_spr = th_zref - th_zref_pr    ! th change at zref due to feedback [K], (+) = warming
    dtref_spr  = t_zref  - t_zref_pr    ! t change at zref due to feedback [K], (+) = warming
    dqref_spr  = q_zref  - q_zref_pr    ! q change at zref due to feedback [kg kg-1], (+) = moistening
    dsref_spr  = s_zref  - s_zref_pr    ! s change at zref due to feedback [-], (+) = incr s
    thstar = -H_S1*kappa/G_S    ! theta flux scale with spray (based on total SHF) [K]
    qstar  = -H_L1*kappa/G_L    ! q flux scale with spray (based on total LHF) [kg kg-1]
    thvstar_pr = thstar_pr*(1. + 0.61*q_1) + 0.61*t_1*qstar_pr    ! thetav flux scale without spray [K]
    thvstar    =    thstar*(1. + 0.61*q_1) + 0.61*t_1*qstar    ! thetav flux scale with spray [K]
    dthvstar_spr = thvstar - thvstar_pr    ! Change to thvstar due to spray [K]

END IF

END SUBROUTINE sprayHFs

!===============================================================================

SUBROUTINE update_HFs(H_Sspr,H_Rspr,H_Lspr,r0,delta_r0,zT,zR,tauf,tauT,Fp,&
            dmdr0,p_0,gammaWB,y0,t_0,rho_a,Dv_a,xs,delspr,z0t,z0q,L,th_0,q_0,&
            G_S,G_L,H_S0pr,H_L0pr,gamma_S,gamma_L,Lv,cp_a,rho_sw,nu,Phi_s,Mw,&
            Ms,cp_sw,g,H_Tspr)

!-------------------------------------------------------------------------------
!
! Subroutine to update spray heat fluxes based on previous values.
!
! Updated (INOUT) Fields:
!     H_Sspr [W m-2], H_Rspr [W m-2], H_Lspr [W m-2]
! Inputs:
!     r0 [m], delta_r0 [m], zT [m], zR [m], tauf [s], tauT [s], Fp [-], 
!     dmdr0 [kg m-2 s-1 m-1], p_0 [Pa], gammaWB [K-1], y0 [-], t_0 [K], 
!     rho_a [kg m-3], Dv_a [m2 s-1], xs [-], delspr [m], z0t [m], z0q [m], 
!     L [m], th_0 [K], q_0 [kg kg-1], G_S [W m-2 K-1], G_L [W m-2], 
!     H_S0pr [W m-2], H_L0pr [W m-2], gamma_S [-], gamma_L [-], Lv [J kg-1], 
!     cp_a [J kg-1 K-1], rho_sw [kg m-3], nu [-], Phi_s [-], Mw [g mol-1], 
!     Ms [g mol-1], cp_sw [J kg-1 K-1], g [m s-2]
! Outputs:
!     H_Tspr [W m-2]
!
!-------------------------------------------------------------------------------

REAL,INTENT(INOUT) :: H_Sspr,H_Rspr,H_Lspr
REAL,INTENT(IN) :: zR,p_0,gammaWB,y0,t_0,rho_a,Dv_a,xs,delspr,z0t,z0q,L,th_0,&
        q_0,G_S,G_L,H_S0pr,H_L0pr,gamma_S,gamma_L,Lv,cp_a,rho_sw,nu,Phi_s,Mw,&
        Ms,cp_sw,g
REAL,DIMENSION(25),INTENT(IN) :: r0,delta_r0,zT,tauf,tauT,Fp,dmdr0
REAL,INTENT(OUT) :: H_Tspr

INTEGER :: i
REAL :: H_S0,H_L0,p_zR,t_zR,q_zR,qsat0_zR,betaWB_zR,s_zR
REAL,DIMENSION(25) :: p_zT,t_zT,q_zT,qsat0_zT,betaWB_zT,s_zT,wetdep_zT,tWB_zT,&
        tdropf,tauR,req,rf,stabIntH_z0tPzT,stabIntH_z0qPzT,stabIntSprayH_z0tPzT,&
        stabIntSprayH_z0qPzT

! 1. Calculate surface heat fluxes (interfacial with feedback) using previous values of spray heat fluxes.
!    If spray heat fluxes are zero, this gives the no-spray interfacial fluxes.
H_S0 = H_S0pr - (1. - gamma_S)*(H_Sspr - H_Rspr)    ! Surface SHF (interfacial with feedback) [W m-2]
H_L0 = H_L0pr - (1. - gamma_L)*H_Lspr    ! Surface LHF (interfacial with feedback) [W m-2]

! 2. Get thermodynamic parameters for calculating H_Tspr
DO i = 1,25
    stabIntH_z0tPzT(i) = stabIntH((z0t+zT(i))/L)
    stabIntH_z0qPzT(i) = stabIntH((z0q+zT(i))/L)
    stabIntSprayH_z0tPzT(i) = stabIntSprayH((z0t+zT(i))/L)
    stabIntSprayH_z0qPzT(i) = stabIntSprayH((z0q+zT(i))/L)
END DO
p_zT = p_0 - rho_a*g*zT    ! Pressure at zT [Pa]
t_zT = (th_0 - 1./G_S*(H_S0*(LOG((z0t+zT)/z0t) - stabIntH_z0tPzT) &
        + zT/delspr*(1. - stabIntSprayH_z0tPzT)*(H_Sspr - H_Rspr)))*(p_zT/1.e5)**0.286    ! Temp at zT [K]
q_zT =   q_0 - 1./G_L*(H_L0*(LOG((z0q+zT)/z0q) - stabIntH_z0qPzT) &
        + zT/delspr*(1. - stabIntSprayH_z0qPzT)*H_Lspr)    ! Spec hum at zT [kg kg-1]
DO i = 1,25
    qsat0_zT(i) = qsat0(t_zT(i),p_zT(i))    ! Saturation specific humidity at zT [kg kg-1]
    s_zT(i) = satratio(t_zT(i),p_zT(i),q_zT(i),0.99999)    ! Saturation ratio at zT [-]
END DO
betaWB_zT = 1./(1. + Lv*gammaWB*(1. + y0)/cp_a*qsat0_zT)    ! Wetbulb coefficient at zT [-]
wetdep_zT = (1. - s_zT/(1. + y0))*(1. - betaWB_zT)/gammaWB    ! Wetbulb depression at zT [-]
tWB_zT = t_zT - wetdep_zT    ! Wetbulb temperature at zT [K]
tdropf = tWB_zT + (t_0 - tWB_zT)*EXP(-tauf/tauT)    ! Final droplet temperature [K]

! 3. Get thermodynamic parameters for calculating H_Rspr
p_zR = p_0 - rho_a*g*zR    ! Pressure at zR [Pa]
t_zR = (th_0 - 1./G_S*(H_S0*(LOG((z0t+zR)/z0t) - stabIntH((z0t+zR)/L)) &
        + zR/delspr*(1. - stabIntSprayH((z0t+zR)/L))*(H_Sspr - H_Rspr)))*(p_zR/1.e5)**0.286    ! Temp at zR [K]
q_zR =   q_0 - 1./G_L*(H_L0*(LOG((z0q+zR)/z0q) - stabIntH((z0q+zR)/L)) &
        + zR/delspr*(1. - stabIntSprayH((z0q+zR)/L))*H_Lspr)    ! Spec hum at zR [kg kg-1]
qsat0_zR = qsat0(t_zR,p_zR)    ! Saturation specific humidity at zR [kg kg-1]
betaWB_zR = 1./(1. + Lv*gammaWB*(1. + y0)/cp_a*qsat0_zR)    ! Wetbulb coefficient at zR [-]
s_zR = satratio(t_zR,p_zR,q_zR,0.99999)    ! Saturation ratio at zR [-]
tauR = rho_sw*r0**2/(rho_a*Dv_a*Fp*qsat0_zR*betaWB_zR*ABS(1. + y0 - s_zR))    ! Char timescale for evap [s]
req = r0*(xs*(1. + nu*Phi_s*Mw/Ms/(1. - s_zR)))**0.33333    ! Equilibrium radius at zR [m]
rf = req + (r0 - req)*EXP(-tauf/tauR)    ! Final droplet radius [m]
IF (ABS(1. + y0 - s_zR) < 1.e-3) THEN
    rf = r0    ! Assume no change if s_zR ~ 1 + y0
END IF

! 4. Calculate updated spray heat fluxes using new thermodynamic parameters
H_Tspr = DOT_PRODUCT(cp_sw*(t_0 - tdropf)*dmdr0,delta_r0)    ! Spray HF due to temp change [W m-2]
H_Sspr = DOT_PRODUCT(cp_sw*SIGN(MIN(ABS(t_0 - tdropf),ABS(t_0 - t_zT)),&
        t_0 - tWB_zT)*dmdr0,delta_r0)    ! Spray SHF [W m-2]
H_Rspr = DOT_PRODUCT(Lv*(1. - (rf/r0)**3)*dmdr0,delta_r0)    ! Spray HF due to size change [W m-2]
H_Lspr = H_Rspr + H_Tspr - H_Sspr    ! Spray LHF [W m-2]

END SUBROUTINE update_HFs

!===============================================================================

SUBROUTINE ssgf_BCF23(eps,swh,dcp,mss,ustar,z0,L,gf,r0,delta_r0,v_g,M_spr,dmdr0)

!-------------------------------------------------------------------------------
!
! Subroutine to calculate seastate-dependent SSGF per BCF23.  Note that the
! gustiness connected to the MO wind profile and the gust height used for spray 
! generation are not related.
!
! Inputs:
!     eps - wave energy dissipation flux [kg s-3]
!     swh - significant wave height [m]
!     dcp - dominant wave phase velocity [m s-1]
!     mss - mean squared waveslope [-]
!     ustar - local turbulence intensity [m s-1], including gustiness
!     z0 - momentum roughness length [m]
!     L - Obukhov stability length [m]
!     gf - gust factor [-]
!     r0 - SSGF radius vector [m]
!     delta_r0 - SSGF bin width vector [m]
!     v_g - droplet settling velocity vector [m s-1]
! Outputs:
!     M_spr - total spray mass flux [kg m-2 s-1]
!     dmdr0 - mass spectrum of ejected droplets [kg m-2 s-1 m-1]
!
!-------------------------------------------------------------------------------

REAL,INTENT(IN) :: eps,swh,dcp,mss,ustar,z0,L,gf
REAL,DIMENSION(25),INTENT(IN) :: r0,delta_r0,v_g
REAL,INTENT(OUT) :: M_spr
REAL,DIMENSION(25),INTENT(OUT) :: dmdr0

REAL,PARAMETER :: fs = 2.2    ! Model coefficient scaling SSGF magnitude [-].  This is the traditional 'sourcestrength'.
REAL,PARAMETER :: C1 = 1.35    ! Additional model coefficient scaling SSGF magnitude [-]
REAL,PARAMETER :: C2 = 0.1116    ! Model coefficient inside exponential [-]
REAL,PARAMETER :: C3 = 0.719    ! Model coefficient on mss [-]
REAL,PARAMETER :: C4 = 2.17    ! Model coefficient on sigma_h [-]
REAL,PARAMETER :: C5 = 0.852    ! Model coefficient in erf [-]
REAL,PARAMETER :: Ceps = 100.    ! Nondimensional mean volumetric dissipation [-], from Sutherland and Melville (2015)
REAL,PARAMETER :: Cbr = 0.8    ! Scale factor on dcp to determine breaking wave crest speed [-], from Banner et al. (2014)
REAL,PARAMETER :: alpha_k = 1.5    ! Kolmogorov constant [-]
REAL,PARAMETER :: kappa = 0.4    ! von Karman constant [-]
REAL,PARAMETER :: rho_sw = 1030.    ! Density of seawater [kg m-3]
REAL,PARAMETER :: nu_w = 0.90e-6    ! Kinematic viscosity of water [m2 s-1]
REAL,PARAMETER :: sigma_surf = 7.4e-5    ! Ratio of surface tension to water density [m3 s-2]
REAL,PARAMETER :: PI = 3.14159

REAL :: Wa,eps_KV_mean,eps_KV,eta_k,h_gust,U_h,U_10,U_crest,sigma_h
REAL,DIMENSION(25) :: dmdr0_form,ejecprob

! Background calculations
Wa = whitecapActive_DLM17(dcp,ustar/SQRT(gf),swh)    ! Actively breaking whitecap fraction [-]
eps_KV_mean = Ceps*eps/swh/rho_sw    ! Mean volumetric kinematic dissipation at surface [m2 s-3], per Sutherland and Melville 2015
eps_KV = eps_KV_mean/Wa    ! Volumetric kinematic dissipation at surface under actively breaking whitecaps [m2 s-3]
eta_k = (nu_w**3/eps_KV)**0.25    ! Kolmogorov microscale [m]
h_gust = 200.*z0    ! Gust height [m]
U_h  = ustar/kappa/gf*(LOG(h_gust/z0) - stabIntM(h_gust/L))    ! Windspeed at gust height [m s-1], gustiness removed
U_10 = ustar/kappa/gf*(LOG(10./z0)    - stabIntM(10./L))    ! 10-m windspeed [m s-1], gustiness removed
U_crest = Cbr*dcp    ! Speed of crest of breaking wave [m s-1]
sigma_h = C4*U_10    ! Standard deviation of windspeed at h_gust [m s-1]

! Calculate SSGF
dmdr0_form = (fs*C1*rho_sw*eps_KV*r0*Wa)/(3.*sigma_surf)*EXP(-1.5*alpha_k*C2*(PI*eta_k/r0)**1.33333)    ! Formation spectrum [kg m-2 s-1 m-1]
ejecprob = 0.5*(1. + ERF((U_h - U_crest - v_g/C3/mss)/sigma_h - C5))    ! Droplet ejection probability [-]
dmdr0 = dmdr0_form*ejecprob    ! Ejected droplet spectrum [kg m-2 s-1 m-1]

! Calculate total spray mass flux
M_spr = DOT_PRODUCT(dmdr0,delta_r0)    ! [kg m-2 s-1]

END SUBROUTINE ssgf_BCF23

!===============================================================================

SUBROUTINE ssgf_F94(U_10,Wform,r0,delta_r0,M_spr,dmdr0)

!-------------------------------------------------------------------------------
!
! Subroutine to calculate wind-dependent Fairall et al. (1994) SSGF.  Universal
! droplet size distribution implemented per Mueller and Veron (2014), and 
! whitecap fraction can be per either Monahan and O'Muircheartaigh (1980) or 
! BCF23.
!
! Inputs:
!     U_10 - 10-m windspeed [m s-1]
!     Wform - string naming which whitecap fraction to use.  Options are:
!         'MOM80': original F94 wind-based whitecap fraction per Monahan and 
!             O'Muircheartaigh (1980).
!         'BCF23': updated wind-based whitcap fraction per BCF23 Eq. A2.
!     r0 - SSGF radius vector [m]
!     delta_r0 - SSGF bin width vector [m]
! Outputs:
!     M_spr - total spray mass flux [kg m-2 s-1]
!     dmdr0 - mass spectrum of ejected droplets [kg m-2 s-1 m-1]
!
!-------------------------------------------------------------------------------

REAL,INTENT(IN) :: U_10
CHARACTER(LEN=100),INTENT(IN) :: Wform
REAL,DIMENSION(25),INTENT(IN) :: r0,delta_r0
REAL,INTENT(OUT) :: M_spr
REAL,DIMENSION(25),INTENT(OUT) :: dmdr0

REAL,PARAMETER :: rho_w = 1030.    ! Density of seawater [kg m-3]
REAL,PARAMETER :: PI = 3.14159
! Coefficients for A92 source function at 11 m/s, based on r80
REAL,PARAMETER :: B0 = 4.405
REAL,PARAMETER :: B1 = -2.646
REAL,PARAMETER :: B2 = -3.156
REAL,PARAMETER :: B3 = 8.902
REAL,PARAMETER :: B4 = -4.482
REAL,PARAMETER :: C1 = 1.02e4
REAL,PARAMETER :: C2 = 6.95e6
REAL,PARAMETER :: C3 = 1.75e17

INTEGER :: i
REAL :: WC_A92_11ms,WC,fs
REAL,DIMENSION(25) :: r0_micrometers,r80,dFdr80,dFdr0_11ms,dFdr0_perWC,&
        dVdr0_perWC

r0_micrometers = r0*1.e6    ! [um]
r80 = 0.518*r0_micrometers**0.976    ! Equilibrium radius at 80% RH [um], per Fitzgerald (1975)

! Define the A92 number source function at 11 m s-1, based on r80 [m-2 s-1 um-1]
DO i = 1,25
    IF (r80(i) < 0.8) THEN
        dFdr80(i) = 0.
    ELSE IF ((r80(i) >= 0.8) .AND. (r80(i) < 15.)) THEN
        dFdr80(i) = 10.**(B0 + B1*LOG10(r80(i)) &
                             + B2*LOG10(r80(i))**2 &
                             + B3*LOG10(r80(i))**3 &
                             + B4*LOG10(r80(i))**4)
    ELSE IF ((r80(i) >= 15.) .AND. (r80(i) < 37.5)) THEN
        dFdr80(i) = C1*r80(i)**(-1)
    ELSE IF ((r80(i) >= 37.5) .AND. (r80(i) < 100.)) THEN
        dFdr80(i) = C2*r80(i)**(-2.8)
    ELSE IF ((r80(i) >= 100.) .AND. (r80(i) < 250.)) THEN
        dFdr80(i) = C3*r80(i)**(-8)
    ELSE IF (r80(i) >= 250.) THEN
        dFdr80(i) = 0.
    END IF
END DO

! Convert to mass SSGF using the selected whitecap fraction.
! dFdr0_11ms used to have a bug where r0 was mistakenly used instead of r0_micrometers -- BWB 240612
dFdr0_11ms = dFdr80*0.506*r0_micrometers**(-0.024)    ! A92 source function at 11 m/s, based on r0 [m-2 s-1 um-1]
IF (Wform == 'MOM80') THEN
    WC_A92_11ms = whitecap_MOM80(11.)    ! Whitecap fraction of A92 source function at 11 m/s [-]
    WC = whitecap_MOM80(U_10)    ! Whitecap fraction at input U_10 [-]
    fs = 0.56    ! SSGF scaling factor [-].  0.56 used instead of 0.4 to correct earlier bug in dFdr0_11ms.
ELSE IF (Wform == 'BCF23') THEN
    WC_A92_11ms = whitecap_BCF23wind(11.)    ! Whitecap fraction of A92 source function at 11 m/s [-]
    WC = whitecap_BCF23wind(U_10)    ! Whitecap fraction at input U_10 [-]
    fs = 3.1    ! SSGF scaling factor [-].  3.1 used instead of 2.2 to correct earlier bug in dFdr0_11ms.
END IF
dFdr0_perWC = dFdr0_11ms/WC_A92_11ms*1.e6    ! F94 number source function per unit whitecap area [m-2 s-1 m-1]
dVdr0_perWC = dFdr0_perWC*1.33333*PI*r0**3    ! F94 volume source function per unit whitecap area [m3 m-2 s-1 m-1]
dmdr0 = fs*rho_w*dVdr0_perWC*WC    ! F94 mass source function [kg m-2 s-1 m-1]

! Calculate total spray mass flux
M_spr = DOT_PRODUCT(dmdr0,delta_r0)    ! [kg m-2 s-1]

END SUBROUTINE ssgf_F94

!===============================================================================

FUNCTION qsat0(t,p) RESULT(q_sat0)

!-------------------------------------------------------------------------------
!
! Function to calculate saturation specific humidity over a plane surface of
! pure water, with saturation vapor pressure per Buck (1981).
! 
! Inputs:
!     t - temperature [K]
!     p - pressure [Pa]
! Outputs:
!     q_sat0 - saturation specific humidity over a plane surface of 
!         pure water [kg kg-1]
!
!-------------------------------------------------------------------------------

REAL,INTENT(IN) :: t,p
REAL :: e_sat0,q_sat0

e_sat0 = 6.1121*EXP(17.502*(t - 273.15)/(t - 273.15 + 240.97))*(1.0007 + 3.46e-8*p)*1.e2    ! Sat vap press [Pa]
q_sat0 = e_sat0*0.622/(p - 0.378*e_sat0)    ! Saturation specific humidity [kg kg-1]

END FUNCTION

!===============================================================================

FUNCTION satratio(t,p,q,max_satratio) RESULT(s)

!-------------------------------------------------------------------------------
!
! Function to calculate saturation ratio, using function qsat0() that defines
! e_sat0 per Buck (1981).  We assume that s = q/qsat0, rather than w/wsat0
! (i.e., we use specific humidity rather than mixing ratio), so that air with 
! q = qsat0 will be exactly saturated.  This is done for convenience and to ease
! interpretation.
!
! Inputs:
!     t - temperature [K]
!     p - pressure [Pa]
!     q - specific humidity [kg kg-1]
!     max_satratio - maximum allowable saturation ratio [-]
! Outputs:
!     s - saturation ratio [-]
!
!-------------------------------------------------------------------------------

REAL,INTENT(IN) :: t,p,q,max_satratio
REAL :: s

s = MIN(q/qsat0(t,p),max_satratio)    ! Saturation ratio [-]

END FUNCTION

!===============================================================================

FUNCTION stabIntM(zeta) RESULT(psiM)

!-------------------------------------------------------------------------------
!
! Function to calculate the integrated stability function for momentum,
! evaluated at zeta.  Implemented per COARE 3.6 Matlab code.
!
! Inputs:
!     zeta - stability parameter [-]
! Outputs:
!     psiM - integrated stability function at zeta [-]
!
!-------------------------------------------------------------------------------

REAL,INTENT(IN) :: zeta
REAL :: psiM,Xk,psik,Xc,psic,f,am,bm,BBm,X

IF (zeta == 0.) THEN    ! Neutral
    psiM = 0.
ELSE IF (zeta < 0.) THEN    ! Unstable
    Xk = (1. - 16.*zeta)**0.25    ! For small negative zeta
    psik = 2.*LOG((1. + Xk)/2.) + LOG((1. + Xk**2)/2.) - 2.*ATAN(Xk) + 2.*ATAN(1.)
    Xc = (1. - 10.15*zeta)**0.33333    ! For large negative zeta
    psic = 1.5*LOG((1. + Xc + Xc**2)/3.) - SQRT(3.)*ATAN((1. + 2.*Xc)/SQRT(3.)) + 4.*ATAN(1.)/SQRT(3.)
    f = zeta**2/(1. + zeta**2)
    psiM = (1. - f)*psik + f*psic
ELSE IF (zeta > 0.) THEN    ! Stable (Grachev 2007, Eq. 12)
    am = 5.
    bm = am/6.5
    BBm = ((1.-bm)/bm)**0.33333
    X = (1. + zeta)**0.33333
    psiM = -(3.*am/bm)*(X-1.) + ((am*BBm)/(2.*bm))*(2.*LOG((BBm+X)/(BBm+1.)) - &
        LOG((BBm**2-BBm*X+X**2)/(BBm**2-BBm+1.)) + 2.*SQRT(3.)*ATAN((2.*X-BBm)/(BBm*SQRT(3.))) - &
        2.*SQRT(3.)*ATAN((2.-BBm)/(BBm*SQRT(3.))))
END IF

END FUNCTION

!===============================================================================

FUNCTION stabIntH(zeta) RESULT(psiH)

!-------------------------------------------------------------------------------
!
! Function to calculate the integrated stability function for heat, evaluated at
! zeta.  Implemented per COARE 3.6 Matlab code.
!
! Inputs:
!     zeta - stability parameter [-]
! Outputs:
!     psiH - integrated stability function at zeta [-]
!
!-------------------------------------------------------------------------------

REAL,INTENT(IN) :: zeta
REAL :: psiH,Yk,psik,Yc,psic,f,a,b,c,BB

IF (zeta == 0.) THEN    ! Neutral
    psiH = 0.
ELSE IF (zeta < 0.) THEN    ! Unstable
    Yk = (1. - 16.*zeta)**0.5    ! For small negative zeta
    psik = 2.*LOG((1. + Yk)/2.)
    Yc = (1. - 34.15*zeta)**0.33333    ! For large negative zeta
    psic = 1.5*LOG((1. + Yc + Yc**2)/3.) - SQRT(3.)*ATAN((1. + 2.*Yc)/SQRT(3.)) + 4.*ATAN(1.)/SQRT(3.)
    f = zeta**2/(1. + zeta**2)
    psiH = (1. - f)*psik + f*psic
ELSE IF (zeta > 0.) THEN    ! Stable (Grachev 2007, Eq. 13)
    a = 5.
    b = 5.
    c = 3.
    BB = SQRT(c**2 - 4.)
    psiH = -(b/2.)*LOG(1.+c*zeta+zeta**2) + \
        (((b*c)/(2.*BB))-(a/BB))*(LOG((2.*zeta+c-BB)/(2.*zeta+c+BB))-LOG((c-BB)/(c+BB)))
END IF

END FUNCTION

!===============================================================================

FUNCTION stabIntSprayH(zeta) RESULT(phisprH)

!-------------------------------------------------------------------------------
!
! Function to calculate the integrated stability function for heat within the
! spray layer, evaluated at zeta.  Implemented per BCF23, using the classical 
! values for the stability functions (gamma = 16, beta = 5, e.g., Dyer 1974).
!
! Inputs:
!     zeta - stability parameter [-]
! Outputs:
!     phisprH - integrated stability function within spray layer at zeta [-]
!
!-------------------------------------------------------------------------------

REAL,INTENT(IN) :: zeta
REAL :: phisprH,Y

IF (zeta == 0.) THEN    ! Neutral
    phisprH = 0.
ELSE IF (zeta > 0.) THEN    ! Stable
    phisprH = -2.5*zeta
ELSE IF (zeta < 0.) THEN    ! Unstable
    Y = (1. - 16.*zeta)**0.5
    phisprH = -(Y - 1.)**2/16./zeta
END IF

END FUNCTION

!===============================================================================

FUNCTION whitecapActive_DLM17(dcp,ustar,swh) RESULT(Wa)

!-------------------------------------------------------------------------------
!
! Function to calculate active breaking whitecap fraction per Deike, Lenain, and
! Melville (2017).  This is a parameterization based on volume of entrained air.
!
! Inputs:
!     dcp - dominant phase speed [m s-1]
!     ustar - bulk friction velocity [m s-1]
!     swh - significant wave height [m]
! Outputs:
!     Wa - actively breaking whitecap fraction [-]
!
!-------------------------------------------------------------------------------

REAL,INTENT(IN) :: dcp,ustar,swh
REAL,PARAMETER :: g = 9.81    ! Acceleration due to gravity [m s-2]
REAL :: Wa

Wa = MIN(0.018*dcp*ustar**2/g/swh,1.0)

END FUNCTION

!===============================================================================

FUNCTION whitecap_MOM80(U_10) RESULT(W)

!-------------------------------------------------------------------------------
!
! Function to calculate whitecap fraction per Monahan and O'Muircheartaigh 
! (1980).
!
! Inputs:
!     U_10 - 10-m windspeed [m s-1]
! Outputs:
!     W - whitecap fraction [-]
!
!-------------------------------------------------------------------------------

REAL,INTENT(IN) :: U_10
REAL :: W

W = MIN(3.8e-6*U_10**3.4,1.0)

END FUNCTION

!===============================================================================

FUNCTION whitecap_BCF23wind(U_10) RESULT(W)

!-------------------------------------------------------------------------------
!
! Function to calculate wind-based whitecap fraction given as BCF23 Eq A2.
!
! Inputs:
!     U_10 - 10-m windspeed [m s-1]
! Outputs:
!     W - whitecap fraction [-]
!
!-------------------------------------------------------------------------------

REAL,INTENT(IN) :: U_10
REAL :: W

W = MIN(6.5e-4*MAX(U_10 - 2.,0.)**1.5,1.0)

END FUNCTION

!===============================================================================

SUBROUTINE waveProps(U_10,ustarb,rho_a,eps,dcp,swh,mss)

!-------------------------------------------------------------------------------
!
! This is a temporary hack to roughly estimate wave properties based on 
! atmospheric variables.  It will be replaced by directly parameterizing the 
! seastate-based SSGF from winds.
!
! Inputs:
!     U_10 - 10-m windspeed [m s-1]
!     ustarb - bulk friction velocity [m s-1]
!     rho_a - air density [kg m-3]
! Outputs:
!     eps - wave energy dissipation flux [W m-2]
!     dcp - dominant wave phase speed [m s-1]
!     swh - significant wave height [m]
!     mss - mean squared waveslope [-]
!
!-------------------------------------------------------------------------------

REAL,INTENT(IN) :: U_10,ustarb,rho_a
REAL,INTENT(INOUT) :: eps,dcp,swh,mss
REAL,PARAMETER :: g = 9.81    ! Acceleration due to gravity [m s-2]
REAL,PARAMETER :: fudge1 = 0.2    ! Fudge factor on eps
REAL,PARAMETER :: fudge2 = 0.7    ! Fudge factor on mss
REAL,PARAMETER :: PI = 3.14159
REAL :: swh_W17,swh_tail

! 1. Significant wave height -- Wang et al. (2017) with hacked tail
swh_W17 = 0.0143*U_10**2 + 0.9626    ! swh per Wang et al. 2017 [m]
swh_tail = 7. + 0.14*U_10    ! Hacked linear tail [m]
swh = (1. - 0.5*(TANH((U_10 - 27.)/7.) + 1.))*swh_W17 &
         + (0.5*(TANH((U_10 - 29.)/7.) + 1.))*swh_tail    ! Smooth merge between two params

! 2. Dominant phase speed -- Wang et al. (2017)
dcp = (g*swh/(0.0628*(2.*PI)**1.5*ustarb**0.5))**0.666667    ! [m s-1]

! 3. Wave energy dissipation flux -- Terray et al. (1996) with fudge factor
eps = rho_a*ustarb**2*0.5*dcp*fudge1    ! [W m-2], factor of 0.5 from T96 Fig 6

! 4. Mean squared waveslope -- Davis et al. (2023) with fudge factor
mss = 0.109*TANH(0.057*U_10)*fudge2    ! [-]

END SUBROUTINE waveProps

!===============================================================================

FUNCTION fall_velocity_PK97(r) RESULT(v_fall)

!-------------------------------------------------------------------------------
!
! Function to calculate terminal fall velocity of spherical droplets, based on 
! Pruppacher and Klett (1997) section 10.3.6.
!
! Inputs:
!     r - droplet radius [m]
! Outputs:
!     v_fall - terminal settling velocity [m s-1]
!
!-------------------------------------------------------------------------------

REAL,DIMENSION(25),INTENT(IN) :: r
REAL,DIMENSION(25) :: v_fall
REAL,PARAMETER :: nu_a = 1.5e-5    ! Kinematic viscosity of air [m2 s-1]
REAL,PARAMETER :: rho_a = 1.25    ! Density of air [kg m-3]
REAL,PARAMETER :: rho_w = 1030.    ! Density of seawater [kg m-3]
REAL,PARAMETER :: g = 9.81    ! Acceleration due to gravity [m s-2]
REAL,PARAMETER :: sigma_aw = 7.4e-2    ! Surface tension of air-water interface [N m-1]
REAL,PARAMETER :: lamda_a0 = 6.6e-8    ! Mean free path at 1013.25 mb, 293.15 K [m]
REAL,DIMENSION(7) :: B1 = (/ -0.318657e1, 0.992696, -0.153193e-2, -0.987059e-3, &
        -0.578878e-3, 0.855176e-4, -0.327815e-5 /)    ! Polynomial coeffs for 10 <= r <= 535 um [-]
REAL,DIMENSION(6) :: B2 = (/ -0.500015e1, 0.523778e1, -0.204914e1, 0.475294, &
        -0.542819e-1, 0.238449e-2 /)    ! Polynomial coefficients for r > 535 um [-]
INTEGER :: i
REAL :: v_stokes,f_slip,CdNRe2,X,Y,NRe,NBo,NP,NBoNP16

DO i = 1,25
    IF (r(i) < 10.e-6) THEN    ! For r < 10 micrometers
        v_stokes = 2.*r(i)**2*g*(rho_w - rho_a)/9./(rho_a*nu_a)    ! Stokes velocity [m s-1]
        f_slip = 1. + 1.26*lamda_a0/r(i)    ! Slip flow Cunningham correction factor [-]
        v_fall(i) = f_slip*v_stokes    ! Settling velocity [m s-1]
    ELSE IF ((r(i) >= 10.e-6) .AND. (r(i) <= 535.e-6)) THEN    ! For r from 10 to 535 micrometers
        CdNRe2 = 32.*r(i)**3*(rho_w - rho_a)/rho_a/nu_a**2*g/3.    ! Product of Cd and Re**2 [-]
        X = LOG(CdNRe2)    ! X in polynomial curve fit [-]
        Y = B1(1) + B1(2)*X + B1(3)*X**2 + B1(4)*X**3 + B1(5)*X**4 + &
            B1(6)*X**5 + B1(7)*X**6    ! Y in polynomial curve fit [-]
        NRe = EXP(Y)    ! Reynolds number [-]
        v_fall(i) = nu_a*NRe/2./r(i)    ! Settling velocity [m s-1]
    ELSE IF (r(i) > 535.e-6) THEN    ! For r > 535 micrometers
        NBo = g*(rho_w - rho_a)*r(i)**2/sigma_aw    ! Bond number [-]
        NP = sigma_aw**3/rho_a**2/nu_a**4/g/(rho_w - rho_a)    ! Physical property number [-]
        NBoNP16 = NBo*NP**0.16666    ! Product of Bond number and physical property number to the 1/6 power [-]
        X = LOG(16./3.*NBoNP16)    ! X in polynomial curve fit [-]
        Y = B2(1) + B2(2)*X + B2(3)*X**2 + B2(4)*X**3 + B2(5)*X**4 + &
            B2(6)*X**5    ! Y in polynomial curve fit [-]
        NRe = NP**0.16666*EXP(Y)    ! Reynolds number [-]
        v_fall(i) = nu_a*NRe/2./r(i)    ! Settling velocity [m s-1]
    END IF
END DO
    
END FUNCTION

!===============================================================================

END MODULE module_sprayHFs
