!=================================================================================================================
 module module_mp_cldfra_cb

 implicit none
 private
 public:: calc_cldfraction_monan


 contains


!=================================================================================================================
      subroutine calc_cldfraction_monan(CLDFRA, qv, qc, qi, qs,         &
     &                 p,t,rho, XLAND, gridkm,zgrid,kpbl,gf_mfx,        &
     &                 ids,ide, jds,jde, kds,kde,                       &
     &                 ims,ime, jms,jme, kms,kme,                       &
     &                 its,ite, jts,jte, kts,kte                        )

      !!    PURPOSE
      !!    -------
      !!**  Routine to diagnose cloud fraction and liquid and ice condensate mixing ratios
      !!**  METHOD
      !!    ------
      !!    Based on the large-scale fields of temperature, water vapor, and possibly
      !!    liquid and solid condensate, the conserved quantities r_t and h_l are constructed
      !!    and then fractional cloudiness, liquid and solid condensate is diagnosed.
      !!
      !!    The total variance is parameterized as the sum of  stratiform/turbulent variance
      !!    and a convective variance.
      !!    The turbulent variance is parameterized as a function of first-order moments, and
      !!    the convective variance is modelled as a function of the convective mass flux (units kg/s m^2)
      !!    as provided by a mass flux convection scheme.
      !!
      !!    Nota: if the host model does not use prognostic values for liquid and solid condensate
      !!    or does not provide a convective mass flux, put all these values to zero.
      !!    Also, it is supposed that vertical model levels are numbered from
      !!    1 to MZP, where 1 is the first model level above the surface
      !!
      !!    ------------------
      !!    REFERENCE
      !!    ---------
      !!      Chaboureau J.P. and P. Bechtold (J. Atmos. Sci. 2002)
      !!      Chaboureau J.P. and P. Bechtold (JGR/AGU 2005)
      !!
      !!    AUTHOR
      !!    ------
      !!      P. BECHTOLD       * Laboratoire d'Aerologie *
      !!
      !!    MODIFICATIONS
      !!    -------------
      !!      Original    13/06/2001
      !!      modified    20/03/2002 : add convective Sigma_s and improve turbulent
      !!                               length-scale in boundary-layer and near tropopause
      !!      adapted     09/12/2016 : adapted to GEOS-5 by Saulo Freitas
      !!      adapted     24/04/2024 : adapted to MPAS/MONAN by Saulo Freitas 
      !-------------------------------------------------------------------------------
      !*       0.    DECLARATIONS
      !              ------------
      implicit none
      !
      !-------------------------------------------------------------------------------
      !
      !*       1.    Set the fundamental thermodynamical constants
      !              these have the same values (not names) as in ARPEGE IFS
      !              -------------------------------------------------------
      real, parameter :: XP00   = 1.e5        ! reference pressure
      real, parameter :: XPI    = 3.141592654 ! C_pi
      real, parameter ::  XG    = 9.80665     ! gravity constant
      real, parameter :: XMD    = 28.9644e-3  ! molecular weight of dry air
      real, parameter :: XMV    = 18.0153e-3  ! molecular weight of water vapor
      real, parameter :: XRD    = 287.05967   ! gaz constant for dry air
      real, parameter :: XRV    = 461.524993  ! gaz constant for water vapor
      real, parameter :: XCPD   = 1004.708845 ! specific heat of dry air
      real, parameter :: XCPV   = 1846.1      ! specific heat of water vapor
      real, parameter :: XRHOLW = 1000.       ! density of liquid water
      real, parameter :: XCL    = 4218.       ! specific heat of liquid water
      real, parameter :: XCI    = 2106.       ! specific heat of ice
      real, parameter :: XTT    = 273.16      ! triple point temperature
      real, parameter :: C_ALVLTT  = 2.5008e6 ! latent heat of vaporisation at XTT
      real, parameter :: XLSTT  = 2.8345e6    ! latent heat of sublimation at XTT
      real, parameter :: XLMTT  = 0.3337e6    ! latent heat of melting at XTT
      real, parameter :: XESTT  = 611.14      ! saturation pressure at XTT
      real, parameter :: XALPW  = 60.22416    ! constants in saturation pressure over liquid water
      real, parameter :: XBETAW = 6822.459384
      real, parameter :: XGAMW  = 5.13948
      real, parameter :: XALPI  = 32.62116    ! constants in saturation pressure over ice
      real, parameter :: XBETAI = 6295.421
      real, parameter :: XGAMI  = 0.56313
      logical, parameter :: LUSERI = .true. ! logical switch to compute both
                                            ! liquid and solid condensate (LUSERI=.TRUE.)
                                            ! or only liquid condensate (LUSERI=.FALSE.)
      !
      !*       0.1   Declarations of input variables :
      !
      !
      integer, intent(in):: ids,ide, jds,jde, kds,kde,                  &
     &                      ims,ime, jms,jme, kms,kme,                  &
     &                      its,ite, jts,jte, kts,kte
      integer, dimension(ims:ime,jms:jme),      intent(in):: kpbl  !index of PBL top
      real, dimension(ims:ime,kms:kme,jms:jme), intent(in):: qv,p,t,rho
      real, dimension(ims:ime,kms:kme,jms:jme), intent(in):: qc,qi,qs,zgrid,gf_mfx
      real, dimension(ims:ime,jms:jme),         intent(in):: xland
      real, dimension(ims:ime,jms:jme),         intent(in):: gridkm
      ! fractional cloudiness (between 0 and 1)
      real, dimension(ims:ime,kms:kme,jms:jme), intent(inout):: cldfra  

      !
      !*       0.1   Declarations of local variables :
      !
      !
      real, dimension(its:ite,kts:kte) :: PPABS  ! pressure (Pa)
      real, dimension(its:ite,kts:kte) :: PZZ    ! height of model levels (m)
      real, dimension(its:ite,kts:kte) :: PT     ! grid scale T  (K)
      real, dimension(its:ite,kts:kte) :: PRV    ! grid scale water vapor mixing ratio (kg/kg)
      real, dimension(its:ite,kts:kte) :: PMFLX  ! convective mass flux (kg/(s m^2))
      real, dimension(its:ite,kts:kte) :: PRC    ! grid scale r_c mixing ratio (kg/kg)
      real, dimension(its:ite,kts:kte) :: PRI    ! grid scale r_i (kg/kg)
      !
      !
      integer  ::  JKT, JKP, JKM,K,i,j     ! loop index
      real,    dimension(its:ite,kts:kte) :: ZTLK, ZRT       ! work arrays for T_l, r_t
      real,    dimension(its:ite,kts:kte) :: ZL              ! length-scale
      integer, dimension(its:ite)   :: ITPL    ! top levels of tropopause/highest inversion
      real,    dimension(its:ite)   :: ZTMIN   ! min Temp. related to ITPL
      !
      real :: ZTEMP, ZLV, ZLS, ZTL, ZPV, ZQSL, ZPIV, ZQSI, ZFRAC, ZCOND, ZCPD ! thermodynamics
      real :: ZLL, DZZ, ZZZ ! length scales
      real :: ZAH, ZA, ZB, ZSBAR, ZQ1, ZSIGMA, ZDRW, ZDTL ! related to computation of Sig_s
      real :: ZSIG_CONV,  ZSIGMA_NOCONV,  ZQ1_NOCONV      ! convective part of Sig_s
      real, parameter :: zq1_tuning = 0.0 ! -0.25 ! 1.0
      !
      !*       0.2  Definition of constants :
      !
      !-------------------------------------------------------------------------------
      !
      real :: ZL0     = 600.        ! tropospheric length scale
                                    ! changed to 600 m instead of 900 m to give a consistent
                                    ! value (linear increase) in general 500 m deep oceanic
                                    ! mixed layer - but could be put back to 900 m if wished
      real :: ZCSIGMA = 0.2         ! constant in sigma_s parameterization
      real :: ZCSIG_CONV = 0.30e-2  ! scaling factor for ZSIG_CONV as function of mass flux
      !
      !
      !logical :: ONLY_CONVECTIVE_CLOUD_FRACTION=.true. ! set .false. for the total cloud fraction
      !-------------------------------------------------------------------------------
      !
      do j = jts,jte
      
      do k = kts,kte
      do i = its,ite
         cldfra(i,k,j) = 0.0
         prv   (i,k)   = qv(i,k,j)
         prc   (i,k)   = qc(i,k,j)
         pri   (i,k)   = qi(i,k,j)
         pt    (i,k)   =  t(i,k,j)
         ppabs (i,k)   =  p(i,k,j)
         pzz   (i,k)   =  zgrid(i,k,j)-zgrid(i,1,j)
         pmflx (i,k)   =  gf_mfx(i,k,j)
      enddo
      enddo

     
      !JKT = KTE+1-KTS
      !-will limit the model vertical column to 60 hPa
      ! do K=KTE,KTS,-1
      !    if(PPABS(i,k) > 60.*100.) then
      !       JKT = k
      !       exit
      !    endif
      ! enddo

      DO k = kts,kte
      DO i = its,ite

         ZTEMP  = PT(i,k)
         !latent heat of vaporisation/sublimation
         ZLV    = C_ALVLTT + ( XCPV - XCL ) * ( ZTEMP - XTT )
         ZLS    = XLSTT    + ( XCPV - XCI ) * ( ZTEMP - XTT )

         !store temperature at saturation and total water mixing ratio
         ZRT(i,k)   = PRV(i,k) + PRC(i,k) + PRI(i,k)
         ZCPD       = XCPD  + XCPV*PRV(i,k) + XCL*PRC(i,k) + XCI*PRI(i,k)
         ZTLK(i,k)  = ZTEMP - ZLV*PRC(i,k)/ZCPD - ZLS*PRI(i,k)/ZCPD
      enddo
      enddo

      !-------------------------------------------------------------------------------
      ! Determine tropopause/inversion  height from minimum temperature

      ITPL (:) = KTS
      ZTMIN(:) = 400.
      do k = KTS+1,KTE-1
      do i = its,ite
         if ( PT(i,k) < ZTMIN(i) ) then
            ZTMIN(i) = PT(i,k)
            ITPL(i)  = K
         endif
      enddo
      enddo

      ! Set the mixing length scale - used for computing the "turbulent part" of Sigma_s

      ZL(:,kts) = 20.
      do k = KTS+1,KTE
          do i = its,ite

         ! free troposphere
         ZL(i,k) = ZL0
         JKP   = ITPL(i)
         ZZZ   = PZZ(i,k) -  PZZ(i,KTS)
            ! approximate length for boundary-layer : linear increase
         if ( ZL0 > ZZZ )  ZL(i,k) = ZZZ
            ! gradual decrease of length-scale near and above tropopause/top inversion
         if ( ZZZ > 0.9*(PZZ(i,JKP)-PZZ(i,KTS)) ) &
            ZL(i,k) = .6 * ZL(i,K-1)
       enddo
      enddo
     !-------------------------------------------------------------------------------
      do k = KTS+1,KTE-1
         JKP=k+1
         JKM=k-1
         do i = its,ite
            !--only allowed above the boundary layer
            if(k <= kpbl(i,j)) cycle
            
            ZTEMP  = PT(i,k)

            !latent heat of vaporisation/sublimation
            ZLV    = C_ALVLTT + ( XCPV - XCL ) * ( ZTEMP - XTT )
            ZLS    = XLSTT    + ( XCPV - XCI ) * ( ZTEMP - XTT )

            ZCPD   = XCPD + XCPV*PRV(i,k) + XCL*PRC(i,k) + XCI*PRI(i,k)
            !temperature at saturation
            ZTL    = ZTEMP - ZLV*PRC(i,k)/ZCPD - ZLS*PRI(i,k)/ZCPD

            !saturated water vapor mixing ratio over liquid water
            ZPV    = MIN(EXP( XALPW - XBETAW / ZTL - XGAMW * LOG( ZTL ) ),0.99*PPABS(i,k))
            ZQSL   = XRD / XRV * ZPV / ( PPABS(i,k) - ZPV )

            !saturated water vapor mixing ratio over ice
            ZPIV   = MIN(EXP( XALPI - XBETAI / ZTL - XGAMI * LOG( ZTL ) ),0.99*PPABS(i,k))
            ZQSI   = XRD / XRV * ZPIV / ( PPABS(i,k) - ZPIV )

            !interpolate between liquid and solid as function of temperature
            ! glaciation interval is specified here to 20 K
            ZFRAC = ( ZTL  - 250.16 ) / ( XTT - 250.16 )  ! liquid/solid fraction
            ZFRAC = MAX( 0., MIN(1., ZFRAC ) )

            if(.not. LUSERI) ZFRAC=1.
            ZQSL = ( 1. - ZFRAC ) * ZQSI + ZFRAC * ZQSL
            ZLV  = ( 1. - ZFRAC ) * ZLS  + ZFRAC * ZLV

            !coefficients a and b
            ZAH  = ZLV * ZQSL / ( XRV * ZTL**2 ) * (XRV * ZQSL / XRD + 1.)
            !orig  ZAH  = ZLV * ZQSL / ( XRV * ZTL**2 )

            ZA  = 1. / ( 1. + ZLV/ZCPD * ZAH )
            ZB  = ZAH * ZA

            !- parameterize Sigma_s with first_order closure
            DZZ  =  PZZ (i,JKP)  - PZZ (i,JKM)
            ZDRW =  ZRT (i,JKP)  - ZRT (i,JKM)
            ZDTL =  ZTLK(i,JKP)  - ZTLK(i,JKM) + XG/ZCPD * DZZ
            ZLL  =  ZL(i,k)

            !- standard deviation due to convection
            ZSIG_CONV = ZCSIG_CONV * PMFLX(i,k) / ZA

            !- turb + conv
            ZSIGMA = SQRT( MAX( 1.e-25, ZCSIGMA*ZCSIGMA* ZLL*ZLL/(DZZ*DZZ) * ( &
                           ZA*ZA*ZDRW*ZDRW - 2.*ZA*ZB*ZDRW*ZDTL                &
                         + ZB*ZB*ZDTL*ZDTL                                   ) &
                         + ZSIG_CONV * ZSIG_CONV )                           )

            !- zsigma should be of order 4.e-4 in lowest 5 km of atmosphere
            ZSIGMA = MAX( ZSIGMA, 1.e-10 )

            !- normalized saturation deficit
            ZSBAR = ZA * ( ZRT (i,k) - ZQSL )
            !- "Q1" parameter
            ZQ1   = ZSBAR / ZSIGMA - zq1_tuning

            !- total cloud fraction
            CLDFRA(I,K,J) = MAX( 0., MIN(1.,0.5+0.36*ATAN(1.55*ZQ1)) )

      !    if(ONLY_CONVECTIVE_CLOUD_FRACTION) then
      !       !- get cloud fraction associated with ONLY the sub-grid scale convective part
      !       !- this sigma does not include the sub-grid scale convective part
      !       ZSIGMA_NOCONV = SQRT( MAX( 1.e-25, ZCSIGMA*ZCSIGMA* ZLL*ZLL/(DZZ*DZZ) * ( &
      !          ZA*ZA*ZDRW*ZDRW - 2.*ZA*ZB*ZDRW*ZDTL   &
      !          + ZB*ZB*ZDTL*ZDTL  )))
      !       !- zsigma should be of order 4.e-4 in lowest 5 km of atmosphere
      !       ZSIGMA_NOCONV = MAX( ZSIGMA_NOCONV, 1.e-10 )
      !       ZQ1_NOCONV = ZSBAR / ZSIGMA_NOCONV
      !
      !       !- cloud fraction associated with ONLY convective part ("total-turb")
      !       PCLDFR(i,k) = 0.36*(ATAN(1.55*ZQ1)-ATAN(1.55*ZQ1_NOCONV))
      ! 
      !       PCLDFR(i,k) = MAX( 0., MIN(1.,PCLDFR(i,k)) )
      !
      !    endif
      !
      !    cycle
      !    !total condensate diagnostic (not being used)
      !    if (ZQ1 > 0. .and. ZQ1 <= 2. ) then
      !       !orig   ZCOND =     EXP(-1.)+.66*ZQ1+.086*ZQ1*ZQ1
      !       ZCOND = MIN(EXP(-1.)+.66*ZQ1+.086*ZQ1**2, 2.) ! We use the MIN function for continuity
      !    else if (ZQ1 > 2.) then
      !       ZCOND = ZQ1
      !    else
      !       ZCOND = EXP( 1.2*ZQ1-1. )
      !    end if
      !    ZCOND = ZCOND * ZSIGMA
      !
      !    if ( zcond < 1.e-12) then
      !       zcond = 0.
      !       pcldfr(i,k) = 0.
      !    end if
      !    if ( pcldfr(i,k) == 0.) then
      !       zcond = 0.
      !    end if
      !
      !    PRC(i,k) = ZFRAC * ZCOND ! liquid condensate
      !    if (LUSERI) then
      !       PRI(i,k) = (1.-ZFRAC) * ZCOND   ! solid condensate
      !    end if
      !
      ! !---
      ! ! compute s'rl'/Sigs^2
      ! ! used in w'rl'= w's' * s'rl'/Sigs^2
      ! !  PSIGRC(i,k) = PCLDFR(i,k)   ! Gaussian relation
      ! !
      ! ! s r_c/ sig_s^2
      ! !    PSIGRC(JI,JJ,JK) = PCLDFR(JI,JJ,JK)  ! use simple Gaussian relation
      ! !
      ! !    multiply PSRCS by the lambda3 coefficient
      ! !
      ! !      PSIGRC(JI,JJ,JK) = 2.*PCLDFR(JI,JJ,JK) * MIN( 3. , MAX(1.,1.-ZQ1) )
      ! ! in the 3D case lambda_3 = 1.
      ! !      INQ1 = MIN( MAX(-22,FLOOR(2*ZQ1) ), 10)
      ! !      ZINC = 2.*ZQ1 - INQ1
      ! !
      ! !      PSIGRC(i,k) =  MIN(1.,(1.-ZINC)*ZSRC_1D(INQ1)+ZINC*ZSRC_1D(INQ1+1))
      ! !
      ! !      PSIGRC(i,k) = PSIGRC(i,k)* MIN( 3. , MAX(1.,1.-ZQ1) )
      ! !---
       ENDDO
       ENDDO
      ENDDO
    end subroutine calc_cldfraction_monan
!=================================================================================================================
  end module module_mp_cldfra_cb
!=================================================================================================================
