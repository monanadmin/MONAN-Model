MODULE VerticalInterpolation
    use mpas_kind_types

  IMPLICIT NONE
   PRIVATE
                               ! Selecting Kinds:
                               ! r4 : Kind for 32-bits Real Numbers
                               ! i4 : Kind for 32-bits Integer Numbers
                               ! RKIND : Kind for 64-bits Real Numbers
   !INTEGER, PARAMETER, PUBLIC :: r4RKIND  = SELECTED_REAL_KIND(6)
   !INTEGER, PARAMETER, PUBLIC :: i4 = SELECTED_INT_KIND(9)
   !INTEGER, PARAMETER, PUBLIC :: RKIND = SELECTED_REAL_KIND(15)
   !INTEGER, PARAMETER, PUBLIC :: r16 = SELECTED_REAL_KIND(15)
   INTEGER, PARAMETER          :: nferr=0

   LOGICAL, PUBLIC :: ExtrapoAdiabatica=.true. ! .TRUE. = extrapolated adiabatically to the pressure levels 
   
   REAL (KIND=RKIND), PARAMETER :: a =17.2693882_RKIND  ! Constant Used in Tetens Formulae for Es
   REAL (KIND=RKIND), PARAMETER :: b =35.86_RKIND       ! Constant Used in Tetens Formulae for Es
   REAL (KIND=RKIND), PARAMETER :: To =273.16_RKIND     ! Constant Used in Tetens Formulae for Es (K)
   REAL (KIND=RKIND), PARAMETER :: Eo =6.1078_RKIND     ! Constant Used in Tetens Formulae for Es (Mb)

   REAL (KIND=RKIND), PARAMETER :: Po = 1000.0_RKIND            ! Reference Surface Pressure Value (Mb)
   REAL (KIND=RKIND), PARAMETER :: Pt =   0.1_RKIND             ! Reference Top Pressure Value (Mb)
   REAL (KIND=RKIND), PARAMETER :: Avogrado=6022.52_RKIND       ! Avogrado's Constant (10**-20/kmol)
   REAL (KIND=RKIND), PARAMETER :: Boltzmann=1.38054_RKIND      ! Boltzmann's Constant (J/K)
   REAL (KIND=RKIND), PARAMETER :: REstar=Avogrado*Boltzmann ! Universal Gas Constant (10**23 J/K/kmol)

   REAL (KIND=RKIND), PARAMETER :: MWWater=18.016_RKIND    ! Molecular Weight of Water (kg/kmol)
   REAL (KIND=RKIND), PARAMETER :: MWN2    =28.0134_RKIND   ! Molecular Weight of Nitrogen (kg/kmol)
   REAL (KIND=RKIND), PARAMETER :: PN2     =0.7809_RKIND    ! Atmosphere Nitrogen Percetage
   REAL (KIND=RKIND), PARAMETER :: MWO2    =31.9988_RKIND   ! Molecular Weight of Oxigen (kg/kmol)
   REAL (KIND=RKIND), PARAMETER :: PO2     =0.2095_RKIND    ! Atmosphere Oxigen Percetage
   REAL (KIND=RKIND), PARAMETER :: MWAr    =39.9480_RKIND   ! Molecular Weight of Argon (kg/kmol)
   REAL (KIND=RKIND), PARAMETER :: PAr     =0.0093_RKIND    ! Atmosphere Argon Percetage
   REAL (KIND=RKIND), PARAMETER :: MWCO2   =44.0103_RKIND   ! Molecular Weight of Carbon Gas (kg/kmol)
   REAL (KIND=RKIND), PARAMETER :: PCO2    =0.0003_RKIND    ! Atmosphere Carbon Gas Percetage
   REAL (KIND=RKIND), PARAMETER :: MWDAir  =PN2*MWN2+PO2*MWO2+PAr*MWAr+PCO2*MWCO2 ! Mean Molecular Weight of Dry Air (kg/kmol)
   REAL (KIND=RKIND), PARAMETER :: Rv      =REstar/MWWater     ! Water Vapor Gas Constant (m2/s2/K)
   REAL (KIND=RKIND), PARAMETER :: Rd      =REstar/MWDAir ! Dry Air Gas Constant (m2/s2/K)
   REAL (KIND=RKIND), PARAMETER :: Cp      =3.5_RKIND*Rd ! Dry Air Gas Specific Heat Capcity at Pressure Constant (m2/s2/K)
   REAL (KIND=RKIND), PARAMETER :: RdByCp  =Rd/Cp       ! Rd by Cp
   REAL (KIND=RKIND), PARAMETER :: TVVTa   =0.0065_RKIND ! Atmosphere Lapse Rate a
   REAL (KIND=RKIND), PARAMETER :: TVVTb   =0.0050_RKIND ! Atmosphere Lapse Rate b
   REAL (KIND=RKIND), PARAMETER :: Tref    =290.66_RKIND  ! Temperature Reference Value (K)
   REAL (KIND=RKIND), PARAMETER :: Zref    =75.0_RKIND    ! Geopotential Height Reference Value (m)
   REAL (KIND=RKIND), PARAMETER :: Grav    =9.80665_RKIND  ! Earth Gravity Aceleration (m/s2)
   REAL (KIND=RKIND), PARAMETER :: RdByGrav=Rd/Grav ! Rd by Grav
   REAL (KIND=RKIND), PARAMETER :: CpByRd  =Cp/Rd ! Cp by Rd
   REAL (KIND=RKIND), PARAMETER :: Eps     =Rd/Rv    ! Rd by Rv
   REAL (KIND=RKIND), PARAMETER :: Eps1    =1.0_RKIND-Eps ! 1-Eps
   REAL (KIND=RKIND), PARAMETER :: PRHcut  =2.0_RKIND ! Pressure Level Cut for Specific Humidity (Mb)
   REAL (KIND=RKIND), PARAMETER :: RHmin   =0.00000001_RKIND ! Minimum Value for Relative Humidity
   REAL (KIND=RKIND), PARAMETER :: RHmax   =0.99999999_RKIND ! Maximum Value for Relative Humidity
   REAL (KIND=RKIND), PARAMETER :: SHmin   =1.0e-7_RKIND     ! Minimum Value for Specific Humidity (kg/kg)


   PUBLIC :: VerticalInterpolation_ZGEO
   PUBLIC :: VerticalInterpolation_Thermal
   PUBLIC :: VerticalInterpolation_Humid
   PUBLIC :: VerticalInterpolation_Thermal2
   PUBLIC :: VerticalInterpolation_Wind
   PUBLIC :: VerticalInterpolation_Single

CONTAINS

  SUBROUTINE VerticalInterpolation_Thermal2(nCells,nVertLevels,Lmax,pressure,temperature, tg, pmand)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nCells
    INTEGER, INTENT(IN) :: nVertLevels
    INTEGER, INTENT(IN) :: Lmax

    REAL (KIND=RKIND), INTENT(IN ) :: pressure   (nVertLevels,nCells)
    REAL (KIND=RKIND), INTENT(IN ) :: temperature(nVertLevels,nCells)

    REAL (KIND=RKIND), INTENT(OUT) :: tg         (Lmax,nCells)
    REAL (KIND=RKIND), INTENT(IN ) :: pmand (Lmax)

    INTEGER :: i, j, k
    INTEGER :: kp
    REAL (KIND=RKIND) :: alnpmd(Lmax)
    REAL (KIND=RKIND) :: pref(Lmax)

    INTEGER           :: ks   (Lmax,nCells)

    REAL (KIND=RKIND) :: ddp  (Lmax,nCells)
 
    REAL (KIND=RKIND) :: pim  (nVertLevels  ,nCells)
    REAL (KIND=RKIND) :: alnpm(nVertLevels  ,nCells)
    REAL (KIND=RKIND) :: alnp (nVertLevels+1,nCells)
    REAL (KIND=RKIND) :: pi   (nVertLevels+1,nCells)
    REAL (KIND=RKIND) :: p    (nVertLevels+1,nCells)
    REAL (KIND=RKIND) :: rdelp(nVertLevels-1,nCells)
    REAL (KIND=RKIND) :: bt   (nVertLevels-1,nCells)

    REAL (KIND=RKIND) :: wtl

    wtl=LOG(Po)
    
    ! compute pressure on each sigma surface 

       DO k=1,nVertLevels
          DO i=1,nCells
            ! transform from Pa to mbar
             p(k,i)= 1.e-2_RKIND * pressure(k,i)
          END DO
       END DO    
       ! set the top of the model at Pt (See Constants).
       p(nVertLevels+1,:)=Pt


       DO kp=1,Lmax
            ! transform from Pa to mbar
          pref(kp) = 1.e-2_RKIND * pmand (kp)
       END DO

    ! log of pressure and exner function on each sigma surface

       DO k=1,nVertLevels+1
          DO i=1,nCells
             alnp(k,i)=LOG(MAX(p(k,i),0.00000001_RKIND))
             pi  (k,i)=EXP(RdByCp*(alnp(k,i)-wtl))
          END DO
       END DO

       DO kp=1,Lmax
          alnpmd(kp) = LOG(MAX(pref(kp),0.00000001_RKIND))
       END DO
    ! mean value of the exner function, log of pim and
    ! the log of the mean value of pressure in each layer

       DO k=1,nVertLevels
          DO i=1,nCells
             pim(k,i)=(pi(k+1,i)-pi(k,i))/ &
                  (RdByCp*(alnp(k+1,i)-alnp(k,i)))
             alnpm(k,i)=wtl+CpByRd*LOG(pim(k,i))
          END DO
       END DO

    ! return interpolated t absolute. ts replaced by temperature

       DO k=1,nVertLevels-1
          DO i=1,nCells
             rdelp(k,i)=1.0_RKIND/(alnpm(k+1,i)-alnpm(k,i))
             bt(k,i)=(temperature(k+1,i)-temperature(k,i))*rdelp(k,i)
          END DO
       END DO

    ! heights to pressure surfaces  -  hydrostatic interpolation
    ! winds and temps by linear interpolation with ln(p)
    ! search for middle of sigma layer above kp

       ks(:,:)=nVertLevels
       DO kp=1,Lmax
          DO k=1,nVertLevels-1
             DO i=1,nCells
                IF ((ks(kp,i) == nVertLevels) .AND. (alnpmd(kp) > alnpm(k,i))) THEN
                   ks(kp,i)=k
                END IF
             END DO
          END DO
       END DO

    ! find values and slopes for upward or downward extrapolation
    ! as well as for interpolation away from tropopause
    ! ts to temperature for absolute t

       DO kp=1,Lmax
          DO i=1,nCells
             IF (ks(kp,i) /= 1) THEN
                ddp(kp,i)=alnpmd(kp)-alnpm(ks(kp,i),i)
                tg (kp,i)=temperature(ks(kp,i),i)+bt(ks(kp,i)-1,i)*ddp(kp,i)
             ELSE
                IF(.NOT.ExtrapoAdiabatica)THEN
                   tg(kp,i)=temperature(1,i)
                ELSE
                   tg(kp,i)=temperature (1,i)*EXP(0.5_RKIND*RdByCp*ABS(alnpm(ks(kp,i),i)-alnpmd(kp)))  
                END IF
             END IF
          END DO
       END DO

  END SUBROUTINE VerticalInterpolation_Thermal2

  SUBROUTINE VerticalInterpolation_Humid(nCells,nVertLevels,Lmax,pressure, gsh,rq, pmand)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nCells
    INTEGER, INTENT(IN) :: nVertLevels
    INTEGER, INTENT(IN) :: Lmax

    REAL (KIND=RKIND), INTENT(IN ) :: pressure   (nVertLevels,nCells)
    REAL (KIND=RKIND), INTENT(IN ) :: gsh        (nVertLevels,nCells)
    REAL (KIND=RKIND), INTENT(OUT) :: rq         (Lmax,nCells)
    REAL (KIND=RKIND), INTENT(IN ) :: pmand (Lmax)

    INTEGER :: i, j, k
    INTEGER :: kp
    REAL (KIND=RKIND) :: alnpmd(Lmax)
    REAL (KIND=RKIND) :: pref(Lmax)

    INTEGER           :: ks   (Lmax,nCells)

    REAL (KIND=RKIND) :: ddp  (Lmax,nCells)
 
    REAL (KIND=RKIND) :: pim  (nVertLevels  ,nCells)
    REAL (KIND=RKIND) :: alnpm(nVertLevels  ,nCells)
    REAL (KIND=RKIND) :: alnp (nVertLevels+1,nCells)
    REAL (KIND=RKIND) :: pi   (nVertLevels+1,nCells)
    REAL (KIND=RKIND) :: p    (nVertLevels+1,nCells)
    REAL (KIND=RKIND) :: rdelp(nVertLevels-1,nCells)
    REAL (KIND=RKIND) :: bqq  (nVertLevels-1,nCells)

    REAL (KIND=RKIND) :: wtl

    wtl=LOG(Po)
    
    ! compute pressure on each sigma surface 

       DO k=1,nVertLevels
          DO i=1,nCells
            ! transform from Pa to mbar
             p(k,i)= 1.e-2_RKIND * pressure(k,i)
          END DO
       END DO    
       ! set the top of the model at Pt (See Constants).
       p(nVertLevels+1,:)=Pt


       DO kp=1,Lmax
            ! transform from Pa to mbar
          pref(kp) = 1.e-2_RKIND * pmand (kp)
       END DO

    ! log of pressure and exner function on each sigma surface

       DO k=1,nVertLevels+1
          DO i=1,nCells
             alnp(k,i)=LOG(MAX(p(k,i),0.00000001_RKIND))
             pi  (k,i)=EXP(RdByCp*(alnp(k,i)-wtl))
          END DO
       END DO

       DO kp=1,Lmax
          alnpmd(kp) = LOG(MAX(pref(kp),0.00000001_RKIND))
       END DO
    ! mean value of the exner function, log of pim and
    ! the log of the mean value of pressure in each layer

       DO k=1,nVertLevels
          DO i=1,nCells
             pim(k,i)=(pi(k+1,i)-pi(k,i))/ &
                  (RdByCp*(alnp(k+1,i)-alnp(k,i)))
             alnpm(k,i)=wtl+CpByRd*LOG(pim(k,i))
          END DO
       END DO

    ! return interpolated t absolute. ts replaced by temperature

       DO k=1,nVertLevels-1
          DO i=1,nCells
             rdelp(k,i)=1.0_RKIND/(alnpm(k+1,i)-alnpm(k,i))
          END DO
       END DO

       DO k=1,nVertLevels-1
          DO i=1,nCells
             bqq(k,i)=(gsh(k+1,i)-gsh(k,i))*rdelp(k,i)
          END DO
       END DO

    ! heights to pressure surfaces  -  hydrostatic interpolation
    ! winds and temps by linear interpolation with ln(p)
    ! search for middle of sigma layer above kp

       ks(:,:)=nVertLevels
       DO kp=1,Lmax
          DO k=1,nVertLevels-1
             DO i=1,nCells
                IF ((ks(kp,i) == nVertLevels) .AND. (alnpmd(kp) > alnpm(k,i))) THEN
                   ks(kp,i)=k
                END IF
             END DO
          END DO
       END DO

    ! find values and slopes for upward or downward extrapolation
    ! as well as for interpolation away from tropopause
    ! ts to temperature for absolute t

       DO kp=1,Lmax
          DO i=1,nCells
             IF (ks(kp,i) /= 1) THEN
                IF (pref(kp) < PRHcut) THEN
                   rq(kp,i)=1.0e-12_RKIND
                ELSE
                   IF (ks(kp,i) <= nVertLevels) THEN
                      rq(kp,i)=gsh(ks(kp,i),i)+bqq(ks(kp,i)-1,i)*(alnpmd(kp)-alnpm(ks(kp,i),i))
                   ELSE
                      rq(kp,i)=gsh(nVertLevels,i)
                   END IF
                END IF
                ddp(kp,i)=alnpmd(kp)-alnpm(ks(kp,i),i)
             ELSE
                IF(.NOT.ExtrapoAdiabatica)THEN
                   rq(kp,i)=gsh(1,i)
                ELSE
                   rq(kp,i)=gsh(1,i)*EXP(0.5_RKIND*RdByCp*ABS(alnpm(ks(kp,i),i)-alnpmd(kp)))
                END IF
             END IF
          END DO
       END DO

  END SUBROUTINE VerticalInterpolation_Humid

  SUBROUTINE VerticalInterpolation_Thermal(nCells,nVertLevels,Lmax,pressure,temperature, gsh, gss, tg, rg,rq, pmand)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nCells !voronoi cel number
    INTEGER, INTENT(IN) :: nVertLevels !vertical layers number
    INTEGER, INTENT(IN) :: Lmax !isobaric layers to interpolate

    REAL (KIND=RKIND), INTENT(IN ) :: pressure   (nVertLevels,nCells) ! pressao da camada atmosferica do modelo [Pa]
    REAL (KIND=RKIND), INTENT(IN ) :: temperature(nVertLevels,nCells) ! campo de temperatura do modelo (voronoi) [K]
    REAL (KIND=RKIND), INTENT(IN ) :: gsh        (nVertLevels,nCells) ! umidade especifica do modelo [kg/kg]

    REAL (KIND=RKIND), INTENT(OUT) :: gss        (nVertLevels,nCells) ! umidade relativa na grade do modelo (no nivel do modelo) 
                                                                      ! nao interpolada [frac. #0~1]
    REAL (KIND=RKIND), INTENT(OUT) :: tg         (Lmax,nCells) !interpolated temperature [k]
    REAL (KIND=RKIND), INTENT(OUT) :: rg         (Lmax,nCells) !interpolated rel. Hum. [frac. #0~1]
    REAL (KIND=RKIND), INTENT(OUT) :: rq         (Lmax,nCells) !interpolated espc. hum. [kg/kg]
    REAL (KIND=RKIND), INTENT(IN ) :: pmand (Lmax) ! isobaric levels to interpolate [Pa] 

    INTEGER :: i, j, k
    INTEGER :: kp
    REAL (KIND=RKIND) :: alnpmd(Lmax)
    REAL (KIND=RKIND) :: pref(Lmax)

    INTEGER           :: ks   (Lmax,nCells)

    REAL (KIND=RKIND) :: ddp  (Lmax,nCells)
 
    REAL (KIND=RKIND) :: pim  (nVertLevels  ,nCells)
    REAL (KIND=RKIND) :: alnpm(nVertLevels  ,nCells)
    REAL (KIND=RKIND) :: alnp (nVertLevels+1,nCells)
    REAL (KIND=RKIND) :: pi   (nVertLevels+1,nCells)
    REAL (KIND=RKIND) :: p    (nVertLevels+1,nCells)
    REAL (KIND=RKIND) :: rdelp(nVertLevels-1,nCells)
    REAL (KIND=RKIND) :: brh  (nVertLevels-1,nCells)
    REAL (KIND=RKIND) :: bqq  (nVertLevels-1,nCells)
    REAL (KIND=RKIND) :: bt   (nVertLevels-1,nCells)

    REAL (KIND=RKIND) :: es
    REAL (KIND=RKIND) :: ee
    REAL (KIND=RKIND) :: wtl

    wtl=LOG(Po)
    
    ! compute pressure on each sigma surface 

       DO k=1,nVertLevels
          DO i=1,nCells
            ! transform from Pa to mbar
             p(k,i)= 1.e-2_RKIND * pressure(k,i)
          END DO
       END DO    
       ! set the top of the model at Pt (See Constants).
       p(nVertLevels+1,:)=Pt


       DO kp=1,Lmax
            ! transform from Pa to mbar
          pref(kp) = 1.e-2_RKIND * pmand (kp)
       END DO

    ! vapour pressure form specific humidity

       DO k=1,nVertLevels
          DO i=1,nCells
             gss(k,i)=MAX(gsh(k,i),SHmin)
             es=Eo*EXP(a*(temperature(k,i)-To)/(temperature(k,i)-b))
             ee=0.5_RKIND*(p(k,i)+p(k+1,i))*gss(k,i)/(eps+eps1*gss(k,i))
             gss(k,i)=ee/es
             gss(k,i)=MIN(gss(k,i),1.0_RKIND)
          END DO
       END DO

    ! log of pressure and exner function on each sigma surface

       DO k=1,nVertLevels+1
          DO i=1,nCells
             alnp(k,i)=LOG(MAX(p(k,i),0.00000001_RKIND))
             pi  (k,i)=EXP(RdByCp*(alnp(k,i)-wtl))
          END DO
       END DO

       DO kp=1,Lmax
          alnpmd(kp) = LOG(MAX(pref(kp),0.00000001_RKIND))
       END DO
    ! mean value of the exner function, log of pim and
    ! the log of the mean value of pressure in each layer

       DO k=1,nVertLevels
          DO i=1,nCells
             pim(k,i)=(pi(k+1,i)-pi(k,i))/ &
                  (RdByCp*(alnp(k+1,i)-alnp(k,i)))
             alnpm(k,i)=wtl+CpByRd*LOG(pim(k,i))
          END DO
       END DO

    ! return interpolated t absolute. ts replaced by temperature

       DO k=1,nVertLevels-1
          DO i=1,nCells
             rdelp(k,i)=1.0_RKIND/(alnpm(k+1,i)-alnpm(k,i))
             bt(k,i)=(temperature(k+1,i)-temperature(k,i))*rdelp(k,i)
          END DO
       END DO

       DO k=1,nVertLevels-1
          DO i=1,nCells
             brh(k,i)=(gss(k+1,i)-gss(k,i))*rdelp(k,i)
             bqq(k,i)=(gsh(k+1,i)-gsh(k,i))*rdelp(k,i)
          END DO
       END DO

    ! heights to pressure surfaces  -  hydrostatic interpolation
    ! winds and temps by linear interpolation with ln(p)
    ! search for middle of sigma layer above kp

       ks(:,:)=nVertLevels
       DO kp=1,Lmax
          DO k=1,nVertLevels-1
             DO i=1,nCells
                IF ((ks(kp,i) == nVertLevels) .AND. (alnpmd(kp) > alnpm(k,i))) THEN
                   ks(kp,i)=k
                END IF
             END DO
          END DO
       END DO

    ! find values and slopes for upward or downward extrapolation
    ! as well as for interpolation away from tropopause
    ! ts to temperature for absolute t

       DO kp=1,Lmax
          DO i=1,nCells
             IF (ks(kp,i) /= 1) THEN
                IF (pref(kp) < PRHcut) THEN
                   rg(kp,i)=0.0_RKIND
                   rq(kp,i)=1.0e-12_RKIND
                ELSE
                   IF (ks(kp,i) <= nVertLevels) THEN
                      rg(kp,i)=gss(ks(kp,i),i)+brh(ks(kp,i)-1,i)*(alnpmd(kp)-alnpm(ks(kp,i),i))
                      rq(kp,i)=gsh(ks(kp,i),i)+bqq(ks(kp,i)-1,i)*(alnpmd(kp)-alnpm(ks(kp,i),i))
                   ELSE
                      rg(kp,i)=gss(nVertLevels,i)
                      rq(kp,i)=gsh(nVertLevels,i)
                   END IF
                END IF
                ddp(kp,i)=alnpmd(kp)-alnpm(ks(kp,i),i)
                tg (kp,i)=temperature(ks(kp,i),i)+bt(ks(kp,i)-1,i)*ddp(kp,i)
             ELSE
                IF(.NOT.ExtrapoAdiabatica)THEN
                   rg(kp,i)=gss(1,i)
                   rq(kp,i)=gsh(1,i)
                   tg(kp,i)=temperature(1,i)
                ELSE
                   rg(kp,i)=gss(1,i)
                   rq(kp,i)=gsh(1,i)*EXP(0.5_RKIND*RdByCp*ABS(alnpm(ks(kp,i),i)-alnpmd(kp)))
                   tg(kp,i)=temperature (1,i)*EXP(0.5_RKIND*RdByCp*ABS(alnpm(ks(kp,i),i)-alnpmd(kp)))  
                END IF
             END IF
             rg(kp,i)=MAX(rg(kp,i),RHmin)
             rg(kp,i)=MIN(rg(kp,i),RHmax)
          END DO
       END DO

  END SUBROUTINE VerticalInterpolation_Thermal



  SUBROUTINE VerticalInterpolation_Wind(nCells,nVertLevels,Lmax,pressure, pref, bfi1, bfo1, bfi2, bfo2, bfi3, bfo3)
    ! fgu, ug, fgv, vg, fgomega, og
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nCells
    INTEGER, INTENT(IN) :: nVertLevels
    INTEGER, INTENT(IN) :: Lmax
    
    REAL (KIND=RKIND), INTENT(IN ) :: pressure(nVertLevels,nCells)
    REAL (KIND=RKIND), INTENT(IN ) :: pref (Lmax)

    REAL (KIND=RKIND), INTENT(IN )           :: bfi1(nVertLevels,nCells)
    REAL (KIND=RKIND), INTENT(OUT)           :: bfo1(Lmax       ,nCells)

    REAL (KIND=RKIND), INTENT(IN ), OPTIONAL :: bfi2(nVertLevels,nCells)
    REAL (KIND=RKIND), INTENT(OUT), OPTIONAL :: bfo2(Lmax,nCells)

    REAL (KIND=RKIND), INTENT(IN ), OPTIONAL :: bfi3(nVertLevels,nCells)
    REAL (KIND=RKIND), INTENT(OUT), OPTIONAL :: bfo3(Lmax       ,nCells)

    REAL (KIND=RKIND) :: wtl

    INTEGER :: nb, i, j, k, kp
    REAL (KIND=RKIND) :: alnpmd (Lmax)
    REAL (KIND=RKIND) :: pmand   (Lmax)
    
    INTEGER           :: ks    (Lmax         ,nCells)

    REAL (KIND=RKIND) :: ddp   (Lmax         ,nCells)
    REAL (KIND=RKIND) :: pim   (nVertLevels  ,nCells)
    REAL (KIND=RKIND) :: alnpm (nVertLevels  ,nCells)
    REAL (KIND=RKIND) :: alnp  (nVertLevels+1,nCells)
    REAL (KIND=RKIND) :: pi    (nVertLevels+1,nCells)
    REAL (KIND=RKIND) :: p     (nVertLevels+1,nCells)
    REAL (KIND=RKIND) :: rdelp (nVertLevels-1,nCells)
    REAL (KIND=RKIND) :: b1    (nVertLevels-1,nCells)
    REAL (KIND=RKIND) :: b2    (nVertLevels-1,nCells)
    REAL (KIND=RKIND) :: b3    (nVertLevels-1,nCells)

    ! initialization

    nb=1
    IF (PRESENT(bfi2)) THEN
       IF (PRESENT(bfo2)) THEN
          nb=2
       ELSE
          WRITE (UNIT=nferr, FMT='(A)') ' sig2po: bfo2 required '
          STOP 9001
       END IF
       IF (PRESENT(bfi3)) THEN
          IF (PRESENT(bfo3)) THEN
             nb=3
          ELSE
             WRITE (UNIT=nferr, FMT='(A)') ' sig2po: bfo3 required '
             STOP 9001
          END IF
       END IF
    END IF


    wtl=LOG(Po)

    ! compute pressure on each sigma_z surface 
    ! set the highest to Pt (See Constants)
    DO k=1,nVertLevels
       DO i=1,nCells
          ! transform from Pa to mbar
          p(k,i)= 1.e-2_RKIND * pressure(k,i)
       END DO
    END DO
    p(nVertLevels+1,:)=Pt
    DO kp=1,Lmax
      ! transform from Pa to mbar
       pmand(kp)= 1.e-2_RKIND *pref(kp)
    END DO
    ! log of pressure and exner function on each sigma surface.

    DO k=1,nVertLevels+1
       DO i=1,nCells
          alnp(k,i)=LOG(MAX(p(k,i),0.00000001_RKIND))
          pi  (k,i)=EXP(RdByCp*(alnp(k,i)-wtl))
       END DO
    END DO

    DO kp=1,Lmax
       alnpmd(kp)= LOG(MAX(pmand(kp),0.00000001_RKIND))
    END DO

    ! mean value of the exner function,
    ! log of pim and 
    ! log of the mean value of pressure 
    ! in each layer

    DO k=1,nVertLevels
       DO i=1,nCells
             pim(k,i) = (pi  (k+1,i)-pi  (k,i))/ &
                (RdByCp*(alnp(k+1,i)-alnp(k,i)))
             alnpm(k,i)=wtl+CpByRd*LOG(pim(k,i))
       END DO
    END DO

    ! interpolated t absolute. ts replaced by tm

    IF (nb == 3) THEN
       DO k=1,nVertLevels-1
          DO i=1,nCells
             rdelp(k,i)=1.0_RKIND/(alnpm(k+1,i)-alnpm(k,i))
             b1(k,i)=(bfi1(k+1,i)-bfi1(k,i))*rdelp(k,i)
             b2(k,i)=(bfi2(k+1,i)-bfi2(k,i))*rdelp(k,i)
             b3(k,i)=(bfi3(k+1,i)-bfi3(k,i))*rdelp(k,i)
          END DO
       END DO
    ELSE IF (nb == 2) THEN
       DO k=1,nVertLevels-1
          DO i=1,nCells
             rdelp(k,i)=1.0_RKIND/(alnpm(k+1,i)-alnpm(k,i))
             b1(k,i)=(bfi1(k+1,i)-bfi1(k,i))*rdelp(k,i)
             b2(k,i)=(bfi2(k+1,i)-bfi2(k,i))*rdelp(k,i)
          END DO
       END DO
    ELSE 
       DO k=1,nVertLevels-1
          DO i=1,nCells
             rdelp(k,i)=1.0_RKIND/(alnpm(k+1,i)-alnpm(k,i))
             b1(k,i)=(bfi1(k+1,i)-bfi1(k,i))*rdelp(k,i)
          END DO
       END DO
    END IF

    ! heights to pressure surfaces  -  hydrostatic interpolation
    ! winds and temps by linear interpolation with ln(p)
    ! search for middle of sigma layer above kp

    ks(:,:)=nVertLevels
    DO kp=1,Lmax
       DO k=1,nVertLevels-1
          DO i=1,nCells
             IF ((ks(kp,i) == nVertLevels) .AND. (alnpmd(kp) > alnpm(k,i))) THEN
                   ks(kp,i)=k
             END IF
          END DO
       END DO
    END DO

    ! find values and slopes for upward or downward extrapolation
    ! as well as for interpolation away from tropopause

    IF (nb == 3) THEN
       DO kp=1,Lmax
          DO i=1,nCells
             IF (ks(kp,i) /= 1) THEN
                ddp(kp,i)=alnpmd(kp)-alnpm(ks(kp,i),i)
                bfo1(kp,i)=bfi1(ks(kp,i),i)+b1(ks(kp,i)-1,i)*ddp(kp,i)
                bfo2(kp,i)=bfi2(ks(kp,i),i)+b2(ks(kp,i)-1,i)*ddp(kp,i)
                bfo3(kp,i)=bfi3(ks(kp,i),i)+b3(ks(kp,i)-1,i)*ddp(kp,i)
             ELSE
                bfo1(kp,i)=bfi1(1,i)
                bfo2(kp,i)=bfi2(1,i)
                bfo3(kp,i)=bfi3(1,i)
             END IF
          END DO
       END DO
    ELSE IF (nb == 2) THEN
       DO kp=1,Lmax
          DO i=1,nCells
             IF (ks(kp,i) /= 1) THEN
                ddp(kp,i)=alnpmd(kp)-alnpm(ks(kp,i),i)
                bfo1(kp,i)=bfi1(ks(kp,i),i)+b1(ks(kp,i)-1,i)*ddp(kp,i)
                bfo2(kp,i)=bfi2(ks(kp,i),i)+b2(ks(kp,i)-1,i)*ddp(kp,i)
             ELSE
                bfo1(kp,i)=bfi1(1,i)
                bfo2(kp,i)=bfi2(1,i)
             END IF
          END DO
       END DO
    ELSE
       DO kp=1,Lmax
          DO i=1,nCells
             IF (ks(kp,i) /= 1) THEN
                ddp(kp,i)=alnpmd(kp)-alnpm(ks(kp,i),i)
                bfo1(kp,i)=bfi1(ks(kp,i),i)+b1(ks(kp,i)-1,i)*ddp(kp,i)
             ELSE
                bfo1(kp,i)=bfi1(1,i)
             END IF
          END DO
       END DO
    END IF

  END SUBROUTINE VerticalInterpolation_Wind






  SUBROUTINE VerticalInterpolation_ZGEO (nVertLevels,nCells,temperature, pressure, top,&
                     lmax, zp, pmand)

    IMPLICIT NONE
    INTEGER       , INTENT(IN ) :: lmax
    INTEGER       , INTENT(IN ) :: nVertLevels
    INTEGER       , INTENT(IN ) :: nCells
    REAL (KIND=RKIND), INTENT(IN ) :: temperature(nVertLevels,nCells)
    REAL (KIND=RKIND), INTENT(IN ) :: pressure   (nVertLevels,nCells)
    REAL (KIND=RKIND), INTENT(IN ) :: top(nCells)
    REAL (KIND=RKIND), INTENT(OUT) :: zp    (Lmax,nCells)
    REAL (KIND=RKIND), INTENT(IN ) :: pmand (Lmax)

    REAL (KIND=RKIND) :: wtl ,part,hmh   
    INTEGER :: i, j, k, kp
    REAL (KIND=RKIND) :: alnpmd(Lmax)
    REAL (KIND=RKIND) :: pref   (Lmax)
    INTEGER       , DIMENSION(nVertLevels,nCells) :: ks
    REAL (KIND=RKIND), DIMENSION(nCells) :: gamma, tmsl, tstr, zlay
    REAL (KIND=RKIND), DIMENSION(nVertLevels+1,nCells) :: alnp, pi, p, z, tlev, balnp
    REAL (KIND=RKIND), DIMENSION(nVertLevels  ,nCells) :: pim, alnpm, dh, bz, zt

    ! constants

    wtl=LOG(Po)

    ! compute pressure on each hybrid sigma_z surface and
    ! set the top of the model at Pt (See Constants)
     
       DO k=1,nVertLevels
         ! call mpas_log_write('--- end interpolate vertical velocity:')
         ! print*,k,MAXVAL(pressure(k,:)),MINVAL(pressure(k,:))
          DO i=1,nCells
             ! transform from Pa to mbar

             p(k,i) = 1.e-2_RKIND * pressure(k,i)
          END DO
       END DO
       ! Reference Top Pressure Value (Mb)
       p(nVertLevels+1,:)=Pt
       DO kp=1,Lmax
         ! transform from Pa to mbar
          pref(kp)= 1.e-2_RKIND *pmand(kp)
       END DO

    ! compute the log of pressure and exner function on each hybrid sigma_z surface

       DO k=1,nVertLevels+1
          DO i=1,nCells
             alnp(k,i)=LOG(p(k,i))
             pi  (k,i)=EXP(RdByCp*(alnp(k,i)-wtl))
          END DO
       END DO

       DO kp=1,Lmax
          alnpmd(kp)= LOG(MAX(pref(kp),0.001_RKIND))
       END DO


    ! compute the mean value of the exner function and
    ! the log of the mean value of pressure in each layer.

       DO k=1,nVertLevels
          DO i=1,nCells
             IF((alnp(k+1,i)-alnp(k,i)) == 0.0_RKIND)THEN
                pim(k,i) = pi(k,i)
                alnpm(k,i)=wtl+CpByRd*LOG(pim(k,i))

	     ELSE
	     
                pim(k,i) = (pi(k+1,i)-pi(k,i))/ &
                          (RdByCp*(alnp(k+1,i)-alnp(k,i))   )
                alnpm(k,i)=wtl+CpByRd*LOG(pim(k,i))
             END IF 
          END DO
       END DO


    ! compute the height of the sigma_z surfaces.

       DO i=1,nCells
          !
          ! height of the sigma_z - topography
          z(1,i)=top(i)
          !
       END DO

       DO k=1,nVertLevels
          DO i=1,nCells
             z(k+1,i)=z(k,i) + RdByGrav*temperature(k,i) * MAX((alnp(k,i)-alnp(k+1,i)),0.00000001_RKIND)
          END DO
       END DO

    ! compute the temperature on sigma_z surfaces. model the temp
    ! profile as linear in log of pressure. extrapolate to determine
    ! the value of temp at the highest and lowest sigma_z surface.

       DO k=2,nVertLevels
          DO i=1,nCells
             IF((alnpm(k-1,i)-alnpm(k,i)) ==0.0_RKIND)THEN
                tlev(k,i)=temperature(k  ,i)
	     ELSE
                tlev(k,i)=((alnp (k  ,i)-alnpm(k,i))*temperature(k-1,i)+ &
                           (alnpm(k-1,i)-alnp (k,i))*temperature(k  ,i))/ &
                           (alnpm(k-1,i)-alnpm(k,i))
             END IF
          END DO
       END DO

       DO i=1,nCells
         !
         ! Surface
         !
          IF((alnpm(1,i)-alnpm(2,i)) ==0.0_RKIND)THEN
             tlev(1,i)=temperature(1,i)
          ELSE  
             tlev(1,i)=((alnp (1,i)-alnpm(2,i))*temperature(1,i)+ &
                        (alnpm(1,i)-alnp (1,i))*temperature(2,i))/ &
                        (alnpm(1,i)-alnpm(2,i))
          END IF
         !
         ! top atmosphere
         !
          IF((alnpm(nVertLevels-1,i)-alnpm(nVertLevels  ,i)) ==0.0_RKIND)THEN
          tlev(nVertLevels+1,i)=temperature(nVertLevels  ,i)
          ELSE
          tlev(nVertLevels+1,i)= &
               ((alnp (nVertLevels+1,i)-alnpm(nVertLevels  ,i))*temperature(nVertLevels-1,i)+ &
                (alnpm(nVertLevels-1,i)-alnp (nVertLevels+1,i))*temperature(nVertLevels  ,i))/ &
                (alnpm(nVertLevels-1,i)-alnpm(nVertLevels  ,i))
          END IF
       END DO


       DO i=1,nCells
          zlay(i)=0.5_RKIND*(z(1,i)+z(2,i))
          tstr(i)=temperature(1,i)+TVVTa*(zlay(i)-z(1,i))
          tmsl(i)=temperature(1,i)+TVVTa*zlay(i)
          gamma(i)=0.0_RKIND
          IF (tmsl(i) > Tref) THEN
             tmsl(i)=Tref
             IF (tstr(i) > Tref) tmsl(i)=Tref - TVVTb*(tstr(i)-Tref)*(tstr(i)-Tref)
          END IF
 
 
          IF (z(1,i) > Zref) gamma(i)=(tstr(i)-tmsl(i))/z(1,i)
       END DO


       DO k=1,nVertLevels
          DO i=1,nCells
             IF((alnp(k+1,i)-alnp(k,i)) == 0.0_RKIND) THEN
                 dh(k,i)=0.0_RKIND
                 bz(k,i)=0.0_RKIND
	     ELSE
                dh(k,i)=alnp(k+1,i)-alnp(k,i)              
                bz(k,i)=RdByGrav*(tlev(k+1,i)-tlev(k,i))/dh(k,i)
             END IF  
             zt(k,i)=0.5_RKIND*(z(k+1,i)+z(k,i))+ 0.125_RKIND*bz(k,i)*dh(k,i)*dh(k,i)

             balnp(k,i)=0.5_RKIND*(alnp(k+1,i)+alnp(k,i))

          END DO
       END DO

    ! start with another search thru the sigma levels
    ! to find the level (ks)  above the desires mandatory level (kp)
    ! heights to pressure surfaces  -  hydrostatic interpolation
    ! winds and temps by linear interpolation with ln(p)
    ! search for middle of sigma layer above kp


       ks(:,:) = nVertLevels+1
       DO kp=1,Lmax
          DO k=1,nVertLevels
             DO i=1,nCells
                IF ((ks(kp,i) == nVertLevels+1) .AND. (alnpmd(kp) > alnp(k,i))) THEN
                   ks(kp,i) = k
                END IF
             END DO
          END DO
       END DO

    ! find values and slopes for upward or downward extrapolation
    ! as well as for interpolation away from tropopause
    ! ts to tm for absolute t

       DO kp=1,Lmax
          DO i=1,nCells
             IF (ks(kp,i) == 1) THEN
                ! zp(i,kp) is below ground.
                part=RdByGrav*(alnp(1,i)-alnpmd(kp))
                
                zp(kp,i)=z(1,i)+tstr(i)*part/(1.0_RKIND-0.5_RKIND*part*gamma(i))
             ELSE
                ! zp references a pressure level in the free air
                hmh=alnpmd(kp)-balnp(ks(kp,i)-1,i)
                zp(kp,i)=zt(ks(kp,i)-1,i)- (RdByGrav*temperature(ks(kp,i)-1,i)+ &
                     0.5_RKIND*bz(ks(kp,i)-1,i)*hmh)*hmh
             END IF
          END DO
       END DO
  END SUBROUTINE VerticalInterpolation_ZGEO
  
  
    SUBROUTINE VerticalInterpolation_Single(nCells,nVertLevels,nlevs, nlevr,Lmax, gausin, gauout, pressure, pmand_in)
    IMPLICIT NONE
    ! vertical interpolation of gaussian grid fields 
    INTEGER, INTENT(IN) :: nCells
    INTEGER, INTENT(IN) :: nVertLevels
    INTEGER, INTENT(IN) :: nlevs
    INTEGER, INTENT(IN) :: nlevr
    INTEGER, INTENT(IN) :: Lmax
    REAL (KIND=RKIND), INTENT(INOUT) :: gausin  (nVertLevels,nCells)
    REAL (KIND=RKIND), INTENT(OUT  ) :: gauout  (nlevr,nCells)
    REAL (KIND=RKIND), INTENT(IN   ) :: pressure(nVertLevels,nCells)
    REAL (KIND=RKIND), INTENT(IN   ) :: pmand_in   (Lmax)
    INTEGER :: k
    INTEGER :: j
    INTEGER :: i
    INTEGER :: lat
    INTEGER :: int
    REAL (KIND=RKIND) :: pmand   (Lmax)

    REAL (KIND=RKIND) :: deltap
    REAL (KIND=RKIND) :: df
    REAL (KIND=RKIND) :: dp
    REAL (KIND=RKIND) :: work(Lmax,nCells)
    REAL (KIND=RKIND) :: pin (nVertLevels,nCells)

    LOGICAL :: Above(Lmax,nCells)

    gauout=0.0_RKIND
    IF (nlevs > 1 .AND. nlevs < nVertLevels) THEN
       DO k=nlevs+1,nVertLevels
          DO i=1,nCells
             gausin(k,i)=0.0_RKIND
          END DO
       END DO
    END IF

       ! single level inputs

       IF (nlevs ==  1) THEN
          DO i=1,nCells
             work(1,i)=gausin(1,i)
          END DO

          ! vertical interpolation

       ELSE
          DO k=1,nVertLevels
             DO i=1,nCells
                ! transform from Pa to mbar
                pin(k,i)= 1.e-2_RKIND * pressure(k,i)
             END DO
          END DO
          DO k=1,Lmax
            ! transform from Pa to mbar
             pmand(k)= 1.e-2_RKIND * pmand_in(k)
          END DO
  
          Above = .FALSE.
          DO k=1,Lmax
             DO i=1,nCells
                IF (pmand(k) > pin(1,i)) THEN
                   work(k,i)=gausin(1,i)
                   Above(k,i)=.TRUE.
                END IF
             END DO
          END DO

          DO k=1,Lmax
             DO i=1,nCells
                IF (.NOT. Above(k,i) .AND. pmand(k) <= pin(nVertLevels,i)) THEN
                   work(k,i)=gausin(nVertLevels,i)
                   Above(k,i)=.TRUE.
                END IF
             END DO
          END DO

          DO k=1,Lmax
             DO int=1,nVertLevels-1
                DO i=1,nCells
                   IF (.NOT. Above(k,i) .AND. &
                        pin(int  ,i)  >= pmand(k) .AND. &
                        pin(int+1,i)  < pmand(k)) THEN

                      ! interpolation linear in p

                      deltap=pmand(k)-pin(int,i)
                      df       =gausin(int+1,i) - gausin(int,i)
                      dp       =pin(int+1   ,i) - pin(int,i)
                      work(k,i)=gausin(int,i)+(df/dp)*deltap
                   END IF
                END DO
             END DO
          END DO
       END IF

       DO k=1,nlevr
          DO i=1,nCells
             gauout(k,i)=work(k,i)
          END DO
       END DO

  END SUBROUTINE VerticalInterpolation_Single
END MODULE  VerticalInterpolation
