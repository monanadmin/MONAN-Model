module BalanceErrorCheckGlacierMod

!!! Check glacier water and energy balance and report error

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use NoahmpFatalErrorMod,only: Noahmp_error_fatal
  use mpas_log

  implicit none

contains

!!!! Water balance check initialization
  subroutine BalanceWaterInitGlacier(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: None (embedded in NOAHMP_GLACIER)
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! --------------------------------------------------------------------
    associate(                                                             &
              SnowWaterEquiv     => noahmp%water%state%SnowWaterEquiv     ,& ! in,  snow water equivalent [mm]
              WaterStorageTotBeg => noahmp%water%state%WaterStorageTotBeg  & ! out, total water storage [mm] at the beginning
             )
! ----------------------------------------------------------------------

    ! compute total glacier water storage before NoahMP processes
    ! need more work on including glacier ice mass underneath snow
    WaterStorageTotBeg = SnowWaterEquiv

    end associate

  end subroutine BalanceWaterInitGlacier


!!!! Water balance check and report error
  subroutine BalanceWaterCheckGlacier(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: ERROR_GLACIER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! --------------------------------------------------------------------
    associate(                                                             &
              GridIndexI         => noahmp%config%domain%GridIndexI       ,& ! in,  grid index in x-direction
              GridIndexJ         => noahmp%config%domain%GridIndexJ       ,& ! in,  grid index in y-direction
              MainTimeStep       => noahmp%config%domain%MainTimeStep     ,& ! in,  main noahmp timestep [s]
              SnowWaterEquiv     => noahmp%water%state%SnowWaterEquiv     ,& ! in,  snow water equivalent [mm]
              WaterStorageTotBeg => noahmp%water%state%WaterStorageTotBeg ,& ! in,  total water storage [mm] at the beginning
              PrecipTotRefHeight => noahmp%water%flux%PrecipTotRefHeight  ,& ! in,  total precipitation [mm/s] at reference height
              EvapGroundNet      => noahmp%water%flux%EvapGroundNet       ,& ! in,  net ground evaporation [mm/s]
              RunoffSurface      => noahmp%water%flux%RunoffSurface       ,& ! in,  surface runoff [mm/s]
              RunoffSubsurface   => noahmp%water%flux%RunoffSubsurface    ,& ! in,  subsurface runoff [mm/s]
              WaterStorageTotEnd => noahmp%water%state%WaterStorageTotEnd ,& ! out, total water storage [mm] at the end
              WaterBalanceError  => noahmp%water%state%WaterBalanceError   & ! out, water balance error [mm] per time step
             )
! ----------------------------------------------------------------------

    ! Error in water balance should be < 0.1 mm
    ! compute total glacier water storage before NoahMP processes
    ! need more work on including glacier ice mass underneath snow
    WaterStorageTotEnd = SnowWaterEquiv
    WaterBalanceError  = WaterStorageTotEnd - WaterStorageTotBeg - &
                         (PrecipTotRefHeight - EvapGroundNet - RunoffSurface - RunoffSubsurface) * MainTimeStep

#ifndef WRF_HYDRO
    if ( abs(WaterBalanceError) > 0.1 ) then
       call mpas_log_write(' ')
       call mpas_log_write('---~---')
       call mpas_log_write('   Noah-MP water budget conservation error (glacier):')
       call mpas_log_write('---~---')
       call mpas_log_write('   GridIndexI                  = $i ', intArgs  = (/ GridIndexI                      /) )
       call mpas_log_write('   GridIndexJ                  = $i ', intArgs  = (/ GridIndexJ                      /) )
       call mpas_log_write('   WaterStorageTotBeg          = $r ', realArgs = (/ WaterStorageTotBeg              /) )
       call mpas_log_write('   WaterStorageTotEnd          = $r ', realArgs = (/ WaterStorageTotEnd              /) )
       call mpas_log_write('   WaterBalanceError           = $r ', realArgs = (/ WaterBalanceError               /) )
       call mpas_log_write('                                 (positive value above means water gain).'              )
       call mpas_log_write('   PrecipTotRefHeight          = $r ', realArgs = (/ PrecipTotRefHeight*MainTimeStep /) )
       call mpas_log_write('   EvapGroundNet               = $r ', realArgs = (/ EvapGroundNet*MainTimeStep      /) )
       call mpas_log_write('   RunoffSurface               = $r ', realArgs = (/ RunoffSurface*MainTimeStep      /) )
       call mpas_log_write('   RunoffSubsurface            = $r ', realArgs = (/ RunoffSubsurface*MainTimeStep   /) )
       call mpas_log_write(' ')
       call Noahmp_error_fatal("Error: Water budget problem in NoahMP LSM (glacier)")
    endif
#endif

    end associate

  end subroutine BalanceWaterCheckGlacier


!!!! Energy balance check and error report
  subroutine BalanceEnergyCheckGlacier(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: ERROR_GLACIER
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! --------------------------------------------------------------------
    associate(                                                              &
              GridIndexI         => noahmp%config%domain%GridIndexI        ,& ! in,  grid index in x-direction
              GridIndexJ         => noahmp%config%domain%GridIndexJ        ,& ! in,  grid index in y-direction
              RadSwDownRefHeight => noahmp%forcing%RadSwDownRefHeight      ,& ! in,  downward shortwave radiation [W/m2] at reference height
              RadSwAbsSfc        => noahmp%energy%flux%RadSwAbsSfc         ,& ! in,  total absorbed solar radiation [W/m2]
              RadSwReflSfc       => noahmp%energy%flux%RadSwReflSfc        ,& ! in,  total reflected solar radiation [W/m2]
              RadLwNetSfc        => noahmp%energy%flux%RadLwNetSfc         ,& ! in,  total net longwave rad [W/m2] (+ to atm)
              HeatSensibleSfc    => noahmp%energy%flux%HeatSensibleSfc     ,& ! in,  total sensible heat [W/m2] (+ to atm)
              HeatLatentGrd      => noahmp%energy%flux%HeatLatentGrd       ,& ! in,  total ground latent heat [W/m2] (+ to atm)
              HeatGroundTot      => noahmp%energy%flux%HeatGroundTot       ,& ! in,  total ground heat flux [W/m2] (+ to soil/snow)
              RadSwAbsGrd        => noahmp%energy%flux%RadSwAbsGrd         ,& ! in,  solar radiation absorbed by ground [W/m2]
              HeatPrecipAdvSfc   => noahmp%energy%flux%HeatPrecipAdvSfc    ,& ! in,  precipitation advected heat - total [W/m2]
              EnergyBalanceError => noahmp%energy%state%EnergyBalanceError ,& ! out, error in surface energy balance [W/m2]
              RadSwBalanceError  => noahmp%energy%state%RadSwBalanceError   & ! out, error in shortwave radiation balance [W/m2]
             )
! ----------------------------------------------------------------------

    ! error in shortwave radiation balance should be <0.01 W/m2
    RadSwBalanceError = RadSwDownRefHeight - (RadSwAbsSfc + RadSwReflSfc)
    ! print out diagnostics when error is large
    if ( abs(RadSwBalanceError) > 0.01 ) then
       call mpas_log_write(' ')
       call mpas_log_write('---~---')
       call mpas_log_write('   Noah-MP solar radiation budget conservation error (glacier):')
       call mpas_log_write('---~---')
       call mpas_log_write('   GridIndexI                          = $i ', intArgs  = (/ GridIndexI           /) )
       call mpas_log_write('   GridIndexJ                          = $i ', intArgs  = (/ GridIndexJ           /) )
       call mpas_log_write('   RadSwBalanceError                   = $r ', realArgs = (/ RadSwBalanceError    /) )
       call mpas_log_write('                                         (positive value above means energy gain).'  )
       call mpas_log_write('   RadSwBalanceError                   = $r ', realArgs = (/ RadSwBalanceError    /) )
       call mpas_log_write('   RadSwDownRefHeight                  = $r ', realArgs = (/ RadSwDownRefHeight   /) )
       call mpas_log_write('   RadSwReflSfc                        = $r ', realArgs = (/ RadSwReflSfc         /) )
       call mpas_log_write('   RadSwAbsGrd                         = $r ', realArgs = (/ RadSwAbsGrd          /) )
       call mpas_log_write('   RadSwAbsSfc                         = $r ', realArgs = (/ RadSwAbsSfc          /) )
       call mpas_log_write('---~---')
       call mpas_log_write(' ')
       call Noahmp_error_fatal("Error: Solar radiation budget problem in NoahMP LSM (glacier)")
    endif

    ! error in surface energy balance should be <0.01 W/m2
    EnergyBalanceError = RadSwAbsGrd + HeatPrecipAdvSfc - (RadLwNetSfc + HeatSensibleSfc + HeatLatentGrd + HeatGroundTot)
    ! print out diagnostics when error is large
    if ( abs(EnergyBalanceError) > 0.01 ) then
       call mpas_log_write(' ')
       call mpas_log_write('---~---')
       call mpas_log_write('   Noah-MP energy budget conservation error (glacier):')
       call mpas_log_write('---~---')
       call mpas_log_write('   GridIndexI                  = $i ', intArgs  = (/ GridIndexI           /) )
       call mpas_log_write('   GridIndexJ                  = $i ', intArgs  = (/ GridIndexJ           /) )
       call mpas_log_write('   EnergyBalanceError          = $r ', realArgs = (/ EnergyBalanceError   /) )
       call mpas_log_write('                                 (positive value above energy gain).'        )
       call mpas_log_write('   Net longwave                = $r ', realArgs = (/ RadLwNetSfc          /) )
       call mpas_log_write('   Total sensible              = $r ', realArgs = (/ HeatSensibleSfc      /) )
       call mpas_log_write('   Ground evap                 = $r ', realArgs = (/ HeatLatentGrd        /) )
       call mpas_log_write('   Total ground                = $r ', realArgs = (/ HeatGroundTot        /) )
       call mpas_log_write('   Precip advected             = $r ', realArgs = (/ HeatPrecipAdvSfc     /) )
       call mpas_log_write('   Absorbed shortwave          = $r ', realArgs = (/ RadSwAbsGrd          /) )
       call mpas_log_write('---~---')
       call mpas_log_write(' ')
       call Noahmp_error_fatal("Error: Energy budget problem in NoahMP LSM (glacier)")
    endif

    end associate

  end subroutine BalanceEnergyCheckGlacier

end module BalanceErrorCheckGlacierMod
