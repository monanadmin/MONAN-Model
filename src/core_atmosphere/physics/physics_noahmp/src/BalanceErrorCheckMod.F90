module BalanceErrorCheckMod

!!! Check water and energy balance and report error

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use NoahmpFatalErrorMod,only: Noahmp_error_fatal
  use mpas_log

  implicit none

contains

!!!! Water balance check initialization
  subroutine BalanceWaterInit(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: None (embedded in NOAHMP_SFLX)
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: LoopInd      ! loop index

! --------------------------------------------------------------------
    associate(                                                                       &
              NumSoilLayer           => noahmp%config%domain%NumSoilLayer           ,& ! in,  number of soil layers
              SurfaceType            => noahmp%config%domain%SurfaceType            ,& ! in,  surface type 1-soil; 2-lake
              ThicknessSnowSoilLayer => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in,  thickness of snow/soil layers [m]
              CanopyLiqWater         => noahmp%water%state%CanopyLiqWater           ,& ! in,  canopy intercepted liquid water [mm]
              CanopyIce              => noahmp%water%state%CanopyIce                ,& ! in,  canopy intercepted ice [mm]
              SnowWaterEquiv         => noahmp%water%state%SnowWaterEquiv           ,& ! in,  snow water equivalent [mm]
              SoilMoisture           => noahmp%water%state%SoilMoisture             ,& ! in,  total soil moisture [m3/m3]
              WaterStorageAquifer    => noahmp%water%state%WaterStorageAquifer      ,& ! in,  water storage in aquifer [mm]
              WaterStorageTotBeg     => noahmp%water%state%WaterStorageTotBeg        & ! out, total water storage [mm] at the beginning
             )
! ----------------------------------------------------------------------

    ! compute total water storage before NoahMP processes
    if ( SurfaceType == 1 ) then  ! soil
       WaterStorageTotBeg = CanopyLiqWater + CanopyIce + SnowWaterEquiv + WaterStorageAquifer
       do LoopInd = 1, NumSoilLayer
          WaterStorageTotBeg = WaterStorageTotBeg + SoilMoisture(LoopInd) * ThicknessSnowSoilLayer(LoopInd) * 1000.0
       enddo
    endif

    end associate

  end subroutine BalanceWaterInit


!!!! Water balance check and report error
  subroutine BalanceWaterCheck(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: ERROR
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    integer                          :: LoopInd      ! loop index
    real(kind=kind_noahmp)           :: SoilWaterStorage

! --------------------------------------------------------------------
    associate(                                                                        &
              NumSoilLayer            => noahmp%config%domain%NumSoilLayer           ,& ! in,    number of soil layers
              SurfaceType             => noahmp%config%domain%SurfaceType            ,& ! in,    surface type 1-soil; 2-lake
              GridIndexI              => noahmp%config%domain%GridIndexI             ,& ! in,    grid index in x-direction
              GridIndexJ              => noahmp%config%domain%GridIndexJ             ,& ! in,    grid index in y-direction
              ThicknessSnowSoilLayer  => noahmp%config%domain%ThicknessSnowSoilLayer ,& ! in,    thickness of snow/soil layers [m]
              MainTimeStep            => noahmp%config%domain%MainTimeStep           ,& ! in,    main noahmp timestep [s]
              FlagCropland            => noahmp%config%domain%FlagCropland           ,& ! in,    flag to identify croplands
              FlagSoilProcess         => noahmp%config%domain%FlagSoilProcess        ,& ! in,    flag to calculate soil process
              IrriFracThreshold       => noahmp%water%param%IrriFracThreshold        ,& ! in,    irrigation fraction parameter
              IrrigationFracGrid      => noahmp%water%state%IrrigationFracGrid       ,& ! in,    total input irrigation fraction
              WaterTableDepth         => noahmp%water%state%WaterTableDepth          ,& ! in,    water table depth [m]
              CanopyLiqWater          => noahmp%water%state%CanopyLiqWater           ,& ! in,    canopy intercepted liquid water [mm]
              CanopyIce               => noahmp%water%state%CanopyIce                ,& ! in,    canopy intercepted ice [mm]
              SnowWaterEquiv          => noahmp%water%state%SnowWaterEquiv           ,& ! in,    snow water equivalent [mm]
              SoilMoisture            => noahmp%water%state%SoilMoisture             ,& ! in,    total soil moisture [m3/m3]
              WaterStorageAquifer     => noahmp%water%state%WaterStorageAquifer      ,& ! in,    water storage in aquifer [mm]
              WaterStorageTotBeg      => noahmp%water%state%WaterStorageTotBeg       ,& ! in,    total water storage [mm] at the beginning
              PrecipTotRefHeight      => noahmp%water%flux%PrecipTotRefHeight        ,& ! in,    total precipitation [mm/s] at reference height
              EvapCanopyNet           => noahmp%water%flux%EvapCanopyNet             ,& ! in,    evaporation of intercepted water [mm/s]
              Transpiration           => noahmp%water%flux%Transpiration             ,& ! in,    transpiration rate [mm/s]
              EvapGroundNet           => noahmp%water%flux%EvapGroundNet             ,& ! in,    net ground (soil/snow) evaporation [mm/s]
              RunoffSurface           => noahmp%water%flux%RunoffSurface             ,& ! in,    surface runoff [mm/dt_soil] per soil timestep
              RunoffSubsurface        => noahmp%water%flux%RunoffSubsurface          ,& ! in,    subsurface runoff [mm/dt_soil] per soil timestep
              TileDrain               => noahmp%water%flux%TileDrain                 ,& ! in,    tile drainage [mm/dt_soil] per soil timestep
              IrrigationRateSprinkler => noahmp%water%flux%IrrigationRateSprinkler   ,& ! in,    rate of irrigation by sprinkler [m/timestep]
              IrrigationRateMicro     => noahmp%water%flux%IrrigationRateMicro       ,& ! in,    micro irrigation water rate [m/timestep]
              IrrigationRateFlood     => noahmp%water%flux%IrrigationRateFlood       ,& ! in,    flood irrigation water rate [m/timestep]
              SfcWaterTotChgAcc       => noahmp%water%flux%SfcWaterTotChgAcc         ,& ! inout, accumulated snow,soil,canopy water change per soil timestep [mm]
              PrecipTotAcc            => noahmp%water%flux%PrecipTotAcc              ,& ! inout, accumulated precipitation per soil timestep [mm]
              EvapCanopyNetAcc        => noahmp%water%flux%EvapCanopyNetAcc          ,& ! inout, accumulated net canopy evaporation per soil timestep [mm]
              TranspirationAcc        => noahmp%water%flux%TranspirationAcc          ,& ! inout, accumulated transpiration per soil timestep [mm]
              EvapGroundNetAcc        => noahmp%water%flux%EvapGroundNetAcc          ,& ! inout, accumulated net ground evaporation per soil timestep [mm]
              WaterStorageTotEnd      => noahmp%water%state%WaterStorageTotEnd       ,& ! out,   total water storage [mm] at the end
              WaterBalanceError       => noahmp%water%state%WaterBalanceError         & ! out,   water balance error [mm] per time step
             )
! ----------------------------------------------------------------------

    ! before water balance check, add irrigation water to precipitation
    if ( (FlagCropland .eqv. .true.) .and. (IrrigationFracGrid >= IrriFracThreshold) ) then
       PrecipTotRefHeight = PrecipTotRefHeight + IrrigationRateSprinkler * 1000.0 / MainTimeStep  ! irrigation
    endif

    ! only water balance check for every soil timestep
    ! Error in water balance should be < 0.1 mm
    if ( SurfaceType == 1 ) then   ! soil
       SoilWaterStorage = 0.0
       do LoopInd = 1, NumSoilLayer
          SoilWaterStorage = SoilWaterStorage + SoilMoisture(LoopInd) * ThicknessSnowSoilLayer(LoopInd) * 1000.0
       enddo
       WaterStorageTotEnd = SoilWaterStorage + CanopyLiqWater + CanopyIce + SnowWaterEquiv + WaterStorageAquifer

       ! accumualted water change (only for canopy and snow during non-soil timestep)
       SfcWaterTotChgAcc = SfcWaterTotChgAcc + (WaterStorageTotEnd - WaterStorageTotBeg)  ! snow, canopy, and soil water change
       PrecipTotAcc      = PrecipTotAcc      + PrecipTotRefHeight * MainTimeStep          ! accumulated precip 
       EvapCanopyNetAcc  = EvapCanopyNetAcc  + EvapCanopyNet      * MainTimeStep          ! accumulated canopy evapo
       TranspirationAcc  = TranspirationAcc  + Transpiration      * MainTimeStep          ! accumulated transpiration
       EvapGroundNetAcc  = EvapGroundNetAcc  + EvapGroundNet      * MainTimeStep          ! accumulated soil evapo

       ! check water balance at soil timestep
       if ( FlagSoilProcess .eqv. .true. ) then
          WaterBalanceError = SfcWaterTotChgAcc - (PrecipTotAcc + IrrigationRateMicro*1000.0 + IrrigationRateFlood*1000.0 - &
                              EvapCanopyNetAcc - TranspirationAcc - EvapGroundNetAcc - RunoffSurface - RunoffSubsurface -   &
                              TileDrain)
#ifndef WRF_HYDRO
          if ( abs(WaterBalanceError) > 0.1 ) then
             call mpas_log_write(' ')
             call mpas_log_write('---~---')
             call mpas_log_write('   Noah-MP water budget conservation error (land):')
             call mpas_log_write('---~---')
             call mpas_log_write('   GridIndexI                  = $i ', intArgs  = (/ GridIndexI                 /) )
             call mpas_log_write('   GridIndexJ                  = $i ', intArgs  = (/ GridIndexJ                 /) )
             call mpas_log_write('   WaterStorageTotBeg          = $r ', realArgs = (/ WaterStorageTotBeg         /) )
             call mpas_log_write('   WaterStorageTotEnd          = $r ', realArgs = (/ WaterStorageTotEnd         /) )
             call mpas_log_write('   WaterBalanceError           = $r ', realArgs = (/ WaterBalanceError          /) )
             call mpas_log_write('                                 (positive value above means water gain).'         )
             call mpas_log_write('--- State:')
             call mpas_log_write('   SoilWaterStorage            = $r ', realArgs = (/ SoilWaterStorage           /) )
             call mpas_log_write('   CanopyLiqWater              = $r ', realArgs = (/ CanopyLiqWater             /) )
             call mpas_log_write('   CanopyIce                   = $r ', realArgs = (/ CanopyIce                  /) )
             call mpas_log_write('   SnowWaterEquiv              = $r ', realArgs = (/ SnowWaterEquiv             /) )
             call mpas_log_write('   WaterStorageAquifer         = $r ', realArgs = (/ WaterStorageAquifer        /) )
             call mpas_log_write('--- Fluxes:')
             call mpas_log_write('   SfcWaterTotChgAcc           = $r ', realArgs = (/ SfcWaterTotChgAcc          /) )
             call mpas_log_write('   PrecipTotAcc                = $r ', realArgs = (/ PrecipTotAcc               /) )
             call mpas_log_write('   IrrigationRateMicro*1000.0  = $r ', realArgs = (/ IrrigationRateMicro*1000.0 /) )
             call mpas_log_write('   IrrigationRateFlood*1000.0  = $r ', realArgs = (/ IrrigationRateFlood*1000.0 /) )
             call mpas_log_write('   EvapCanopyNetAcc            = $r ', realArgs = (/ EvapCanopyNetAcc           /) )
             call mpas_log_write('   EvapGroundNetAcc            = $r ', realArgs = (/ EvapGroundNetAcc           /) )
             call mpas_log_write('   TranspirationAcc            = $r ', realArgs = (/ TranspirationAcc           /) )
             call mpas_log_write('   RunoffSurface               = $r ', realArgs = (/ RunoffSurface              /) )
             call mpas_log_write('   RunoffSubsurface            = $r ', realArgs = (/ RunoffSubsurface           /) )
             call mpas_log_write('   WaterTableDepth             = $r ', realArgs = (/ WaterTableDepth            /) )
             call mpas_log_write('   TileDrain                   = $r ', realArgs = (/ TileDrain                  /) )
             call mpas_log_write('--- Soil state:')
             call mpas_log_write('   NumSoilLayer = $i ',intArgs=(/ NumSoilLayer /)) 
             do LoopInd=1,NumSoilLayer
                call mpas_log_write('   Layer $i - Thickness = $r ; SoilMoisture = $r'                       &
                                   , intArgs  = (/ LoopInd /)                                                &
                                   , realArgs = (/ ThicknessSnowSoilLayer(LoopInd), SoilMoisture(LoopInd) /) )
             end do
             call mpas_log_write('---~---')
             call mpas_log_write(' ')
             call Noahmp_error_fatal("Error: Water budget problem in NoahMP LSM (land)")
          endif
#endif
       endif ! FlagSoilProcess

    else ! water point
       WaterBalanceError = 0.0
    endif

    end associate

  end subroutine BalanceWaterCheck


!!!! Energy balance check and error report
  subroutine BalanceEnergyCheck(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: ERROR
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! --------------------------------------------------------------------
    associate(                                                                 &
              GridIndexI           => noahmp%config%domain%GridIndexI         ,& ! in,  grid index in x-direction
              GridIndexJ           => noahmp%config%domain%GridIndexJ         ,& ! in,  grid index in y-direction
              RadSwDownRefHeight   => noahmp%forcing%RadSwDownRefHeight       ,& ! in,  downward shortwave radiation [W/m2] at reference height
              VegFrac              => noahmp%energy%state%VegFrac             ,& ! in,  greeness vegetation fraction
              RadSwAbsSfc          => noahmp%energy%flux%RadSwAbsSfc          ,& ! in,  total absorbed solar radiation [W/m2]
              RadSwReflSfc         => noahmp%energy%flux%RadSwReflSfc         ,& ! in,  total reflected solar radiation [W/m2]
              RadSwReflVeg         => noahmp%energy%flux%RadSwReflVeg         ,& ! in,  reflected solar radiation by vegetation [W/m2]
              RadSwReflGrd         => noahmp%energy%flux%RadSwReflGrd         ,& ! in,  reflected solar radiation by ground [W/m2]
              RadLwNetSfc          => noahmp%energy%flux%RadLwNetSfc          ,& ! in,  total net longwave rad [W/m2] (+ to atm)
              HeatSensibleSfc      => noahmp%energy%flux%HeatSensibleSfc      ,& ! in,  total sensible heat [W/m2] (+ to atm)
              HeatLatentCanopy     => noahmp%energy%flux%HeatLatentCanopy     ,& ! in,  canopy latent heat flux [W/m2] (+ to atm)
              HeatLatentGrd        => noahmp%energy%flux%HeatLatentGrd        ,& ! in,  total ground latent heat [W/m2] (+ to atm)
              HeatLatentTransp     => noahmp%energy%flux%HeatLatentTransp     ,& ! in,  latent heat flux from transpiration [W/m2] (+ to atm)
              HeatGroundTot        => noahmp%energy%flux%HeatGroundTot        ,& ! in,  total ground heat flux [W/m2] (+ to soil/snow)
              RadSwAbsVeg          => noahmp%energy%flux%RadSwAbsVeg          ,& ! in,  solar radiation absorbed by vegetation [W/m2]
              RadSwAbsGrd          => noahmp%energy%flux%RadSwAbsGrd          ,& ! in,  solar radiation absorbed by ground [W/m2]
              HeatPrecipAdvSfc     => noahmp%energy%flux%HeatPrecipAdvSfc     ,& ! in,  precipitation advected heat - total [W/m2]
              HeatPrecipAdvBareGrd => noahmp%energy%flux%HeatPrecipAdvBareGrd ,& ! in,  precipitation advected heat - bare ground net [W/m2]
              HeatPrecipAdvVegGrd  => noahmp%energy%flux%HeatPrecipAdvVegGrd  ,& ! in,  precipitation advected heat - under canopy net [W/m2]
              HeatPrecipAdvCanopy  => noahmp%energy%flux%HeatPrecipAdvCanopy  ,& ! in,  precipitation advected heat - vegetation net [W/m2]
              HeatLatentIrriEvap   => noahmp%energy%flux%HeatLatentIrriEvap   ,& ! in,  latent heating due to sprinkler evaporation [W/m2]
              HeatCanStorageChg    => noahmp%energy%flux%HeatCanStorageChg    ,& ! in,  canopy heat storage change [W/m2]
              EnergyBalanceError   => noahmp%energy%state%EnergyBalanceError  ,& ! out, error in surface energy balance [W/m2]
              RadSwBalanceError    => noahmp%energy%state%RadSwBalanceError    & ! out, error in shortwave radiation balance [W/m2]
             )
! ----------------------------------------------------------------------

    ! error in shortwave radiation balance should be <0.01 W/m2
    RadSwBalanceError = RadSwDownRefHeight - (RadSwAbsSfc + RadSwReflSfc)
    ! print out diagnostics when error is large
    if ( abs(RadSwBalanceError) > 0.01 ) then
       call mpas_log_write(' ')
       call mpas_log_write('---~---')
       call mpas_log_write('   Noah-MP solar radiation budget conservation error (land):')
       call mpas_log_write('---~---')
       call mpas_log_write('   GridIndexI                          = $i ', intArgs  = (/ GridIndexI           /) )
       call mpas_log_write('   GridIndexJ                          = $i ', intArgs  = (/ GridIndexJ           /) )
       call mpas_log_write('   RadSwBalanceError                   = $r ', realArgs = (/ RadSwBalanceError    /) )
       call mpas_log_write('                                         (positive value above means energy gain).'  )
       call mpas_log_write('--- Vegetation ')
       call mpas_log_write('   RadSwDownRefHeight * VegFrac        = $r ', realArgs = (/ RadSwDownRefHeight*VegFrac        /) )
       call mpas_log_write('   VegFrac*RadSwAbsVeg + RadSwAbsGrd   = $r ', realArgs = (/ VegFrac*RadSwAbsVeg+RadSwAbsGrd   /) )
       call mpas_log_write('   VegFrac*RadSwReflVeg + RadSwReflGrd = $r ', realArgs = (/ VegFrac*RadSwReflVeg+RadSwReflGrd /) )
       call mpas_log_write('--- Ground ')
       call mpas_log_write('   (1 - VegFrac) * RadSwDownRefHeight  = $r ', realArgs = (/ (1.0-VegFrac)*RadSwDownRefHeight  /) )
       call mpas_log_write('   (1 - VegFrac) * RadSwAbsGrd         = $r ', realArgs = (/ (1.0-VegFrac)*RadSwAbsGrd         /) )
       call mpas_log_write('   (1 - VegFrac) * RadSwReflGrd        = $r ', realArgs = (/ (1.0-VegFrac)*RadSwReflGrd        /) )
       call mpas_log_write('   RadSwReflVeg                        = $r ', realArgs = (/ RadSwReflVeg                      /) )
       call mpas_log_write('   RadSwReflGrd                        = $r ', realArgs = (/ RadSwReflGrd                      /) )
       call mpas_log_write('   RadSwReflSfc                        = $r ', realArgs = (/ RadSwReflSfc                      /) )
       call mpas_log_write('   RadSwAbsVeg                         = $r ', realArgs = (/ RadSwAbsVeg                       /) )
       call mpas_log_write('   RadSwAbsGrd                         = $r ', realArgs = (/ RadSwAbsGrd                       /) )
       call mpas_log_write('   RadSwAbsSfc                         = $r ', realArgs = (/ RadSwAbsSfc                       /) )
       call mpas_log_write('---~---')
       call mpas_log_write(' ')
       call Noahmp_error_fatal("Error: Solar radiation budget problem in NoahMP LSM (land)")
    endif

    ! error in surface energy balance should be <0.01 W/m2
    EnergyBalanceError = RadSwAbsVeg + RadSwAbsGrd + HeatPrecipAdvSfc -                     &
                        (RadLwNetSfc + HeatSensibleSfc + HeatLatentCanopy + HeatLatentGrd + &
                         HeatLatentTransp + HeatGroundTot + HeatLatentIrriEvap + HeatCanStorageChg)
    ! print out diagnostics when error is large
    if ( abs(EnergyBalanceError) > 0.01 ) then
       call mpas_log_write(' ')
       call mpas_log_write('---~---')
       call mpas_log_write('   Noah-MP energy budget conservation error (land):')
       call mpas_log_write('---~---')
       call mpas_log_write('   GridIndexI                  = $i ', intArgs  = (/ GridIndexI           /) )
       call mpas_log_write('   GridIndexJ                  = $i ', intArgs  = (/ GridIndexJ           /) )
       call mpas_log_write('   EnergyBalanceError          = $r ', realArgs = (/ EnergyBalanceError   /) )
       call mpas_log_write('                                 (positive value above energy gain).'        )
       call mpas_log_write('   Net solar                   = $r ', realArgs = (/ RadSwAbsSfc          /) )
       call mpas_log_write('   Net longwave                = $r ', realArgs = (/ RadLwNetSfc          /) )
       call mpas_log_write('   Total sensible              = $r ', realArgs = (/ HeatSensibleSfc      /) )
       call mpas_log_write('   Canopy evap                 = $r ', realArgs = (/ HeatLatentCanopy     /) )
       call mpas_log_write('   Ground evap                 = $r ', realArgs = (/ HeatLatentGrd        /) )
       call mpas_log_write('   Transpiration               = $r ', realArgs = (/ HeatLatentTransp     /) )
       call mpas_log_write('   Total ground                = $r ', realArgs = (/ HeatGroundTot        /) )
       call mpas_log_write('   Sprinkler                   = $r ', realArgs = (/ HeatLatentIrriEvap   /) )
       call mpas_log_write('   Canopy heat storage change  = $r ', realArgs = (/ HeatCanStorageChg    /) )

       call mpas_log_write('   Precip advected (surface)   = $r ', realArgs = (/ HeatPrecipAdvSfc     /) )
       call mpas_log_write('   Precip advected (canopy)    = $r ', realArgs = (/ HeatPrecipAdvCanopy  /) )
       call mpas_log_write('   Precip advected (veg. grnd) = $r ', realArgs = (/ HeatPrecipAdvVegGrd  /) )
       call mpas_log_write('   Precip advected (bare grnd) = $r ', realArgs = (/ HeatPrecipAdvBareGrd /) )
       call mpas_log_write('   Veg fraction                = $r ', realArgs = (/ VegFrac              /) )
       call mpas_log_write('---~---')
       call mpas_log_write(' ')
       call Noahmp_error_fatal("Error: Energy budget problem in NoahMP LSM (land)")
    endif

    end associate

  end subroutine BalanceEnergyCheck

end module BalanceErrorCheckMod
