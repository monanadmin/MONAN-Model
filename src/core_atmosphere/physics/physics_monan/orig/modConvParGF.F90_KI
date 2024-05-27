module modConvParGF
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! ## Grell-Freitas Convective parameterization
   !!
   !! ![](https://i.ibb.co/LNqGy3S/logo-Monan-Color-75x75.png)
   !! ## MONAN
   !!
   !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
   !!
   !! E-mail: <mailto:saulo.freitas@inpe.br>
   !!
   !! Date: 2014
   !!
   !! #####Version: 0.1.0
   !!
   !! ---
   !! **Full description**:
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!- This convective parameterization is build to attempt                      !
   !!  a smooth transition to cloud resolving scales as proposed                 !
   !!  by Arakawa et al (2011, ACP). The scheme is  described                    !
   !!  in the paper Grell and Freitas (ACP, 2014).                               ! 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!- Implemented in GEOS5 GCM by Saulo Freitas (July 2016)                     !
   !!- Use the following references for this implementation:                     !
   !!- Freitas et al (2018, JAMES/AGU, https://doi.org/10.1029/2017MS001251)     !
   !!- Freitas et al (2021, GMD/EGU,   https://doi.org/10.5194/gmd-14-5393-2021) !
   !!- Please, contact Saulo Freitas (saulo.r.de.freitas@gmail.com) for comments !
   !!- questions, bugs, etc.                                                     !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! ** History**:
   !! -Adapted for BRAMS 6.0 by Saulo Freitas (November 2021)                    
   !! -Refactoring by Luiz Flavio Rodrigues at 20 December 2021 (Monday)         
   !!  Keywords using ; are separeted, some loops receives exit instead goto,    
   !!  The identation was fixed and all keywords are lowercase                   
   !! -Refactoring by GCC (INPE) at 20 January 2023 using fprettify and manual
   !!  changes according MONAN rules code patterns DTN 01
   !!
   !! --- 
   !! ** Licence **:
   !!
   !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
   !!
   !! This program is free software: you can redistribute it and/or modify
   !! it under the terms of the GNU General Public License as published by
   !! the  Free  Software  Foundation, either version 3 of the License, or
   !! (at your option) any later version.
   !!
   !! This program is distributed in the hope that it  will be useful, but
   !! ** WITHOUT  ANY  WARRANTY **;  without  even  the   implied   warranty  of
   !! **MERCHANTABILITY** or **FITNESS FOR A  PARTICULAR PURPOSE**.  See  the, GNU
   !! GNU General Public License for more details.
   !!
   !! You should have received a copy  of the GNU General  Public  License
   !! along with this program.  
   !! If not, see [GNU Public License](https://www.gnu.org/licenses/gpl-3.0.html).
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !
   use modGate, only: cupout, rundata, p_nvar_grads, jl, p_use_gate, ppres, ptemp, pq, pu &
                  ,   runlabel, runname, pv, pvervel, pgeo, zqr, zadvq, zadvt

   use modConstants, only: c_rgas, c_cp, c_rm, c_p00, c_tcrit, c_grav, c_cpor, c_alvl, c_pi &    
                        ,  c_akmin, c_ccnclean, c_t00, c_t_ice, c_xlf, p_max_qsat           &
                        ,  p_xmbmaxshal, p_mintracer, p_smaller_qv, c_t01, c_t100, c_temp0i &
                        ,  c_rgas_atm, c_hplus, c_r2es, c_r3les, c_r4ies, c_r4les, c_retv   &
                        ,  c_rticecu, c_rtwat_rticecu_r, c_r3ies, c_r5alscp, c_r5alvcp      &
                        ,  c_ralsdcp, c_ralvdcp, c_rtice, c_rtwat_rtice_r, p_tkmin, c_tice  &
                        ,  p_ccnclean,  c_rtt, c_rtwat
   !
   use modHenrysLawCts, only: getHenrysLawCts

   !-- for future use
   ! use Lib_modConvParGF, only : gen_random,get_kbmax_kdet_k22,get_edt &
   !                                , get_lambdaU,reset_1d,get_capmax,reset_2d&
   !                                , get_lcl, get_cloud_bc, set_entr_detr_rates&
   !                                , set_scale_dep_factor, get_cloud_top_by_inv_layers&
   !                                , get_1st_guess_MSE_profile, get_updraft_profile  &
   !                                , get_delmix,get_downdraft_profile, check_mass_conserv&
   !                                , get_dellas,smooth_tend,get_env_change &
   !                                , tridiag, fct1d3,get_pr_ens,apply_sub_microphys &
   !                                , get_max_dd_height
   !

   !use node_mod, only: mynum
   use modVector, only: get_num_elements, get_data_value, init, insert_range, remove &
                       ,free_memory, vector_t, insert_unique

   !use module_mp_wsm3, only: WSM3_GF
   
   implicit none

   private
   public  maxiens, icumulus_gf, closure_choice, deep, shal, mid &
      ,use_scale_dep,dicycle,tau_deep,tau_mid,hcts                       &
      ,use_tracer_transp, use_tracer_scaven,use_memory,convection_tracer &
      ,use_flux_form,use_tracer_evap,downdraft,use_fct                   &
      ,use_rebcb, vert_discr, clev_grid, apply_sub_mp, alp1              &
      ,sgs_w_timescale, lightning_diag, tau_ocea_cp,  tau_land_cp        &
      ,autoconv, overshoot,use_wetbulb                                   &
      ,c0_deep, qrc_crit,lambau_deep,lambau_shdn,c0_mid               &
      ,cum_max_edt_land  ,cum_max_edt_ocean, cum_hei_down_land           &
      ,cum_hei_down_ocean,cum_hei_updf_land, cum_hei_updf_ocean          &
      ,use_momentum_transp,cum_entr_rate                                 &
      ,nmp, lsmp, cnmp,moist_trigger,frac_modis,max_tq_tend              &
      ,cum_fadj_massflx, cum_use_excess, cum_ave_layer, adv_trigger      &
      ,use_smooth_prof, output_sound,use_cloud_dissipation               &
      ,use_smooth_tend,beta_sh,c0_shal                                   &
      ,use_linear_subcl_mf,cap_maxs,liq_ice_number_conc                  &
      ,sig_factor,lcl_trigger, rh_dicycle, add_coldpool_prop             &
      ,add_coldpool_clos,add_coldpool_trig,mx_buoy1, mx_buoy2, cum_t_star&
      ,add_coldpool_diff,n_cldrop,use_gustiness, use_random_num          &
      ,dcape_threshold,modConvParGF_initialized,use_pass_cloudvol        &
      ,use_lcl_ctrl_entr,use_rhu_ctrl_entr

   public convParGFDriver,makeDropletNumber ,makeIceNumber,FractLiqF    &
      ,coldPoolStart, readGFConvParNML, initModConvParGF    
   !
   !-
   !- plume spectral size
   integer, parameter  :: maxiens = 3, deep=1 ,shal=2 , mid = 3
   character(len=10),parameter,dimension(maxiens)  :: cumulus_type = (/ &
       'deep      ' &
      ,'shallow   ' &
      ,'mid       ' &
      /)
   !- number of microphysics schemes in the host model
   !integer, parameter  :: nmp = 2, lsmp = 1, cnmp = 2 ! --- for GEOS-5
    integer, parameter  :: nmp = 2, lsmp = 1, cnmp = 2 ! --- for BRAMS
    integer, parameter ::               &
      maxens  = 1,  & ! 1  ensemble one on cap_max
      maxens2 = 1,  & ! 1  ensemble two on precip efficiency
      maxens3 = 16, & !16 ensemble three done in cup_forcing_ens16 for G3d
      ensdim  = maxens*maxens2*maxens3,&
      ens4    = 1
    integer, parameter :: shall_closures = 12 

   !------------------- namelist variables -----------------------------------

   !-- plume to be activated (1 true, 0 false): deep, shallow, congestus
   integer, dimension(maxiens) :: icumulus_gf   ! = (/1,1,1/)
   !
   !-- choice for the closures:
   !--  deep   : 0 ensemble (all)          , 1 GR, 4 ll omega, 7 moist conv, 10 PB
   !--  shallow: 0 ensemble (all)          , 1 Wstar, 4 heat-engine, 7 BLQE, 10 TKE-based
   !--  mid    : 0 ensemble (Wstar/BLQE/PB), 1 Wstar, 2 BLQE, 3 PB, 4 PB_BL

   integer, dimension(maxiens) :: closure_choice  != (/10,  10,  3/) ! deep, shallow, congestus

   !-- gross entrainment rate: deep, shallow, congestus
   real,    dimension(maxiens) :: cum_entr_rate != (/&
                                                ! 6.3e-4  & !deep
                                                !,1.0e-3  & !shallow
                                                !,5.0e-4  & !mid
                                                !/)
   integer :: use_lcl_ctrl_entr != 0/1     - default 0

   integer :: use_rhu_ctrl_entr != 0/1     - default 1

   integer :: use_pass_cloudvol != 0/1     - default 1

   integer :: use_tracer_transp != 0/1     - default 1

   integer :: use_tracer_scaven != 0/1/2/3 - default 2

   integer :: use_flux_form     != 1/2/3   - default 1

   integer :: use_fct           != 0/1     - default 1 (only for use_flux_form     = 2)

   integer :: use_tracer_evap   != 0/1     - default 1 (only for use_tracer_scaven > 0)

   integer :: convection_tracer != 0/1:  turn on/off the "convection" tracer

   integer :: use_memory        != -1/0/1/2 .../10    !-

   integer :: add_coldpool_prop != -1,0,1,2,3 add coldpool propagation

   integer :: add_coldpool_clos ! add the the mass flux associated to the W @ leading of the gust front
   
   integer :: add_coldpool_trig ! add triggering criteria based on cold pool presence

   integer :: add_coldpool_diff ! add vert/horizontal diffusion to the cold pool propaga

   integer :: use_scale_dep     != 0/1:  scale dependence flag, default = 1

   integer :: dicycle           != 0/1/2:  diurnal cycle closure, default = 1
                                != 2 uses qadv closure (becker et al 2021)
   integer :: rh_dicycle        ! controls of rh on the diurnal cycle (see tian et al 2022 grl)

   integer :: clev_grid         != 0/1/2: interpolation method to define environ state at the
                                != cloud levels (at face layer), default = 0
                                != clev_grid = 0 default method
                                != clev_grid = 1 interpolation method based on tiedtke (1989)
                                != clev_grid = 2 for gate soundings only

   integer :: use_rebcb         != 0/1: turn on/off rainfall evap below cloud base, default = 0

   integer :: vert_discr        != 0/1: 1=new vert discretization, default = 0

   integer :: sgs_w_timescale   != 0/1: vertical velocity for tau_ecmwf, default = 0

   integer :: lightning_diag    != 0/1: do lightning_diaggnostics based on lopez (2016, mwr)

   integer :: apply_sub_mp      != 0/1: subsidence transport applied the to grid-scale/anvil ice/liq mix
                                !=  ratio and cloud fraction

   real    :: alp1              != 0/0.5/1: apply subsidence transport of ls/anvil cloud fraction using
                                !=      time implicit discretization

   integer :: use_wetbulb       != 0/1

   real    :: overshoot         != 0, 1

   integer :: autoconv          != 1, 3 or 4 autoconversion formulation: (1) kessler,
                                !  (3) kessler with temp dependence, (4) sundvisqt
   real    ::  c0_deep          != default= 3.e-3   conversion rate (cloud to rain, m-1) - for deep      plume
   real    ::  c0_mid           != default= 2.e-3   conversion rate (cloud to rain, m-1) - for congestus plume
   real    ::  c0_shal          != default= 0.e-3   conversion rate (cloud to rain, m-1) - for shallow   plume
   real    ::  qrc_crit         != default= 2.e-4   kg/kg

   integer :: use_momentum_transp != 0/1:  turn on/off conv transp of momentum
   real    ::  lambau_deep        != default= 2.0 lambda parameter for deep/congestus convection momentum transp
   real    ::  lambau_shdn        != default= 2.0 lambda parameter for shallow/downdraft convection momentum transp

   integer :: downdraft           != 0/1:  turn on/off downdrafts, default = 1


   real    ::  tau_deep            != deep      convective timescale
   real    ::  tau_mid             != congestus convective timescale

   real    :: max_tq_tend          != max t,q tendency allowed (100 k/day)

   integer :: use_smooth_prof      != 1 makes the normalized mass flux, entr and detrainment profiles smoother

   integer :: use_smooth_tend      != 0 => off, > 0 produces smoother tendencies (e.g.: for 1=> makes average between k-1,k,k+1)
   !---                                      deep, shallow, congestus
   integer :: irainevap           

   real,   dimension(maxiens) :: cum_hei_down_land != [0.2,0.8] height of the max z downdraft , default = 0.50
   real,   dimension(maxiens) :: cum_hei_down_ocean!= [0.2,0.8] height of the max z downdraft , default = 0.50
   real,   dimension(maxiens) :: cum_hei_updf_land != [0.2,0.8] height of the max z updraft   , default = 0.35
   real,   dimension(maxiens) :: cum_hei_updf_ocean!= [0.2,0.8] height of the max z updraft   , default = 0.35
   real,   dimension(maxiens) :: cum_max_edt_land  != maximum evap fraction allowed over the land  ,default= 0.9
   real,   dimension(maxiens) :: cum_max_edt_ocean != maximum evap fraction allowed over the ocean ,default= 0.9

   real,   dimension(maxiens) :: cum_fadj_massflx  != multiplicative factor for tunning the mass flux at cloud base
                                                   != default = 1.0
   real,   dimension(maxiens) :: cum_ave_layer     != layer depth for average the properties
                                                   != of source air parcels (mbar)
   real,   dimension(maxiens) :: cum_t_star        != scale temperature for the diurnal cycle closure

   integer,dimension(maxiens) :: cum_use_excess    != use t,q excess sub-grid scale variability

   integer :: moist_trigger  != relative humidity effects on the cap_max trigger function
   integer :: frac_modis     != use fraction liq/ice content derived from modis/calipo sensors
   integer :: adv_trigger    != dcape trigger
   real    :: dcape_threshold!= cape time rate threshold for adv_trigger = 1 (j kg^-1 hr^-1)
                             != typical range is [-200,200] j/kg/hr, wu et all (2007) recomends ~ 70 j/kg/hr
                             != 55 J/kg/hr is indicated for the Amazon basin (Song&Zhang 2017)
   integer :: lcl_trigger    != greater than zero, activates the lcl trigger which requires the lcl height
                             != be lower than the pbl height, only for shallow convection
   
   integer :: output_sound   != outputs a vertical profile for the gf stand alone model

   real    :: tau_ocea_cp    != cold pool lifetime over the ocean
   real    :: tau_land_cp    != cold pool lifetime over land
   real    :: mx_buoy1       !=   250.5 J/kg
   real    :: mx_buoy2       != 20004.0 J/kg: temp exc=10 K, q deficit=4 g/kg (=> mx_buoy ~ 20 kJ/kg)

   real    :: use_cloud_dissipation != to acccount for the cloud dissipation at the decayment phase
   integer :: use_gustiness         != not in use
   real    :: use_random_num        != stochastic pertubation for the height of maximum Zu

   real    :: beta_sh               != only for shallow plume
   integer :: use_linear_subcl_mf   != only for shallow plume
   real    :: cap_maxs              != max distance (hPa) the air parcel is allowed to go up looking for the LFC
   integer :: liq_ice_number_conc   != include drop/ice number mixing ratio convective tendencies
   real    :: sig_factor            != exponential factor for the sigma determination (orig = 0.1)
   real    :: N_cldrop              != cloud drop number concentration (cm\u02c6-3)


   !------------------- internal variables  --------------------------------------

   integer, parameter :: ON = 1, OFF = 0 !=  ON/OFF integer paremeters

   real    ::  HEI_DOWN_LAND     != [0.2,0.8] height of the max Z Downdraft , default = 0.50
   real    ::  HEI_DOWN_OCEAN    != [0.2,0.8] height of the max Z Downdraft , default = 0.50
   real    ::  HEI_UPDF_LAND     != [0.2,0.8] height of the max Z Updraft   , default = 0.35
   real    ::  HEI_UPDF_OCEAN    != [0.2,0.8] height of the max Z Updraft   , default = 0.35
   real    ::  MAX_EDT_LAND      != default= 0.9 - maximum evap fraction allowed over the land
   real    ::  MAX_EDT_OCEAN     != default= 0.9 - maximum evap fraction allowed over the ocean
   real    ::  FADJ_MASSFLX      != default= 1.0 - multiplicative factor for the mass flux at cloud base
   real    ::  T_star            != Scale Temperature for the DC closure
   real    ::  AVE_LAYER         != layer depth for average the properties of source air parcels (mbar)
   real    ::  c0                != autoconversion constant
   integer ::  USE_EXCESS        != default= 1   - use T,Q excess sub-grid scale variability

   !-- General internal controls for the diverse options in GF

   logical, parameter :: coupl_mphysics = .true.  != coupling with cloud microphysics (do not change to false)

   logical, parameter :: melt_glac      = .true.  != turn on/off ice phase/melting

   logical            :: first_guess_w  = .false. != use it to calculate a 1st guess of the updraft vert velocity

   integer, parameter :: aeroevap = 1             !=rainfall evaporation (1) orig  - (2) mix orig+new - (3) new

   logical :: ave_from_surface = .false.          != 
   !
   !- proportionality constant to estimate pressure
   !- gradient of updraft (Zhang and Wu, 2003, JAS) => REAL, PARAMETER ::    pgcon=-0.55
   real, parameter ::    pgcon= 0.0

   integer, parameter :: MAX_NSPEC=200
   integer           ,dimension(MAX_NSPEC)    :: ind_chem
   character(len=100),dimension(MAX_NSPEC)    ::  CHEM_NAME
   integer           ,dimension(MAX_NSPEC)    ::  CHEM_NAME_MASK,CHEM_NAME_MASK_EVAP
   real              ,dimension(MAX_NSPEC)    ::  CHEM_ADJ_AUTOC
   integer :: ispc_CO
   type Hcts_vars
      real :: hstar,dhr,ak0,dak
   end type Hcts_vars
   type (Hcts_vars), allocatable :: Hcts(:)

   integer :: whoami_all, JCOL,itime1_in
   real    :: time_in
   logical :: wrtgrads = .false.
   integer :: nrec = 0, ntimes = 0
   real    :: int_time = 0.

   integer :: vec_max_size
   !! max size control loop vector can assume
   type(vector_t) :: vec_ok
   !! vector where loop will execute
   type(vector_t) :: vec_removed
   !! vector of removed indexes
   logical :: is_removed, is_inserted

   logical :: modConvParGF_initialized

contains


   !---------------------------------------------------------------------------------------------------

   subroutine convParGFDriver(mxp,myp,mzp,mtp,nmp, time, itime1 &
                           ,ims,ime, jms,jme, kms,kme               &
                           ,its,ite, jts,jte, kts,kte               &
                           ,flip                                    &
                           ,FSCAV                                   &
                           ,mynum                 &
                           ,dt                    &
                           ,dx2d                  &
                           ,stochastic_sig        &
                           ,zm                    &
                           ,zt                    &
                           ,dm                    &

                           ,lons                  &
                           ,lats                  &
                           ,aot500                &
                           ,temp2m                &
                           ,sflux_r               &
                           ,sflux_t               &
                           ,qexcp                 &
                           ,hexcp                 &
                           ,wlpool                &      
                           ,topt                  &
                           ,xland                 &
                           ,sfc_press             &
                           ,kpbl                  &
                           ,tke_pbl               &
                           ,turb_len_scale        & !NEW 26FEB2024

                           ,col_sat               &
                           ,u                     &
                           ,v                     &
                           ,w                     &
                           ,temp                  &
                           ,press                 &
                           ,rvap                  &
                           ,mp_ice                &
                           ,mp_liq                &
                           ,mp_cf                 &
                           ,curr_rvap             &
                           ,TRACER                &!-note: uses GEOS-5 data structure
                           ,cnvcf                 & ! conv cloud fraction
                           !---- forcings---
                           ,buoy_exc              &
                           ,rthften               &! gsf_t
                           ,rqvften               &! gsf_q
                           ,rth_advten            &!advf_t
                           ,rthblten              &!sgsf_t
                           ,rqvblten              &!sgsf_q
                           !---- output ----
                           ,conprr                &
                           ,lightn_dens           &
                           ,rh_dicycle_fct        &
                           ,rthcuten              &
                           ,rqvcuten              &
                           ,rqccuten              &
                           ,rnlcuten              &
                           ,rnicuten              &
                           ,rucuten               &
                           ,rvcuten               &
                           ,sub_mpqi              &
                           ,sub_mpql              &
                           ,sub_mpcf              &
                           ,rbuoycuten            &
                           ,rchemcuten            &
                           ,revsu_gf              &
                           ,prfil_gf              &
                           !
                           ,do_this_column        &
                           ,ierr4d                &
                           ,jmin4d                &
                           ,klcl4d                &
                           ,k224d                 &
                           ,kbcon4d               &
                           ,ktop4d                &
                           ,kstabi4d              &
                           ,kstabm4d              &
                           ,cprr4d                &
                           ,xmb4d                 &
                           ,edt4d                 &
                           ,pwav4d                &
                           ,sigma4d               &
                           ,pcup5d                &
                           ,up_massentr5d         &
                           ,up_massdetr5d         &
                           ,dd_massentr5d         &
                           ,dd_massdetr5d         &
                           ,zup5d                 &
                           ,zdn5d                 &
                           ,prup5d                &
                           ,prdn5d                &
                           ,clwup5d               &
                           ,tup5d                 &
                           ,conv_cld_fr5d         &
                           !-- for debug/diagnostic
                           ,AA0,AA1,AA1_ADV,AA1_RADPBL,AA1_BL,AA2,AA3,AA1_CIN,TAU_BL,TAU_EC &
                           ,VAR2d,VAR3d_aGF,VAR3d_bGF,VAR3d_cGF,VAR3d_dGF)

      implicit none
      !include "mpif.h"
      !------------------------------------------------------------------------
      integer, intent(in) :: ims,ime, jms,jme, kms,kme,    &
         its,ite, jts,jte, kts,kte,    &
         mynum,mzp,mxp,myp,mtp,nmp, itime1

      real,    intent(in) :: DT, time

      integer, intent(in), dimension(mzp) :: flip

      real:: FSCAV(mtp)

      real,    dimension(kts:kte,its:ite,jts:jte), intent(in)  ::       &
         zm,        &
         zt,        &
         u,         &
         v,         &
         w,         &
         rvap,      &
         temp,      &
         press,     &
         dm,        &
         curr_rvap, &
         buoy_exc,  &
         cnvcf,     &
         qexcp,     &
         hexcp 

      integer, dimension(its:ite,jts:jte), intent(in) :: kpbl

      !-- intent (in)
      real,    dimension(its:ite,jts:jte) :: topt ,aot500 ,temp2m ,sfc_press &
                                            ,sflux_r ,sflux_t,xland,lons,lats &
                                            ,dx2d,col_sat,stochastic_sig,tke_pbl,wlpool

      real,    dimension(kts:kte,its:ite,jts:jte), intent(in) :: rthften    &
         ,rqvften    &
         ,rth_advten &
         ,rthblten   &
         ,rqvblten   &
         ,turb_len_scale

      real,    dimension(its:ite,jts:jte),         intent(out  ) ::   CONPRR,LIGHTN_DENS
      real,    dimension(its:ite,jts:jte),         intent(inout) ::   rh_dicycle_fct

      real,    dimension(kts:kte,its:ite,jts:jte), intent(out) ::     rthcuten   &
         ,rqvcuten   &
         ,rqccuten   &
         ,rnlcuten   &
         ,rnicuten   &
         ,rucuten    &
         ,rvcuten    &
         ,rbuoycuten &
         ,revsu_gf   &
         ,prfil_gf   &
         ,var3d_agf  &
         ,var3d_bgf  &
         ,var3d_cgf  &
         ,var3d_dgf


      real,    dimension(nmp,kts:kte,its:ite,jts:jte), intent(in)  :: mp_ice  &
         ,mp_liq  &
         ,mp_cf

      real,    dimension(nmp,kts:kte,its:ite,jts:jte), intent(out) :: sub_mpqi   &
         ,sub_mpql   &
         ,sub_mpcf

      !-***** TRACER has different data structure   (i,j,k,ispc) *********
      real,    dimension(its:ite,jts:jte,kts:kte,mtp), intent(in )  :: TRACER
      !-***** rchemcuten uses the GF data structure (ispc,k,i,j) *********
      real,    dimension(mtp,kts:kte,its:ite,jts:jte), intent(out)  :: rchemcuten

      integer, dimension(its:ite,jts:jte), intent(inout) :: do_this_column

      !- for convective transport and cloud/radiation (OUT)
      integer,dimension(mxp,myp,maxiens)::  &
         !  integer, dimension(its:ite,jts:jte,maxiens) , INTENT(OUT) ::    &
          ierr4d                    &
         ,jmin4d                    &
         ,klcl4d                    &
         ,k224d                     &
         ,kbcon4d                   &
         ,ktop4d                    &
         ,kstabi4d                  &
         ,kstabm4d

      real,dimension(mxp,myp,maxiens)::     &
         !   real,dimension(its:ite,jts:jte,maxiens)    , INTENT(OUT) ::    &
          cprr4d                    &
         ,xmb4d                     &
         ,edt4d                     &
         ,pwav4d                    &
         ,sigma4d
      real,dimension(mzp,mxp,myp,maxiens):: &
         !   real,dimension(its:ite,jts:jte,kts:kte,maxiens), INTENT(OUT) ::    &
          pcup5d                    &
         ,up_massentr5d             &
         ,up_massdetr5d             &
         ,dd_massentr5d             &
         ,dd_massdetr5d             &
         ,zup5d                     &
         ,zdn5d                     &
         ,prup5d                    &
         ,prdn5d                    &
         ,clwup5d                   &
         ,tup5d                     &
         ,conv_cld_fr5d
      !--for debug
      real   ,dimension(mxp,myp)  ,intent(inout)  :: aa0,aa1,aa1_adv,aa1_radpbl    &
                                                    ,aa2,aa3,aa1_bl,aa1_cin,tau_bl &
                                                    ,tau_ec,var2d

      !----------------------------------------------------------------------
      ! LOCAL VARS

      ! basic environmental input includes
      ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
      ! convection for this call only and at that particular gridpoint
      ! 
      real,   dimension (its:ite,jts:jte) ::  rtgt

      real,   dimension (kts:kte,its:ite) ::  &
                                              zo,temp_old,qv_old,po,us,vs,rhoi,phil            &
                                             ,temp_new_dp,qv_new_dp,temp_new_sh,qv_new_sh,z2d  &
                                             ,tkeg,rcpg,dhdt,temp_new_md,qv_new_md             &
                                             ,temp_new_bl,qv_new_bl,dm2d,temp_tendqv,qv_curr   &
                                             ,buoy_exc2d,revsu_gf_2d,prfil_gf_2d,var3d_agf_2d  &
                                             ,var3d_bgf_2d,temp_new,qv_new,tpert_2d            &
                                             ,temp_new_adv,qv_new_adv,cnvcf2d,turb_len_scale2d

      real,   dimension (kts:kte,its:ite,maxiens) ::  outt,outq,outqc,outu,outv,outbuoy &
                                                     ,outnliq,outnice

      real,   dimension (mtp,kts:kte,its:ite)         :: se_chem
      real,   dimension (mtp,kts:kte,its:ite,maxiens) :: out_chem

      real,   dimension (nmp,kts:kte,its:ite)         :: mpqi,mpql,mpcf
      real,   dimension (nmp,kts:kte,its:ite,maxiens) :: outmpqi,outmpql,outmpcf

      real,   dimension (its:ite)   :: ter11, xlandi,pbl,zws,ccn,psur             &
                                      ,ztexec,zqexec,h_sfc_flux,le_sfc_flux,tsur  &
                                      ,xlons,xlats,fixout_qv,cum_ztexec,cum_zqexec&
                                      ,zlcl_sfc,plcl_sfc,tlcl_sfc

      real,   dimension (kts:kte,its:ite,1:ens4)      ::  omeg

      real,   dimension (kts:kte) :: min_tend,distance
      integer,dimension (its:ite) :: kpbli

      integer :: i,j,k,kr,n,itf,jtf,ktf,ispc,zmax,status

      real :: dp,dq,exner, dtdt,pten,pqen,paph,zrho,pahfs,pqhfl,zkhvfl,pgeoh
      real :: fixouts,dt_inv

      real,   dimension(mxp,myp,-1:5) :: dummy_precip
      integer :: imemory,irun,jlx,kk,kss,plume,ii_plume

      !----------------------------------------------------------------------
      !-do not change this
      itf=ite
      ktf=kte-1
      jtf=jte
      int_time   = int_time + dt
      WHOAMI_ALL = mynum
      time_in    = time
      itime1_in  = itime1
      !----------------------------------------------------------------------

      !-- big loop over j dimension
      do j = jts,jtf
         JCOL = J
          conprr     (:,j) = 0.0
          lightn_dens(:,j) = 0.0
          var2d      (:,j) = 0.0

         !-- initialization
         rtgt(:,:) = 1.0
         ztexec   (:) = 0.0
         zqexec   (:) = 0.0
         fixout_qv(:) = 1.0
            !
            !--- (k,i)
         revsu_gf_2d (:,:) = 0.0
         prfil_gf_2d (:,:) = 0.0
         var3d_agf_2d(:,:) = 0.0
         var3d_bgf_2d(:,:) = 0.0
         Tpert_2d    (:,:) = 0.0
         temp_tendqv (:,:) = 0.0
            !
         omeg   (:,:,:) = 0.0
         !- tendencies (w/ maxiens)
         outt   (:,:,:) = 0.0
         outu   (:,:,:) = 0.0
         outv   (:,:,:) = 0.0
         outq   (:,:,:) = 0.0
         outqc  (:,:,:) = 0.0
         outnice(:,:,:) = 0.0
         outnliq(:,:,:) = 0.0
         outbuoy(:,:,:) = 0.0

         if(APPLY_SUB_MP == 1) then
            !- tendencies (w/ nmp and maxiens)
               outmpqi(:,:,:,:) = 0.0
               outmpql(:,:,:,:) = 0.0
               outmpcf(:,:,:,:) = 0.0
         endif

         if(USE_TRACER_TRANSP==1) then
            out_chem(:,:,:,:) = 0.0
         endif
         !
         if(autoconv == 2) then
            do i= its,itf
               ccn(i) = max( 100., ( 370.37*(0.01+MAX(0.,aot500(i,j))))**1.555 )
            enddo
         else
            do i= its,itf
               ccn(i) = 100.
            enddo
         endif

         do i=its,itf

            xlandi(i) = xland(i,j)!flag < 1 para land
                                  !flag  =1 para water
            psur  (i) = sfc_press(i,j)*1.e-2 ! mbar
            tsur  (i) = temp2m(i,j)
            ter11 (i) = max(0.,topt(i,j))
            kpbli (i) = kpbl(i,j)
            xlons (i) = lons(i,j)*180./3.14159
            xlats (i) = lats(i,j)*180./3.14159
         enddo

         do i=its,itf
            do k=kts,ktf
               kr=k   !+1   !<<<< only kr=k (the input was already converted to the BRAMS vertical grid,
                            !                see cup_grell3.f90 routine)

               !- heigths, current pressure, temp and water vapor mix ratio
               zo      (k,i)  = zt(kr,i,j)*rtgt(i,j)+topt(i,j)
               po      (k,i)  = press(kr,i,j)*1.e-2 !mbar
               temp_old(k,i)  = temp(kr,i,j)

               qv_old  (k,i)  = rvap     (kr,i,j) ! @ begin of the timestep
               qv_curr (k,i)  = curr_rvap(kr,i,j) ! current (after dynamics + physical processes called before GF)

               !- air density, TKE and cloud liq water mixing ratio
               rhoi    (k,i)  = 1.e2*po (k,i)/( 287.04*temp_old(k,i)*(1.+0.608*qv_old(k,i)))
               tkeg    (k,i)  = 1.e-5
               rcpg    (k,i)  = 0.
              
               !- cloud fraction 
               cnvcf2d  (k,i)   =  cnvcf(kr,i,j)

               !- wind velocities
               us      (k,i)  =  u (kr,i,j)
               vs      (k,i)  =  v (kr,i,j)
               dm2d    (k,i)  =  dm(kr,i,j)
               omeg    (k,i,:)= -c_grav*rhoi(k,i)*w(kr,i,j)
               !-buoyancy excess
               buoy_exc2d(k,i)= buoy_exc(kr,i,j)
               !- temp/water vapor modified only by advection
               temp_new_ADV(k,i)= temp_old(k,i)  +  (rth_advten(kr,i,j) )*dt
               qv_new_ADV  (k,i)=   qv_old(k,i)  +  (rqvften   (kr,i,j) )*dt
               !-- turb length scale
               turb_len_scale2d (k,i)  = turb_len_scale (kr,i,j)
            enddo
         enddo

         if(APPLY_SUB_MP == 1) then
            do i=its,itf
               do k=kts,ktf
                  kr=k   !+1   !<<<< only kr=k
                  !- microphysics ice and liq mixing ratio, and cloud fraction of the host model
                  !- (only subsidence is applied)
                  mpqi   (:,k,i) = mp_ice  (:,kr,i,j) ! kg/kg
                  mpql   (:,k,i) = mp_liq  (:,kr,i,j) ! kg/kg
                  mpcf   (:,k,i) = mp_cf   (:,kr,i,j) ! 1
               enddo
            enddo
         endif
         if(USE_TRACER_TRANSP==1) then
            do i=its,itf
               do k=kts,kte
                  kr=k !+1
                  !- atmos composition
                  do ispc=1,mtp
                     se_chem(ispc,k,i) = max(p_mintracer, TRACER(i,j,flip(kr),ispc))
                  enddo
               enddo
            enddo
         endif
         !- pbl  (i) = depth of pbl layer (m)
         !- kpbli(i) = index of zo(k,i)
         do i=its,itf
            pbl  (i)  = zo(kpbli(i),i) - topt(i,j)
         enddo

         !- begin: for GATE soundings-------------------------------------------
         !- this section is intended for model developments only and must
         !- not be used for normal runs.
         if(P_USE_GATE) then
            if(CLEV_GRID == 0) stop "use_gate requires CLEV_GRID 1 or 2"
            if(USE_TRACER_TRANSP==1) then
               ispc_CO=1
               if( .not. allocated(Hcts)) allocate(Hcts(mtp))
               CHEM_NAME_MASK (:) = 1
               !--- dummy initization FSCAV
               do i=1,mtp
                  !FSCAV(i) = 0.1  !km^-1

                  FSCAV(i) = 1.e-5  !km^-1
                  Hcts(i)%hstar  = 0.0 !8.300e+4! 2.4E+3 !59.
                  Hcts(i)%dhr    = 0.0 !7400.   !5000.  !4200.
                  Hcts(i)%ak0    = 0.0
                  Hcts(i)%dak    = 0.0
                   ! H2O2      0.00000      8.300e+4    7400.00000       0.00000       0.00000
                   ! HNO3      0.00000      2.100e+5    8700.00000       0.00000       0.00000
                   ! NH3       0.00000      59.00000    4200.00000       0.00000       0.00000
                   ! SO2       0.00000      2.400e+3    5000.00000       0.00000       0.00000
               enddo
               do i=its,itf
                  se_chem(1:mtp,kts:kpbli(i)-1,i) = 1.+1.e-6
                  do k=kpbli(i),kte
                     se_chem(1:mtp,k,i) = 1.*exp(-(max(0.,0.9*float(k-kpbli(i)))/float(kpbli(i))))+1.e-6
                  enddo
                  do k=kts+1,kte-1
                     se_chem(1:mtp,k,i) = 1./3. *( se_chem(1:mtp,k,i) + se_chem(1:mtp,k-1,i) + se_chem(1:mtp,k+1,i))
                  enddo
               enddo
            endif

            !--- only for GATE soundingg
            if(trim(RUNDATA) == "GATE.dat") then
               jlx= jl
               !jlx= 10 ! to run with only one soundings
               !buoy_exc2d(:,:) = float(jl)*30.
               !print*,"GATE",jl,jlx,buoy_exc2d(1,1)
               !jlx= 42 ! to run with only one soundings

               do i=its,itf
                  do k=kts,kte
                     po       (k,i) = 0.5*(ppres(jlx,k)+ppres(jlx,min(kte,k+1)))
                     temp_old (k,i) = ptemp(jlx,k)+273.15
                     qv_old   (k,i) = pq(jlx,k)/1000.
                     us       (k,i) = pu(jlx,k)
                     vs       (k,i) = pv(jlx,k)
                     omeg     (k,i,:)=pvervel(jlx,k)
                     phil     (k,i) = pgeo(jlx,k)*c_grav   !geo
                     rhoi     (k,i) = 1.e2*po(k,i)/(c_rgas*temp_old(k,i))
                  enddo

                  do k=kts,kte
                     mpql     (:,k,i) = 0.
                     mpql     (:,k,i) = 0.
                     mpcf     (:,k,i) = 0.
                     if(po(k,i) > 900. .or. po(k,i)<300.) cycle
                     pqen  =  exp((-3.e-5*(po(k,i)-550.)**2))
                     pten  =  min(1., (max(0.,(temp_old(k,i)-c_Tice))/(c_T00-c_Tice))**2)
                     mpql  (:,k,i) =3.*pqen* pten
                     mpqi  (:,k,i) =3.*pqen*(1.- pten)
                     mpcf  (:,k,i) = (mpqi  (:,k,i)+mpql  (:,k,i))*100.
                  enddo

                  do k=kts,kte
                     zo       (k,i) = 0.5*(phil(k,i)+phil(min(kte,k+1),i))/c_grav    !meters
                  enddo
                  ter11(i)  = phil(1,i)/c_grav  ! phil is given in g*h.
                  psur (i)  = ppres(jlx,1)
                  tsur (i)  = temp2m(i,j) !temp_old(i,1)
                  kpbli(i)  = 5
                  pbl  (i)  = zo(kpbli(i),i)
                  zws  (i)  = 1.0 ! wstar
                  do k=kts,ktf
                     temp_new(k,i) = temp_old(k,i) + dt *(zadvt(jlx,k)+zqr(jlx,k))/86400.
                     qv_new  (k,i) = qv_old  (k,i) + dt * zadvq(jlx,k)

                     temp_new_dp (k,i) = temp_old(k,i) + dt *(zadvt(jlx,k)+zqr(jlx,k))/86400.
                     qv_new_dp   (k,i) = qv_old  (k,i) + dt * zadvq(jlx,k)

                     temp_new_md (k,i) = temp_new_dp(k,i)
                     qv_new_md   (k,i) = qv_new_dp  (k,i)
                     temp_new_bl (k,i) = temp_new_dp(k,i)
                     qv_new_bl   (k,i) = qv_new_dp  (k,i)
                     temp_new_adv(k,i) = temp_old   (k,i) + dt * zadvt(jlx,k)/86400.
                     qv_new_adv  (k,i) = qv_old     (k,i) + dt * zadvq(jlx,k)
                  enddo
               enddo
            endif
         endif !- end:   for GATE soundings-------------------------------------------
         !
         !- get excess T and Q for source air parcels
         do i=its,itf
            pten = temp_old(1,i)
            pqen = qv_old  (1,i)
            paph = 100.*psur(i)
            zrho = paph/(287.04*(temp_old(1,i)*(1.+0.608*qv_old(1,i))))
            !- sensible and latent sfc fluxes for the heat-engine closure
            h_sfc_flux (i)=zrho*c_cp  *sflux_t(i,j)!W/m^2
            le_sfc_flux(i)=zrho*c_alvl*sflux_r(i,j)!W/m^2
            !
            !- local le and h fluxes for W*
            pahfs=-sflux_t(i,j) *zrho*1004.64  !W/m^2
            pqhfl=-sflux_r(i,j)                !kg/m^2/s
            !- buoyancy flux (h+le)
            zkhvfl= (pahfs/1004.64+0.608*pten*pqhfl)/zrho ! K m s-1
            !- depth of 1st model layer
            !- (zo(1)-top is ~ 1/2 of the depth of 1st model layer, => mult by 2)
            pgeoh =  2.*( zo(1,i)-topt(i,j) )*c_grav ! m+2 s-2
            !-convective-scale velocity w*
            !- in the future, change 0.001 by ustar^3
            zws(i) = max(0.,0.001-1.5*0.41*zkhvfl*pgeoh/pten) ! m+3 s-3
            !
            !-- get LCL properties for parcels from surface level
            call calc_lcl(pten,paph,pqen,tlcl_sfc(i),plcl_sfc(i),zlcl_sfc(i))
            zlcl_sfc(i) = max(zlcl_sfc(i), 0.)
            !print*,'lcl',minval(zlcl_sfc(its:itf)),maxval(zlcl_sfc(its:itf))
            !
            if(zws(i) > tiny(pgeoh)) then
               !-convective-scale velocity w*
               zws(i) = 1.2*zws(i)**.3333
               !- temperature excess
               ztexec(i)     = max(0.,-1.5*pahfs/(zrho*zws(i)*1004.64)) ! K
               !print*,"exce1=",pahfs,zrho,ztexec(i),zws(i),pgeoh,zo(1,i),topt(i,j)
               !call flush(6)
               !- moisture  excess
               zqexec(i)     = max(0.,-1.5*pqhfl/(zrho*zws(i)))        !kg kg-1
            endif   ! zws > 0
            !
            !- zws for shallow convection closure (Grant 2001)
            !- depth of the pbl
            pgeoh = pbl(i)*c_grav
            !-convective-scale velocity W* (m/s)
            zws(i) = max(0.,0.001-1.5*0.41*zkhvfl*pgeoh/pten)
            zws(i) = 1.2*zws(i)**.3333
         enddo
         
         !
         !------ CALL CUMULUS PARAMETERIZATION
         !

         do ii_plume = 1, maxiens

            if(ii_plume == 1) then
               plume = shal
               c0 = c0_shal
            endif
            if(ii_plume == 2) then
               plume = deep
               c0 = c0_deep
            endif
            if(ii_plume == 3) then
               plume = mid
               c0 = c0_mid
            endif

            if(icumulus_gf(plume) /= ON ) cycle

            hei_down_land  =  cum_hei_down_land  (plume)
            hei_down_ocean =  cum_hei_down_ocean (plume)
            hei_updf_land  =  cum_hei_updf_land  (plume)
            hei_updf_ocean =  cum_hei_updf_ocean (plume)
            max_edt_land   =  cum_max_edt_land   (plume)
            max_edt_ocean  =  cum_max_edt_ocean  (plume)
            fadj_massflx   =  cum_fadj_massflx   (plume)
            use_excess     =  cum_use_excess     (plume)
            ave_layer      =  cum_ave_layer      (plume)
            T_star         =  cum_t_star         (plume)
            !print*,"plume=",plume,shal,mid,deep

            !-- set minimum/max for excess of T and Q
            if(use_excess == 0) then
               cum_ztexec(:)= 0.
               cum_zqexec(:)= 0.
            elseif (use_excess == 1) then
               cum_ztexec(:)= ztexec(:)
               cum_zqexec(:)= zqexec(:)
            elseif (use_excess == 2) then
               do i=its,itf
                  cum_zqexec(i)=min(5.e-4, max(1.e-4,zqexec(i)))! kg kg^-1
                  cum_ztexec(i)=min(0.5,   max(0.2  ,ztexec(i)))! Kelvin
               enddo
            else
               do i=its,itf
                  if(xlandi(i) > 0.98) then ! ocean
                     cum_zqexec(i)=min(8.e-4, max(5.e-4,zqexec(i)))! kg kg^-1
                     cum_ztexec(i)=min(1.,    max(0.5  ,ztexec(i)))! Kelvin
                  else                      ! land
                     cum_ztexec(i)= ztexec(i)
                     cum_zqexec(i)= zqexec(i)
                  endif
               enddo
            endif
            !
            !--- replace 'q' and 't' excess in case of use of the cold pool scheme
            !
            if(convection_tracer == 1 .and. plume == deep) then           
                 if(use_gustiness == 1 ) then
                    cum_ztexec(:)= ztexec(:)
                    cum_zqexec(:)= zqexec(:)
                  endif

            !   if(use_gustiness == 1 ) then           
            !      k = 2 ! surface in brams
            !      do i=its,itf
            !       cum_ztexec(i)= (hexcp(k,i,j) - qexcp(k,i,j)*c_alvl)/cp
            !       cum_zqexec(i)= qexcp(k,i,j)
            !      enddo
            !   else 
            !       cum_ztexec(:)= 0.
            !       cum_zqexec(:)= 0.
            !   endif
            endif
            !
            !
            !-- shallow convection
            !
            if(plume == shal) then
               do i=its,itf
                  do k=kts,ktf
                     kr=k!+1 <<<<
                     if(p_use_gate) then
                        dhdt    (k,i)= c_cp*(temp_new_dp(k,i)-temp_old(k,i))+c_alvl*(qv_new_dp(k,i)-qv_old(k,i))
                        temp_new(k,i)= temp_new_dp(k,i)
                        qv_new  (k,i)= qv_new_dp  (k,i)
                     else

                        temp_new(k,i)=temp_old(k,i) + (rthblten(kr,i,j)+rthften(kr,i,j))*dt
                        qv_new  (k,i)=  qv_old(k,i) + (rqvblten(kr,i,j)+rqvften(kr,i,j))*dt
                        qv_new  (k,i)= max(p_smaller_qv,qv_new  (k,i))

                        !- only pbl forcing changes moist static energy
                        dhdt(k,i)= c_cp   *(rthblten(kr,i,j)) +  c_alvl  *(rqvblten(kr,i,j))

                        !- all forcings change moist static energy
                        dhdt(k,i)=dhdt(k,i) + c_cp*rthften(kr,i,j) + c_alvl*rqvften(kr,i,j)

                     endif
                  enddo
               enddo
            endif
            !
            !--- deep convection
            if(plume == deep) then

               if(p_use_gate) then
                  do k=kts,ktf
                     do i=its,itf
                        temp_new(k,i) = temp_new_dp(k,i)
                        qv_new  (k,i) = qv_new_dp  (k,i)
                     enddo
                  enddo
               else
                  do k=kts,ktf
                     do i=its,itf
                        kr=k!+1 <<<<
                        temp_new    (k,i)= temp_old(k,i)  +  (rthblten(kr,i,j) + rthften (kr,i,j))*dt
                        qv_new      (k,i)=   qv_old(k,i)  +  (rqvblten(kr,i,j) + rqvften (kr,i,j))*dt

                        temp_new_BL (k,i)= temp_old(k,i)  +  (rthblten(kr,i,j) )*dt
                        qv_new_BL   (k,i)=   qv_old(k,i)  +  (rqvblten(kr,i,j) )*dt
                     enddo
                  enddo
               endif
            endif
            !
            !--- mid/congestus type convection
            if(plume == mid) then

               if(p_use_gate) then
                  do k=kts,ktf
                     do i=its,itf
                        temp_new(k,i) = temp_new_dp(k,i)
                        qv_new  (k,i) = qv_new_dp  (k,i)
                     enddo
                  enddo
               else
                  do i=its,itf
                     do k=kts,ktf
                        kr=k!+1 <<<<

                        temp_new(k,i)=temp_old(k,i) + (rthblten(kr,i,j)+rthften(kr,i,j))*dt
                        qv_new  (k,i)=  qv_old(k,i) + (rqvblten(kr,i,j)+rqvften(kr,i,j))*dt
                        qv_new  (k,i)= max(p_smaller_qv,qv_new  (k,i))

                        !- only pbl forcing changes moist static energy
                        dhdt(k,i)= c_cp   *(rthblten(kr,i,j)) +  c_alvl  *(rqvblten(kr,i,j))

                        !- all forcings change moist static energy
                        dhdt(k,i)=dhdt(k,i) + c_cp*rthften(kr,i,j) + c_alvl*rqvften(kr,i,j)

                        !- temp/water vapor modified only by bl processes
                        temp_new_BL(k,i)= temp_old(k,i)  +  (rthblten(kr,i,j) )*dt
                        qv_new_BL  (k,i)= qv_old  (k,i)  +  (rqvblten(kr,i,j) )*dt

                     enddo
                  enddo
               endif
            endif
            !

            call CUP_GF(its,ite,kts,kte, itf,ktf, mtp, nmp, FSCAV  &
                        ,cumulus_type  (plume)            &
                        ,closure_choice(plume)            &
                        ,cum_entr_rate (plume)            &
                        ,cum_use_excess(plume)            &
                        !- input data
                        ,dx2d          (:,j)              &
                        ,stochastic_sig(:,j)              &
                        ,col_sat       (:,j)              &
                        ,tke_pbl       (:,j)              &
                        ,rh_dicycle_fct(:,j)              &
                        ,wlpool        (:,j)              &
                        ,dt                               &
                        ,kpbli                            &
                        ,cum_ztexec                       &
                        ,cum_zqexec                       &
                        ,ccn                              &
                        ,rhoi                             &
                        ,omeg                             &
                        ,temp_old                         &
                        ,qv_old                           &
                        ,ter11                            &
                        , h_sfc_flux                      &
                        ,le_sfc_flux                      &
                        ,zlcl_sfc                         &
                        ,xlons                            &
                        ,xlats                            &
                        ,xlandi                           &
                        ,temp_new                         &
                        ,qv_new                           &
                        ,temp_new_BL                      &
                        ,qv_new_BL                        &
                        ,temp_new_ADV                     &
                        ,qv_new_ADV                       &
                        ,zo                               &
                        ,po                               &
                        ,tsur                             &
                        ,psur                             &
                        ,us                               &
                        ,vs                               &
                        ,dm2d                             &
                        ,se_chem                          &
                        ,zws                              &
                        ,dhdt                             &
                        ,buoy_exc2d                       &
                        ,cnvcf2d                          &
                        ,turb_len_scale2d                 &
                        ,mpqi                             &
                        ,mpql                             &
                        ,mpcf                             &
                        !output data
                        ,outt                 (:,:,plume) &
                        ,outq                 (:,:,plume) &
                        ,outqc                (:,:,plume) &
                        ,outu                 (:,:,plume) &
                        ,outv                 (:,:,plume) &
                        ,outnliq              (:,:,plume) &
                        ,outnice              (:,:,plume) &
                        ,outbuoy              (:,:,plume) &
                        ,outmpqi            (:,:,:,plume) &
                        ,outmpql            (:,:,:,plume) &
                        ,outmpcf            (:,:,:,plume) &
                        ,out_chem           (:,:,:,plume) &
                        !- for convective transport
                        ,ierr4d               (:,j,plume) &
                        ,jmin4d               (:,j,plume) &
                        ,klcl4d               (:,j,plume) &
                        ,k224d                (:,j,plume) &
                        ,kbcon4d              (:,j,plume) &
                        ,ktop4d               (:,j,plume) &
                        ,kstabi4d             (:,j,plume) &
                        ,kstabm4d             (:,j,plume) &
                        ,cprr4d               (:,j,plume) &
                        ,xmb4d                (:,j,plume) &
                        ,edt4d                (:,j,plume) &
                        ,pwav4d               (:,j,plume) &
                        ,sigma4d              (:,j,plume) &
                        ,pcup5d             (:,:,j,plume) &
                        ,up_massentr5d      (:,:,j,plume) &
                        ,up_massdetr5d      (:,:,j,plume) &
                        ,dd_massentr5d      (:,:,j,plume) &
                        ,dd_massdetr5d      (:,:,j,plume) &
                        ,zup5d              (:,:,j,plume) &
                        ,zdn5d              (:,:,j,plume) &
                        ,prup5d             (:,:,j,plume) &
                        ,prdn5d             (:,:,j,plume) &
                        ,clwup5d            (:,:,j,plume) &
                        ,tup5d              (:,:,j,plume) &
                        ,conv_cld_fr5d      (:,:,j,plume) &
                        !-- for debug/diag
                        ,AA0(:,j),AA1(:,j),AA1_ADV(:,j),AA1_RADPBL(:,j),AA2(:,j),AA3(:,j) &
                        ,AA1_BL(:,j),AA1_CIN(:,j),TAU_BL(:,j),TAU_EC(:,j) &
                        !-- for diag
                        ,lightn_dens  (:,j)               &
                        ,var2d        (:,j)               &
                        ,revsu_gf_2d                      &
                        ,prfil_gf_2d                      &
                        ,var3d_agf_2d                     &
                        ,var3d_bgf_2d                     &
                        ,Tpert_2d                         &
                        )
         enddo !- plume

         !--- reset ierr4d to value different of zero in case the correspondent
         !--- plume (shalllow, congestus, deep) was not actually used
         do n=1,maxiens
            if(icumulus_gf(n) == OFF ) ierr4d (:,j,n) = -99
         enddo

         do i=its,itf
            do_this_column(i,j) = 0
            loop1:  do n=1,maxiens
               if(ierr4d (i,j,n) == 0 ) then
                  do_this_column(i,j) = 1
                  exit loop1
               endif
            enddo loop1
         enddo
         !----------- check for negative water vapor mix ratio
         do i=its,itf
            if(do_this_column(i,j) == 0) cycle
            do k = kts,ktf
               temp_tendqv(k,i)= outq (k,i,shal) + outq (k,i,deep) + outq (k,i,mid )
            enddo

            do k = kts,ktf
               distance(k)= qv_curr(k,i) + temp_tendqv(k,i) * dt
            enddo

            if(minval(distance(kts:ktf)) < 0.) then
               zmax   =  MINLOC(distance(kts:ktf),1)

               if( abs(temp_tendqv(zmax,i) * dt) <  p_mintracer) then
                  fixout_qv(i)= 0.999999
                !fixout_qv(i)= 0.
               else
                  fixout_qv(i)= ( (p_smaller_qv - qv_curr(zmax,i))) / (temp_tendqv(zmax,i) *dt)
               endif
               fixout_qv(i)=max(0.,min(fixout_qv(i),1.))
            endif
         enddo
         !------------ feedback
         !-- deep convection
         do i=its,itf
            if(do_this_column(i,j) == 0) cycle
            cprr4d(i,j,deep) =  cprr4d(i,j,deep)* fixout_qv(i)
            cprr4d(i,j,mid)  =  cprr4d(i,j,mid) * fixout_qv(i)
            cprr4d(i,j,shal) =  cprr4d(i,j,shal)* fixout_qv(i)
            conprr(i,j)= (cprr4d(i,j,deep) + cprr4d(i,j,mid) + cprr4d(i,j,shal))
            conprr(i,j)= max(0.,conprr(i,j))
         enddo

         !-- deep + shallow + mid convection
         do i = its,itf
            if(do_this_column(i,j) == 0) cycle
            do k = kts,kte
               kr=k!+1
               !- feedback the tendencies from convection
               RTHCUTEN (kr,i,j)= (outt (k,i,shal) + outt (k,i,deep) + outt (k,i,mid)) *fixout_qv(i)

               RQVCUTEN (kr,i,j)= (outq (k,i,shal) + outq (k,i,deep) + outq (k,i,mid)) *fixout_qv(i)

               RQCCUTEN (kr,i,j)= (outqc(k,i,shal) + outqc(k,i,deep) + outqc(k,i,mid)) *fixout_qv(i)

               REVSU_GF (kr,i,j)= revsu_gf_2d(k,i)*fixout_qv(i) !-- already contains deep and mid amounts.

               !---these arrays are only for the deep plume mode
               PRFIL_GF (kr,i,j)= prfil_gf_2d (k,i)*fixout_qv(i) !-- ice/liq prec flux of the deep plume
               !VAR3d_aGF(kr,i,j)= var3d_gf_2d(k,i)              !-- vertical velocity of the deep plume
               VAR3d_aGF(kr,i,j)= outt (k,i,mid)*fixout_qv(i)   !--
               VAR3d_bGF(kr,i,j)= outq (k,i,mid)*fixout_qv(i)   !--

               if(icumulus_gf(shal) == OFF) then
                  VAR3d_cGF(kr,i,j)= outqc (k,i,deep)*fixout_qv(i)  !--
                  VAR3d_dGF(kr,i,j)= outqc (k,i,mid) *fixout_qv(i)  !--
               else
                  VAR3d_cGF(kr,i,j)= outt (k,i,shal)*fixout_qv(i)   !--
                  VAR3d_dGF(kr,i,j)= outq (k,i,shal)*fixout_qv(i)   !--
               endif

            enddo
         enddo
         if(USE_MOMENTUM_TRANSP > 0) then
            do i = its,itf
               if(do_this_column(i,j) == 0) cycle
               do k = kts,kte
                  kr=k!+1
                  RUCUTEN (kr,i,j)= (outU(k,i,deep)+outU(k,i,mid)+outU(k,i,shal)) *fixout_qv(i)
                  RVCUTEN (kr,i,j)= (outV(k,i,deep)+outV(k,i,mid)+outV(k,i,shal)) *fixout_qv(i)
               enddo
            enddo
         endif

         if(APPLY_SUB_MP == 1) then
            do i = its,itf
               if(do_this_column(i,j) == 0) cycle
               do k = kts,kte
                  kr=k!+1
                  SUB_MPQL (:,kr,i,j)= (outmpql(:,k,i,deep)+outmpql(:,k,i,mid)+outmpql(:,k,i,shal)) *fixout_qv(i)
                  SUB_MPQI (:,kr,i,j)= (outmpqi(:,k,i,deep)+outmpqi(:,k,i,mid)+outmpqi(:,k,i,shal)) *fixout_qv(i)
                  SUB_MPCF (:,kr,i,j)= (outmpcf(:,k,i,deep)+outmpcf(:,k,i,mid)+outmpcf(:,k,i,shal)) *fixout_qv(i)
               enddo
            enddo
         endif

         if(LIQ_ICE_NUMBER_CONC == 1) then
            do i = its,itf
               if(do_this_column(i,j) == 0) cycle
               do k = kts,kte
                  kr=k!+1
                  RNICUTEN (kr,i,j)= (outnice(k,i,shal) + outnice(k,i,deep) + outnice(k,i,mid)) *fixout_qv(i)
                  RNLCUTEN (kr,i,j)= (outnliq(k,i,shal) + outnliq(k,i,deep) + outnliq(k,i,mid)) *fixout_qv(i)
               enddo
            enddo
         endif

         if(USE_TRACER_TRANSP==1) then
            do i = its,itf
               if(do_this_column(i,j) == 0) cycle
               do k = kts,kte
                  kr=k!+1
                  RCHEMCUTEN (:,kr,i,j)= (out_CHEM(:,k,i,deep) +out_CHEM(:,k,i,mid)+out_CHEM(:,k,i,shal)) *fixout_qv(i)
               enddo
            enddo

            !- constrain positivity for tracers
            do i = its,itf
               if(do_this_column(i,j) == 0) cycle

               do ispc=1,mtp
                  if(CHEM_NAME_MASK (ispc) == 0 ) cycle

                  do k=kts,ktf
                     distance(k)= se_chem(ispc,k,i) + RCHEMCUTEN(ispc,k,i,j)* dt
                  enddo

                  !-- fixer for mass of tracer
                  if(minval(distance(kts:ktf)) < 0.) then
                     zmax   =  MINLOC(distance(kts:ktf),1)

                     if( abs(RCHEMCUTEN(ispc,zmax,i,j)*dt) <  p_mintracer) then
                        fixouts= 0.999999
                      !fixouts= 0.
                     else
                        fixouts=  ( (p_mintracer - se_chem(ispc,i,zmax))) / (RCHEMCUTEN(ispc,zmax,i,j)*dt)
                     endif
                     if(fixouts > 1. .or. fixouts <0.)fixouts=0.

                     RCHEMCUTEN(ispc,kts:ktf,i,j)=fixouts*RCHEMCUTEN(ispc,kts:ktf,i,j)
                  endif
               enddo
            enddo
         endif

         if(CONVECTION_TRACER==1) then
            do i = its,itf
               if(do_this_column(i,j) == 0) cycle
               do k = kts,kte
                  kr=k!+1
                  RBUOYCUTEN (kr,i,j)= (outbuoy(k,i,deep)+outbuoy(k,i,mid)+outbuoy(k,i,shal)) *fixout_qv(i)
               enddo
            enddo
         endif

  !do i = its,itf
  ! AA1_ADV(i,j)=zt(ktop4d(i,j,deep),i,j)*rtgt(i,j)
  ! !AA0    (i,j)=zlcl_sfc(i)
  !enddo

      enddo



   end subroutine convParGFDriver
   !---------------------------------------------------------------------------------------------------

   subroutine CUP_GF (its,ite,kts,kte ,itf,ktf, mtp, nmp &
                     ,fscav             &
                     ,cumulus           &
                     ,ichoice           &
                     ,entr_rate_input   &
                     ,use_excess        &
                     !input data
                     ,dx                &
                     ,stochastic_sig    &
                     ,col_sat           &
                     ,tke_pbl           &
                     ,rh_dicycle_fct    &
                     ,wlpool            &
                     ,dtime             &
                     ,kpbl              &
                     ,ztexec            &
                     ,zqexec            &
                     ,ccn               &
                     ,rho               &
                     ,omeg              &
                     ,t                 &
                     ,q                 &
                     ,z1                &
                     , h_sfc_flux       &
                     ,le_sfc_flux       &
                     ,zlcl_sfc          &
                     ,xlons             &
                     ,xlats             &
                     ,xland             &
                     ,tn                &
                     ,qo                &
                     ,tn_bl             &
                     ,qo_bl             &
                     ,tn_adv            &
                     ,qo_adv            &
                     ,zo                &
                     ,po                &
                     ,tsur              &
                     ,psur              &
                     ,us                &
                     ,vs                &
                     ,dm2d              &
                     ,se_chem           &
                     ,zws               &
                     ,dhdt              &
                     ,buoy_exc          &
                     ,cnvcf             &
                     ,turb_len_scale    &
                     ,mpqi              &
                     ,mpql              &
                     ,mpcf              &
                     !output data
                     ,outt              &
                     ,outq              &
                     ,outqc             &
                     ,outu              &
                     ,outv              &
                     ,outnliq           &
                     ,outnice           &
                     ,outbuoy           &
                     ,outmpqi           &
                     ,outmpql           &
                     ,outmpcf           &
                     ,out_chem          &
                     !- for convective transport
                     ,ierr              &
                     ,jmin              &
                     ,klcl              &
                     ,k22               &
                     ,kbcon             &
                     ,ktop              &
                     ,kstabi            &
                     ,kstabm            &
                     ,pre               &
                     ,xmb               &
                     ,edto              &
                     ,pwavo             &
                     ,sig               &
                     ,po_cup            &
                     ,up_massentro      &
                     ,up_massdetro      &
                     ,dd_massentro      &
                     ,dd_massdetro      &
                     ,zuo               &
                     ,zdo               &
                     ,pwo               &
                     ,pwdo              &
                     ,qrco              &
                     ,tup               &
                     ,clfrac            &
                     !- for convective transport-end
                     !- for debug/diag
                     ,AA0_,AA1_,AA1_ADV_,AA1_RADPBL_,AA2_,AA3_,AA1_BL_,AA1_CIN_,TAU_BL_,TAU_EC_   &
                     ,lightn_dens       &
                     ,var2d             &
                     ,revsu_gf          &
                     ,prfil_gf          &
                     ,var3d_agf         &
                     ,var3d_bgf         &
                     ,Tpert             &
                     )
      implicit none

      !-local settings
      logical, parameter:: USE_INV_LAYERS=.true. 

      character*(*),intent(in) :: cumulus
      integer      ,intent(in) :: itf,ktf,its,ite,kts,kte,ichoice,use_excess,mtp, nmp
      integer      ,intent(inout),  dimension (:) ::   kpbl
      !
      ! outtem = output temp tendency (per s)
      ! outq   = output q tendency (per s)
      ! outqc  = output qc tendency (per s)
      ! pre    = output precip
      real,    dimension (:,:)  ,intent (inout)  :: outu,outv,outt,outq,outqc&
                                                   ,outbuoy,outnliq ,outnice &
                                                   ,revsu_gf,prfil_gf,var3d_agf ,var3d_bgf

      real,    dimension (:)    ,intent (out  )  :: pre,sig,lightn_dens,var2d
      !
      ! basic environmental input
      !
      real,    dimension (:,:)  ,intent (inout)  :: dhdt,rho,t,po,us,vs,tn,dm2d &
                                                   ,buoy_exc,tn_bl,tn_adv,cnvcf,turb_len_scale
      
      real,    dimension (:,:,:),intent (inout)  :: omeg
      real,    dimension (:,:)  ,intent (inout)  :: q,qo,Tpert,qo_bl,qo_adv
      real,    dimension (:)    ,intent (inout)  :: ccn,z1,psur,xland,xlons,xlats  &
                                                   ,h_sfc_flux,le_sfc_flux,tsur,dx &
                                                   ,zlcl_sfc

      real,    dimension (:)    ,intent (in   )  :: col_sat,stochastic_sig,tke_pbl
      real,    dimension (:)    ,intent (inout)  :: zws,ztexec,zqexec, rh_dicycle_fct,wlpool
      real                      ,intent (in   )  :: dtime,entr_rate_input
      real,    dimension (:,:,:),intent (inout)  :: mpqi,mpql,mpcf
      real,    dimension (:,:,:),intent (inout)  :: outmpqi,outmpql,outmpcf
      real,    dimension (:)    ,intent (inout)  :: aa0_,aa1_,aa2_,aa3_,aa1_bl_&
                                                   ,aa1_cin_,tau_bl_,tau_ec_   &
                                                   ,aa1_radpbl_,aa1_adv_
      integer, dimension (:)     ,intent (inout)  :: &
          ierr              &
         ,jmin              &
         ,klcl              &
         ,k22               &
         ,kbcon             &
         ,ktop              &
         ,kstabi            &
         ,kstabm

      real,    dimension (:)     ,intent (inout)  :: &
          xmb               &
         ,edto              &
         ,pwavo
      real,    dimension (:,:)   ,intent (inout)  :: &
          po_cup            &
         ,up_massentro      &
         ,up_massdetro      &
         ,dd_massentro      &
         ,dd_massdetro      &
         ,zuo               &
         ,zdo               &
         ,pwo               &
         ,pwdo              &
         ,qrco              &
         ,tup               &
         ,clfrac

      !------------------------------------------------------
      ! local ensemble dependent variables in this routine
      real,    dimension (its:ite,1:maxens2) ::   edtc
      real,    dimension (its:ite,1:ensdim)  ::   xf_ens,pr_ens
      !
      !*******the following are your basic environmental
      !          variables. They carry a "_cup" if they are
      !          on model cloud levels (staggered). They carry
      !          an "o"-ending (z becomes zo), if they are the forced
      !          variables. They are preceded by x (z becomes xz)
      !          to indicate modification by some typ of cloud
      !
      ! z           = heights of model levels
      ! q           = environmental mixing ratio
      ! qes         = environmental saturation mixing ratio
      ! t           = environmental temp
      ! p           = environmental pressure
      ! he          = environmental moist static energy
      ! hes         = environmental saturation moist static energy
      ! z_cup       = heights of model cloud levels
      ! q_cup       = environmental q on model cloud levels
      ! qes_cup     = saturation q on model cloud levels
      ! t_cup       = temperature (Kelvin) on model cloud levels
      ! p_cup       = environmental pressure
      ! he_cup = moist static energy on model cloud levels
      ! hes_cup = saturation moist static energy on model cloud levels
      ! gamma_cup = gamma on model cloud levels
      !
      !
      ! hcd = moist static energy in downdraft
      ! zd normalized downdraft mass flux
      ! dby = buoancy term
      ! entr = entrainment rate
      ! zd   = downdraft normalized mass flux
      ! entr= entrainment rate
      ! hcd = h in model cloud
      ! bu = buoancy term
      ! zd = normalized downdraft mass flux
      ! gamma_cup = gamma on model cloud levels
      ! qcd = cloud q (including liquid water) after entrainment
      ! qrch = saturation q in cloud
      ! pwd = evaporate at that level
      ! pwev = total normalized integrated evaoprate (I2)
      ! entr= entrainment rate
      ! z1 = terrain elevation
      ! entr = downdraft entrainment rate
      ! jmin = downdraft originating level
      ! kdet = level above ground where downdraft start detraining
      ! psur    = surface pressure
      ! z1      = terrain elevation
      ! pr_ens  = precipitation ensemble
      ! xf_ens  = mass flux ensembles
      ! massfln = downdraft mass flux ensembles used in next timestep
      ! omeg    = omega from large scale model
      ! mconv   = moisture convergence from large scale model
      ! zd      = downdraft normalized mass flux
      ! zu      = updraft normalized mass flux
      ! dir     = "storm motion"
      ! mbdt    = arbitrary numerical parameter
      ! dtime   = dt over which forcing is applied
      ! iact_gr_old = flag to tell where convection was active
      ! kbcon       = LFC of parcel from k22
      ! k22         = updraft originating level
      ! ichoice     = flag if only want one closure (usually set to zero!)
      ! dby  = buoancy term
      ! ktop = cloud top (output)
      ! xmb  = total base mass flux
      ! hc   = cloud moist static energy
      ! hkb  = moist static energy at originating level
      ! cd  = detrainment function for updraft
      ! cdd = detrainment function for downdraft
      ! dellat = change of temperature per unit mass flux of cloud ensemble
      ! dellaq = change of q per unit mass flux of cloud ensemble
      ! dellaqc = change of qc per unit mass flux of cloud ensemble

      character*128 :: ierrc(its:ite)      

      real :: dsubh_aver,dellah_aver,x_add,cap_max_inc
      real :: day,dz,dzo,radius,entrd_rate,zcutdown,depth_min,zkbmax,z_detr,zktop       &
        ,massfld,dh,trash,frh,xlamdd,radiusd,frhd,effec_entrain,detdo1,detdo2,entdo     &
        ,dp,subin,detdo,entup,detup,subdown,entdoj,entupk,detupk,totmas,min_entr_rate   &
        ,tot_time_hr,beta,env_mf,env_mf_p,env_mf_m,dts,denom,denomU,umean,xh_env_eff

      integer :: iversion
      integer :: ipr=0,jpr=0
      integer :: k,i,iedt,nens,nens3
      integer :: vtp_index
      
      integer, dimension (its:ite) :: kzdown,kdet,kb, kbmax,start_level
      real,    dimension (its:ite,1:maxens3) ::  xff_mid
      real,    dimension (its:ite,shall_closures) :: xff_shal

      ! aa0 cloud work function for downdraft
      ! edt = epsilon
      ! aa0     = cloud work function without forcing effects
      ! aa1     = cloud work function with forcing effects
      ! xaa0    = cloud work function with cloud effects
      ! edt     = epsilon

      real,    dimension (its:ite) ::                                                      &
         edt,aa1,aa0,xaa0,hkb,hkbo,xhkb,qkb, pwevo,bu,bud,cap_max,xland1,vshear            &
        ,cap_max_increment,psum,psumh,sigd,mconv,rescale_entrain,entr_rate,mentrd_rate     &
        ,v_ratio, aa0_bl,aa1_bl,tau_bl,tau_ecmwf,wmean,aa1_fa,aa1_tmp    &
        ,aa2,aa3,cin0,cin1,edtmax,edtmin,aa_tmp,aa_ini,aa_adv   &
        ,daa_adv_dt,wlpool_bcon,xk_x,xf_dicycle,mbdt,xf_coldpool, ke_gustfront& 
        ,vvel1d, x_add_buoy,lambau_dn,lambau_dp,q_wetbulb,t_wetbulb,col_sat_adv &
        ,Q_adv,alpha_adv,aa1_radpbl,aa1_adv,p_cwv_ave,cape, depth_neg_buoy,frh_bcon &
        ,check_sig,random,rh_entr_factor,rntot,delqev,delq2,qevap,rn,qcond,rainevap

      real,    dimension (kts:kte,its:ite) ::                                         &
          entr_rate_2d,mentrd_rate_2d, he, hes, qes, z,heo,heso,qeso,zo               &
         ,xhe,xhes,xqes,xz,xt,xq, qes_cup, q_cup, he_cup, hes_cup, z_cup, p_cup       &
         ,gamma_cup, t_cup,qeso_cup,qo_cup,heo_cup,heso_cup,zo_cup,gammao_cup,tn_cup  &
         ,xqes_cup,xq_cup,xhe_cup,xhes_cup,xz_cup,xt_cup,hcot,evap_bcb,dby,hc,clw_all &                   
         ,dbyo,qco,qrcdo,hcdo,qcdo,dbydo,hco,xdby,xhc,cupclw,pwo_eff                  &
         ,cd,cdd,dellah,dellaq,dellat,dellaqc,dsubq,dsubh,dellabuoy,u_cup,v_cup       &
         ,uc,vc,ucd,vcd,dellu,dellv,subten_H,subten_Q,subten_T,vvel2d,tempco,tempcdo  &
         ,p_liq_ice,melting_layer,melting ,up_massentru,up_massdetru,dd_massentru     &
         ,dd_massdetru, prec_flx,evap_flx,qrs,zenv

      real, dimension (nmp,kts:kte,its:ite) ::  dellampqi,dellampql,dellampcf

      !-- atmos composition arrays
      real, dimension (: )   ,intent (in)    ::   fscav   !(mtp )  
      real, dimension (:,:,:),intent (inout) ::   se_chem !(mtp,kts:kte,its:ite)
      real, dimension (:,:,:),intent (inout) ::   out_chem!(mtp,kts:kte,its:ite)

      integer :: ispc,kmp,istep,lstep
      real, dimension (mtp,kts:kte,its:ite) ::  se_cup_chem,sc_up_chem,sc_dn_chem,pw_up_chem,pw_dn_chem
      real, dimension (mtp,its:ite)         ::  tot_pw_up_chem,tot_pw_dn_chem
      real                                  ::  evap_(mtp),wetdep_(mtp),trash_(mtp),trash2_(mtp) &
                                               ,massi,massf,dtime_max,evap,wetdep,residu_(mtp)   &
                                               ,s1,s2,q1,q2,rzenv,factor,CWV,entr_threshold      &
                                               ,resten_H,resten_Q,resten_T,alp0,beta1,beta2 
      !----------------------------------------------------------------------
      !-- only for debug (atmos composition)
      real, allocatable, dimension (:,:,:),save    ::   se_chem_update

      !--only for debug
      if(p_use_gate) then
         if( .not. allocated(se_chem_update)) allocate(se_chem_update(3,kts:kte,its:ite))
         if(jl==1) then
         !    se_chem_update(1,:,:) = mpql (lsmp,:,:)
         !    se_chem_update(2,:,:) = mpqi (lsmp,:,:)
         !    se_chem_update(3,:,:) = mpcf (lsmp,:,:)
         else
         !    mpql (lsmp,:,:)= se_chem_update(1,:,:)
         !    mpqi (lsmp,:,:)= se_chem_update(2,:,:)
         !    mpcf (lsmp,:,:)= se_chem_update(3,:,:)
         endif
      endif
      !----------------------------------------------------------------------
      !
      !-- init the vector vec_ok with the all indexes to process
      vec_max_size = ite - its + 1
      call init(vec_ok, vec_max_size)
      call insert_range(vec_ok, its, ite)
      !-- vec removed will be inserted when removing
      call init(vec_removed, vec_max_size)

      !
      !--- maximum depth (mb) of capping inversion (larger cap = no convection)
      call get_capmax(cumulus,itf,ktf,its,ite,kts,kte,cap_max_inc &
                     ,cap_max_increment,cap_max,cap_maxs,moist_trigger)
      !     
      !--- lambda_U parameter for momentum transport
      call get_lambdaU(cumulus,itf,ktf,its,ite,kts,kte,lambau_dp,lambau_dn        &
                      ,lambau_deep,lambau_shdn,pgcon)
      !
      !-- init/reset 1-d and 2-d local vars
      call reset_1d(its,ite,ierrc,xland,xland1,aa0,aa1,aa2,aa3                    &
                   ,aa1_bl,aa1_fa,aa0_bl,q_adv,aa1_radpbl,aa1_adv,alpha_adv,cin1  &
                   ,xk_x,edt,edto,tau_bl,q_wetbulb,t_wetbulb,tau_ecmwf,xf_dicycle &
                   ,x_add_buoy,xf_coldpool,wlpool_bcon,ke_gustfront,random,mbdt)

      call reset_2d(its,ite,kts,kte,zo,z,xz,hcdo,cupclw,qrcdo,hcot,xf_ens,pr_ens  &
                   ,evap_bcb,uc,vc,hc,hco, zuo,zdo,zenv)
      !  
      !---  create a real random number in the interval [-use_random_num, +use_random_num]
      if( cumulus == 'deep' .and. use_random_num > 1.e-6) &
         call gen_random(its,ite,use_random_num,random)
      !
      !-- define limits of evaporation by the downdrafts
      call get_edt(cumulus,itf,ktf,its,ite,kts,kte,xland,edtmin,edtmax,max_edt_ocean,max_edt_land &
                  ,c0_mid)  
      !
      !--- environmental conditions, FIRST HEIGHTS
      !--- calculate moist static energy, heights, qes
      !
      call cup_env(z ,qes ,he ,hes ,t ,q ,po,z1 ,psur,ierr,-1,itf,ktf,its,ite, kts,kte)
      call cup_env(zo,qeso,heo,heso,tn,qo,po,z1, psur,ierr,-1,itf,ktf,its,ite, kts,kte)

      !
      !--- outputs a model sounding for the stand-alone code (part 1)
      !
      if(OUTPUT_SOUND == 1) then
         call SOUND(1,cumulus, int_time, dtime, ens4, itf, its, ite, kts, kte, xlats, xlons, jcol, whoami_all &
                    , z, qes, he, hes, t, q, po, z1, psur, zo, qeso, heo, heso, tn, qo, us, vs, omeg, xz &
                    , h_sfc_flux, le_sfc_flux, tsur, dx, stochastic_sig, zws, ztexec, zqexec, xland &
                    , kpbl, k22, klcl, kbcon, ktop, aa0, aa1, sig, xaa0, hkb, xmb, pre, edto &
                    , zo_cup, dhdt, rho, zuo, zdo, up_massentro, up_massdetro, outt, outq, outqc, outu, outv)
      endif 
      !
      !--- environmental values on cloud levels
      !
      call cup_env_clev(t,qes,q,he,hes,z,po,qes_cup,q_cup,he_cup,us,vs,u_cup,v_cup &
                       , hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur,tsur,ierr,z1     &
                       ,itf,ktf,its,ite, kts,kte)

      call cup_env_clev(tn,qeso,qo,heo,heso,zo,po,qeso_cup,qo_cup, heo_cup,us,vs,u_cup,v_cup   &
                       ,heso_cup,zo_cup,po_cup,gammao_cup ,tn_cup,psur,tsur,ierr,z1            &
                       ,itf,ktf,its,ite, kts,kte)
      !
      !-- get the pickup of ensemble ave prec, following Neelin et al 2009.
      !
      call precip_cwv_factor(itf,ktf,its,ite,kts,kte,ierr,tn,po,qo,po_cup,cumulus,p_cwv_ave)
      !
      !
      !-- partition between liq/ice cloud contents
      !
      call get_partition_liq_ice(ierr,tn,z1,zo_cup,po_cup,p_liq_ice,melting_layer &
                                ,itf,ktf,its,ite,kts,kte,cumulus)
      !
      !-- get several indexes
      !
      call get_kbmax_kdet_k22(cumulus,itf,ktf,its,ite,kts,kte,ierr,ierrc,z1,zo_cup,heo_cup&
                             ,depth_min,z_detr,zkbmax,kbmax,kdet,k22,kstabm)
      !
      !-- define entrainment/detrainment profiles for updrafts
      !
      !- initial entrainment/detrainment
      entr_rate   (:  ) = entr_rate_input
      min_entr_rate     = entr_rate_input * 0.1
      !
      !-- controls of LCL and PBL turbulence on the entrainment rate
      if(use_lcl_ctrl_entr > 0 ) then 
         call get_lcl2(cumulus,ave_layer,its,ite,itf,kts,kte,ktf,ierr,zqexec,ztexec,xland &
                      ,po,t_cup,p_cup,z_cup,q_cup,k22,klcl,kpbl,psur,zlcl_sfc)
         call LCL_and_PBL_ctrl_on_entrainment(cumulus,its,ite,itf,kts,kte,ktf,min_entr_rate &
                                             ,entr_rate,zlcl_sfc,turb_len_scale,AA0_,AA1_,AA1_BL_)
      endif
      !
      !--- determine the entrainment dependent on environmental moist (here relative humidity)
      !--- also the controls of RH on the diurnal cycle (see Tian et al 2022 GRL)
      if(cumulus == 'deep') &
          call rh_controls(whoami_all,itf,ktf,its,ite,kts,kte,ierr,tn,po,qo,qeso,po_cup,cumulus,rh_entr_factor, &
                           rh_dicycle_fct,entr_rate_input, entr_rate ,xlons,dtime)         
      !
      !-- cold pool parameterization and convective memory
      !
      if (convection_tracer == 1 .and. trim(cumulus) == 'deep') then
         call coldPoolConvMem(cumulus,its, itf, kts, kte, ktf, aa2_, entr_rate_input, ztexec, zqexec&
                            , xland, po, buoy_exc, ierr, cap_max, wlpool, aa3_, x_add_buoy          &
                            , min_entr_rate, entr_rate)
      end if
      !
      !-- determine LCL for the air parcels around K22
      !
      call get_lcl(cumulus,convection_tracer,use_memory,ave_layer,its,ite,itf,kts,kte,ktf,ierr &
                  ,x_add_buoy,zqexec,ztexec,xland,po,t_cup,p_cup,z_cup,q_cup,k22,klcl)

      !--- start_level
      !
      start_level(:)=  KLCL(:) !start_level(:)=  KTS
      !
      !-- check if LCL height is below PBL height to allow shallow convection
      !
      if(lcl_trigger > 0 .and. cumulus == 'shallow')then
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            if(klcl(i) > max(1,kpbl(i)-lcl_trigger)) then
                  ierr(i)=21 ; is_removed  = remove(vec_ok, i)
                  ierrc(i)='for shallow convection:  LCL height < PBL height'
            endif
         enddo
         !print*,"LCL",maxval(klcl),minval(klcl),maxval(kpbl),minval(kpbl)
      endif 
      !
      !
      !
      !--- determine the vertical entrainment/detrainment rates, the level of convective cloud base -kbcon-
      !--- and the scale dependence factor (sig).
      !
      call set_entr_detr_rates(cumulus,its,ite,itf,kts,kte,ktf,ierr,klcl,min_entr_rate,entr_rate &
                              ,entr_rate_2d,cd,mentrd_rate,cdd,qo_cup, qeso_cup,cnvcf)
      !
      !
      !--- determine the moist static energy of air parcels at source level
      !
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)

         x_add = (c_alvl*zqexec(i)+c_cp*ztexec(i)) +  x_add_buoy(i)
         call get_cloud_bc(cumulus,ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i),he_cup (kts:kte,i),hkb (i)&
                          ,k22(i),x_add,Tpert(kts:kte,i))
         call get_cloud_bc(cumulus,ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i),heo_cup(kts:kte,i),hkbo(i)&
                          ,k22(i),x_add,Tpert(kts:kte,i))
      enddo
      !
      !--- determine the level of convective cloud base  - kbcon
      !
      call cup_cloud_limits(cumulus,ierrc,ierr,cap_max_increment,cap_max,heo_cup,heso_cup,qo_cup   &
                           ,qeso_cup,po,po_cup,zo_cup,heo,hkbo,qo,qeso,entr_rate_2d,hcot,k22,kbmax &
                           ,klcl,kbcon,ktop,depth_neg_buoy,frh_bcon,Tpert,start_level              &
                           ,use_excess,zqexec,ztexec,x_add_buoy,xland,itf,ktf,its,ite, kts,kte)
      !
      !--- scale dependence factor (sig)
      !
      call set_scale_dep_factor(cumulus,USE_SCALE_DEP,DOWNDRAFT,its,ite,itf,kts,kte,ktf   &
                               ,ierr,sig_factor,stochastic_sig,dx,sig,sigd)
      !
      !--- increase detrainment in stable layers
      !
      call cup_minimi(HEso_cup,Kbcon,kstabm,kstabi,ierr,itf,ktf,its,ite, kts,kte)
      !
      !--- option for using the inversion layers as a barrier for the convection development
      !
      call get_cloud_top_by_inv_layers(cumulus,USE_INV_LAYERS,its,ite,itf,kts,kte,ktf &
                                      ,ierr,ierrc,psur,po_cup,tn_cup,zo_cup,mid,shal,kbcon,ktop)
      !
      !--- determine the normalized mass flux profile for updraft
      !
      do vtp_index = get_num_elements(vec_ok),1,-1
          i = get_data_value(vec_ok,vtp_index)
         call get_zu_zd_pdf(vtp_index,trim(cumulus),trim(cumulus)//"_up",ierr(i),k22(i),ktop(i),zuo(kts:kte,i),kts,kte,ktf  &
                           ,kpbl(i),k22(i),kbcon(i),klcl(i),po_cup(kts:kte,i),psur(i),xland(i),random(i),x_add_buoy(i))
      enddo
      !
      !-- calculate mass entrainment and detrainment
      !
      call get_lateral_massflux(itf,ktf, its,ite, kts,kte,min_entr_rate,ierr,ktop,zo_cup,zuo,cd  &
                               ,entr_rate_2d,po_cup,up_massentro, up_massdetro,cumulus,kbcon,k22,kpbl &
                               ,up_massentru,up_massdetru,lambau_dp)

      !
      !-- 1st guess for moist static energy and dbyo (not including ice phase)
      !
      call get_1st_guess_MSE_profile(cumulus,ave_layer,its,ite,itf,kts,kte,ktf,ierr,start_level,k22    &
                                    ,ktop,hkbo,zqexec,ztexec,x_add_buoy,xland,zuo,zuo,po,heo_cup,tpert &
                                    ,up_massentro,up_massdetro,heso_cup,heo,hco,cnvcf,heso)
      !-- Get buoyancy of updrafts
      !
      call get_buoyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop,hco,heo_cup,heso_cup,dbyo,zo_cup)

      if(first_guess_w .or. autoconv == 4 .or. autoconv == 3) then
         call cup_up_moisture_light(cumulus,start_level,klcl,ierr,ierrc,zo_cup,qco,qrco,pwo,pwavo,hco,tempco,xland   &
                                   ,po,p_cup,kbcon,ktop,cd,dbyo,clw_all,t_cup,qo,gammao_cup,zuo          &
                                   ,qeso_cup,k22,qo_cup,zqexec,use_excess,rho,up_massentro,up_massdetro    &
                                   ,psum,psumh,x_add_buoy,1,itf,ktf,its,ite, kts,kte         )

         call cup_up_vvel(vvel2d,vvel1d,zws,entr_rate_2d,cd,zo,zo_cup,zuo,dbyo,GAMMAo_CUP,tn_cup &
                         ,tempco,qco,qrco,qo,start_level,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte&
                         ,wlpool,wlpool_bcon)
      endif

      !
      !--- calculate moisture properties of updraft
      !
      call cup_up_moisture(cumulus,start_level,klcl,ierr,ierrc,zo_cup,qco,qrco,pwo,pwavo,hco,tempco,xland   &
                          ,po,p_cup,kbcon,ktop,cd,dbyo,clw_all,t_cup,qo,gammao_cup,zuo,qeso_cup   &
                          ,k22,qo_cup,zqexec,use_excess,ccn,rho,up_massentro,up_massdetro,psum    &
                          ,psumh,x_add_buoy,vvel2d,vvel1d,zws,entr_rate_2d                        &
                          ,1,itf,ktf,its,ite, kts,kte,z1,dtime,cnvcf,qes                              )

      !
      !-- get melting profile
      !
      call get_melting_profile(ierr,tn_cup,po_cup, p_liq_ice,melting_layer,qrco    &
                              ,pwo,edto,pwdo,melting,itf,ktf,its,ite, kts,kte, cumulus)

      !
      !-- updraft moist static energy + momentum budget
      !
       call get_updraft_profile(cumulus,ave_layer,use_linear_subcl_mf,its,ite,itf,kts,kte,ktf,ierr,start_level&
                               ,k22,ktop,hkb,hkbo,zqexec,ztexec,x_add_buoy,xland,pgcon,po,u_cup,v_cup,he,heo  &
                               ,he_cup,heo_cup,heso_cup ,hes_cup,p_liq_ice,qrco,Tpert ,us,vs,zuo,up_massdetru &
                               ,up_massentru,up_massdetro,up_massentro,hc,uc,vc,hco,cnvcf,hes,heso)
      
      !
      !--- Get buoyancy of updrafts
      !
      call get_buoyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop,hc, he_cup, hes_cup,  dby, z_cup)
      call get_buoyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop,hco,heo_cup,heso_cup,dbyo,zo_cup)
      !
      !
      !--- calculate workfunctions for updrafts
      !
      call cup_up_aa0(aa0,z_cup ,zuo,dby  ,GAMMA_CUP  ,t_cup ,k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte)
      call cup_up_aa0(aa1,zo_cup,zuo,dbyo ,GAMMAo_CUP ,tn_cup,k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte)

      ! do vtp_index = get_num_elements(vec_ok),1,-1
      !    i = get_data_value(vec_ok,vtp_index)
      !    if(aa1(i).eq.0.)then
      !       ierr(i)=17 ; is_removed = remove(vec_ok, i)
      !       ierrc(i)="cloud work function zero"
      !    endif
      ! enddo
      !---
      !
      if(.not. FIRST_GUESS_W .or. use_pass_cloudvol == 3) then
         !--- calculate in-cloud/updraft air temperature for vertical velocity
         !
         tempco (:,:) = tn_cup(:,:)
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            do k=kts,ktf
               tempco (k,i) = (1./c_cp)*(hco (k,i)-c_grav*zo_cup(k,i)-c_alvl*qco (k,i))
            enddo
         enddo
      endif
      !
      !---gust front impact on removal instabil time-scale and trigger function based on KE > CIN
      !
      if(convection_tracer == 1) then 
            wlpool_bcon  (:) = wlpool(:)
            ke_gustfront (:) = 0.5 * max(wlpool_bcon(:)**2,  zws(:)**2) + 1.e-6
            if(cumulus == 'deep')then
               call cup_up_aa0(cin1,zo_cup,zuo,dbyo ,GAMMAo_CUP   ,tn_cup  ,k22,klcl,kbcon,ktop &
                              ,ierr,itf,ktf,its,ite, kts,kte,'CIN')
            
               if( add_coldpool_trig == 2)then
                  do vtp_index = get_num_elements(vec_ok),1,-1
                      i = get_data_value(vec_ok,vtp_index)
                      if (ke_gustfront(i) < abs( min(cin1(i), 0.))) then
                        ierr(i)=500; is_removed = remove(vec_ok, i)
                        ierrc(i)="ke_gustfront less than -cin1"
                      endif 
                  enddo
               endif
            endif 
      endif
      !
      !--- vertical velocity
      !
      if(.not. FIRST_GUESS_W .or. use_pass_cloudvol == 3) then
           call cup_up_vvel(vvel2d,vvel1d,zws,entr_rate_2d,cd,zo,zo_cup,zuo,dbyo,GAMMAo_CUP,tn_cup      &
                           ,tempco,qco,qrco,qo,start_level,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte&
                           ,wlpool,wlpool_bcon)
      endif

      !---- new rain (for future use, keep this)
      !
      !--- calculate rain mixing ratio in updrafts
      !
      !       call cup_up_rain(cumulus,klcl,kbcon,ktop,k22,ierr,xland         &
      !                       ,zo_cup,qco,qrco,pwo,pwavo,po,p_cup,t_cup,tempco&
      !                       ,zuo,up_massentro,up_massdetro,vvel2d,rho         &
      !                       ,qrs                                            &
      !                       ,itf,ktf,its,ite, kts,kte)

      !---- new rain
      !
      !-----------------------------------------------------------------------------------------
      !--- downdraft section
      !-----------------------------------------------------------------------------------------
      !
      !--- downdraft max origin level (kzdown)
      !
      call get_max_dd_height(cumulus,zcutdown,itf,ktf,its,ite,kts,kte,ktop,kzdown,ierr,z1,zo_cup)
      !                        
      !--- downdraft originating level - jmin
      !
      call cup_minimi(heso_cup,k22,kzdown,jmin,ierr,itf,ktf,its,ite, kts,kte)

      call get_jmin(cumulus,itf,ktf,its,ite, kts,kte,ierr,kdet,ktop,kbcon,jmin,ierrc  &
                   ,beta,depth_min,heso_cup,zo_cup,melting_layer)
      !
      !-- get downdraft normalized mass flux
      !
      do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            call get_zu_zd_pdf(vtp_index,trim(cumulus),"DOWN",ierr(i),kdet(i),jmin(i),zdo(:,i),kts,kte,ktf &
                              ,kpbl(i),k22(i),kbcon(i),klcl(i),po_cup(kts:kte,i),psur(i),xland(i)&
                              ,random(i),x_add_buoy(i))
      enddo
      !
      !--- get lateral mass fluxes associated with downdrafts
      !
      call get_lateral_massflux_down(trim(cumulus),itf,ktf, its,ite, kts,kte,ierr,jmin       &
                                    ,zo_cup,zdo,cdd,mentrd_rate_2d,dd_massentro,dd_massdetro &
                                    ,cumulus,mentrd_rate,dd_massentru,dd_massdetru,lambau_dn)
      !
      !---  get wet bulb temperature and moisture at jmin
      !
      if(USE_WETBULB == 1 .and. cumulus /= 'shallow' ) then
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            k = jmin(i)
            call get_wetbulb(jmin(i),qo_cup(k,i),t_cup(k,i),po_cup(k,i),q_wetbulb(i),t_wetbulb(i))
         enddo
      endif
      !
      !
      !--- downdraft moist static energy + momentum budget
      !
      call get_downdraft_profile(cumulus,ave_layer,pgcon,use_wetbulb,itf,ktf, its,ite, kts,kte        &
                                ,ierrc,ierr,jmin,t_wetbulb,q_wetbulb,hc,heo,heso_cup,u_cup,v_cup,us,vs&
                                ,zdo,dd_massdetro,dd_massentro,zo_cup,dd_massdetru,dd_massentru       &
                                ,bud,hcdo,ucd,vcd,dbydo)
      !
      !--- calculate moisture properties of downdraft
      !
      call cup_dd_moisture(cumulus,ierrc,zdo,hcdo,heso_cup,qcdo,qeso_cup,pwdo,qo_cup,zo_cup &
                          ,dd_massentro,dd_massdetro,jmin,ierr,gammao_cup,pwevo,bu,qrcdo,qo &
                          ,heo,tn_cup,t_wetbulb,q_wetbulb,qco,pwavo,tempcdo,itf,ktf,its,ite, kts,kte)
      !
      !--- determine downdraft strength in terms of windshear
      !
      call cup_dd_edt(cumulus,ierr,us,vs,zo,ktop,kbcon,edt,po,pwavo,pwo,ccn,pwevo,edtmax,edtmin,maxens2&
                     ,edtc,psum,psumh,rho,aeroevap,itf,ktf,its,ite, kts,kte,vshear)

      !
      !-- determine epsilon =  dd_mass_flux / up_mass_flux 
      !
      if( cumulus /= 'shallow' ) then
         do iedt=1,maxens2
            do vtp_index = get_num_elements(vec_ok),1,-1
             i = get_data_value(vec_ok,vtp_index)
             edto(i)=sigd(i)*edtc(i,iedt)
           end do
         end do
      !
      endif  
      !-----------------------------------------------------------------------------------------



      !
      !--- Implements Becker et al (2021) closure, part 1
      !
      if (DICYCLE == 2 .and. trim(cumulus) == 'deep') then
         call BeckerEtAlClosure(cumulus,ave_layer,its, itf, ite, kts, ktf, kte, k22, start_level, ktop, klcl, kbcon &
                              , t, tn, tn_adv, q, qo, qo_adv, psur, tsur,us, vs, zqexec, ztexec                     &
                              , x_add_buoy, xland, tpert, zuo, up_massdetro , up_massentro, zuo, p_liq_ice          &
                              , qrco, ierr, zo, po, z1,aa1_radpbl, aa1_adv)
      end if
  
      !
      !--- Bechtold et al 2008 time-scale of cape removal
      !
      if(cumulus=='deep') then
               tau_ecmwf(:)=tau_deep;   wmean(:) = 3. !  mean vertical velocity m/s
      else
               tau_ecmwf(:)=tau_mid ;   wmean(:) = 3. 
      endif
 
      if(sgs_w_timescale == 1 .and. cumulus=='deep') then
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            !- mean vertical velocity based on integration of vertical veloc equation
            wmean(i) = min(max(vvel1d(i),3.),20.)

            !- time-scale cape removal from Bechtold et al. 2008
            tau_ecmwf(i)=( zo_cup(ktop(i),i)- zo_cup(kbcon(i),i) ) / wmean(i)
            tau_ecmwf(i)= min(10800., max(600.,tau_ecmwf(i)))
         enddo
      endif

      ! 
      !--- diurnal cycle section
      !
      !--- Implements the Bechtold et al (2014)
      !
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         umean= 2.0+sqrt(0.5*(us(1,i)**2+vs(1,i)**2+us(kbcon(i),i)**2+vs(kbcon(i),i)**2))
         !--                    - over land -            -          over ocean       s   -
         tau_bl(i)=(1.-xland(i))*tau_ecmwf(i) + xland(i)*(zo_cup(kbcon(i),i)- z1(i)) /umean
      enddo
      !
      !-- calculate "pcape" or equivalent cloud work function from the BL forcing only
      !
      if( dicycle <= 2 .and. cumulus == 'deep') then
         iversion=0
         call cup_up_aa1bl(iversion,aa1_bl,aa1_fa,aa1,t,tn,q,qo,dtime,po_cup,zo_cup,zuo,dbyo        &
                          ,gammao_cup,tn_cup,rho,klcl,kpbl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte & 
                          ,xland,ztexec,xlons,xlats, h_sfc_flux,le_sfc_flux,tau_bl,tau_ecmwf        &
                          ,t_star,cumulus,tn_bl,qo_bl  )

         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            aa1_bl(i) = (aa1_bl(i)/T_star) * tau_bl(i) ! units J/kg
           !aa1_bl(i) = (aa1_bl(i)/T_star) * tau_bl(i) - cin1(i)
            aa1_bl(i) = min(2000., abs(aa1_bl(i)))*sign(1.,aa1_bl(i))
         enddo
         !
         !--- Adds Becker et al (2021) closure, part 2
         !
         if(dicycle==2) then
            call get_Qadv(cumulus,itf,ktf,its,ite,kts,kte,ierr,dtime,q,qo,qo_adv,po,po_cup &
                         ,qeso, q_adv,col_sat_adv,alpha_adv,tau_bl,zo_cup,kbcon,ktop)
         endif
      endif
      !
      !
      !--- get the environmental normalized mass flux
      !
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         zenv(:,i) = zuo(:,i)-edto(i)*zdo(:,i)
      enddo
      !
      !--- check mass conservation
      !
      call check_mass_conserv(itf,ktf, its,ite, kts,kte,ierr,ktop,edto,zo,zuo,zdo,dd_massdetro  &
                             ,dd_massentro,up_massentro,up_massdetro)
      !
      !--- change per unit mass that a model cloud would modify the environment
      !
      call get_dellas(vert_discr,use_fct,alp1,dtime,itf,ktf, its,ite, kts              &
                     ,kte,ierr,ktop,edto,po_cup,hco,heo,heo_cup,hcdo,qo,qco,qrco       &
                     ,qrcdo,qo_cup,qcdo,pwo,pwdo,p_liq_ice,melting,dbydo,us,vs         &
                     ,uc,vc,u_cup,v_cup,ucd,vcd,zuo,zdo,zenv,dd_massdetro,dd_massentro &
                     ,up_massentro,up_massdetro,dellu, dellv, dellah, dellat, dellaq   &
                     ,dellaqc, dellabuoy, subten_H, subten_Q, subten_T)

      !
      !--- apply environmental subsidence on grid-scale ice and liq water contents, and cloud fraction (Upwind scheme)
      !
      if(APPLY_SUB_MP == 1) &
         call apply_sub_microphys(cumulus,alp1,itf,ktf,its,ite,kts,kte,ktop,nmp,ierr,po_cup,zenv &
                                 ,mpql,mpqi,mpcf,dellampqi,dellampql,dellampcf)
      !
      !
      !--- make the smoothness procedure
      !
      if(USE_SMOOTH_TEND > 0) &
         call smooth_tend(use_smooth_tend,itf,ktf, its,ite, kts,kte,ierr,ktop &
                         ,po_cup,dellu, dellv, dellah, dellat, dellaq,dellaqc)
      !
      !
      !--- using dellas, calculate changed environmental profiles
      !
      call get_env_change(coupl_mphysics,itf,ktf,its,ite,kts,kte,ierr,heo,qo &
                         ,tn,p_liq_ice,mbdt,dellu,dellv,dellah,dellat,dellaq &
                         ,dellaqc,subten_h,subten_q,subten_t,xhe,xq,xt)
      !
      !--- calculate moist static energy, heights, qes
      !
      call cup_env(xz,xqes,xhe,xhes,xt,xq,po,z1,psur,ierr,-1,itf,ktf,its,ite,kts,kte)
      !
      !--- environmental values on cloud levels
      !
      call cup_env_clev(xt,xqes,xq,xhe,xhes,xz,po,xqes_cup,xq_cup, xhe_cup,us,vs,u_cup,v_cup&
                       ,xhes_cup,xz_cup,po_cup,gamma_cup,xt_cup,psur,tsur,ierr,z1,itf,ktf,its,ite, kts,kte)
      !
      !--- static control
      !--- moist static energy inside cloud - part 1
      do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            xhc(:,i) = xhes_cup(:,i)
            !-note that hkbo already contains the contribuition from ztexec, zqexec and x_add_buoy
            call get_cloud_bc(cumulus,ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i),dellah (kts:kte,i)&
                             ,dellah_aver,k22(i))
            xhkb(i) = dellah_aver*mbdt(i) + hkbo(i)
            xhc (kts:start_level(i),i) = xhkb(i)           
      enddo
      !
      !--- option to produce linear fluxes in the sub-cloud layer.
      if(cumulus == 'shallow' .and. use_linear_subcl_mf == 1) then
            do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               call get_delmix(cumulus,kts,kte,ktf,xland(i),start_level(i),po(kts:kte,i) &
                              ,xhe_cup (kts:kte,i), xhc (kts:kte,i))
           enddo
      endif
      !--- moist static energy inside cloud - part 2
      do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
         
            do k=start_level(i) + 1,ktop(i) + 1 
              denom= (zuo(k-1,i)-.5*up_massdetro(k-1,i)+up_massentro(k-1,i)) +1.e-12
              
              if(use_pass_cloudvol == 2 .or. use_pass_cloudvol == 12 .or. use_pass_cloudvol == 3) then 
                  xh_env_eff = (1.-cnvcf(k-1,i))*xhe(k-1,i) + cnvcf(k-1,i)*xhes(k-1,i)
                  xhc(k,i)=(xhc(k-1,i)*zuo(k-1,i)-.5*up_massdetro(k-1,i)*xhc(k-1,i)+ &
                                                     up_massentro(k-1,i)*xh_env_eff) / denom

              else
                  xhc(k,i)=(xhc(k-1,i)*zuo(k-1,i)-.5*up_massdetro(k-1,i)*xhc(k-1,i)+ &
                                                     up_massentro(k-1,i)*xhe(k-1,i)) / denom
              endif

              if(k==start_level(i)+1)  then
                     x_add = (c_alvl*zqexec(i)+c_cp*ztexec(i)) +  x_add_buoy(i)
                     xhc(k,i)= xhc(k,i) + x_add*up_massentro(k-1,i)/denom
              endif
              !                          ------ ice content --------
              xhc (k,i)= xhc (k,i)+ c_xlf*(1.-p_liq_ice(k,i))*qrco(k,i)
            enddo
      enddo
      !
      !--- buoyancy and workfunctions for updraft
      call get_buoyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop,xhc,xhe_cup,xhes_cup,xdby,xz_cup)
      call cup_up_aa0  (xaa0,xz_cup,zuo,xdby,GAMMA_CUP,xt_cup, k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte)

      !
      !--- LARGE SCALE FORCING

      !--- precip ensemble
      !
      call get_pr_ens(cumulus,c0_mid,itf,ktf,its,ite,kts,kte,maxens3,ktop,ierr,ierrc,pwo,pwdo,edto,pr_ens)

      !--- calculate cloud base mass flux
      !

      if(cumulus == 'deep') &
            call cup_forcing_ens_deep(itf,ktf,its,ite, kts,kte,ens4,ensdim,ichoice,maxens,maxens2,maxens3 &
                                     ,ierr,k22,kbcon,ktop,xland1,aa0,aa1,xaa0,mbdt,dtime,xf_ens,mconv,qo  &
                                     ,po_cup,omeg,zdo,zuo,pr_ens,edto,tau_ecmwf,aa1_bl,xf_dicycle, xk_x   &
                                     ,alpha_adv,Q_adv,aa1_radpbl,aa1_adv,wlpool_bcon,xf_coldpool,cin1     )

      if(cumulus == 'mid') &
            call cup_forcing_ens_mid(aa0,aa1,xaa0,mbdt,dtime,ierr ,po_cup,ktop,k22,kbcon,kpbl,ichoice    &
                                    ,maxens,maxens3,itf,ktf,its,ite, kts,kte,tau_ecmwf,aa1_bl,xf_dicycle &
                                    ,dhdt,xff_mid,zws,hc,hco,he_cup,heo_cup,wlpool_bcon,xf_coldpool)

      if(cumulus == 'shallow') then
            call cup_up_cape_cin(cape,zo_cup,zuo,dby,gamma_cup,t_cup,klcl,k22,kbcon,ktop,ierr    &
                                ,tempco,qco,qrco,qo_cup,itf,ktf,its,ite, kts,kte, 'CAPE')

            call cup_forcing_ens_shal(itf,ktf,its,ite,kts,kte,dtime,ichoice,ierrc,ierr,klcl,kpbl,kbcon,k22,ktop  &
                                     ,xmb,tsur,cape,h_sfc_flux,le_sfc_flux,zws,po, hco, heo_cup,po_cup,t_cup     &
                                     ,dhdt,rho,xff_shal,xf_dicycle,tke_pbl,wlpool_bcon,xf_coldpool,ke_gustfront)
      endif
      !
      !
      !--- include kinetic energy dissipation converted to heating
      !
      call ke_to_heating(itf,ktf,its,ite, kts,kte,ktop,ierr,po_cup,us,vs,dellu,dellv,dellat)
      !
      !--- get the feedback to the 'large-scale' fields
      !
      call cup_output_ens(cumulus,xff_shal,xff_mid,xf_ens,ierr,edto,dellat,dellaq     &
                         ,dellaqc,outt, outq,outqc,zuo,pre,pwo,pwdo,pwo_eff,xmb       &
                         ,ktop,maxens2,maxens,pr_ens,maxens3,ensdim,sig,xland1        &
                         ,ichoice,itf,ktf,its,ite, kts,kte                    &
                         ,xf_dicycle,outu,outv,dellu,dellv,dtime,po_cup,kbcon         &
                         ,dellabuoy,outbuoy,dellampqi,outmpqi,dellampql,outmpql       &
                         ,dellampcf,outmpcf,nmp,rh_dicycle_fct,xf_coldpool,wlpool_bcon)
      !
      !--- get the net precipitation flux (after downdraft evaporation)
      !
      call get_precip_fluxes(cumulus,klcl,kbcon,ktop,k22,ierr,xland,pre,xmb,pwo,pwavo,edto &
                            ,pwevo,pwdo,t_cup,tempco,prec_flx,evap_flx,itf,ktf,its,ite, kts,kte)

      !
      !--- includes rainfall evap below cloud base
      !
      if(use_rebcb == 1)                                                                 &
         call rain_evap_below_cloudbase(cumulus,itf,ktf, its,ite, kts,kte,ierr,kbcon,ktop&
                                       ,xmb,psur,xland,qo_cup,t_cup,po_cup,qes_cup,pwavo &
                                       ,edto,pwevo,pwo,pwdo,pre,prec_flx,evap_flx        &
                                       ,outt,outq,outbuoy,evap_bcb)
      !
      !--- includes effects of the remained cloud dissipation into the enviroment
      !
      if(use_cloud_dissipation >= 0.)                                                &
         call cloud_dissipation(cumulus,itf,ktf, its,ite, kts,kte,ierr,kbcon,ktop    &
                               ,dtime,xmb,xland,qo_cup,qeso_cup,po_cup,outt,outq     &
                               ,outqc,zuo,vvel2d,qrco,sig,tempco,qco,tn_cup,heso_cup,zo,zo_cup)
      !
      !--- includes rain evaporation as in SAS
      !
      if(irainevap == 1)                                                         &
          call sas_rainevap(cumulus,itf,ktf, its,ite, kts,kte,ierr,kbcon,ktop    &
                           ,dtime,sig,xmb,pre,xland,edto,pwo,pwdo,po_cup,qo,tn   &
                           ,outt, outq,outbuoy,rntot,delqev,delq2,qevap,rn,qcond,rainevap,qeso)
      !
      !
      !--- get lightning flashes density (parameterization from Lopez 2016, MWR)
      !
      if( lightning_diag == 1 .and. trim(cumulus) == 'deep') then
         call cup_up_cape_cin(cape,zo_cup,zuo,dby,gamma_cup,t_cup,klcl,k22,kbcon,ktop,ierr &
                             ,tempco,qco,qrco,qo_cup,itf,ktf,its,ite, kts,kte, 'CAPE')

         call cup_up_lightning(itf,ktf,its,ite, kts,kte, ierr, kbcon,ktop,xland,cape &
                              ,zo,zo_cup,t_cup,t,tempco,qrco,po_cup,rho,prec_flx,lightn_dens)
      endif
      !
      !--- outputs the cloud fraction source for the passive cloud volume parameterization
      if(use_pass_cloudvol == 3) then
         do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)

               !--- estimation of the fractional area of updraft that will remains 
               !--- in the ambient
               do k=kbcon(i),ktop(i)
                  !--- units 1/sec
                  clfrac(k,i) = ((xmb(i)/sig(i))*zuo(k,i)/(rho(k,i)*max(2.,vvel2d(k,i))))/dtime
               enddo
         enddo
         !print*,'clfrac',trim(cumulus),maxval(clfrac),minval(clfrac)
      endif

      !
      !--- convert mass fluxes, etc...
      !
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         pwavo       (i)   = xmb(i)*pwavo       (i)
         pwevo       (i)   = xmb(i)*pwevo       (i)
         zuo         (:,i) = xmb(i)*zuo         (:,i)
         zdo         (:,i) = xmb(i)*zdo         (:,i)
         pwo         (:,i) = xmb(i)*pwo         (:,i)
         pwdo        (:,i) = xmb(i)*pwdo        (:,i)
         up_massentro(:,i) = xmb(i)*up_massentro(:,i)
         up_massdetro(:,i) = xmb(i)*up_massdetro(:,i)
         dd_massentro(:,i) = xmb(i)*dd_massentro(:,i)
         dd_massdetro(:,i) = xmb(i)*dd_massdetro(:,i)
         zenv        (:,i) = xmb(i)*zenv        (:,i)
         subten_Q    (:,i) = xmb(i)*subten_Q    (:,i)
         subten_H    (:,i) = xmb(i)*subten_H    (:,i)
         subten_T    (:,i) = xmb(i)*subten_T    (:,i)
      end do

      !
      !--- outputs a model sounding for the stand-alone code (part 2)
      !
      if(OUTPUT_SOUND == 1) then
         call SOUND(2,cumulus, int_time, dtime, ens4, itf, its, ite, kts, kte, xlats, xlons, jcol, whoami_all &
                    , z, qes, he, hes, t, q, po, z1, psur, zo, qeso, heo, heso, tn, qo, us, vs, omeg, xz &
                    , h_sfc_flux, le_sfc_flux, tsur, dx, stochastic_sig, zws, ztexec, zqexec, xland &
                    , kpbl, k22, klcl, kbcon, ktop, aa0, aa1, sig, xaa0, hkb, xmb, pre, edto &
                    , zo_cup, dhdt, rho, zuo, zdo, up_massentro, up_massdetro, outt, outq, outqc, outu, outv)
      endif
      
      !
      !-- get tendencies for ice/liq drops number concentration for cloud microphysics
      !
      if(liq_ice_number_conc == 1) then
         call get_liq_ice_number_conc(itf,ktf,its,ite, kts,kte,ierr,ktop &
                                     ,dtime,rho,outqc,tempco,outnliq,outnice)
      endif
      !
      !
      !--------------------------------------------------------------------------------------------!
      !-- section for atmospheric composition convective tracer and wet removal
      if(use_tracer_transp==1)  then  
        call tracerConvTransGF(its, itf, jl, jmin, ite, kts, kte, ktf, mtp, chem_name_mask, ave_layer, ierr, k22, ktop &
                              ,   start_level, dd_massdetro, dd_massentro, dtime, edto, fscav, po, po_cup, pw_up_chem &
                              ,   pwavo, pwevo, pwdo, pwo, qrco, sc_up_chem, tempco, tot_pw_up_chem &
                              ,   up_massdetro, up_massentro, vvel2d, xland, zdo, zo_cup, zuo, cumulus &
                              ,   se_chem, massf, out_chem, pw_dn_chem, sc_dn_chem, se_cup_chem, tot_pw_dn_chem, zenv)
      
      endif !- end of section for atmospheric composition
      !---------------------------------------------------------------------------------------------!
      !
      !-- GATE soundings
      if(p_use_gate .or. wrtgrads) then
       call gateSoundings(its, itf, jl, kte, kts, ktf, ierr, jmin, k22, kbcon, klcl, ktop, aa1, cd, clfrac, dby, dd_massentro &
                     , dd_massdetro, dellah, dellaq, dellaqc, dtime, edto, entr_rate_2d, evap_bcb, hc, hco, he_cup, HKB &
                     , massf, massi, melting_layer, mpcf, mpqi, mpql, out_chem, outmpcf, outmpqi, outmpql, outnice, outnliq &
                     , outq, outqc, outt, outu, outv, p_liq_ice, po, po_cup, pre, prec_flx, pwdo, pwo, q_cup, q, qcdo, qco &
                     , qeso_cup, qo_cup, qo , qrco, qrs, sc_dn_chem, sc_up_chem, se_chem, se_cup_chem, subten_h, subten_q &
                     , subten_t , t_cup, t, tn, tot_pw_dn_chem, tot_pw_up_chem, up_massentro, up_massdetro, us, vvel1d &
                     , vvel2d, xmb, z_cup, zdo, zenv, zo, zqexec, zuo, zws,  heso_cup, heo_cup,cumulus)
      endif
      ! 
      !-- miscelaneous procedures for output and/or debug
      !
      !-- get the total (deep+congestus) evaporation flux for output (units kg/kg/s)
      !
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         do k=kts,ktop(i)
            dp=100.*(po_cup(k,i)-po_cup(k+1,i))
             !--- add congestus and deep plumes, and convert to kg/kg/s
            revsu_gf(k,i) = revsu_gf(k,i) + evap_flx(k,i)*c_grav/dp
         enddo
      enddo
      !
      !--- for tracer convective transport / outputs
      !
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         do k=kts,ktf
            tup   (k,i) = tempco(k,i) !in-updraft temp
            cupclw(k,i) = qrco  (k,i) !in-updraft condensed water
         end do
         tup (kte,i) = t_cup(kte,i)
      end do
      !--- for outputs (only deep plume)
      !
      if( trim(cumulus) == 'deep') then
      do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            var2d(i) = p_cwv_ave(i)
            do k=kts,ktop(i)+1
               prfil_gf  (k,i) = prec_flx(k,i)
               var3d_agf (k,i) = vvel2d  (k,i)
            enddo
         enddo
      endif
      !--- for output/debug/diag only 
      if( trim(cumulus) == 'deep') then
            !AA1_       (:) = AA1       (:)
            !AA0_       (:) = AA0       (:)
            AA1_RADPBL_(:) = AA1_RADPBL(:)
               !aa0_    (:)  = aa0      (:)
               !aa1_    (:)  = aa1      (:)
               !aa1_bl_ (:)  = zlcl_sfc (:)! aa1_bl   (:)
               !tau_bl_ (:)  = tau_bl   (:)
               tau_ec_ (:)  = tau_ecmwf(:)
            do i=its,itf
               if(ierr(i) == 0) then
                 aa1_adv_ (i) = z_cup(ktop(i),i)*1.e-3
                 !aa0_     (i) = entr_rate(i)*1.e+3
                 !aa1_bl_  (i) = zlcl_sfc (i)!*1.e-3
               else
                 aa1_adv_ (i) = -0.9990000E+34
                 !aa0_     (i) = -0.9990000E+34
                 !aa1_bl_  (i) = -0.9990000E+34
               endif
            enddo  

            !if(dicycle == 2) then 
            !  aa1_adv_ (:) = q_adv     (:)
            !else
            !  aa1_adv_ (:) = vshear (:)
            !endif
          
            do i=its,itf
               if(ierr(i) == 0) cycle
               kbcon(i)=1
               ktop(i)=1
               klcl(i)=1
               jmin(i)=1
               k22(i)=1
            enddo
      endif
      !
      !
      !
      !
      call free_memory(vec_ok)
      call free_memory(vec_removed)

   end subroutine CUP_GF
   !------------------------------------------------------------------------------------
   subroutine cup_dd_edt(cumulus,ierr,us,vs,z,ktop,kbcon,edt,p,pwav,pw,ccn,pwev,edtmax  &
                        ,edtmin,maxens2,edtc,psum2,psumh,rho,aeroevap          &
                        ,itf,ktf,its,ite, kts,kte,vshear )

      implicit none

      character *(*)  ,intent (in) :: cumulus
      integer         ,intent (in) :: aeroevap,itf,ktf,its,ite, kts,kte
      integer         ,intent (in) :: maxens2
      !
      ! ierr error value, maybe modified in this routine
      !
      real,    dimension (:,:),intent (in   ) :: rho,us,vs,z,p,pw
      real,    dimension (:)  ,intent (in   ) :: pwav,pwev,ccn,psum2,psumh,edtmax,edtmin
      integer, dimension (:)  ,intent (in   ) :: ktop,kbcon
      integer, dimension (:)  ,intent (inout) :: ierr
      real,    dimension (:,:),intent (out  ) :: edtc
      real,    dimension (:)  ,intent (out  ) :: edt,vshear
      !
      !  local variables in this routine
      !
      integer i,k,kk,vtp_index
      real    einc,pef,pefb,prezk,zkbc
      real,    dimension (its:ite)         ::   vws,sdp
      real :: pefc,aeroadd,rhoc,dp,prop_c
      real, parameter ::  alpha3 = 1.9 ,beta3  = -1.13

      !
      !--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
      !
      ! */ calculate an average wind shear over the depth of the cloud
      !
      edt   =0.
      vws   =0.
      sdp   =0.
      vshear=0.
      edtc  =0.

      if(cumulus=='shallow') return

      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         do kk = kbcon(i),ktop(i)
            dp = p(kk,i) - p(kk+1,i)
            vws(i) = vws(i) + (abs((us(kk+1,i)-us(kk,i))/(z(kk+1,i)-z(kk,i)))  + &
                               abs((vs(kk+1,i)-vs(kk,i))/(z(kk+1,i)-z(kk,i)))) * dp
            sdp(i) = sdp(i) + dp
         enddo
         vshear(i) = 1.e3 * vws(i) / sdp(i)
      end do

      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         pef = (1.591-0.639*vshear(i)+0.0953*(vshear(i)**2) -0.00496*(vshear(i)**3))

         !print*,"shear=",vshear(i),pef,1-max(min(pef,0.9),0.1)
         pef = min(pef,0.9)
         pef = max(pef,0.1)
         edt(i) = 1.-pef

         !
         !--- cloud base precip efficiency
         !
         if(use_rebcb == 0) then 
           zkbc  = z(kbcon(i),i)*3.281e-3
           prezk = 0.02
           if(zkbc > 3.0) prezk=0.96729352+zkbc*(-0.70034167+zkbc* &
                        (0.162179896+zkbc*(- 1.2569798e-2+zkbc*(4.2772e-4-zkbc*5.44e-6))))
           if(zkbc > 25.) prezk=2.4
           pefb = 1./(1.+prezk)
           pefb = min(pefb,0.9)
           pefb = max(pefb,0.1)
           edt(i) = 1.-0.5*(pefb+pef)
         endif


         if(aeroevap.gt.1)then
            aeroadd=(p_ccnclean**beta3)*((psumh(i))**(alpha3-1)) !*1.e6
            !if(i.eq.ipr)write(0,*)'edt',p_ccnclean,psumh(i),aeroadd
            !prop_c=.9/aeroadd
            prop_c=.5*(pefb+pef)/aeroadd
            aeroadd=(ccn(i)**beta3)*((psum2(i))**(alpha3-1)) !*1.e6
            !if(i.eq.ipr)write(0,*)'edt',ccn(i),psum2(i),aeroadd,prop_c
            aeroadd=prop_c*aeroadd
            pefc=aeroadd
            if(pefc.gt.0.9)pefc=0.9
            if(pefc.lt.0.1)pefc=0.1
            EDT(I)=1.-pefc
            if(aeroevap.eq.2)EDT(I)=1.-.25*(pefb+pef+2.*pefc)
         endif

      enddo
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         edtc(i,1) = -edt(i)*pwav(i)/pwev(i)
         edtc(i,1) = min(edtmax(i), edtc(i,1))
         edtc(i,1) = max(edtmin(i), edtc(i,1))
      enddo

   end subroutine cup_dd_edt
   !------------------------------------------------------------------------------------
   subroutine cup_dd_moisture(cumulus,ierrc,zd,hcd,hes_cup,qcd,qes_cup                &
                             ,pwd,q_cup,z_cup,dd_massentr,dd_massdetr,jmin,ierr       &
                             ,gamma_cup,pwev,bu,qrcd, q,he,t_cup,t_wetbulb            &
                             ,q_wetbulb,qco,pwavo,tempcdo,itf,ktf,its,ite, kts,kte            )

      implicit none

      character(len=*), intent(in) :: cumulus
      integer         , intent(in) :: itf,ktf,its,ite, kts,kte
      ! q       = environmental q on model levels
      ! q_cup   = environmental q on model cloud levels
      ! qes_cup = saturation q on model cloud levels
      ! hes_cup = saturation h on model cloud levels
      ! hcd = h in model cloud
      ! bu = buoancy term
      ! zd = normalized downdraft mass flux
      ! gamma_cup = gamma on model cloud levels
      ! mentr_rate = entrainment rate
      ! qcd  = cloud q (including liquid water) after entrainment
      ! qrch = saturation q in cloud
      ! pwd  = evaporate at that level
      ! pwev = total normalized integrated evaoprate (I2)
      ! entr = entrainment rate
      ! cdd  = detrainment function
      !
      real,    dimension (:) ,intent (in  )          ::           &
         t_wetbulb,q_wetbulb, pwavo
      real,    dimension (:,:) ,intent (in   ) ::           &
         zd,t_cup,hes_cup,hcd,qes_cup,q_cup,z_cup,                      &
         dd_massentr,dd_massdetr,gamma_cup,q,he,qco
      integer, dimension (:) ,intent (in   )         ::           &
         jmin
      integer, dimension (:) ,intent (inout)         ::           &
         ierr
      real,    dimension (:,:) ,intent (out  ) ::           &
         qcd,qrcd,pwd,tempcdo
      real,    dimension (:) ,intent (out  )         ::           &
         pwev,bu
      character*(*), intent (inout), dimension(:)  :: ierrc
      !
      !  local variables in this routine
      !
      integer   ::     i,k,vtp_index
      real      ::     dh,dz,dq_eva,denom,fix_evap
      !
      bu  =0.  !-- buoyancy
      qcd =0.  !-- in-downdradt water vapor mixing ratio
      qrcd=0.  !-- saturation water vapor mixing ratio
      pwev=0.  !-- column integrated rain evaporation (normalized)
      pwd =0.  !-- rain evaporation at layer k
      tempcdo = t_cup  !-- in-cloud downdraft air temperature

      if(cumulus == 'shallow') return
      !
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)

         !-- boundary condition in jmin ('level of free sinking')
         k = jmin(i)
         dz= z_cup(k+1,i)-z_cup(k,i)

         qcd(k,i)=q_cup(k,i)

         if(use_wetbulb==1) then
            !--option 1
            !qcd(k,i)=q_wetbulb(i)
            !--option 2
            qcd(k,i)=0.5*(q_wetbulb(i)+qco(k,i)) ! mixture 50% env air + updraft
         endif

         dh=hcd(k,i)-hes_cup(k,i)

         if(dh < 0.)then
            qrcd(k,i)=(qes_cup(k,i)+(1./c_alvl)*(gamma_cup(k,i) /(1.+gamma_cup(k,i)))*dh)
         else
            qrcd(k,i)=qes_cup(k,i)
         endif

         pwd (k,i) = zd(k,i)*min(0.,qcd(k,i)-qrcd(k,i))
         qcd (k,i) = qrcd(k,i)
         pwev(i)   = pwev(i)+pwd(k,i)
         bu  (i)   = dz*dh

         do k=jmin(i)-1,kts,-1

            dz=z_cup(k+1,i)-z_cup(k,i)

            !-- downward transport + mixing
            denom   = (zd(k+1,i)-0.5*dd_massdetr(k,i)+dd_massentr(k,i)) + 1.e-12
            qcd(k,i)= (qcd(k+1,i)*zd(k+1,i) -0.5*dd_massdetr(k,i)*qcd(k+1,i)+ &
                                                 dd_massentr(k,i)*q  (k,i)    )/ denom
            !
            !--- to be negatively buoyant, hcd should be smaller than hes!
            !--- ideally, dh should be negative till dd hits ground, but that is not always
            !--- the case
            !
            dh    = hcd(k,i)-hes_cup(k,i)
            bu(i) = bu(i)+dz*dh
            qrcd(k,i)=qes_cup(k,i)+(1./c_alvl)*(gamma_cup(k,i) /(1.+gamma_cup(k,i)))*dh

            !-- rain water evaporation amount at layer k
            dq_eva=qcd(k,i)-qrcd(k,i)

            if(dq_eva > 0.)then
               dq_eva=0.
               qrcd(k,i)=qcd(k,i)
            endif
            !-- amount of the evaporated rain water
            pwd(k,i)=zd(k,i)*dq_eva  ! kg[water vapor]/kg[air]

            !-- source term for in-downdraft water vapor mixing ratio
            qcd(k,i)=qrcd(k,i)     ! => equiv to qcd = qcd - dq_eva !( -dq_eva >0 => source term for qcd)

            !-- total evaporated rain water
            pwev(i)=pwev(i)+pwd(k,i)

            !-- for GEOS diagnostic
            ! evap(k,i) = - edt * xmb * zd * dq_eva = - edt * xmb * pwd (k,i)
            ! downdfrat temp = (hcd(k,i)-qcd(k,i)*c_alvl-c_grav*z_cup(k,i))/c_cp - 273.15

            !--- calculate in-cloud downdraft air temperature 
            tempcdo(k,i) = (1./c_cp)*(hcd(k,i)-c_grav*z_cup(k,i)-c_alvl*qcd(k,i))
         
         end do

         if(pwev(i) >= 0.)then
            ierr(i)=70 ; is_removed = remove(vec_ok, i)
            ierrc(i)=" pwev >= 0 in cup_dd_moisture"
         endif
         if(bu(i) >= 0.)then
            ierr(i)=73 ; is_removed = remove(vec_ok, i)
            ierrc(i)=" bu >= 0   in cup_dd_moisture"
         endif

         !-- fix evap, in case of not conservation
         if(abs(pwev(i)) > pwavo(i) .and. ierr(i) == 0)then
               fix_evap = pwavo(i)/(1.e-16+abs(pwev(i)))
               pwev(i)  = 0.
               do k=jmin(i),kts,-1
                  pwd(k,i) = pwd (k,i)*fix_evap
                  pwev(i)  = pwev(i) + pwd(k,i)
                  dq_eva   = pwd (k,i)/(1.e-16+zd(k,i))
                  qcd(k,i) = qrcd(k,i) + dq_eva
               end do
               if(pwev(i) >= 0.)then
                  ierr(i)=75 ; is_removed = remove(vec_ok, i)
                  ierrc(i)=" pwev >= 0 in cup_dd_moisture"
               endif
         endif

      end do
   end subroutine cup_dd_moisture

   !------------------------------------------------------------------------------------

   subroutine cup_env(z,qes,he,hes,t,q,p,z1,psur,ierr,itest,itf,ktf,its,ite, kts,kte  )

      implicit none
      integer ,intent (in) :: itf,ktf,its,ite, kts,kte
      !
      ! ierr error value, maybe modified in this routine
      ! q           = environmental mixing ratio
      ! qes         = environmental saturation mixing ratio
      ! t           = environmental temp
      ! tv          = environmental virtual temp
      ! p           = environmental pressure
      ! z           = environmental heights
      ! he          = environmental moist static energy
      ! hes         = environmental saturation moist static energy
      ! psur        = surface pressure
      ! z1          = terrain elevation
      !
      real,    dimension (:,:),intent (in   )  :: p,t,q
      real,    dimension (:,:),intent (out  )  :: he,hes,qes
      real,    dimension (:,:),intent (in   )  :: z
      real,    dimension (:)  ,intent (in   )  :: psur,z1
      integer, dimension (:)  ,intent (in   )  :: ierr
      integer                 ,intent (in   )  :: itest
      !
      !  local variables in this routine
      !
      integer ::  i,k, vtp_index
      real    ::  pqsat

      !he  =0.0
      !hes =0.0
      !qes =0.0

      !--- better formulation for the mixed phase regime
      do k=kts,ktf
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            pqsat=satur_spec_hum(t(k,i),p(k,i))
            qes(k,i)=pqsat

            qes(k,i) = min(p_max_qsat, max(1.e-08,qes(k,i)))
            qes(k,i) = max(qes(k,i), q(k,i))
               
         end do
      end do
      !
      !--- calculate moist static energy - HE
      !    saturated moist static energy - HES
      !
      do k=kts,ktf
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            he (k,i)=c_grav*z(k,i)+c_cp*t(k,i)+c_alvl*q  (k,i)
            hes(k,i)=c_grav*z(k,i)+c_cp*t(k,i)+c_alvl*qes(k,i)
            he (k,i)=min(hes(k,i), he(k,i))
         end do
      end do

   end subroutine cup_env
   !------------------------------------------------------------------------------------
   subroutine cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,he_cup,us, vs,u_cup,v_cup &
                          ,hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur,tsur              &
                          ,ierr,z1,itf,ktf,its,ite, kts,kte                           )
      implicit none
      integer ,intent (in   )    :: itf,ktf, its,ite, kts,kte
      !
      ! ierr error value, maybe modified in this routine
      ! q           = environmental mixing ratio
      ! q_cup       = environmental mixing ratio on cloud levels
      ! qes         = environmental saturation mixing ratio
      ! qes_cup     = environmental saturation mixing ratio on cloud levels
      ! t           = environmental temp
      ! t_cup       = environmental temp on cloud levels
      ! p           = environmental pressure
      ! p_cup       = environmental pressure on cloud levels
      ! z           = environmental heights
      ! z_cup       = environmental heights on cloud levels
      ! he          = environmental moist static energy
      ! he_cup      = environmental moist static energy on cloud levels
      ! hes         = environmental saturation moist static energy
      ! hes_cup     = environmental saturation moist static energy on cloud levels
      ! gamma_cup   = gamma on cloud levels
      ! psur        = surface pressure
      ! z1          = terrain elevation
      !
      real,    dimension (:,:) ,intent (in )  :: qes,q,he,hes,z,p,t,us, vs
      real,    dimension (:,:) ,intent (out)  :: qes_cup,q_cup,he_cup,hes_cup,z_cup &
                                                ,p_cup,gamma_cup,t_cup,u_cup,v_cup
      real,    dimension (:)   ,intent (in )  :: psur,z1,tsur
      integer, dimension (:)   ,intent (in )  :: ierr
      !
      !  local variables in this routine
      !
      integer                              :: i,k, vtp_index 
      real                                 :: p1,p2,ct1,ct2,rho

      ! qes_cup  =0.
      ! q_cup    =0.
      ! hes_cup  =0.
      ! he_cup   =0.
      ! z_cup    =0.
      ! p_cup    =0.
      ! t_cup    =0.
      ! gamma_cup=0.
      ! u_cup    =0.
      ! v_cup    =0.

      if( clev_grid == 2 ) then
         !--original formulation
         do k=kts+1,ktf
            do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               qes_cup(k,i)=.5*(qes(k-1,i)+qes(k,i))
               q_cup  (k,i)=.5*(q  (k-1,i)+q  (k,i))
               hes_cup(k,i)=.5*(hes(k-1,i)+hes(k,i))
               he_cup (k,i)=.5*(he (k-1,i)+he (k,i))
               if(he_cup(k,i).gt.hes_cup(k,i))he_cup(k,i)=hes_cup(k,i)
               z_cup  (k,i)=.5*(z(k-1,i)+z(k,i))
               p_cup  (k,i)=.5*(p(k-1,i)+p(k,i))
               t_cup  (k,i)=.5*(t(k-1,i)+t(k,i))
               gamma_cup(k,i)=(c_alvl/c_cp)*(c_alvl/(c_rm*t_cup(k,i) &
                              *t_cup(k,i)))*qes_cup(k,i)
               u_cup  (k,i)=.5*(us(k-1,i)+us(k,i))
               v_cup  (k,i)=.5*(vs(k-1,i)+vs(k,i))

            end do
         end do
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            qes_cup(1,i)=qes(1,i)
            q_cup(1,i)=q(1,i)
            !hes_cup(1,i)=hes(1,i)
            !he_cup(1,i)=he(1,i)
            hes_cup(1,i)=c_grav*z1(i)+c_cp*t(1,i)+c_alvl*qes(1,i)
            he_cup (1,i)=c_grav*z1(i)+c_cp*t(1,i)+c_alvl*q  (1,i)
            !z_cup(1,i)=.5*(z(1,i)+z1(i))
            !p_cup(1,i)=.5*(p(1,i)+psur(i))
            z_cup(1,i)=z1(i)
            p_cup(1,i)=psur(i)
            t_cup(1,i)=t(1,i)
            gamma_cup(1,i)=c_alvl/c_cp*(c_alvl/(c_rm*t_cup(1,i) &
                           *t_cup(1,i)))*qes_cup(1,i)
            u_cup(1,i)=us(1,i)
            v_cup(1,i)=vs(1,i)
         end do

       !do k=kts,ktf
       ! i=1
       !        print*,"air_dens=",k,z_cup(k,i),p_cup(k,i),(p_cup(k,i)-p_cup(k+1,i))/(z_cup(k+1,i)-z_cup(k,i))/c_grav
       !enddo

      elseif( clev_grid == 0) then
         !--- weigthed mean
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)

            p_cup  (1,i)=psur(i)
            z_cup  (1,i)=z1(i)
            do k=kts,ktf-1
               p_cup (k+1,i) = 2.0*p(k,i) - p_cup(k,i)
               z_cup (k+1,i) = 2.0*z(k,i) - z_cup(k,i)
            end do

            ! ----------- p,T          k+1
            !p1
            ! ----------- p_cup,T_cup  k+1
            !p2
            ! ----------- p,T          k
            !
            ! ----------- p_cup,T_cup  k

            do k=kts,ktf-1
               p1=abs((p    (k+1,i) - p_cup(k+1,i))/(p(k+1,i)-p(k,i)))
               p2=abs((p_cup(k+1,i) - p    (k,i  ))/(p(k+1,i)-p(k,i)))

               t_cup  (k+1,i) = p1*t  (k,i) + p2*t  (k+1,i)

               u_cup  (k+1,i) = p1*us (k,i) + p2*us (k+1,i)
               v_cup  (k+1,i) = p1*vs (k,i) + p2*vs (k+1,i)
               q_cup  (k+1,i) = p1*q  (k,i) + p2*q  (k+1,i)
               he_cup (k+1,i) = p1*he (k,i) + p2*he (k+1,i)

               qes_cup(k+1,i) = p1*qes(k,i) + p2*qes(k+1,i)
               hes_cup(k+1,i) = p1*hes(k,i) + p2*hes(k+1,i)

               if(he_cup(k+1,i).gt.hes_cup(k+1,i))he_cup(k+1,i)=hes_cup(k+1,i)

               gamma_cup(k+1,i)=(c_alvl/c_cp)*(c_alvl/(c_rm*t_cup(k+1,i) &
                                *t_cup(k+1,i)))*qes_cup(k+1,i)

            end do
            !--- surface level from X(kts) and X_cup(kts+1) determine X_cup(kts)
            k=kts
            p1=abs(p    (k,i  )-p_cup(k,i))
            p2=abs(p_cup(k+1,i)-p_cup(k,i))

            ct1=(p1+p2)/p2
            ct2=p1/p2

            t_cup  (k,i) = ct1*t  (k,i) - ct2*t_cup(k+1,i)
            q_cup  (k,i) = ct1*q  (k,i) - ct2*q_cup(k+1,i)

            u_cup  (k,i) = ct1*us (k,i) - ct2*u_cup(k+1,i)
            v_cup  (k,i) = ct1*vs (k,i) - ct2*v_cup(k+1,i)
            qes_cup(k,i) = ct1*qes(k,i) - ct2*qes_cup(k+1,i)

            hes_cup(k,i)=c_grav*z_cup(k,i)+c_cp*t_cup(k,i)+c_alvl*qes_cup(k,i)
            he_cup (k,i)=c_grav*z_cup(k,i)+c_cp*t_cup(k,i)+c_alvl*q_cup  (k,i)

            if(he_cup(k,i).gt.hes_cup(k,i))he_cup(k,i)=hes_cup(k,i)

            gamma_cup(k,i)=c_alvl/c_cp*(c_alvl/(c_rm*t_cup(k,i)*t_cup(k,i)))*qes_cup(k,i)
         end do

      elseif(clev_grid == 1) then
         !--- based on Tiedke (1989)
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            do k=ktf, kts+1,-1

               qes_cup(k,i) = qes(k,i)
               q_cup  (k,i) = q  (k,i)
               p_cup  (k,i) = 0.5*(p(k-1,i)+p(k,i))
               z_cup  (k,i) = 0.5*(z(k-1,i)+z(k,i))
               t_cup  (k,i) = (max(c_cp*t(k-1,i)+c_grav*z(k-1,i),c_cp*t(k,i)+c_grav*z(k,i)) - c_grav*z_cup(k,i))/c_cp

               if(qes(k,i) < p_max_qsat) &
                  call get_interp(qes_cup(k,i),t_cup(k,i),p_cup(k,i),qes_cup(k,i),t_cup(k,i))

               q_cup  (k,i) = min(q(k,i),qes(k,i)) + qes_cup(k,i) - qes(k,i)
               q_cup  (k,i) = max(q_cup  (k,i) ,0.0)

            end do
            !---level kts
            qes_cup(1,i)= qes (1,i)
            q_cup  (1,i)= q   (1,i)
            z_cup  (1,i)= z1  (i)
            p_cup  (1,i)= psur(i)

            t_cup  (1,i)= (c_cp*t(1,i)+c_grav*z(1,i) - c_grav*z_cup(1,i))/c_cp

            hes_cup(1,i)=c_grav*z_cup(1,i)+c_cp*t_cup(1,i)+c_alvl*qes_cup(1,i)
            he_cup (1,i)=c_grav*z_cup(1,i)+c_cp*t_cup(1,i)+c_alvl*q_cup  (1,i)

            gamma_cup(1,i)=c_alvl/c_cp*(c_alvl/(c_rm*t_cup(1,i)*t_cup(1,i)))*qes_cup(1,i)
            u_cup(1,i)=us(1,i)
            v_cup(1,i)=vs(1,i)

            do k=ktf,kts+1,-1
               p1=max(c_cp*t_cup(k,i)+c_grav*z_cup(k,i), c_cp*t_cup(k-1,i)+c_grav*z_cup(k-1,i))
               t_cup(k,i) = (p1-c_grav*z_cup(k,i))/c_cp

               hes_cup(k,i)=c_cp*t_cup(k,i)+c_alvl*qes_cup(k,i)+c_grav*z_cup  (k,i)
               he_cup (k,i)=c_cp*t_cup(k,i)+c_alvl*q_cup  (k,i)+c_grav*z_cup  (k,i)
               he_cup (k,i)=min(hes_cup(k,i),he_cup(k,i))

               gamma_cup(k,i)=(c_alvl/c_cp)*(c_alvl/(c_rm*t_cup(k,i)*t_cup(k,i)))*qes_cup(k,i)
               u_cup    (k,i)=us(k,i)
               v_cup    (k,i)=vs(k,i)
            end do
         end do
      else
         stop "cup_env_clev"
      endif

      return
      !-- for checking only
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            do k=kts,kte-1
               rho=100*(p_cup(k,i)-p_cup(k+1,i))/(z_cup(k+1,i)-z_cup(k,i))/c_grav ! air dens by hidrostatic balance (kg/m3)
               write(23,101) i,k,z_cup(k,i),p_cup(k,i),t_cup(k,i),q_cup(k,i)*1000.,he_cup(k,i),u_cup(k,i),v_cup(k,i),rho

               rho=100*(p(k,i)-p(k+1,i))/(z(k+1,i)-z(k,i))/c_grav
               write(25,101) i,k,z    (k,i),p    (k,i),t    (k,i),q    (k,i)*1000.,he   (k,i),us   (k,i),vs   (k,i),rho

101            format(2i3,8F15.5)
            end do
            exit
      end do
     

   end subroutine cup_env_clev
   !------------------------------------------------------------------------------------
   subroutine cup_forcing_ens_mid(aa0,aa1,xaa0,mbdt,dtime,ierr,po_cup,ktop,k22,kbcon,kpbl &
                                 ,ichoice, maxens,maxens3,itf,ktf,its,ite, kts,kte        &
                                 ,tau_ecmwf,aa1_bl,xf_dicycle, dhdt,xff_mid,zws,hc,hco    &
                                 ,he_cup,heo_cup,wlpool,xf_coldpool)

      implicit none
      ! aa0     = cloud work function without forcing effects
      ! aa1     = cloud work function with forcing effects
      ! xaa0    = cloud work function with cloud effects
      ! mbdt    = arbitrary numerical parameter
      ! dtime   = dt over which forcing is applied
      ! kbcon   = LFC of parcel from k22
      ! k22     = updraft originating level
      ! ichoice = flag if only want one closure
      ! name    = deep,mid or shallow convection flag
      !
      integer  ,intent (in   )    :: itf,ktf,its,ite, kts,kte,maxens,maxens3
      integer, dimension (:)   ,intent (in) :: k22,kbcon,ktop,kpbl
      real,    dimension (:,:) ,intent (in) :: po_cup,dhdt,hc,hco,he_cup,heo_cup
      real,    dimension (:)   ,intent (in) :: xaa0
      real,    dimension (:)   ,intent (in) :: aa1,zws,mbdt,  aa0
      real                     ,intent (in) :: dtime
      integer                  ,intent (in) :: ichoice
      real,    dimension (:)   ,intent (in) :: aa1_bl,tau_ecmwf,wlpool
      integer, dimension (:)   ,intent (inout) :: ierr
      real,    dimension (:)   ,intent (inout) :: xf_dicycle,xf_coldpool
      real,    dimension (:,:) ,intent (out)   :: xff_mid
      !
      !  local variables in this routine
      !
      real,    dimension (1:maxens) :: xk
      integer                       :: i,k,vtp_index
      real                          :: xff_dicycle, trash, blqe,xff_ens1,mf_ens1

      !-initialization
      xff_mid     (:,:)= 0.
      xf_dicycle  (:)  = 0.

      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)

         xk(1)=(xaa0(i)-(aa1(i)))/mbdt(i)

         if(xk(1).le.0.and.xk(1).gt.-0.1*mbdt(i)) xk(1)=-0.1*mbdt(i)
         if(xk(1).gt.0.and.xk(1).lt.1.e-2       ) xk(1)=1.e-2

         !- closure 3 for mid
         if(xk(1) < 0.) xff_mid(i,3)=max(0., -(aa1(i)/tau_ecmwf(i))/xk(1))
      enddo

      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         !- Boundary layer quasi-equilibrium (Raymond 1995)
         if(k22(i).lt.kpbl(i)+1)then
            blqe=0.
            do k=kts,kbcon(i) !- orig formulation
               !do k=kts,kpbl(i)
               blqe = blqe+100.*dhdt(k,i)*(po_cup(k,i)-po_cup(k+1,i))/c_grav
            enddo
            !trash = max((hc (kbcon(i),i)-he_cup (kbcon(i),i)),1.e1)!- orig formulation
            trash = max((hco(kbcon(i),i)-heo_cup(kbcon(i),i)),1.e1)
            xff_mid(i,2) = max(0.,blqe/trash)
         endif

         !- W* closure (Grant,2001)
         xff_mid(i,1)=0.03*zws(i)
      enddo
      
   end subroutine cup_forcing_ens_mid
   !------------------------------------------------------------------------------------
   subroutine cup_minimi(array,ks,kend,kt,ierr,itf,ktf,its,ite, kts,kte  )

      implicit none
      integer   ,intent (in   )               :: itf,ktf,its,ite, kts,kte
      ! array input array
      ! x output array with return values
      ! kt output array of levels
      ! ks,kend  check-range
      real,    dimension (:,:),intent (in   ) :: array
      integer, dimension (:)  ,intent (in   ) :: ierr,ks,kend
      integer, dimension (:)  ,intent (out  ) :: kt
      
      !-- local vars
      real,    dimension (its:ite)            :: x
      integer                                 :: i,k,kstop,vtp_index

      kt(:) = ks(:)
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         x(i)=array(ks(i),i)
         kstop=max(ks(i)+1,kend(i))
         !
         do k=ks(i)+1,kstop
            if(array(k,i) < x(i)) then
                  x (i) = array(k,i)
                  kt(i) = k
            endif
         enddo
      enddo
   end subroutine cup_minimi
   !------------------------------------------------------------------------------------
   subroutine cup_up_aa0(aa0,z_cup,zu,dby,GAMMA_CUP,t_cup,k22,klcl,kbcon,ktop  &
                        ,ierr,itf,ktf,its,ite, kts,kte,integ_interval          )
      implicit none
      ! aa0 cloud work function
      ! gamma_cup = gamma on model cloud levels
      ! t_cup = temperature (Kelvin) on model cloud levels
      ! dby = buoancy term
      ! zu= normalized updraft mass flux
      ! z = heights of model levels
      ! ierr error value, maybe modified in this routine
      ! on input
      character(len=*),optional, intent(in) :: integ_interval
      integer                  , intent(in) :: itf,ktf,its,ite, kts,kte
      integer, dimension (:)   , intent (in):: k22,klcl,kbcon,ktop

      real,    dimension (:,:) , intent (in):: z_cup,zu,gamma_cup,t_cup,dby
      !
      ! input and output
      integer, dimension (:) ,intent (in   )  :: ierr
      real,    dimension (:) ,intent (out  )  :: aa0
      !
      !  local variables in this routine
      integer                             ::  i,k,vtp_index
      real                                ::  dz,da,aa_2,aa_1
      integer, dimension (its:ite)        ::  kbeg,kend
      !
      !
      !  initialize array to zero.
      aa0(:)=0.

      !--  set domain of integration
      if(present(integ_interval)) then
         if(integ_interval == 'BL') then
            kbeg(:) = kts
            kend(:) = kbcon(:)-1
         elseif(integ_interval == 'CIN') then
            kbeg(:) = k22(:) 
            kend(:) = kbcon(:)-1
         else
            stop "unknown range in cup_up_aa0"
         endif
      else
         kbeg(:) = kbcon(:)
         kend(:) = ktop (:)
      endif

      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         do k= kbeg(i),kend(i)

            dz=z_cup(k+1,i)-z_cup(k,i)
            aa_1=zu(k,i  )*(c_grav/(c_cp*t_cup(k,i  )))*dby(k,i  )/(1.+gamma_cup(k,i  ))
            aa_2=zu(k+1,i)*(c_grav/(c_cp*t_cup(k+1,i)))*dby(k+1,i)/(1.+gamma_cup(k+1,i))
            da=0.5*(aa_1+aa_2)*dz

            aa0(i)=aa0(i)+da

           !aa0(i)=aa0(i)+max(0.,da)
         enddo
      enddo

   end subroutine cup_up_aa0
   !------------------------------------------------------------------------------------
   subroutine cup_up_moisture(name,start_level,klcl,ierr,ierrc,z_cup,qc,qrc,pw,pwav,hc,tempc &
                             ,xland,po,p_cup,kbcon,ktop,cd,dby,clw_all,t_cup,q,gamma_cup     &
                             ,zu,qes_cup,k22,qe_cup,zqexec,use_excess,ccn,rho,up_massentr    &
                             ,up_massdetr,psum,psumh,x_add_buoy,vvel2d,vvel1d,zws            &
                             ,entr_rate_2d,itest,itf,ktf,its,ite, kts,kte,z1,dtime,cnvcf,qes )

      implicit none
      real, parameter :: FRACT = 1.
      !
      !  on input
      integer  ,intent (in   ) ::  use_excess,itest,itf,ktf    &
                                  ,its,ite, kts,kte
      real, intent(in) :: dtime
      ! cd= detrainment function
      ! q = environmental q on model levels
      ! qe_cup = environmental q on model cloud levels
      ! qes_cup = saturation q on model cloud levels
      ! dby = buoancy term
      ! cd= detrainment function
      ! zu = normalized updraft mass flux
      ! gamma_cup = gamma on model cloud levels
      !
      character *(*)              ,intent (in) ::  name
      integer, dimension (:)      ,intent (in) ::  kbcon,ktop,k22,klcl,start_level
      real,    dimension (:,:)    ,intent (in) ::  t_cup,p_cup,rho,q,zu,gamma_cup       &
                                                  ,qe_cup,hc,po,up_massentr,up_massdetr &
                                                  ,dby,qes_cup,z_cup,cd,qes

      real,  dimension (:)        ,intent (in) ::  zqexec,xland,x_add_buoy
      real,  dimension (:)        ,intent (in) ::  zws,ccn,z1
      real,  dimension (:,:)      ,intent (in) ::  entr_rate_2d
      real,  dimension (:,:)      ,intent (in) ::  cnvcf
      real,  dimension (:,:)      ,intent (in) ::  vvel2d
      real,  dimension (:)        ,intent (in) ::  vvel1d
      !
      ! input and output
      !
      ! ierr error value, maybe modified in this routine
      integer, dimension (:)  ,intent (inout)  ::  ierr
      ! qc = cloud q (including liquid water) after entrainment
      ! qrch = saturation q in cloud
      ! qrc = liquid water content in cloud after rainout
      ! pw = condensate that will fall out at that level
      ! pwav = totan normalized integrated condensate (I1)
      ! c0 = conversion rate (cloud to rain)

      real,   dimension (:,:)     ,intent (out)   :: qc,qrc,pw,clw_all,tempc
      real,   dimension (:)       ,intent (out)   :: pwav,psum,psumh
      character*128               ,intent (inout) :: ierrc(:)
      !
      !  local variables in this routine
      !
      integer :: iounit,iprop,i,k,k1,k2,n,nsteps,vtp_index
      real :: dp,rhoc,dh,qrch,dz,radius,berryc0,q1,berryc
      real :: qaver,denom,aux,cx0,qrci,step,cbf,qrc_crit_BF,min_liq,qavail
      real :: delt,tem1,qrc_0,cup,hei,qrc_crit_yin,dz1m,vs,Nc,w_upd,Faut
      real :: q_env_eff
      integer, parameter :: n_smooth=1
      real, parameter :: qrc_crit_yin_ref=-500.*log(2000./9500)*1.e-6/0.93
      real, parameter :: peaut     = .55       ! collection efficiency
      real, parameter :: r0        = .8e-5     ! 8 microm  in contrast to 10 micro m
      real, parameter :: xmyu      = 1.718e-5  ! the dynamic viscosity kgm-1s-1
      real, parameter :: rhowater  = 1000.  ! density of liquid water at 0^oc (kg m^-3)
      real, parameter :: rhosnow   = 100.   ! density of snow (kg m^-3)
      real, parameter :: rhoair0   = 1.28   ! density of dry air at 0^oc and 1000mb pressure (kg m^-3)

      pwav    (:)   = 0.0
      psum    (:)   = 0.0
      psumh   (:)   = 0.0
      pw      (:,:) = 0.0
      clw_all (:,:) = 0.0
      qrc     (:,:) = 0.0         !--- liq/ice water
      tempc   (:,:) = t_cup (:,:)
      qc      (:,:) = qe_cup(:,:) !--- total water: liq/ice = vapor water

      !--- get boundary condition for qc
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         call get_cloud_bc(name,ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i),qe_cup (kts:kte,i),qaver,k22(i))
         qc  (kts:start_level(i),i) = qaver + zqexec(i) + FRACT* x_add_buoy(i)/c_alvl
         qrc (kts:start_level(i),i) = 0.
      enddo

      !--- option to produce linear fluxes in the sub-cloud layer.
      if(name == 'shallow' .and. use_linear_subcl_mf == 1) then
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            call get_delmix(name,kts,kte,ktf,xland(i),start_level(i),po(kts:kte,i) &
                           ,qe_cup(kts:kte,i), qc(kts:kte,i))
         enddo
      endif
      
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)

         do k=start_level(i) + 1,ktop(i) + 1
            dz=z_cup(k+1,i)-z_cup(k,i)
            !
            !--- saturation  in cloud, this is what is allowed to be in it
            !
            qrch = qes_cup(k,i)+(1./c_alvl)*(gamma_cup(k,i)/(1.+gamma_cup(k,i)))*dby(k,i)

            !--- steady state plume equation, for what could be in cloud without condensation
            denom =  (zu(k-1,i)-.5*up_massdetr(k-1,i)+up_massentr(k-1,i)) + 1.e-12
            
            if(use_pass_cloudvol == 2 .or. use_pass_cloudvol == 12 .or. use_pass_cloudvol == 3) then 
               q_env_eff = (1.-cnvcf(k-1,i))*q (k-1,i) + cnvcf(k-1,i)*qes(k-1,i)

               qc (k,i)  = (qc (k-1,i)*zu(k-1,i)-.5*up_massdetr(k-1,i)* qc(k-1,i) +   &
                                                    up_massentr(k-1,i)* q_env_eff)/ denom
            else
               qc (k,i)=  (qc (k-1,i)*zu(k-1,i)-.5*up_massdetr(k-1,i)* qc(k-1,i) +   &
                                                   up_massentr(k-1,i)* q (k-1,i))/ denom
            endif
            


            if(k==start_level(i)+1) &
            qc(k,i)= qc(k,i) + (zqexec(i) + FRACT*x_add_buoy(i)/c_alvl) * up_massentr(k-1,i)/denom

            !--- assuming no liq/ice water in the environment
            qrc(k,i)=  (qrc(k-1,i)*zu(k-1,i)-.5*up_massdetr(k-1,i)* qrc(k-1,i))/ denom

            !-- updraft temp
            tempc(k,i) = (1./c_cp)*(hc(k,i)-c_grav*z_cup(k,i)-c_alvl*QRCH)

            !--- total condensed water before rainout
            clw_all(k,i)= max(0.,qc(k,i)-qrch)

            qrc   (k,i) = min(clw_all(k,i),qrc(k,i))

            !--- production term => condensation/diffusional growth
            cup         = max(0.,qc(k,i)-qrch-qrc(k,i))/dz

            if(c0 < 1.e-6)then
               qrc (k,i) = clw_all(k,i)
               qc  (k,i) = qrc(k,i)+min(qc(k,i),qrch)
               pwav(i)   = 0.
               psum(i)   = psum(i)+clw_all(k,i)*zu(k,i) *dz
               cycle
            endif

            if (autoconv == 1 ) then
               min_liq  = qrc_crit * ( xland(i)*1. + (1.-xland(i))*0.7 )
               if(trim(name) == 'mid') min_liq = min_liq*0.5 
               
               cx0     = c0*DZ
               !cx0  = 1.e-1 * float(JL)/5. *DZ * 1.e-3
               qrc(k,i)= clw_all(k,i)/(1.+cx0)
               pw (k,i)= cx0*max(0.,qrc(k,i) - min_liq)! units kg[rain]/kg[air]

               !--- convert pw to normalized pw
               pw (k,i)=pw(k,i)*zu(k,i)

             elseif (autoconv == 2 ) then
               
               min_liq  = qrc_crit * ( xland(i)*1. + (1.-xland(i))*0.7 )
                              
               !-- Yin, J., et al, 2015. J Meteorol Res.
               !-- includes Yin's factor 
               !hei = max(2000.,min(9000.,z1(i)+0.5*(z_cup(k,i)+z_cup(k-1,i))))/9500.
               !qrc_crit_yin = -500.*log(hei)*1.e-6/rho(k,i)  ! units kg/kg
               !min_liq=(qrc_crit_yin/qrc_crit_yin_ref)*min_liq
               !---
               !
               if(trim(name).eq.'mid') min_liq = min_liq*0.5 

               !-- Eq: qrc + pw = clw_all,  if qrc>min_liq
               !--     pw       = cx0 * (qrc-min_liq) * dz
               !-- =>  qrc      = (clw_all + c0*dz*min_liq) / (1 + c0*dz)

               if(clw_all(k,i) <= min_liq) then !=> more heating at upper levels, more detrained ice
                  qrc(k,i) = clw_all(k,i)
                  pw (k,i) = 0.
               else
                 !-- autoconversion factor
                  cx0     = c0*DZ
                 !cx0  = 1.e-2 * float(JL)/1. *DZ * 1.e-3
              
                  qrc(k,i) = (clw_all(k,i)+min_liq*cx0)/(1.+cx0)  ! units kg[cloud/ice]/kg[air]
                  pw (k,i) = cx0*(qrc(k,i) - min_liq)             ! units kg[rain]/kg[air]
                  !--- convert pw to normalized pw
                  pw (k,i) = pw(k,i)*zu(k,i)
          
               endif

            elseif (autoconv == 3 ) then
              !NC=300.e+6 over the land, NC=50.e+6 over the oceans
               Nc = n_cldrop * 1.e+6; Faut= 1.5
               w_upd = min(15.,max(vvel2d(k,i),2.))
               !cx0 = 1350. * (Nc * 1.e-6)**(-1.79) * (dz/w_upd)*Faut
               !pw (k,i)= cx0*(clw_all(k,i))**2.47! units kg[rain]/kg[air]
             
               cx0 = 7.98e+10 * (Nc * 1.e-6)**(-3.01) *Faut
               !-- Yin, J., et al, 2015. J Meteorol Res.
               !-- includes Yin's factor 
               !hei = max(2000.,min(9000.,z1(i)+0.5*(z_cup(k,i)+z_cup(k-1,i))))/9500.
               !qrc_crit_yin = -500.*log(hei)*1.e-6/rho(k,i)  ! units kg/kg
               !cx0 = 1.0*(qrc_crit_yin/qrc_crit_yin_ref)*cx0
               
              !pw (k,i)= cx0*exp(log(clw_all(k,i))*4.22)*(dz/w_upd)! units kg[rain]/kg[air]
               pw (k,i)= cx0*exp(log(clw_all(k,i))*4.22)*dtime     ! units kg[rain]/kg[air]
               
               pw (k,i)= min( pw(k,i),clw_all(k,i) )
               qrc(k,i)=clw_all(k,i)-pw(k,i)     ! units kg[cloud/ice]/kg[air]
               
               if(qrc(k,i) < 0.) stop 'qrc(k,i) < 0. AUTO 3'
               !--- convert pw to normalized pw
               pw (k,i)=pw(k,i)*zu(k,i)

            elseif (autoconv == 4 ) then
              !NC=300.e+6 over the land, NC=50.e+6 over the oceans
               Nc = n_cldrop *1.e+6 
               min_liq  = 4./3.*c_pi*rhowater*r0**3*Nc/rhoair0  ! 0.419e-3 -- .61e-3
              
              if(clw_all(k,i) <= min_liq) then !=> more heating at upper levels, more detrained ice
              
                  qrc(k,i)= clw_all(k,i)
                  pw(k,i) = 0.
              
              else
                  w_upd = min(15.,max(vvel2d(k,i),2.))
                 
                  cx0 = .104*9.8*peaut/(Nc*rhowater)**(1./3.)/xmyu*rhoair0**(4./3.)! 7.03

                 !pw (k,i)= cx0*(clw_all(k,i)**7./3.)* dz/w_upd ! units kg[rain]/kg[air]
                  pw (k,i)= cx0*exp(log(clw_all(k,i))*((7./3.)))* dz/w_upd ! units kg[rain]/kg[air]

                  pw (k,i)= min( pw(k,i),clw_all(k,i) )
                  qrc(k,i)=clw_all(k,i)-pw(k,i)          ! units kg[cloud/ice]/kg[air]
               
                  if(qrc(k,i) < 0.) stop 'qrc(k,i) < 0. AUTO 4'
                  !--- convert pw to normalized pw
                  pw (k,i)=pw(k,i)*zu(k,i)
               endif

            elseif (autoconv == 5 ) then
               !  c0_deep     = 1.5e-3; c0_mid     = 1.5e-3 ; qrc_crit        = 1.e-4 !(kg/kg)
               min_liq  = qrc_crit * ( xland(i)*0.4 + (1.-xland(i))*1. )
               if(clw_all(k,i) <= min_liq) then !=> more heating at upper levels, more detrained ice

                  qrc(k,i)= clw_all(k,i)
                  pw(k,i) = 0.
               else

                  cx0     =  c0*(1.+ 0.33*FractLiqF(tempc(k,i)))
                  !--- v0
                  qrc(k,i)= qrc(k,i)*exp(-cx0*dz) + (cup/cx0)*(1.-exp(-cx0*dz))
                  qrc(k,i)= max(qrc(k,i),min_liq)
                  pw (k,i)= max(0.,clw_all(k,i)-qrc(k,i)) ! units kg[rain]/kg[air]
                  qrc(k,i)= clw_all(k,i)-pw(k,i)
                  !--- v1
                  !  qrc_0   = qrc(k,i)
                  !  qrc(k,i)= (qrc_0-min_liq)*exp(-cx0*dz) + (cup/cx0)*(1.-exp(-cx0*dz))+min_liq
                  !  qrc(k,i)= max(qrc(k,i),min_liq)
                  !  pw (k,i)= max(0.,clw_all(k,i)-qrc(k,i)) ! units kg[rain]/kg[air]
                  !  qrc(k,i)= clw_all(k,i)-pw(k,i)

                  !  qrc(k,i)= (clw_all(k,i)-min_liq)*exp(-cx0*dz)+min_liq
                  !  pw (k,i)= clw_all(k,i)-qrc(k,i) ! units kg[rain]/kg[air]
                  !--- v3
                  !  qrc(k,i)= (clw_all(k,i)-min_liq) / (1.+cx0*dz)+min_liq
                  !  pw (k,i)= cx0*dz*(qrc(k,i)-min_liq) ! units kg[rain]/kg[air]
                  !  print*,"BG=",k,real(cx0*1.e+3,4),real(pw(k,i),4),real(qrc(k,i),4)&
                  !              ,real(clw_all(k,i)-pw(k,i)-qrc(k,i),4) !==> must be zero

                  !--- convert pw to normalized pw
                  pw (k,i)= pw(k,i)*zu(k,i)
               endif
            endif
            !- total water (vapor + condensed) in updraft after the rainout
            qc(k,i)=qrc(k,i)+min(qc(k,i),qrch)

            !--- integrated normalized condensates
            pwav(i)=pwav(i)+pw(k,i)
            psum(i)=psum(i)+clw_all(k,i)*zu(k,i) *dz

         enddo
         if(pwav(i) < 0.) then
               ierr(i)=66 ; is_removed = remove(vec_ok, i)
               ierrc(i)="pwav negative"
         endif

      enddo

      !--- get back water vapor qc
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         do k=start_level(i)+1,kbcon(i) + 1
            vs  = 0.
            dz1m= 0.
            do k1 = max(k-n_smooth,kts),min(k+n_smooth,ktf)
               dz   = z_cup(k1+1,i)-z_cup(k1,i)
               vs   =  vs + dz*qrc(k1,i)
               dz1m = dz1m + dz
            enddo
            qrc(k,i) = vs/(1.e-16+dz1m)
            
         enddo
         do k=kts,ktop(i)+1
            qc(k,i)= qc(k,i)-qrc(k,i)
         enddo
      enddo
   end subroutine cup_up_moisture
   !------------------------------------------------------------------------------------
   subroutine cup_up_moisture_light(name,start_level,klcl,ierr,ierrc,z_cup,qc,qrc,pw,pwav,hc          &
                                   ,tempc,xland,po,p_cup,kbcon,ktop,cd,dby,clw_all,t_cup,q,gamma_cup  &
                                   ,zu,qes_cup,k22,qe_cup,zqexec,use_excess,rho ,up_massentr          &
                                   ,up_massdetr,psum,psumh,x_add_buoy,itest,itf,ktf,its,ite, kts,kte)

      implicit none
      !
      !  on input
      integer  ,intent (in   ) ::  use_excess,itest,itf,ktf    &
                                  ,its,ite, kts,kte
      ! cd= detrainment function
      ! q = environmental q on model levels
      ! qe_cup = environmental q on model cloud levels
      ! qes_cup = saturation q on model cloud levels
      ! dby = buoancy term
      ! cd= detrainment function
      ! zu = normalized updraft mass flux
      ! gamma_cup = gamma on model cloud levels
      !
      character *(*)           ,intent (in) ::  name
      integer, dimension (:)   ,intent (in) ::  kbcon,ktop,k22,klcl,start_level
      real,    dimension (:,:) ,intent (in) ::  t_cup,p_cup,rho,q,zu,gamma_cup       &
                                               ,qe_cup,hc,po,up_massentr,up_massdetr &
                                               ,dby,qes_cup,z_cup,cd

      real,  dimension (:)     ,intent (in) ::  zqexec,xland,x_add_buoy
      !
      ! input and output
      !
      ! ierr error value, maybe modified in this routine
      integer, dimension (:)  ,intent (inout)  ::  ierr
      !
      ! qc = cloud q (including liquid water) after entrainment
      ! qrch = saturation q in cloud
      ! qrc = liquid water content in cloud after rainout
      ! pw = condensate that will fall out at that level
      ! pwav = totan normalized integrated condensate (I1)
      ! c0 = conversion rate (cloud to rain)

      real,   dimension (:,:) ,intent (out)   :: qc,qrc,pw,clw_all,tempc
      real,   dimension (:)   ,intent (out)   :: pwav,psum,psumh
      character*128           ,intent (inout) :: ierrc(:)
      !
      !  local variables in this routine
      !
      real, parameter :: FRACT = 1.
      integer         :: iounit,iprop,i,k,k1,k2,n,nsteps,vtp_index
      real            :: dp,rhoc,dh,qrch,dz,radius,berryc0,q1,berryc            &
                        ,qaver,denom,aux,cx0,qrci,step,cbf,qrc_crit_BF,min_liq  &
                        ,qavail,delt_hc_glac,delt,tem1

      pwav    (:)  =0.0
      psum    (:)  =0.0
      psumh   (:)  =0.0
      pw      (:,:)=0.0
      qrc     (:,:)=0.0
      clw_all (:,:)=0.0
      tempc   (:,:)=t_cup (:,:)
      qc      (:,:)=qe_cup(:,:)

      !--- get boundary condition for qc
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         call get_cloud_bc(name,ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i),qe_cup (kts:kte,i),qaver,k22(i))
         qc (kts:start_level(i),i) = qaver + zqexec(i) + FRACT*x_add_buoy(i)/c_alvl
      enddo

      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)

         do k=start_level(i)+1,ktop(i) + 1

            dz=z_cup(k+1,i)-z_cup(k,i)
            !
            !--- saturation  in cloud, this is what is allowed to be in it
            !
            qrch = qes_cup(k,i)+(1./c_alvl)*(gamma_cup(k,i)/(1.+gamma_cup(k,i)))*dby(k,i)

            !-    1. steady state plume equation, for what could
            !-       be in cloud without condensation
            denom   =  (zu(k-1,i)-.5*up_massdetr(k-1,i)+up_massentr(k-1,i)) + 1.e-12

            qc (k,i)=  ( qc (k-1,i)*zu(k-1,i)-.5*up_massdetr(k-1,i)* qc(k-1,i) +   &
                                                 up_massentr(k-1,i)* q (k-1,i)     )/ denom
            
            if(k==start_level(i)+1) &
              qc(k,i)= qc(k,i) + (zqexec(i) + fract*x_add_buoy(i)/c_alvl)*up_massentr(k-1,i)/denom

            !--- total condensed water before rainout
            clw_all(k,i)=max(0.,qc(k,i)-qrch)
            !--- updraft temp
            tempc(k,i) = (1./c_cp)*(hc(k,i) - c_grav*z_cup(k,i)-c_alvl*qrch)

            !--add glaciation effect on the mse
            if(melt_glac) then
               delt_hc_glac = clw_all(k,i)*(1.- FractLiqF(tempc(k,i)))*c_xlf

               tempc(k,i) = tempc(k,i)+(1./c_cp)*delt_hc_glac
            endif

            cx0     = c0*dz
            if(c0 < 1.e-6) cx0 = 0.

            qrc(k,i)= clw_all(k,i)/(1.+cx0)
            pw (k,i)= cx0*max(0.,qrc(k,i) -qrc_crit)! units kg[rain]/kg[air]
            !--- convert pw to normalized pw
            pw (k,i)= pw(k,i)*zu(k,i)

            !- total water (vapor + condensed) in updraft after the rainout
            qc(k,i)=qrc(k,i)+min(qc(k,i),qrch)

         enddo
      enddo

      !- get back water vapor qc
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         do k=kts,ktop(i)+1
            qc(k,i)= qc(k,i)-qrc(k,i)
         enddo
      enddo
   end subroutine cup_up_moisture_light
   !------------------------------------------------------------------------------------
   subroutine cup_up_aa1bl(version,aa1_bl,aa1_fa,aa1,t,tn,q,qo,dtime,po_cup,z_cup,zu,dby   &
                          ,gamma_cup,t_cup,rho,klcl,kpbl,kbcon,ktop,ierr,itf,ktf,its,ite   &
                          ,kts,kte,xland,ztexec,xlons,xlats, h_sfc_flux,le_sfc_flux,tau_bl &
                          ,tau_ecmwf,t_star,cumulus ,tn_bl,qo_bl)

      implicit none
      character*(*), intent (in)            :: cumulus
      integer      , intent (in)            ::itf,ktf,its,ite, kts,kte,version
      ! aa0 cloud work function
      ! gamma_cup = gamma on model cloud levels
      ! t_cup = temperature (Kelvin) on model cloud levels
      ! dby = buoancy term
      ! zu= normalized updraft mass flux
      ! z = heights of model levels
      ! ierr error value, maybe modified in this routine
      !
      integer, dimension (:)   ,intent (in  ) :: klcl,kbcon,ktop,kpbl
      real,    dimension (:,:) ,intent (in  ) :: z_cup,zu,gamma_cup,t_cup,dby,t,tn,q,qo &
                                                ,po_cup,rho,tn_bl,qo_bl
      real                     ,intent(in   ) :: dtime, t_star

      real,    dimension (:)   ,intent (in  ) :: xland,ztexec,xlons,xlats, h_sfc_flux &
                                                ,le_sfc_flux,aa1,tau_bl,tau_ecmwf
      !
      ! input and output
      integer, dimension (:)   ,intent (inout) :: ierr
      real,    dimension (:)   ,intent (out  ) :: aa1_bl,aa1_fa
      !
      !  local variables in this routine
      !
      integer                               ::    i,k,vtp_index
      real                                  ::    dz,da,aa_1,aa_2,tcup,da_bl,a1_bl
      !
      !
      aa1_bl (:)=0.
      if(version == 0 ) then
      do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            !***       do k=kts,kbcon(i)
            do k=kts,kpbl(i)
               dz = c_grav * (z_cup (k+1,i)-z_cup (k,i))
               da = dz*(tn(k,i)*(1.+0.608*qo(k,i))-t(k,i)*(1.+0.608*q(k,i)))/dtime
               aa1_bl(i)=aa1_bl(i)+da ! Units : J K / (kg seg)
            enddo
         enddo
      elseif(version==1) then
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            do k=kts,kpbl(i)
               dz = (z_cup (k+1,i)-z_cup (k,i))
               aa_1=(c_grav/(c_cp*t_cup(k,i  )))*dby(k,i  )*zu(k,i  )
               aa_2=(c_grav/(c_cp*t_cup(k+1,i)))*dby(k+1,i)*zu(k+1,i)
               da=0.5*(aa_1+aa_2)*dz! Units : J / kg
               aa1_bl(i)=aa1_bl(i)+da
            enddo
         enddo
      else
         stop "unknown version option in routine: cup_up_aa1bl"
      endif

      return
      
      aa1_fa (:)=0.
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         do k= kbcon(i),ktop(i)

            dz=z_cup(k+1,i)-z_cup(k,i)
            aa_1=(c_grav/(c_cp*((t_cup(k,i  )))))*dby(k,i  )/(1.+gamma_cup(k,i  ))*zu(k,i)
            aa_2=(c_grav/(c_cp*((t_cup(k+1,i)))))*dby(k+1,i)/(1.+gamma_cup(k+1,i))*zu(k+1,i)
            da=0.5*(aa_1+aa_2)*dz
            aa1_fa(i)=aa1_fa(i)+da
         enddo
      enddo

   end subroutine cup_up_aa1bl
   !------------------------------------------------------------------------------------
   subroutine get_lateral_massflux(itf,ktf, its,ite, kts,kte,min_entr_rate              &
                                  ,ierr,ktop,zo_cup,zuo,cd,entr_rate_2d,po_cup          &
                                  ,up_massentro, up_massdetro                           &
                                  ,draft,kbcon,k22,kpbl,up_massentru,up_massdetru,lambau)
      implicit none
      character *(*), intent (in) :: draft
      integer, intent(in) :: itf,ktf, its,ite, kts,kte
      real,    intent(in) :: min_entr_rate
      integer, intent(in)   , dimension(:)            :: ierr,ktop,kbcon,k22,kpbl
      real,    intent(in)   , dimension(:), optional  :: lambau
      real,    intent(in)   , dimension(:,:) :: zo_cup,zuo,po_cup
      real,    intent(inout), dimension(:,:) :: cd,entr_rate_2d
      real,    intent(out)  , dimension(:,:) :: up_massentro, up_massdetro
      real,    intent(out)  , dimension(:,:), optional :: up_massentru,up_massdetru
      !
      !-- local vars
      integer :: i,k, turn, ismooth1,ismooth2, vtp_index
      real :: dz, mass1,mass2,dp,rho,zuo_ave
      logical  :: SMOOTH
      integer, parameter :: MASS_U_OPTION = 1
      integer, parameter :: SMOOTH_DEPTH  = 2 ! -- increasing this parameter,
                                              ! -- strongly damps the heat/drying rates, precip ...
      integer ::  incr1=1, incr2=1, nlay, k_ent

      !---
      SMOOTH = .false.
      if(USE_SMOOTH_PROF == 1)  SMOOTH = .true.
      
      up_massentro(:,:)=0.
      up_massdetro(:,:)=0.
      if(present(up_massentru) .and. present(up_massdetru))then
         up_massentru(:,:)=0.
         up_massdetru(:,:)=0.
      endif
      nlay = int(kte/90)

      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         
         !-will not allow detrainment below the location of the maximum zu
         ! if(draft=='shallow'.or.draft == 'mid') cd(i,1:maxloc(zuo(:,i),1)-2)=0.0

          !-will not allow detrainment below cloud base or in the PBL
         if(draft=='shallow') then
            cd(1:max(kbcon(i),kpbl(i))+nlay,i)=0.0

         else
            cd(1:maxloc(zuo(:,i),1)+nlay,i)=0.0
         endif

         !- mass entrainment and detrainment are defined on model levels
         do k=kts,maxloc(zuo(:,i),1)
            !=> below location of maximum value zu -> change entrainment

            dz=zo_cup(k+1,i)-zo_cup(k,i)
            zuo_ave = 0.5*(zuo(k+1,i)+zuo(k,i))

            up_massdetro(k,i)=cd(k,i)*dz*zuo_ave
            up_massentro(k,i)=zuo(k+1,i)-zuo(k,i)+up_massdetro(k,i)
            !-- check limits of allowed entrainment rates 
            up_massentro(k,i)=max(up_massentro(k,i),min_entr_rate*dz*zuo_ave)
            

            !-- check up_massdetro in case of up_massentro has been changed above
            up_massdetro(k,i)=-zuo(k+1,i)+zuo(k,i)+up_massentro(k,i)

            cd          (k,i)=up_massdetro(k,i)/(dz*zuo_ave)
            entr_rate_2d(k,i)=up_massentro(k,i)/(dz*zuo_ave)

         enddo

         !--- limit the effective entrainment rate
         k_ent=maxloc(zuo(:,i),1)
         do k=k_ent+1,ktop(i)-1
            entr_rate_2d(k,i)=entr_rate_2d(k_ent,i)*(min(zo_cup(k_ent,i)/zo_cup(k,i),1.))
            entr_rate_2d(k,i)=max(min_entr_rate, entr_rate_2d(k,i))
         enddo
         entr_rate_2d(ktop(i):kte,i)=0.

         !=================
         if(SMOOTH .and. draft /= 'shallow' ) then
            !---smoothing the transition zone (from maxloc(zu)-1 to maxloc(zu)+1)

            ismooth1 = max(kts+2, maxloc(zuo(:,i),1) - SMOOTH_DEPTH)
            ismooth2 = min(ktf-2, maxloc(zuo(:,i),1) + SMOOTH_DEPTH)
            !if(draft == 'shallow') ismooth1 = max(ismooth1,max(kbcon(i),kpbl(i))+nlay)+1

            do k=ismooth1,ismooth2
               dz=zo_cup(k+1,i)-zo_cup(k,i)

               zuo_ave = 0.5*(zuo(k+1,i)+zuo(k,i))

               up_massentro(k,i)=0.5*(entr_rate_2d(k,i)*dz*zuo_ave+up_massentro(k-1,i))

               up_massdetro(k,i)=zuo(k,i)+up_massentro(k,i)-zuo(k+1,i)

               if(up_massdetro(k,i).lt.0.)then
                  up_massdetro(k,i)=0.
                  up_massentro(k,i)=zuo(k+1,i)-zuo(k,i)
                  entr_rate_2d(k,i)=(up_massentro(k,i))/(dz*zuo_ave)
               endif
               if(zuo_ave > 0.) &
                  cd(k,i)=up_massdetro(k,i)/(dz*zuo_ave)
            enddo

            do k=ismooth1,ismooth2
               dz=zo_cup(k+1,i)-zo_cup(k,i)

               zuo_ave = 0.5*(zuo(k+1,i)+zuo(k,i))

               up_massdetro(k,i)=0.5*(cd(k,i)*dz*zuo_ave+up_massdetro(k-1,i))
               up_massentro(k,i)=zuo(k+1,i)-zuo(k,i)+up_massdetro(k,i)

               if(up_massentro(k,i).lt.0.)then
                  up_massentro(k,i)=0.
                  up_massdetro(k,i)=zuo(k,i)-zuo(k+1,i)
                  cd(k,i)=up_massdetro(k,i)/(dz*zuo_ave)
               endif
               if(zuo_ave > 0.) &
                  entr_rate_2d(k,i)=(up_massentro(k,i))/(dz*zuo_ave)
            enddo
          !-----end of the transition zone
         endif
         !=================

         do k=maxloc(zuo(:,i),1)+incr1 ,ktop(i)
            !=> above location of maximum value zu -> change detrainment
            dz=zo_cup(k+1,i)-zo_cup(k,i)
            zuo_ave = 0.5*(zuo(k+1,i)+zuo(k,i))

            up_massentro(k,i)=entr_rate_2d(k,i)*dz*zuo_ave

            up_massdetro(k,i)=zuo(k,i)+up_massentro(k,i)-zuo(k+1,i)
            up_massdetro(k,i)=max(up_massdetro(k,i),0.0)
            !-- check up_massentro in case of dd_up_massdetro has been changed above
            up_massentro(k,i)=-zuo(k,i)+up_massdetro(k,i)+zuo(k+1,i)

            if(zuo_ave.gt.0.) then
               cd          (k,i)=up_massdetro(k,i)/(dz*zuo_ave)
               entr_rate_2d(k,i)=up_massentro(k,i)/(dz*zuo_ave)
            endif
         enddo

         !do k=kts,kte
         !   up_massentr(k,i)=up_massentro(k,i)
         !   up_massdetr(k,i)=up_massdetro(k,i)
         !enddo
         if(present(up_massentru) .and. present(up_massdetru))then
            if(mass_U_option==1) then
               do k=kts+1,kte
                  !--       for weaker mixing
                  up_massentru(k-1,i)=up_massentro(k-1,i)+lambau(i)*up_massdetro(k-1,i)
                  up_massdetru(k-1,i)=up_massdetro(k-1,i)+lambau(i)*up_massdetro(k-1,i)
               !--       for stronger mixing
               ! up_massentru(k-1,i)=up_massentro(k-1,i)+lambau(i)*up_massentro(k-1,i)
               ! up_massdetru(k-1,i)=up_massdetro(k-1,i)+lambau(i)*up_massentro(k-1,i)
               enddo
            else
               turn=maxloc(zuo(:,i),1)
               do k=kts+1,turn
                  up_massentru(k-1,i)=up_massentro(k-1,i)+lambau(i)*up_massentro(k-1,i)
                  up_massdetru(k-1,i)=up_massdetro(k-1,i)+lambau(i)*up_massentro(k-1,i)
               enddo
               do k=turn+1,kte
                  up_massentru(k-1,i)=up_massentro(k-1,i)+lambau(i)*up_massdetro(k-1,i)
                  up_massdetru(k-1,i)=up_massdetro(k-1,i)+lambau(i)*up_massdetro(k-1,i)
               enddo
            endif
         endif
         do k=ktop(i)+1,kte
            cd          (k,i)=0.
            entr_rate_2d(k,i)=0.
         enddo
      enddo ! i

      return 
      !---- check mass conservation
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         do k=kts+1,kte
            dz =      zo_cup(k,i)-zo_cup(k-1,i)
            dp = 100*(po_cup(k,i)-po_cup(k-1,i))
            rho= -dp/dz/c_grav
            mass1= (zuo(k,i)-zuo(k-1,i)) - up_massentro(k-1,i)+up_massdetro(k-1,i)
            !print*,"masscons=",mass1!,-rho*c_grav*(zuo(k,i)-zuo(k-1,i))/dp, (zuo(k,i)-zuo(k-1,i))/dz,( up_massentro(k-1,i)-up_massdetro(k-1,i))/dz,rho
            mass2= (zuo(k,i)-zuo(k-1,i)) - up_massentru(k-1,i)+up_massdetru(k-1,i)
         enddo
      enddo
   end subroutine get_lateral_massflux
   !------------------------------------------------------------------------------------
   subroutine get_lateral_massflux_down(cumulus,itf,ktf, its,ite, kts,kte ,ierr,jmin,zo_cup,zdo,cdd &
                                       ,mentrd_rate_2d,dd_massentro,dd_massdetro,draft,mentrd_rate  &
                                       ,dd_massentru,dd_massdetru,lambau)

      implicit none
      character *(*), intent (in)            :: draft,cumulus
      integer, intent(in)                    :: itf,ktf, its,ite, kts,kte
      real,    intent(in)   , dimension(:)   :: mentrd_rate
      integer, intent(in)   , dimension(:)   :: ierr,jmin
      real,    intent(in)   , dimension(:)   :: lambau
      real,    intent(in)   , dimension(:,:) :: zo_cup,zdo
      real,    intent(inout), dimension(:,:) :: cdd,mentrd_rate_2d
      real,    intent(  out), dimension(:,:) :: dd_massentro, dd_massdetro
      real,    intent(  out), dimension(:,:), optional :: dd_massentru, dd_massdetru
      integer ::i,ki,vtp_index
      real :: dzo

      cdd          = 0.
      dd_massentro = 0.
      dd_massdetro = 0.
      if(present(dd_massentru).and.present(dd_massdetru))then
         dd_massentru = 0.
         dd_massdetru = 0.
      endif
      if(cumulus == 'shallow') return

      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)

         mentrd_rate_2d(1:jmin(i)  ,i) = mentrd_rate(i)
         cdd           (1:jmin(i)-1,i) = mentrd_rate(i)
         mentrd_rate_2d(1          ,i) = 0.
 
         do ki=jmin(i)   ,maxloc(zdo(:,i),1),-1

            !=> from jmin to maximum value zd -> change entrainment
            dzo=zo_cup(ki+1,i)-zo_cup(ki,i)
            dd_massdetro(ki,i)=cdd(ki,i)*dzo*zdo(ki+1,i)
            !XXX
            dd_massentro(ki,i)= zdo(ki,i)-zdo(ki+1,i)+dd_massdetro(ki,i)
            dd_massentro(ki,i)= MAX(0.,dd_massentro(ki,i))
            !-- check dd_massdetro in case of dd_massentro has been changed above
            dd_massdetro(ki,i)= dd_massentro(ki,i)-zdo(ki,i)+zdo(ki+1,i)

              !~ if(dd_massentro(ki,i).lt.0.)then
              !~ dd_massentro(ki,i)=0.
              !~ dd_massdetro(ki,i)=zdo(ki+1,i)-zdo(ki,i)
              !~ if(zdo(ki+1,i) > 0.0)&
              !~ cdd(ki,i)=dd_massdetro(ki,i)/(dzo*zdo(ki+1,i))
              !~ endif
              !~ if(zdo(ki+1,i) > 0.0)&
              !~ mentrd_rate_2d(ki,i)=dd_massentro(ki,i)/(dzo*zdo(ki+1,i))
         enddo

         do ki=maxloc(zdo(:,i),1)-1,kts,-1
            !=> from maximum value zd to surface -> change detrainment
            dzo=zo_cup(ki+1,i)-zo_cup(ki,i)
            dd_massentro(ki,i)=mentrd_rate_2d(ki,i)*dzo*zdo(ki+1,i)
            !XXX
            dd_massdetro(ki,i) = zdo(ki+1,i)+dd_massentro(ki,i)-zdo(ki,i)
            dd_massdetro(ki,i) = MAX(0.0,dd_massdetro(ki,i))
            !-- check dd_massentro in case of dd_massdetro has been changed above
            dd_massentro(ki,i) = dd_massdetro(ki,i)+zdo(ki,i)-zdo(ki+1,i)


             !~ if(dd_massdetro(ki,i).lt.0.)then
             !~ dd_massdetro(ki,i)=0.
             !~ dd_massentro(ki,i)=zdo(ki,i)-zdo(ki+1,i)
             !~ if(zdo(ki+1,i) > 0.0)&
             !~ mentrd_rate_2d(ki,i)=dd_massentro(ki,i)/(dzo*zdo(ki+1,i))
             !~ endif
             !~ if(zdo(ki+1,i) > 0.0)&
             !~ cdd(ki,i)= dd_massdetro(ki,i)/(dzo*zdo(ki+1,i))
         enddo

         if(present(dd_massentru).and.present(dd_massdetru))then
            do ki=jmin(i),kts,-1
               dd_massentru(ki,i)=dd_massentro(ki,i)+lambau(i)*dd_massdetro(ki,i)
               dd_massdetru(ki,i)=dd_massdetro(ki,i)+lambau(i)*dd_massdetro(ki,i)
            enddo
         endif
      enddo

   end subroutine get_lateral_massflux_down
   !------------------------------------------------------------------------------------
   subroutine get_zu_zd_pdf(vtp_index,cumulus, draft,ierr,kb,kt,zu,kts,kte,ktf,kpbli &
                           ,k22,kbcon,klcl,po_cup,psur,xland,random,x_add_buoy)

      implicit none
      character*(*), intent(in) ::draft,cumulus
      integer, intent(in   ) :: kts,kte,ktf,kpbli,k22,kbcon,kt,kb,klcl,vtp_index
      integer, intent(inout) :: ierr
      real   , intent(in   ) :: po_cup(:),psur,xland,random,x_add_buoy
      real   , intent(inout) :: zu(:)
      !
      !- local var
      integer :: i,k,kb_adj,kpbli_adj,level_max_zu
      real :: zumax,ztop_adj,beta, alpha,kratio,tunning,FZU,krmax,dzudk,hei_updf,hei_down
      real :: zul(kts:kte),  pmaxzu ! pressure height of max zu for deep

      integer:: minzu,maxzul,maxzuh,kstart
      logical :: do_smooth

      integer :: k1
      real :: wgty,dp_layer,slope,adj_cp_entr
     
      DO_SMOOTH = .true.
      !if(USE_SMOOTH_PROF == 1)  DO_SMOOTH = .true.

      !-- fill zu with zeros
      zul=0.0

      adj_cp_entr = 1.
      if(draft == "deep_up" .and. (use_memory == 22 .or. use_memory == 222)) then
           !-- higher hei_max_zu_updraft, lower effec entrainment, lower precip
           ! adj_cp_entr = min(1.3,max(1.0, 2.-coldPoolStart(x_add_buoy)))
           !
           !-- lower hei_max_zu_updraft, larger effec entrainment, larger precip
             adj_cp_entr = min(1.0,max(0.7,    coldPoolStart(x_add_buoy)))

      endif
      if(draft == "deep_up" .or. draft == "mid_up" ) then       !--- land/ocean

         hei_updf=(1.-xland)*hei_updf_LAND+xland*hei_updf_OCEAN
         
         !- add a randomic perturbation
         hei_updf = hei_updf + random
         !
         !- adjust for cloud organization
         hei_updf = hei_updf * adj_cp_entr
         !
         !- sanity check
         hei_updf = max(0.1, min(0.9, hei_updf)) 
         
         !- for gate soundings
         !hei_updf = max(0.2, min(0.8,(float(JL))/100.))
         
         !--hei_updf parameter goes from 0 to 1 = rainfall decreases with hei_updf
         pmaxzu  =  (psur-100.) * (1.- 0.5*hei_updf) + 0.6*( po_cup(kt) ) * 0.5*hei_updf

         !- beta parameter: must be larger than 1, higher makes the profile sharper around the maximum zu
         !beta   = max(1.1, 2.2 - 0.8*hei_updf)
         beta    = min(3.0, 3.5 - 1.8*hei_updf)
         
         !- for gate soundings
         !beta   = 3.1-float(JL)/100. 
         
         kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)
         kb_adj=max(kb,kb_adj)
         kb_adj=min(kb_adj,kt)

         !- this alpha constrains the location of the maximun ZU to be at "kb_adj" vertical level
         alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))

         !
         do k=klcl-1,min(kte,kt)
            kratio= float(k)/float(kt+1)
            zu(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
         enddo

         !-- special treatment below k22/klcl
         do k=klcl,kts+1,-1
            zu(k)=zu(k+1)*0.5
         enddo

        !if(use_linear_subcl_mf == 1) then
        !    zu(kts)=0.
        !    kstart=kbcon
        !    slope=(zu(kstart)-zu(kts))/(po_cup(kstart)-po_cup(kts) + 1.e-6)
        !    do k=kstart-1,kts+1,-1
        !       zu(k) = zu(kstart)-slope*(po_cup(kstart)-po_cup(k))
        !      !print*,"k=",zu(kstart),zu(k),zu(kts)
        !    enddo
        !    go to 333
        ! endif
  
         !-- smooth section
         if(do_smooth) then
            !--from surface
            zul(kts+1)=zu(kts+1)*0.25
            do k=kts+2,maxloc(zu,1)
               zul(k)=zu(k-1)*0.8+zu(k)*0.2
            enddo
            do k=kts+1,maxloc(zu,1)
               zu(k)=(zul(k)+zu(k))*0.5
            enddo

            !--from ktop
            ! zul(kt)=zu(kt)*0.1
            ! do k=kt-1,max( kt-min(maxloc(zu,1),5),kts) ,-1
            !    zul(k)=(zul(k+1)+zu(k))*0.5
            ! enddo
            ! wgty=0.0
            ! do k=kt,max( kt-min(maxloc(zu,1),5), kts ),-1
            !    wgty=wgty+1./(float(min(maxloc(zu,1),5))+1)
            !    zu(k)=zul(k)*(1.-wgty)+ zu(k)*wgty
            ! enddo
         endif
         zu(kts)=0.
      !---------------------------------------------------------
      elseif(draft == "shallow_up") then
         kb_adj   =kts     ! level where mass flux starts
         kpbli_adj=kpbli
         if(kpbli_adj < kb_adj .or. kpbli_adj >= kt ) then
            kpbli_adj = kb_adj + 1
         endif

         !- location of the maximum Zu: dp_layer mbar above PBL height
         !dp_layer     = 10. !mbar
         !level_max_zu = minloc(abs(po_cup(kts:kt+1)-(po_cup(kpbli_adj)-dp_layer)),1)
         !

         k1           = max(kbcon,kpbli_adj)
         !- location of the maximum Zu: dp_layer mbar above k1 height
         hei_updf     =(1.-xland)*hei_updf_LAND+xland*hei_updf_OCEAN

         !hei_updf = (float(JL)-20)/40. ; print*,"JL=",jl,hei_updf

         dp_layer     = hei_updf*(po_cup(k1)-po_cup(kt))

         level_max_zu = minloc(abs(po_cup(kts:kt+1)-(po_cup(k1)-dp_layer)),1)
         level_max_zu = min(level_max_zu,kt -1)
         level_max_zu = max(level_max_zu,kts+1)

         krmax        = float(level_max_zu)/float(kt+1)
         krmax        = min(krmax,0.99)

         beta= beta_sh!smaller => sharper detrainment layer
         !beta= ((1.-xland)*0.43 +xland)*beta_sh

         !beta= 3.0!smaller => sharper detrainment layer
         !beta = 1.+4.*(float(JL))/40.

         !- this alpha imposes the maximum zu at kpbli
         alpha=1.+krmax*(beta-1.)/(1.-krmax)
         !alpha=min(6.,alpha)

         !- to check if dZu/dk = 0 at k=kpbli_adj
         !kratio=krmax
         !dzudk=(alpha-1.)*(kratio**(alpha-2.)) * (1.-kratio)**(beta-1.) - &
         !          (kratio**(alpha-1.))*((1.-kratio)**(beta-2.))*(beta-1.)

         !- Beta PDF
         do k=kts+1,min(kte,kt)
            kratio= float(k)/float(kt+1)
            zu(k)=kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
         enddo
         zu(kts)=0.
         !
         !-- special treatment below kbcon - linear Zu
         if(use_linear_subcl_mf == 1) then
            kstart=kbcon
            slope=(zu(kstart)-zu(kts))/(po_cup(kstart)-po_cup(kts) + 1.e-6)
            do k=kstart-1,kts+1,-1
               zu(k) = zu(kstart)-slope*(po_cup(kstart)-po_cup(k))
              !print*,"k=",zu(kstart),zu(k),zu(kts)
            enddo
         endif
         !-- special treatment below kclcl
         !do k=(klcl-1),kts+1,-1
         !  zu(k)=zu(k+1)*0.5
         !enddo
         !
         !-- smooth section
         !IF( .not. do_smooth) then
         ! zul(kts+1)=zu(kts+1)*0.1
         ! do k=kts+2,maxloc(zu,1)
         !    zul(k)=(zu(k-1)+zu(k))*0.5
         ! enddo
         ! do k=kts+1,maxloc(zu,1)
         !    zu(k)=(zul(k)+zu(k))*0.5
         ! enddo
         !ENDIF
         !zu(kts)=0.

      !---------------------------------------------------------
      elseif(draft == "DOWN" ) then
         if(cumulus == 'shallow') return

         beta =2.5
         hei_down=(1.-xland)*hei_down_LAND+xland*hei_down_OCEAN

         pmaxzu= hei_down * po_cup(kt) + (1.-hei_down)*psur
         kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)

         !- this alpha constrains the location of the maximun ZU to be at "kb_adj" vertical level
         alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))

         do k=kts+1,min(kt+1,ktf)
            kratio= float(k)/float(kt+1)
            zu(k)= kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
         enddo
         !-- smooth section
         if(do_smooth) then
            zul(kts+1)=zu(kts+1)*0.1
            wgty=0.
            do k=kts+2,maxloc(zu,1)
               wgty=wgty+1./(float(max(2,maxloc(zu,1)))-1.)
               wgty=0.5
               !print*,"zD1=",k,zu(k),zul(k-1),wgty,zul(k-1)*(1.-wgty)+ zu(k)*wgty
               zul(k)=zul(k-1)*(1.-wgty)+ zu(k)*wgty
            enddo
            wgty=0.
            do k=kts+1,maxloc(zu,1)
               wgty=wgty+1./(float(max(2,maxloc(zu,1)))-1.)
               wgty=0.5
               !print*,"zD2=",k,zu(k),zul(k),wgty,zul(k)*(1.-wgty)+ zu(k)*wgty
               zu(k)=zul(k)*(1.-wgty)+ zu(k)*wgty
            enddo
         endif
         zu(kts)=0.

      endif

      !---------------------------------------------------------

      if( maxval(zu(kts:min(kte,kt+1)),1) <= 0.0) then
         zu=0.0
         ierr=51  ; is_removed = remove(vec_ok, vtp_index) !ierr(i)=51
      else
         !- normalize ZU
         zu(kts:min(kte,kt+1))= zu(kts:min(kte,kt+1))/ (1.e-9+maxval(zu(kts:min(kte,kt+1)),1))
      endif

   end subroutine get_zu_zd_pdf
   !------------------------------------------------------------------------------------
   subroutine cup_up_cape_cin(aa0,z_cup,zu,dby,GAMMA_CUP,t_cup,klcl ,k22,kbcon,ktop,ierr     &
                             ,tempco,qco,qrco, qo_cup,itf,ktf,its,ite, kts,kte,integ_interval)

      implicit none
      integer ,intent (in   )                   ::   itf,ktf, its,ite, kts,kte

      ! aa0 = dummy array for CAPE (total cape)
      ! gamma_cup = gamma on model cloud levels
      ! t_cup = temperature (Kelvin) on model cloud levels
      ! dby = buoancy term
      ! zu= normalized updraft mass flux
      ! z = heights of model levels
      ! ierr = error value, maybe modified in this routine
      ! tempco = in-cloud temperature (Kelvin) on model cloud levels
      ! qco    = in-cloud water vapor mixing ratio on model cloud levels
      ! qo_cup = environ water vapor mixing ratio on model cloud levels
      ! qrco   = in-cloud liquid water mixing ratio on model cloud levels
      character(len=*), intent(in) :: integ_interval

      real,    dimension (:,:) ,intent (in   ) :: z_cup,zu,gamma_cup,t_cup,dby,tempco &
                                                 ,qco,qrco, qo_cup
      integer, dimension (:)   ,intent (in   ) :: k22,kbcon,ktop,klcl
      !
      ! input and output
      integer, dimension (:)   ,intent (inout) ::  ierr
      real,    dimension (:)   ,intent (out  ) ::  aa0
      !
      !  local variables in this routine
      integer       ::   i,k,vtp_index
      real          ::   dz,daa0,daa1,daa2
      integer, dimension (its:ite) ::  kbeg,kend
      !

      if(integ_interval == 'CIN') then 
         kbeg(:) = k22(:)
         kend(:) = kbcon(:)-1
      else                       ! CAPE
         kbeg(:) = kbcon(:)
         kend(:) = ktop (:)
      endif

      aa0(:)=0.
      do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            do k= kbeg(i),kend(i)
               dz=z_cup(k+1,i)-z_cup(k,i)
               daa1= c_grav*dz*(   (tempco(k,i)  *(1.+0.608*qco   (k,i))  - t_cup(k,i)*(1.+0.608*qo_cup(k,i)))&
                            / (t_cup (k,i)  *(1.+0.608*qo_cup(k,i))) )

               daa2=c_grav*dz*(    (tempco(k+1,i)*(1.+0.608*qco   (k+1,i))  - t_cup(k+1,i)*(1.+0.608*qo_cup(k+1,i)))&
                            / (t_cup (k+1,i)*(1.+0.608*qo_cup(k+1,i))) )
               
               aa0(i)=aa0(i)+0.5*(daa1+daa2)
            enddo
            ! print*,trim(integ_interval),AA0(I),kbcon(i)
      enddo
   end subroutine cup_up_cape_cin
   !------------------------------------------------------------------------------------
   subroutine cup_cloud_limits(name,ierrc,ierr,cap_inc,cap_max_in,heo_cup,heso_cup,qo_cup &
                              ,qeso_cup,po,po_cup,z_cup,heo,hkbo,qo,qeso,entr_rate_2d,hcot&
                              ,k22,kbmax,klcl,kbcon,ktop,depth_neg_buoy,frh,Tpert         &
                              ,start_level_,use_excess,zqexec,ztexec, x_add_buoy,xland    &
                              ,itf,ktf,its,ite, kts,kte)

      implicit none
      character *(*)            ,intent (in   ) ::  name
      character*128             ,intent (inout) :: ierrc(:)
      integer                   ,intent (in   ) :: itf,ktf,its,ite, kts,kte,use_excess
      integer, dimension (:)    ,intent (in   ) :: kbmax,start_level_
      integer, dimension (:)    ,intent (inout) :: kbcon,ierr,ktop,klcl,k22
      real,    dimension (:,:)  ,intent (in   ) :: heo_cup,heso_cup,po_cup,z_cup,heo &
                                                  ,qo_cup,qeso_cup,po,qo,qeso,Tpert
      real,    dimension (:)    ,intent (in   ) :: cap_max_in,cap_inc,xland
      real,    dimension (:)    ,intent (in   ) :: zqexec,ztexec,x_add_buoy
      real,    dimension (:)    ,intent (inout) :: hkbo,depth_neg_buoy,frh
      real,    dimension (:,:)  ,intent (in)    :: entr_rate_2d
      real,    dimension (:,:)  ,intent (inout) :: hcot

      !  local variables in this routine
      real, parameter              :: frh_crit_O=0.7
      real, parameter              :: frh_crit_L=0.7  !--- test 0.5
      real                         :: delz_oversh !--- height of cloud overshoot is 10% higher than the LNB.
                                                  !--- Typically it can 2 - 2.5km higher, but it depends on
                                                  !--- the severity of the thunderstorm.
      real,    dimension (its:ite) ::   cap_max
      integer                      ::   i,k,k1,k2,kfinalzu,vtp_index
      real                         ::   plus,hetest,dz,dbythresh,denom &
                                       ,dzh,del_cap_max,fx,x_add,Z_overshoot,frh_crit
      real   , dimension (kts:kte) ::   dby
      integer, dimension (its:ite) ::   start_level

      delz_oversh = OVERSHOOT
      start_level(:) = 0
      ktop       (:) = 0 !ktf-1
      dby        (:) = 0.0
      cap_max    (:) = cap_max_in(:)
      hcot       (:,:) = 0.0

      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         start_level(i)             = start_level_(i)
         hcot(kts:start_level(i),i) = hkbo(i) ! assumed no entrainment between these layers
      enddo
      !
      !--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
      !
      loop0: do i=its,itf
         kbcon         (i) = kbmax(i)+3
         frh           (i) = 0.
         depth_neg_buoy(i) = 0.

         if(ierr(i) /= 0) cycle
      !
      !loop0: do vtp_index = get_num_elements(vec_ok),1,-1
      !  i = get_data_value(vec_ok,vtp_index)

         loop1:  do while(ierr(i) == 0)

            kbcon(i)=start_level(i)
            do k=start_level(i)+1,KBMAX(i)+3
               dz=z_cup(k,i)-z_cup(k-1,i)
               hcot(k,i)= ( (1.-0.5*entr_rate_2d(k-1,i)*dz)*hcot(k-1,i)     &
                                  + entr_rate_2d(k-1,i)*dz *heo (k-1,i) ) / &
                          (  1.+0.5*entr_rate_2d(k-1,i)*dz)
               if(k==start_level(i)+1) then
                  x_add    = (c_alvl*zqexec(i)+c_cp*ztexec(i)) + x_add_buoy(i)
                  hcot(k,i)= hcot(k,i) +  x_add
               endif
            enddo

            loop2:      do while (hcot(kbcon(i),i) < HESO_cup(kbcon(i),i))
               kbcon(i)=kbcon(i)+1
               if(kbcon(i).gt.kbmax(i)+2) then
                  ierr(i)=3 ; is_removed = remove(vec_ok, i)
                  ierrc(i)="could not find reasonable kbcon in cup_kbcon : above kbmax+2 "
                  exit loop2
               endif
                !print*,"kbcon=",kbcon(i);call flush(6)
            enddo loop2

            if(ierr(i) /= 0) cycle loop0

            !---     cloud base pressure and max moist static energy pressure
            !---     i.e., the depth (in mb) of the layer of negative buoyancy
            depth_neg_buoy(i) = - (po_cup(kbcon(i),i)-po_cup(start_level(i),i))

            if(MOIST_TRIGGER == 1) then
               frh(i)=0.
               dzh = 0
               do k=k22(i),kbcon(i)
                  dz     = z_cup(k,i)-z_cup(max(k-1,kts),i)
                  frh(i) = frh(i) + dz*(qo(k,i)/qeso(k,i))
                  dzh    = dzh + dz
                  !print*,"frh=", k,dz,qo(k,i)/qeso(k,i)
               enddo
               frh(i) = frh(i)/(dzh+1.e-16)
               frh_crit =frh_crit_O*xland(i) + frh_crit_L*(1.-xland(i))

               !fx     = 4.*(frh(i) - frh_crit)* abs(frh(i) - frh_crit) !-quadratic
               fx     = ((2./0.78)*exp(-(frh(i) - frh_crit)**2)*(frh(i) - frh_crit)) !- exponential
               fx     = max(-1.,min(1.,fx))

               del_cap_max = fx* cap_inc(i)
               cap_max(i)  = min(max(cap_max_in(i) + del_cap_max, 10.),150.)
               !print*,"frh=", frh(i),kbcon(i),del_cap_max, cap_max(i)!,  cap_max_in(i)
            endif

            !- test if the air parcel has enough energy to reach the positive buoyant region
            if(cap_max(i) >= depth_neg_buoy(i)) cycle loop0


            !--- use this for just one search (original k22)
            !            if(cap_max(i) < depth_neg_buoy(i)) then
            !                    ierr(i)=3
            !                    ierrc(i)="could not find reasonable kbcon in cup_cloud_limits"
            !            endif
            !            cycle loop0
            !---

            !- if am here -> kbcon not found for air parcels from k22 level
            k22(i)=k22(i)+1
            !--- increase capmax
            !if(USE_MEMORY == 2000) cap_max(i)=cap_max(i)+cap_inc(i)

            !- get new hkbo
            x_add = (c_alvl*zqexec(i)+c_cp*ztexec(i)) +  x_add_buoy(i)
            call get_cloud_bc(name,ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i),heo_cup (kts:kte,i),hkbo (i),k22(i),x_add,Tpert(kts:kte,i))
            !
            start_level(i)=start_level(i)+1
            !
            hcot(start_level(i),i)=hkbo (i)
         enddo loop1
         !--- last check for kbcon
         if(kbcon(i) == kts) then
            ierr(i)=33 ; is_removed = remove(vec_ok, i)
            ierrc(i)="could not find reasonable kbcon in cup_kbcon = kts"
         endif
      enddo loop0
      !
      !
      !--- DETERMINE THE LEVEL OF NEUTRAL BUOYANCY - KTOP
      !
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)

         start_level(i)=kbcon(i)

         do k=start_level(i)+1,ktf-1
            dz=z_cup(k,i)-z_cup(k-1,i)
            denom = 1.+0.5*entr_rate_2d(k-1,i)*dz + 1.e-12
            hcot(k,i)=( (1.-0.5*entr_rate_2d(k-1,i)*dz)*hcot(k-1,i)    &
                               +entr_rate_2d(k-1,i)*dz *heo (k-1,i) )/ denom
         enddo
         do k=start_level(i)+1,ktf-1

            if(hcot(k,i) < heso_cup(k,i) )then
               ktop(i)  =  k - 1
               exit
            endif
         enddo
         if(ktop(i) <= kbcon(i)+1) then
              ierr(i)=41 ; is_removed = remove(vec_ok, i)
         endif 

         !----------------
         if(OVERSHOOT > 1.e-6 .and. ierr(i) == 0) then
            Z_overshoot = (1. + delz_oversh) * z_cup(ktop(i),i)
            do k=ktop(i),ktf-2
               if(Z_overshoot < z_cup(k,i)) then
                  ktop(i) = min(k-1, ktf-2)
                  exit
               endif
            enddo
         endif
      enddo
   end subroutine cup_cloud_limits
   !------------------------------------------------------------------------------------
   subroutine get_buoyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop &
                          ,hc,he_cup,hes_cup,dby,z_cup)

      implicit none
      integer                    ,intent (in)  :: itf,ktf,its,ite, kts,kte
      integer, dimension (:)     ,intent (in)  :: ierr,klcl,kbcon,ktop
      real,    dimension (:,:)   ,intent (in)  :: hc,he_cup,hes_cup,z_cup
      real,    dimension (:,:)   ,intent (out) :: dby
      integer :: i,k,vtp_index

      dby(:,:) = 0.0
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         do k=kts,klcl(i)-1
             dby (k,i)=hc (k,i)-0.5*(he_cup (k,i)+hes_cup (k,i))
         enddo
          do  k=klcl(i),ktop(i)+1
             dby (k,i)=hc (k,i)-hes_cup (k,i)
         enddo
         !do  k=kts,ktop(i)+1
         !   dby (k,i)=hc (k,i)-hes_cup (k,i)
         !enddo
      enddo
   end subroutine get_buoyancy

   !------------------------------------------------------------------------------------
   subroutine cup_up_vvel(vvel2d,vvel1d,zws,entr_rate_2d,cd ,z,z_cup,zu,dby,GAMMA_CUP,t_cup &
                         ,tempco,qco,qrco,qo,start_level,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte &
                         ,wlpool,wlpool_bcon)

      implicit none
      integer                  ,intent (in   )  ::  itf,ktf,its,ite, kts,kte
      real,    dimension (:,:) ,intent (in   )  ::  z,z_cup,zu,gamma_cup,t_cup,dby &
                                                   ,entr_rate_2d,cd,tempco,qco,qrco,qo

      integer, dimension (:)   ,intent (in   )  ::  klcl,kbcon,ktop,start_level
      real,    dimension (:)   ,intent (in   )  ::  zws
      real,    dimension (:)   ,intent (inout)  ::  wlpool

      ! input and output
      integer, dimension (:)   ,intent (inout) ::  ierr
      real,    dimension (:,:) ,intent (out  ) ::  vvel2d
      real,    dimension (:)   ,intent (out  ) ::  vvel1d
      real,    dimension (:)   ,intent (inout) ::  wlpool_bcon

      !
      !  local variables in this routine
      integer            ::  i,k,k1,nvs,vtp_index
      real               ::  dz,BU,dw2,dw1,kx,dz1m,Tv,Tve,vs,ke
      real   , parameter :: f=2., C_d=0.506, gam=0.5, beta=1.875, eps=0.622
      real   , parameter :: ftun1=0.10, ftun2=1.0
      logical, parameter :: smooth=.true.
      integer, parameter :: n_smooth=1
      
      !-- initialize arrays to zero.
      vvel1d(:  ) = 0.0
      vvel2d(:,:) = 0.0

      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)

         vvel2d(kts:kbcon(i),i)= max(1.,max(wlpool_bcon(i)**2,zws(i)**2))

         loop0:  do k= kbcon(i),ktop(i)

            dz=z_cup(k+1,i)-z_cup(k,i)

            Tve= 0.5* ( t_cup (k,i  )*(1.+(qo (k,i  )/eps)/(1.+qo (k,i  ))) + &
                        t_cup (k+1,i)*(1.+(qo (k+1,i)/eps)/(1.+qo (k+1,i))))

            Tv = 0.5* ( tempco(k,i  )*(1.+(qco(k,i  )/eps)/(1.+qco(k,i  ))) + &
                        tempco(k+1,i)*(1.+(qco(k+1,i)/eps)/(1.+qco(k+1,i)) ))

            BU = c_grav*( (Tv-Tve)/Tve -  ftun2*0.50*(qrco(k+1,i)+qrco(k,i) ))

            dw1 = 2./(f*(1.+gam)) * BU * dz
            kx  = (1.+beta*C_d)*max(entr_rate_2d(k,i),cd(k,i))*dz*ftun1
            
            dw2 =  (vvel2d(k,i)) -2.*kx * (vvel2d(k,i))

            vvel2d(k+1,i)=(dw1+dw2)/(1.+kx)

            if( vvel2d(k+1,i)< 0.) then
               vvel2d(k+1,i) = 0.5* vvel2d(k,i)
            endif

         enddo loop0
       enddo
       if(smooth) then
            do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               do k=kts,ktop(i)+1
                  vs =0.
                  dz1m= 0.
                  do k1 = max(k-n_smooth,kts),min(k+n_smooth,ktf)
                     dz   = z_cup(k1+1,i)-z_cup(k1,i)
                     vs   =  vs + dz*vvel2d(k1,i)
                     dz1m = dz1m + dz
                  enddo
                  vvel2d(k,i) = vs/(1.e-16+dz1m)
               !if(k>ktop(i)-3)print*,"v2=",k,ktop(i),sqrt(vvel2d(k,i)),sqrt(vvel2d(ktop(i),i))
               enddo
            enddo
       endif

      !-- convert to vertical velocity
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         vvel2d(:,i)= sqrt(max(0.1,vvel2d(:,i)))

         if(maxval(vvel2d(:,i)) < 1.0) then
            ierr(i)=54 ; is_removed = remove(vec_ok, i)
         !  print*,"ierr=54",maxval(vvel2d(:,i))
         endif


         !-- sanity check
         where(vvel2d(:,i) < 1. ) vvel2d(:,i) = 1.
         where(vvel2d(:,i) > 30.) vvel2d(:,i) = 30.
         vvel2d(ktop(i)+1:kte,i) = 0.1

         !-- get the column average vert velocity
         do k= kbcon(i),ktop(i)
            dz=z_cup(k+1,i)-z_cup(k,i)
            vvel1d(i)=vvel1d(i)+vvel2d(k,i)*dz
         !print*,"w=",k,z_cup(k,i),vvel2d(k,i)
         enddo
         vvel1d(i)=vvel1d(i)/(z_cup(ktop(i)+1,i)-z_cup(kbcon(i),i)+1.e-16)
         vvel1d(i)=max(1.,vvel1d(i))
      enddo

   end subroutine cup_up_vvel
   !------------------------------------------------------------------------------------
   subroutine cup_output_ens(cumulus,xff_shal,xff_mid,xf_ens,ierr,edto,dellat,dellaq,dellaqc&
                            ,outtem,outq,outqc,zu,pre,pwo,pwdo,pwo_eff,xmb,ktop        &
                            ,nx,nx2,pr_ens, maxens3,ensdim,sig,xland1                  &
                            ,ichoice,itf,ktf,its,ite, kts,kte                  &
                            ,xf_dicycle,outu,outv,dellu,dellv,dtime,po_cup,kbcon       &
                            ,dellabuoy,outbuoy, dellampqi,outmpqi,dellampql,outmpql    &
                            ,dellampcf,outmpcf ,nmp, rh_dicycle_fct,xf_coldpool,wlpool_bcon)
      implicit none
      !
      integer   ,intent (in   )            :: ichoice,itf,ktf,its,ite, kts,kte
      integer   ,intent (in   )            :: ensdim,nx,nx2,maxens3,nmp
      ! xf_ens = ensemble mass fluxes
      ! pr_ens = precipitation ensembles
      ! dellat = change of temperature per unit mass flux of cloud ensemble
      ! dellaq = change of q per unit mass flux of cloud ensemble
      ! dellaqc = change of qc per unit mass flux of cloud ensemble
      ! outtem = output temp tendency (per s)
      ! outq   = output q tendency (per s)
      ! outqc  = output qc tendency (per s)
      ! pre    = output precip
      ! xmb    = total base mass flux
      ! xfac1  = correction factor
      ! pw = pw -epsilon*pd (ensemble dependent)
      ! ierr error value, maybe modified in this routine
      !
      character *(*)            ,intent (in   ) :: cumulus
      real,    dimension (:,:)  ,intent (inout) :: xf_ens,pr_ens
      real,    dimension (:,:)  ,intent (out  ) :: outtem,outq,outqc,outu,outv,outbuoy
      real,    dimension (:,:,:),intent (out  ) :: outmpqi,outmpql,outmpcf
      
      real,    dimension (:,:)  ,intent (in   ) :: zu,po_cup,pwo,pwdo
      real,    dimension (:)    ,intent (in   ) :: sig, rh_dicycle_fct,wlpool_bcon,edto
    
      real,    dimension (:,:)  ,intent (in   ) :: xff_mid
      real,    dimension (:)    ,intent (out  ) :: pre,xmb
      
      real,    dimension (:)    ,intent (in   ) :: xland1
      real,    dimension (:,:)  ,intent (inout) :: pwo_eff,dellat,dellaqc,dellaq &
                                                  ,dellu,dellv,dellabuoy
      real,    dimension (:,:,:),intent (in   ) :: dellampqi,dellampql,dellampcf

      integer, dimension (:)    ,intent (in   ) :: ktop,kbcon
      integer, dimension (:)    ,intent (inout) :: ierr
      real,    dimension (:)    ,intent (inout) :: xf_dicycle,xf_coldpool
      real,                      intent (in   ) :: dtime
      real,   dimension (:,:)   ,intent (in   ) :: xff_shal
      !
      !  local variables in this routine
      !
      integer                          :: i,k,n,ncount,zmax,kk,kqmx,ktmx,vtp_index
      real                             :: outtes,ddtes,dtt,dtq,dtqc &
                                         ,dtpw,prerate,fixouts,dp   &
                                         ,xfixQ,xfixT,dtts,dtqs,fsum, rcount
      real,    dimension (its:ite)     :: xmb_ave,xmbmax
      real,    dimension (8)           :: tend1d
      !
      !-- init/reset 
      outtem  = 0.
      outq    = 0.
      outqc   = 0.
      outu    = 0.
      outv    = 0.
      outbuoy = 0.
      pre     = 0.
      xmb     = 0.
      xmb_ave = 0.
      pwo_eff = 0.
      !
      !--- get the net precipitation (per unit of mass)
      !
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         pwo_eff(:,i)=pwo(:,i)+edto(i)*pwdo(:,i)
      enddo
!
      do vtp_index = get_num_elements(vec_ok),1,-1
          i = get_data_value(vec_ok,vtp_index)
          do n=1,maxens3
               if(pr_ens(i,n) <= 0.)then
                  xf_ens(i,n)=0.
               endif
         enddo
      enddo
      !
      !--- calculate ensemble average mass fluxes
      !
      if(cumulus == 'deep') then
      do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            xmb_ave(i)=0.
            do n=1,maxens3
                  xmb_ave(i)=xmb_ave(i)+xf_ens(i,n)
            enddo
            !- 'ensemble' average mass flux
            xmb_ave(i)=xmb_ave(i)/float(maxens3)
         enddo

      !- mid (congestus type) convection
      elseif(cumulus=='mid') then
         if(ichoice .le. 3) then
            do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               if(ichoice == 0) then
                  xmb_ave(i)=0.3333*(xff_mid(i,1)+xff_mid(i,2)+xff_mid(i,3))
               else
                  xmb_ave(i)= xff_mid(i,ichoice)
               endif
            enddo
         else
            stop 'For mid ichoice must be 0,1,2,3'
         endif

      !- shallow  convection
      elseif(cumulus=='shallow') then
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)

            if(ichoice > 0) then
               xmb_ave(i)=xff_shal(i,ichoice)
            else
               fsum=0.
               xmb_ave(i)=0.
               do k=1,shall_closures
                  !- heat engine closure is not working properly
                  !- turning it off for now.
                  if(k.ge.4 .and. k.le.6) cycle                  
                  xmb_ave(i)=xmb_ave(i)+xff_shal(i,k)
                  fsum=fsum+1.
               enddo
               !- ensemble average of mass flux
               xmb_ave(i)=xmb_ave(i)/fsum
            endif
         enddo
      endif
      !- apply the mean tropospheric RH control on diurnal cycle (Tian GRL 2022)
       if(cumulus == 'deep' .and. rh_dicycle == 1) then
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            xf_dicycle(i) =  xf_dicycle(i) * rh_dicycle_fct(i)
         enddo
      endif

     !- set the updraft mass flux, do not allow negative values and apply the diurnal cycle closure
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         !- mass flux of updradt at cloud base
         xmb(i) = xmb_ave(i)

         !- add kinetic energy at the gust front of the cold pools
         xmb(i) = xmb(i) + xf_coldpool(i)

         !- diurnal cycle closure
         xmb(i) = xmb(i) - xf_dicycle(i)
         if(xmb(i) .le. 0.)then
            ierr(i)=13 ; is_removed = remove(vec_ok, i)
            xmb (i)=0.
         endif
      enddo
      !-apply the scale-dependence Arakawa's approach
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         !- scale dependence
         xmb(i)=sig(i)*xmb(i)

         !- apply the adjust factor for tunning
         !xmb(i) = FADJ_MASSFLX * xmb(i)

         if(xmb(i) == 0. ) then 
           ierr(i)=14 ; is_removed = remove(vec_ok, i)
         endif 
         if(xmb(i) > 100.) then 
           ierr(i)=15 ; is_removed = remove(vec_ok, i)
         endif 
      enddo

      !--- sanity check for mass flux
      !
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         xmbmax(i)=100.*(po_cup(kbcon(i),i)-po_cup(kbcon(i)+1,i))/(c_grav*dtime)
         xmb(i) = min(xmb(i),xmbmax(i))
      enddo

      !--- check outtem and and outq for high values
      !--- criteria: if abs (dT/dt or dQ/dt) > 100 K/day => fix xmb
      if( MAX_TQ_TEND < -1.e-2) then
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            fixouts=xmb(i) *86400.*max(maxval(abs(dellat(kts:ktop(i),i))),&
                         (c_alvl/c_cp)*maxval(abs(dellaq(kts:ktop(i),i))) )

            if(fixouts > abs(MAX_TQ_TEND)) then ! K/day
               fixouts=abs(MAX_TQ_TEND)/(fixouts)
               xmb   (i)  = xmb   (i)  *fixouts
               xf_ens(i,:)= xf_ens(i,:)*fixouts
            endif
         enddo
      endif
       !--- criteria: if abs (dT/dt or dQ/dt) > 100 K/day => fix dT/dt, dQ/dt and xmb
      if( MAX_TQ_TEND > 1.e-2) then
         do vtp_index = get_num_elements(vec_ok),1,-1
           i = get_data_value(vec_ok,vtp_index)
           
           tend1d=0.
           do k=kts,ktop(i)
              dp         = (po_cup(k,i)-po_cup(k+1,i))
              tend1d(1)  = tend1d(1)  +  dp*xmb(i) * 86400.*(dellat(k,i))
 
              if(xmb(i) * 86400.*abs(dellat(k,i)) > MAX_TQ_TEND ) & 
                 dellat(k,i)=MAX_TQ_TEND/(xmb(i)*86400)*sign(1., dellat(k,i))

              tend1d(2)  = tend1d(2)  +  dp*xmb(i) * 86400.*(dellat(k,i))         
           enddo
      
           do k=kts,ktop(i)
              dp         = (po_cup(k,i)-po_cup(k+1,i))
              tend1d(3)  = tend1d(3)  +  dp*xmb(i) * 86400.*(dellaq(k,i))*(c_alvl/c_cp)

              if(xmb(i) * 86400.*abs(dellaq(k,i))*(c_alvl/c_cp)  > MAX_TQ_TEND ) &    
                 dellaq(k,i)=MAX_TQ_TEND/(xmb(i)*86400*(c_alvl/c_cp))*sign(1., dellaq(k,i))
            
              tend1d(4)  = tend1d(4)  +  dp*xmb(i) * 86400.*(dellaq(k,i))*(c_alvl/c_cp)
           enddo
           xfixT = tend1d(1)/(1.e-6+tend1d(2))
           xfixQ = tend1d(3)/(1.e-6+tend1d(4))
         
           xmb(i) = xmb(i)/ max(1.,max(xfixQ,xfixT)) 
          !   print*,"tend",
        enddo
      endif 
      !
      !-- now do feedback
      !
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         do k=kts,ktop(i)
            pre (i)  = pre(i) + pwo_eff(k,i)*xmb(i)

            outtem   (k,i)= dellat     (k,i)*xmb(i)
            outq     (k,i)= dellaq     (k,i)*xmb(i)
            outqc    (k,i)= dellaqc    (k,i)*xmb(i)
            outu     (k,i)= dellu      (k,i)*xmb(i)
            outv     (k,i)= dellv      (k,i)*xmb(i)
            outbuoy  (k,i)= dellabuoy  (k,i)*xmb(i)
         enddo
         xf_ens (i,:)= sig(i)*xf_ens(i,:)


         if(APPLY_SUB_MP == 1) then
            do k=kts,ktop(i)
               outmpqi(:,k,i)= dellampqi(:,k,i)*xmb(i)
               outmpql(:,k,i)= dellampql(:,k,i)*xmb(i)
               outmpcf(:,k,i)= dellampcf(:,k,i)*xmb(i)
            enddo
            outmpqi(:,ktop(i)+1:ktf,i)=0.
            outmpql(:,ktop(i)+1:ktf,i)=0.
            outmpcf(:,ktop(i)+1:ktf,i)=0.
         endif
      enddo
      !
   end subroutine cup_output_ens
   !------------------------------------------------------------------------------------
   subroutine cup_forcing_ens_deep(itf,ktf,its,ite, kts,kte,ens4,ensdim,ichoice,maxens         &
                                  ,maxens2,maxens3,ierr,k22,kbcon,ktop,xland,aa0,aa1,xaa0      &
                                  ,mbdt,dtime,xf_ens,mconv,qo,p_cup,omeg,zd,zu,pr_ens,edt      &
                                  ,tau_ecmwf,aa1_bl,xf_dicycle,xk_x,alpha_adv,Q_adv,aa1_radpbl &
                                  ,aa1_adv,wlpool,xf_coldpool,cin1)

      implicit none
      integer  ,intent (in   ) :: itf,ktf,its,ite, kts,kte,ens4
      integer  ,intent (in   ) :: ensdim,maxens,maxens2,maxens3
      !
      ! ierr error value, maybe modified in this routine
      ! pr_ens = precipitation ensemble
      ! xf_ens = mass flux ensembles
      ! massfln = downdraft mass flux ensembles used in next timestep
      ! omeg = omega from large scale model
      ! mconv = moisture convergence from large scale model
      ! zd      = downdraft normalized mass flux
      ! zu      = updraft normalized mass flux
      ! aa0     = cloud work function without forcing effects
      ! aa1     = cloud work function with forcing effects
      ! xaa0    = cloud work function with cloud effects
      ! edt     = epsilon
      ! dir     = "storm motion"
      ! mbdt    = arbitrary numerical parameter
      ! dtime   = dt over which forcing is applied
      ! iact_gr_old = flag to tell where convection was active
      ! kbcon       = LFC of parcel from k22
      ! k22         = updraft originating level
      ! ichoice       = flag if only want one closure (usually set to zero!)
      ! name        = deep or shallow convection flag
      !
      real,    dimension (:,:)    ,intent (inout)  :: pr_ens
      real,    dimension (:,:)    ,intent (out  )  :: xf_ens
      real,    dimension (:,:)    ,intent (in   )  :: zd,zu,p_cup,qo
      real,    dimension (:,:,:)  ,intent (in   )  :: omeg
      real,    dimension (:)      ,intent (in   )  :: xaa0
      real,    dimension (:)      ,intent (in   )  :: aa1,edt,xland
      real,    dimension (:)      ,intent (inout)  :: mconv
      real,    dimension (:)      ,intent (in   )  :: aa0
      real,    dimension (:)      ,intent (in   )  :: mbdt
      real                        ,intent (in   )  :: dtime
      integer, dimension (:)      ,intent (in   )  :: k22,kbcon,ktop
      integer, dimension (:)      ,intent (inout)  :: ierr
      integer                     ,intent (in   )  :: ichoice

      real,    dimension (:)      ,intent(in)      :: aa1_bl,tau_ecmwf,alpha_adv&
                                                     ,Q_adv,aa1_radpbl,aa1_adv,wlpool &
                                                     ,cin1
      real,    dimension (:)      ,intent(inout)   :: xf_dicycle,xk_x,xf_coldpool
     
      !
      !  local variables in this routine
      !
      integer                          :: i,k,kk,nall,n,ne,nens,nens3,vtp_index
      real                             :: xff_dicycle
      real                             :: a1,a_ave,xff0,xomg,KE_gf,W_cb
      real, parameter                  :: c1 = 0.06, c2 = 1., c3 = 0.28,  c4 = 0.0 !0.64 
      real, dimension (1:maxens3)      :: xff_ens3
      real, dimension (its:ite)        :: xk
      real, dimension (its:ite)        :: ens_adj
      !
      ens_adj(:)=1.
      xf_ens(:,1:16)= 0.
      !--- LARGE SCALE FORCING
      !
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)

         xff0 = (aa1(i)-aa0(i))/dtime
         !-- default
         xff_ens3(1) = max(0.,(aa1(i)-aa0(i))/dtime)

         xff_ens3(2) = xff_ens3(1)
         xff_ens3(3) = xff_ens3(1)
         xff_ens3(16)= xff_ens3(1)
         !
         !--- more like brown (1979), or frank-cohen (199?)
         !--- omeg is in pa/s
         xomg=0.
         kk=0
         do k=max(kts,kbcon(i)-1),kbcon(i)+1
            xomg=xomg-omeg(k,i,1)/c_grav
            kk=kk+1
         enddo
         xff_ens3(4) = max(0.0,xomg/float(kk)) ! kg[air]/m^3 * m/s
         xff_ens3(5) = xff_ens3(4)
         xff_ens3(6) = xff_ens3(4)
         xff_ens3(14)= xff_ens3(4)
         !
         !--- more like krishnamurti et al.;
         !
         mconv(i) = 0.
         do k=kbcon(i),ktop(i)
             mconv(i)=mconv(i)+omeg(k,i,1)*(qo(k+1,i)-qo(k,i))/c_grav
         enddo

         mconv(i)  = max(0., mconv(i))
         xff_ens3(7) = mconv(i)
         xff_ens3(8) = xff_ens3(7)
         xff_ens3(9) = xff_ens3(7)
         xff_ens3(15)= xff_ens3(7)
         !
         !---- more like  betchold et al (2014). note that aa1 already includes the forcings tendencies
         xff_ens3(10)= aa1(i)/tau_ecmwf(i)

         xff_ens3(11)= xff_ens3(10)
         xff_ens3(12)= xff_ens3(10)
         xff_ens3(13)= xff_ens3(10)

         !
         if(ichoice == 0)then
            if(xff0 < 0.)then
               xff_ens3( 1)=0.
               xff_ens3( 2)=0.
               xff_ens3( 3)=0.
               xff_ens3(16)=0.

               xff_ens3(10)=0.
               xff_ens3(11)=0.
               xff_ens3(12)=0.
               xff_ens3(13)=0.
            endif
         endif

         xk(i)=(xaa0(i)-(aa1(i)))/mbdt(i)
         if(xk(i).le.0. .and. xk(i).gt.-0.1*mbdt(i)) xk(i)=-0.1*mbdt(i)
         if(xk(i).gt.0. .and. xk(i).lt.1.e-2       ) xk(i)=1.e-2
         !
         !--- special treatment for stability closures
         !
         if(xk(i).lt.0.)then
            if(xff_ens3( 1).gt.0.)xf_ens(i, 1)=max(0.,-xff_ens3( 1)/xk(i))
            if(xff_ens3( 2).gt.0.)xf_ens(i, 2)=max(0.,-xff_ens3( 2)/xk(i))
            if(xff_ens3( 3).gt.0.)xf_ens(i, 3)=max(0.,-xff_ens3( 3)/xk(i))
            if(xff_ens3(16).gt.0.)xf_ens(i,16)=max(0.,-xff_ens3(16)/xk(i))
         else
            xff_ens3(1 )=0.
            xff_ens3(2 )=0.
            xff_ens3(3 )=0.
            xff_ens3(16)=0.
         endif

         xf_ens(i,4) =max(0.,xff_ens3(4) )
         xf_ens(i,5) =max(0.,xff_ens3(5) )
         xf_ens(i,6) =max(0.,xff_ens3(6) )
         xf_ens(i,14)=max(0.,xff_ens3(14))

         a1=max(1.e-3,pr_ens(i,7) )
         xf_ens(i,7) =max(0.,xff_ens3(7)/a1)
         a1=max(1.e-3,pr_ens(i,8) )
         xf_ens(i,8) =max(0.,xff_ens3(8)/a1)
         a1=max(1.e-3,pr_ens(i,9) )
         xf_ens(i,9) =max(0.,xff_ens3(9)/a1)
         a1=max(1.e-3,pr_ens(i,15))
         xf_ens(i,15)=max(0.,xff_ens3(15)/a1)
         if(xk(i).lt.0.)then
            xf_ens(i,10)= max(0.,-xff_ens3(10)/xk(i))
            xf_ens(i,11)= max(0.,-xff_ens3(11)/xk(i))
            xf_ens(i,12)= max(0.,-xff_ens3(12)/xk(i))
            xf_ens(i,13)= max(0.,-xff_ens3(13)/xk(i))
         else
            xf_ens(i,10)= 0.
            xf_ens(i,11)= 0.
            xf_ens(i,12)= 0.
            xf_ens(i,13)= 0.
         endif


         if(ichoice.ge.1)then
            xf_ens(i,1:16) =xf_ens(i,ichoice)
         endif

         !---special combination for 'ensemble closure':
         !---over the land, only applies closures 1 and 10.
         !if(ichoice == 0 .and. xland(i) < 0.1)then
         !  xf_ens(i,1:16) =0.5*(xf_ens(i,10)+xf_ens(i,1))
         !endif
      !------------------------------------
      enddo
      !-
      !-- diurnal cycle mass flux closure
      !-
      xf_dicycle(:)=0.0

      if(DICYCLE==1 .or. DICYCLE==2 )then
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
         
            !if(p_cup(kbcon(i),i)< 950. )cycle

            !--- Bechtold et al (2014)
            xff_dicycle  = (aa1(i)-aa1_bl(i))

            !--- Bechtold et al (2014) + Becker et al (2021)
            !xff_dicycle  = (1.- alpha_adv(i))* AA1(i) +  alpha_adv(i)*AA1_RADPBL(i) &
            !                  + alpha_adv(i)*Q_adv(i) - AA1_BL(i)

            xff_dicycle  = xff_dicycle /tau_ecmwf(i)

            if(xk(i) < 0.) xf_dicycle(i)= max(0.,-xff_dicycle/xk(i))

            xf_dicycle(i) = xf_ens(i,10)-xf_dicycle(i)

         enddo
      endif
      !-----------------------------------------------------
      !-
      !- add the mass flux associated to the vertical velocity
      !- at leading edge of the cold pool gust front as a surplus
      !- for the mass flux already determined.
      !-
      if(convection_tracer == 1 .and. add_coldpool_clos > 0) then 
         if(add_coldpool_clos == 1 )then
            do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               ke_gf = 0.5*wlpool(i)**2 + 1.0e-6
               w_cb = c3*sqrt(ke_gf) + c4
               xf_coldpool(i) = c1 * w_cb * exp (-c2* abs(min(cin1(i),0.))/ke_gf)
            enddo
            !if(maxval(xf_coldpool) > 0.1) print*,'xf_coldpool',maxval(xf_coldpool)

         endif
         if(add_coldpool_clos == 2 )then
            do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               if(xk(i) >= 0.) cycle
               xf_coldpool(i) = -(0.5*wlpool(i)**2/tau_ecmwf(i)) /xk(i)
            enddo
         endif
      endif 
   end subroutine cup_forcing_ens_deep
   !------------------------------------------------------------------------------------
   subroutine get_partition_liq_ice(ierr,tn,z1,zo_cup,po_cup, p_liq_ice,melting_layer  &
                                   ,itf,ktf,its,ite, kts,kte, cumulus          )
      implicit none
      character(len=*), parameter :: procedureName = 'getPartitionLiqIce' ! Subroutine Name
      character *(*), intent (in)              :: cumulus
      integer  ,intent (in   )                 :: itf,ktf, its,ite, kts,kte
      integer  ,intent (in   ), dimension(:)   :: ierr
      real     ,intent (in   ), dimension(:)   :: z1
      real     ,intent (in   ), dimension(:,:) :: tn,po_cup,zo_cup
      real     ,intent (inout), dimension(:,:) :: p_liq_ice,melting_layer
      
      !Local variables:
      integer :: i,k, vtp_index
      real    :: dp, height
      real, dimension(its:ite) :: norm
      real, parameter ::  T1=276.16, Z_meltlayer1=4000.,Z_meltlayer2=6000.,delT=3.

      p_liq_ice    (:,:) = 1.
      melting_layer(:,:) = 0.
      !
      !-- get function of T for partition of total condensate into liq and ice phases.
      if(MELT_GLAC .and. cumulus == 'deep') then
         do k=kts,ktf
            do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               p_liq_ice(k,i) = FractLiqF(tn(k,i))
            enddo
         enddo
         !        go to 650
         !
         !-- define the melting layer (the layer will be between C_t00+1 < TEMP < T_1
         !-- definition em terms of temperatura
         do k=kts,ktf
            do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               if    (tn(k,i) <= C_t00-delT) then
                  melting_layer(k,i) = 0.

               elseif(  tn(k,i) < C_t00+delT .and. tn(k,i) > C_t00-delT) then
                  melting_layer(k,i) =  ((tn(k,i)-(C_t00-delt))/(2.*delT))**2

               else
                  melting_layer(k,i) = 1.
               endif
               melting_layer(k,i) = melting_layer(k,i)*(1.-melting_layer(k,i))
            enddo
         enddo
         !go to 655
         !650 continue
         !        !-- definition em terms of height above local terrain
         !        DO k=kts,ktf
         !          DO i=its,itf
         !             if(ierr(i) /= 0) cycle
         !             height= zo_cup(k,i)+z1(i)
         !             if   (height > Z_meltlayer2 ) then
         !                melting_layer(k,i) = 0.
         !
         !             elseif(height > Z_meltlayer1  .and. height < Z_meltlayer2 ) then
         !
         !                melting_layer(k,i) =  ((height - Z_meltlayer1)/(Z_meltlayer2-Z_meltlayer1))**2.
         !
         !
         !             else
         !                melting_layer(k,i) = 1.
         !             endif
         !             melting_layer(k,i) = melting_layer(k,i)*(1.-melting_layer(k,i))
         !          ENDDO
         !        ENDDO
         !
         !
         !         655 continue
         !-normalize vertical integral of melting_layer to 1
         norm(:)=0.
         do k=kts,ktf-1
            do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               dp = 100.*(po_cup(k,i)-po_cup(k+1,i))
               norm(i) = norm(i) + melting_layer(k,i)*dp/c_grav
            enddo
         enddo
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            melting_layer(:,i)=melting_layer(:,i)/(norm(i)+1.e-6)*(100*(po_cup(kts,i)-po_cup(ktf,i))/c_grav)
          !print*,"i2=",i,maxval(melting_layer(:,i)),minval(melting_layer(:,i)),norm(i)
         enddo
      !--check
      !       norm(:)=0.
      !        DO k=kts,ktf-1
      !          DO i=its,itf
      !             dp = 100.*(po_cup(k,i)-po_cup(k+1,i))
      !             norm(i) = norm(i) + melting_layer(k,i)*dp/c_grav/(100*(po_cup(kts,i)-po_cup(ktf,i))/c_grav)
      !             !print*,"n=",i,k,norm(i)
      !          ENDDO
      !        ENDDO

      !~ ELSE
         !~ p_liq_ice    (:,:) = 1.
         !~ melting_layer(:,:) = 0.
      end if
   end  subroutine get_partition_liq_ice
   !------------------------------------------------------------------------------------
   subroutine get_melting_profile(ierr,tn_cup,po_cup, p_liq_ice,melting_layer,qrco      &
                                 ,pwo,edto,pwdo,melting,itf,ktf,its,ite, kts,kte, cumulus)
      implicit none
      character *(*), intent (in)              :: cumulus
      integer  ,intent (in   )                 :: itf,ktf, its,ite, kts,kte
      integer  ,intent (in   ), dimension(:)   :: ierr
      real     ,intent (in   ), dimension(:)   :: edto
      real     ,intent (in   ), dimension(:,:) :: tn_cup,po_cup,qrco,pwo &
                                                 ,pwdo,p_liq_ice,melting_layer
      real     ,intent (inout), dimension(:,:) :: melting
      !-local vars
      integer :: i,k,vtp_index
      real    :: dp
      real, dimension(its:ite)         :: norm,total_pwo_solid_phase
      real, dimension(kts:kte,its:ite) :: pwo_solid_phase,pwo_eff

      if(MELT_GLAC .and. cumulus == 'deep') then

         norm                  = 0.0
         pwo_solid_phase       = 0.0
         pwo_eff               = 0.0
         melting               = 0.0
         !-- set melting mixing ratio to zero for columns that do not have deep convection
         !do i=its,itf
         !   if(ierr(i) > 0) melting(:,i) = 0.
         !enddo

         !-- now, get it for columns where deep convection is activated
         total_pwo_solid_phase(:)=0.

         do k=kts,ktf-1
            do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               dp = 100.*(po_cup(k,i)-po_cup(k+1,i))

               !-- effective precip (after evaporation by downdraft)
               !-- pwdo is not defined yet
               !pwo_eff(k,i) = 0.5*(pwo(k,i)+pwo(k+1,i) + edto(i)*(pwdo(k,i)+pwdo(k+1,i)))
               pwo_eff(k,i) = 0.5*(pwo(k,i)+pwo(k+1,i))

               !-- precipitation at solid phase(ice/snow)
               pwo_solid_phase(k,i) = (1.-p_liq_ice(k,i))*pwo_eff(k,i)

               !-- integrated precip at solid phase(ice/snow)
               total_pwo_solid_phase(i) = total_pwo_solid_phase(i)+pwo_solid_phase(k,i)*dp/c_grav
            enddo
         enddo

         do k=kts,ktf
            do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               !-- melting profile (kg/kg)
               melting(k,i) = melting_layer(k,i)*(total_pwo_solid_phase(i)/(100*(po_cup(kts,i)-po_cup(ktf,i))/c_grav))
               !print*,"mel=",k,melting(k,i),pwo_solid_phase(k,i),po_cup(k,i)
            enddo
         enddo

      !-- check conservation of total solid phase precip
      !       norm(:)=0.
      !        DO k=kts,ktf-1
      !          DO i=its,itf
      !             dp = 100.*(po_cup(k,i)-po_cup(k+1,i))
      !             norm(i) = norm(i) + melting(k,i)*dp/c_grav
      !          ENDDO
      !        ENDDO
      !
      !       DO i=its,itf
      !         print*,"cons=",i,norm(i),total_pwo_solid_phase(i)
      !        ENDDO
      !--

      else
         !-- no melting allowed in this run
         melting     (:,:) = 0.
      endif
   end  subroutine get_melting_profile
   !------------------------------------------------------------------------------------
   subroutine  ke_to_heating(itf,ktf,its,ite, kts,kte,ktop,ierr,po_cup,us,vs,dellu,dellv,dellat)

      implicit none
      integer                  ,intent (in   ) :: itf,ktf,its,ite, kts,kte
      integer, dimension (:)   ,intent (in   ) :: ierr,ktop
      real   , dimension (:,:) ,intent (in   ) :: po_cup,us,vs,dellu,dellv
      real   , dimension (:,:) ,intent (inout) :: dellat

      real :: dts,fp,dp,fpi
      integer ::i,k,vtp_index

      ! since kinetic energy is being dissipated, add heating accordingly (from ECMWF)
      !
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         dts=0.
         fpi=0.
         do k=kts,ktop(i)
            dp=(po_cup(k,i)-po_cup(k+1,i))*100.
            !total KE dissiptaion estimate
            dts = dts - (dellu(k,i)*us(k,i)+dellv(k,i)*vs(k,i))*dp/c_grav
            !
            ! fpi needed for calcualtion of conversion to pot. energyintegrated
            fpi = fpi + sqrt(dellu(k,i)*dellu(k,i) + dellv(k,i)*dellv(k,i))*dp
         enddo
         if(fpi > 0.)then
            do k=kts,ktop(i)
               fp= sqrt((dellu(k,i)*dellu(k,i)+dellv(k,i)*dellv(k,i)))/fpi
               dellat(k,i) = dellat(k,i) + fp*dts*c_grav/c_cp
            enddo
         endif
      enddo

   end subroutine  ke_to_heating
   !---------------------------------------------------------------------------------------------------------------------------
   subroutine gateSoundings(its, itf, jl, kte, kts, ktf, ierr, jmin, k22, kbcon, klcl, ktop, aa1, cd, clfrac, dby, dd_massentro &
                          , dd_massdetro, dellah, dellaq, dellaqc, dtime, edto, entr_rate_2d, evap_bcb, hc, hco, he_cup, HKB &
                          , massf, massi, melting_layer, mpcf, mpqi, mpql, out_chem, outmpcf, outmpqi, outmpql, outnice, outnliq &
                          , outq, outqc, outt, outu, outv, p_liq_ice, po, po_cup, pre, prec_flx, pwdo, pwo, q_cup, q_in, qcdo, qco &
                          , qeso_cup, qo_cup, qo , qrco, qrr, sc_dn_chem, sc_up_chem, se_chem, se_cup_chem, subten_h, subten_q &
                          , subten_t , t_cup, t_in, tn, tot_pw_dn_chem, tot_pw_up_chem, up_massentro, up_massdetro, us, vvel1d &
                          , vvel2d, xmb, z_cup, zdo, zenv, zo, zqexec, zuo, zws, heso_cup, heo_cup, cumulus)
      !! ## Gate soundings
      !!
      !! Author: Rodrigues, L. F. [LFR]
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 24Fevereiro2023 09:16
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Gate soundings
      !!
      !! ** History**:
      !!
      !! ---
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!

      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_name = 'gateSoundings'
      !! subroutine name

      ! Variables (input, output, inout)
      integer, intent(in) :: its, itf, jl, kte, kts, ktf
      integer, intent(in) :: ierr(:)
      integer, intent(in) :: jmin(:)
      integer, intent(in) :: k22(:)
      integer, intent(in) :: kbcon(:)
      integer, intent(in) :: klcl(:)
      integer, intent(in) :: ktop(:)

      real, intent(in) :: aa1(:)
      real, intent(in) :: cd(:,:)
      real, intent(in) :: clfrac(:,:)
      real, intent(in) :: dby(:,:)
      real, intent(in) :: dd_massentro(:,:)
      real, intent(in) :: dd_massdetro(:,:)
      real, intent(in) :: dellah(:,:)
      real, intent(in) :: dellaq(:,:)
      real, intent(in) :: dellaqc(:,:)
      real, intent(in) :: dtime
      real, intent(in) :: edto(:)
      real, intent(in) :: entr_rate_2d(:,:)
      real, intent(in) :: evap_bcb(:,:)
      real, intent(in) :: hc(:,:)
      real, intent(in) :: hco(:,:)
      real, intent(in) :: he_cup(:,:)
      real, intent(in) :: heso_cup(:,:)
      real, intent(in) :: heo_cup(:,:)
      real, intent(in) :: HKB(:)
      real, intent(in) :: massf
      real, intent(in) :: massi
      real, intent(in) :: melting_layer(:,:)
      real, intent(in) :: mpcf(:,:,:)
      real, intent(in) :: mpqi(:,:,:)
      real, intent(in) :: mpql(:,:,:)
      real, intent(in) :: out_chem(:,:,:)
      real, intent(in) :: outmpcf(:,:,:)
      real, intent(in) :: outmpqi(:,:,:)
      real, intent(in) :: outmpql(:,:,:)
      real, intent(in) :: outnice(:,:)
      real, intent(in) :: outnliq(:,:)
      real, intent(in) :: outq(:,:)
      real, intent(in) :: outqc(:,:)
      real, intent(in) :: outt(:,:)
      real, intent(in) :: outu(:,:)
      real, intent(in) :: outv(:,:)
      real, intent(in) :: p_liq_ice(:,:)
      real, intent(in) :: po(:,:)
      real, intent(in) :: po_cup(:,:)
      real, intent(in) :: pre(:)
      real, intent(in) :: prec_flx(:,:)
      real, intent(in) :: pwdo(:,:)
      real, intent(in) :: pwo(:,:)
      real, intent(in) :: q_cup(:,:)
      real, intent(in) :: q_in(:,:)
      real, intent(in) :: qcdo(:,:)
      real, intent(in) :: qco(:,:)
      real, intent(in) :: qeso_cup(:,:)
      real, intent(in) :: qo_cup(:,:)
      real, intent(in) :: qo(:,:)
      real, intent(in) :: qrco(:,:)
      real, intent(in) :: qrr(:,:)
      real, intent(in) :: sc_dn_chem(:,:,:)
      real, intent(in) :: sc_up_chem(:,:,:)
      real, intent(in) :: se_chem(:,:,:)
      real, intent(in) :: se_cup_chem(:,:,:)
      real, intent(in) :: subten_h(:,:)
      real, intent(in) :: subten_q(:,:)
      real, intent(in) :: subten_t(:,:)
      real, intent(in) :: t_cup(:,:)
      real, intent(in) :: t_in(:,:)
      real, intent(in) :: tn(:,:)
      real, intent(in) :: tot_pw_dn_chem(:,:)
      real, intent(in) :: tot_pw_up_chem(:,:)
      real, intent(in) :: up_massentro(:,:)
      real, intent(in) :: up_massdetro(:,:)
      real, intent(in) :: us(:,:)

      real, intent(in) :: vvel1d(:)
      real, intent(in) :: vvel2d(:,:)
      real, intent(in) :: xmb(:)
      real, intent(in) :: z_cup(:,:)
      real, intent(in) :: zdo(:,:)
      real, intent(in) :: zenv(:,:)
      real, intent(in) :: zo(:,:)
      real, intent(in) :: zqexec(:)
      real, intent(in) :: zuo(:,:)
      real, intent(in) :: zws(:)

      character(len=*), intent(in) :: cumulus

      ! Local variables:
      integer :: i, k, nvar, nvar_begin, kmp
      character :: cty
      real :: dp, e_dn, c_up, trash, trash2, env_mf
      real :: resten_h, resten_q, resten_t

      if (trim(cumulus) == 'deep') then
         cty = '1'
         nvar_begin = 0
      end if
      if (trim(cumulus) == 'shallow') then
         cty = '2'
         nvar_begin = 101
      end if
      if (trim(cumulus) == 'mid') then
         cty = '3'
         nvar_begin = 201
      end if
      do i = its, itf
         !if(ierr(i).eq.0) then
         !- 2-d section
         do k = kts, ktf  !max(1,ktop(i))
            nvar = nvar_begin

            if (trim(cumulus) == 'deep') &
               call setGradsVar(jl, k, nvar, zo(k,i), "zo"//cty, ' height', '3d')
            !        call set_grads_var(jl,k,nvar,po(k,i),"po"//cty ,' press','3d')

            dp = 100.*(po_cup(k,i) - po_cup(k+1,i))
            e_dn = -0.5*(pwdo(k,i) + pwdo(k+1,i))*c_grav/dp*edto(i)*86400.*real(c_alvl)/real(c_cp) &
                 * xmb(i) ! pwdo < 0 and E_dn must > 0
            c_up = dellaqc(k,i) + (zuo(k+1,i)*qrco(k+1,i) - zuo(k,i)*qrco(k,i)) &
                 * c_grav/dp + 0.5*(pwo(k,i) + pwo(k+1,i))*c_grav/dp
            c_up = -c_up*86400.*real(c_alvl)/real(c_cp)*xmb(i)

            trash = -(zuo(k+1,i)*(qco(k+1,i) - qo_cup(k+1,i)) - zuo(k,i) &
                  * (qco(k,i) - qo_cup(k,i)))*c_grav/dp
            trash2 = +(zdo(k+1,i)*(qcdo(k+1,i) - qo_cup(k+1,i)) - zdo(k,i) &
                   * (qcdo(k,i) - qo_cup(k,i)))*c_grav/dp &
                     *edto(i)

            trash = trash*86400.*real(c_alvl)/real(c_cp)*xmb(i)
            trash2 = trash2*86400.*real(c_alvl)/real(c_cp)*xmb(i)

            env_mf = 0.5*((zuo(k+1,i) - zdo(k+1,i)*edto(i)) + (zuo(k,i) - zdo(k,i) &
                   * edto(i)))
            resten_h = dellah(k,i) - subten_h(k,i)
            resten_q = dellaq(k,i) - subten_q(k,i)
            resten_t = (1./real(c_cp))*(resten_h - real(c_alvl)*resten_q)
            !trash2 = qco   (k,i  )! zuo(k+1,i)*(qco (k+1,i)-qo_cup(k+1,i) ) !*g/dp
            !trash  = qo_cup(k,i  )! zuo(k,i  )*(qco (k,i  )-qo_cup(k,i  ) ) !*g/dp
            trash2 = zuo(k+1,i)*(qco(k+1,i) - qo_cup(k+1,i))*1000 !*g/dp
            trash = zuo(k,i)*(qco(k,i) - qo_cup(k,i))*1000  !*g/dp

            call setGradsVar(jl, k, nvar, out_chem(1,k,i)*86400, "outchem"//cty, ' outchem', '3d')
            call setGradsVar(jl, k, nvar, sc_up_chem(1,k,i), "scup"//cty, ' sc_chem', '3d')
            call setGradsVar(jl, k, nvar, sc_dn_chem(1,k,i), "scdn"//cty, ' sc_chem', '3d')
            call setGradsVar(jl, k, nvar, massi, "mi"//cty, ' initial mass', '2d')
            call setGradsVar(jl, k, nvar, massf, "mf"//cty, ' final mass', '2d')
            call setGradsVar(jl, k, nvar, se_chem(1,k,i), "se"//cty, ' se_chem', '3d')
            call setGradsVar(jl, k, nvar, se_cup_chem(1,k,i), "secup"//cty, ' se_cup_chem', '3d')
            !-- only for debug
            !call set_grads_var(jl,k,nvar,se_chem_update(1,i,k),"newse"//cty ,' new se_chem','3d')
            ! if (APPLY_SUB_MP == 1) then
            !    kmp = p_lsmp
            !    call setGradsVar(jl, k, nvar, outmpqi(kmp, i, k)*86400*1000, "outqi"//cty, ' outmpqi', '3d')
            !    call setGradsVar(jl, k, nvar, outmpql(kmp, i, k)*86400*1000, "outql"//cty, ' outmpql', '3d')
            !    call setGradsVar(jl, k, nvar, outmpcf(kmp, i, k)*86400, "outcf"//cty, ' outmpcf', '3d')
            !    call setGradsVar(jl, k, nvar, mpqi(kmp, i, k), "mpqi"//cty, ' mpqi', '3d')
            !    call setGradsVar(jl, k, nvar, mpql(kmp, i, k), "mpql"//cty, ' mpql', '3d')
            !    call setGradsVar(jl, k, nvar, mpcf(kmp, i, k), "mpcf"//cty, ' mpcf', '3d')
            ! end if
            call setGradsVar(jl, k, nvar, env_mf, "sub"//cty, ' sub', '3d')
            if (LIQ_ICE_NUMBER_CONC == 1) then
               call setGradsVar(jl, k, nvar, outnice(k,i)*86400., "outnice"//cty, 'out # ice1/day', '3d')
               call setGradsVar(jl, k, nvar, outnliq(k,i)*86400., "outnliq"//cty, 'out # liq /day', '3d')
            end if
            call setGradsVar(jl, k, nvar, zuo(k,i)/(1.e-12+xmb(i)), "zup"//cty, 'norm m flux up ', '3d')
            call setGradsVar(jl, k, nvar, zdo(k,i)/(1.e-12+xmb(i)), "zdn"//cty, 'norm m flux dn ', '3d')
            call setGradsVar(jl, k, nvar, zenv(k,i)/(1.e-12+xmb(i)), "zenv"//cty, 'norm m flux env ', '3d')
!           call setGradsVar(jl, k, nvar, -edto(i)*xmb(i)*zdo(k,i), "mdn"//cty, 'm flux down (kg/s/m^2)', '3d')
            call setGradsVar(jl, k, nvar, -zdo(k,i), "mdn"//cty, 'm flux down (kg/s/m^2)', '3d')
            call setGradsVar(jl, k, nvar, up_massentro(k,i), "upent"//cty, 'up_massentr(kg/s/m^2)', '3d')
            call setGradsVar(jl, k, nvar, up_massdetro(k,i)/(1.e-12+xmb(i)), "updet"//cty, 'up_massdetr(kg/s/m^2)', '3d')
            call setGradsVar(jl, k, nvar, outt(k,i)*86400., "outt"//cty, 'outt K/day', '3d')
            call setGradsVar(jl, k, nvar, resten_t*86400., "rest"//cty, 'residuo T K/day', '3d')
            call setGradsVar(jl, k, nvar, resten_h*86400./real(c_cp), "resh"//cty, 'residuo H J/kg/day', '3d')
            call setGradsVar(jl, k, nvar, resten_q*86400.*real(c_alvl)/real(c_cp), "resq"//cty, 'residuo q K/day   ', '3d')
            call setGradsVar(jl, k, nvar, subten_t(k,i)*86400., "subt"//cty, 'subT K/day', '3d')
            call setGradsVar(jl, k, nvar, subten_h(k,i)*86400./real(c_cp), "subh"//cty, 'subH J/kg/day', '3d')
            call setGradsVar(jl, k, nvar, subten_q(k,i)*86400.*real(c_alvl)/real(c_cp), "subq"//cty &
                          , 'subq K/day   ', '3d')
            call setGradsVar(jl, k, nvar, outq(k,i)*86400.*real(c_alvl)/real(c_cp), "outq"//cty, 'outq K/s', '3d')
            call setGradsVar(jl, k, nvar, outqc(k,i)*86400.*real(c_alvl)/real(c_cp), "outqc"//cty, 'outqc K/day', '3d')
            call setGradsVar(jl, k, nvar, pre(i)*3600., "precip"//cty, 'precip mm', '2d')
            call setGradsVar(jl, k, nvar, prec_flx(k,i)*3600., "precflx"//cty, 'prec flx mm', '3d')
            call setGradsVar(jl, k, nvar, pwo(k,i), "pwo"//cty, ' xx', '3d')
            call setGradsVar(jl, k, nvar, outu(k,i)*86400., "outu"//cty, 'out_U m/s/day', '3d')
            call setGradsVar(jl, k, nvar, outv(k,i)*86400., "outv"//cty, 'out_V m/s/day', '3d')
            call setGradsVar(jl, k, nvar, xmb(i), "xmb"//cty, 'xmb kg/m2/s', '2d')
            call setGradsVar(jl, k, nvar, vvel2d(k,i), "W2d"//cty, 'W /m/s', '3d')
            call setGradsVar(jl, k, nvar, vvel1d(i), "W1d"//cty, 'W1s /m/s', '2d')
            call setGradsVar(jl, k, nvar, us(k,i), "us"//cty, 'U /m/s', '3d')
            call setGradsVar(jl, k, nvar, outu(k,i)*86400./(1.e-16 + xmb(i)), "delu"//cty, 'dellu', '3d')
            call setGradsVar(jl, k, nvar, evap_bcb(k,i)*1000., "evcb"//cty, 'g/kg', '3d')

            call setGradsVar(jl, k, nvar, tot_pw_up_chem(1, i), "pwup"//cty, 'pwup', '2d')
            call setGradsVar(jl, k, nvar, tot_pw_dn_chem(1, i), "pwdn"//cty, 'pwdn', '2d')
            !----
            !----
            call setGradsVar(jl, k, nvar, xmb(i)*dellah(k,i)*86400./real(c_cp), "delh"//cty, 'dellah K/day', '3d')
            call setGradsVar(jl, k, nvar, xmb(i)*dellaq(k,i)*86400.*real(c_alvl)/real(c_cp), "dellq"//cty &
                         , 'dellaq K/day', '3d')
            call setGradsVar(jl, k, nvar, xmb(i)*dellaqc(k,i)*86400.*real(c_alvl)/real(c_cp), "dellqc"//cty &
                          , 'dellaqc K/day' &
                             , '3d')
            call setGradsVar(jl, k, nvar, xmb(i), "xmb"//cty, 'm flux up (kg/s/m^2)', '2d')
            call setGradsVar(jl, k, nvar, aa1(i), "aa1"//cty, 'AA1 J/kg3)', '2d')
            call setGradsVar(jl, k, nvar, float(ierr(i)), "ierr"//cty, 'ierr #', '2d')
            call setGradsVar(jl, k, nvar, xmb(i)*dd_massentro(k,i), "ddent"//cty, 'dd_massentr(kg/s/m^2)', '3d')
            call setGradsVar(jl, k, nvar, xmb(i)*dd_massdetro(k,i), "dddet"//cty, 'dd_massdetr(kg/s/m^2)', '3d')
                  !!     cycle
            call setGradsVar(jl, k, nvar, hc(k,i)/c_cp, "hc"//cty, ' hc', '3d')
            call setGradsVar(jl, k, nvar, hco(k,i)/c_cp, "hco"//cty, ' hco', '3d')
            call setGradsVar(jl, k, nvar, heso_cup(k,i)/c_cp, "heso"//cty, ' heso_cup', '3d')
            call setGradsVar(jl, k, nvar, heo_cup(k,i)/c_cp, "heo"//cty, ' heo_cup', '3d')
            call setGradsVar(jl, k, nvar, dby(k,i), "dby"//cty, ' dbuo', '3d')
            !call set_grads_var(jl,k,nvar,QCUP(k,i),"qcup"//cty ,'C_UP','3d')
            call setGradsVar(jl, k, nvar, t_cup(k,i) - 273.15, "te"//cty, ' K', '3d')
            call setGradsVar(jl, k, nvar, 1000.*q_cup(k,i), "qe"//cty, ' kg kg-1', '3d')
            call setGradsVar(jl, k, nvar, he_cup(k,i), "he"//cty, ' he', '3d')
            call setGradsVar(jl, k, nvar, HKB(i), "hkb"//cty, ' H', '2d')
            call setGradsVar(jl, k, nvar, HKB(i), "hkb"//cty, ' H', '2d')
            call setGradsVar(jl, k, nvar, 1000.*zqexec(i), "qex"//cty, ' qex', '2d')
            call setGradsVar(jl, k, nvar, z_cup(max(1, k22(i))  ,i), "zs"//cty, ' m', '2d')
            call setGradsVar(jl, k, nvar, z_cup(max(1, kbcon(i)),i), "zbcon"//cty, ' m', '2d')
            call setGradsVar(jl, k, nvar, z_cup(max(1, ktop(i)) ,i), "ztop"//cty, ' m', '2d')
            call setGradsVar(jl, k, nvar, z_cup(max(1, klcl(i)) ,i), "zlcl"//cty, ' m', '2d')
            call setGradsVar(jl, k, nvar, z_cup(max(1, jmin(i)) ,i), "zjmin"//cty, ' m', '2d')
            call setGradsVar(jl, k, nvar, zws(i), "ws"//cty, ' m/s', '2d')
            call setGradsVar(jl, k, nvar, clfrac(k,i), "clfrac"//cty, 'shcf #', '3d')
            call setGradsVar(jl, k, nvar, entr_rate_2d(k,i)*1000., "entr"//cty, ' km-1', '3d')
            call setGradsVar(jl, k, nvar, cd(k,i)*1000., "detr"//cty, ' km-1', '3d')
            call setGradsVar(jl, k, nvar, pwdo(k,i), "pwd"//cty, ' xx', '3d')
            call setGradsVar(jl, k, nvar, edto(i), "edt"//cty, 'edt kg/m2/s', '2d')
            call setGradsVar(jl, k, nvar, e_dn, "EVAP"//cty, ' xx', '3d')
            call setGradsVar(jl, k, nvar, c_up, "CUP"//cty, ' xx', '3d')
            !       call set_grads_var(jl,k,nvar,trash,"TUP"//cty ,' xx','3d')
            !       call set_grads_var(jl,k,nvar,trash2,"TDN"//cty ,' xx','3d')
            call setGradsVar(jl, k, nvar, trash, "F1"//cty, ' F1', '3d')
            call setGradsVar(jl, k, nvar, trash2, "F2"//cty, ' F2', '3d')
            call setGradsVar(jl, k, nvar, p_liq_ice(k,i), "pli"//cty, '#', '3d')
            call setGradsVar(jl, k, nvar, melting_layer(k,i), "cpli"//cty, '#', '3d')
            call setGradsVar(jl, k, nvar, t_in(k,i), "t"//cty, 'temp K', '3d')
            call setGradsVar(jl, k, nvar, tn(k,i), "tn"//cty, 'temp K', '3d')
            call setGradsVar(jl, k, nvar, 1000.*q_in(k,i), "q"//cty, 'q g/kg', '3d')
            call setGradsVar(jl, k, nvar, 1000.*qo(k,i), "qn"//cty, 'q g/kg', '3d')
            call setGradsVar(jl, k, nvar, 1000.*qrco(k,i), "qrc"//cty, 'q g/kg', '3d')
            call setGradsVar(jl, k, nvar, 1000.*(q_in(k,i) + outq(k,i)*dtime), "qnc"//cty, 'q upd conv g/kg' &
                           , '3d')
            call setGradsVar(jl, k, nvar, 1000.*(qo(k,i) + outq(k,i)*dtime), "qnall"//cty, 'q upd all g/kg' &
                          ,  '3d')
            call setGradsVar(jl, k, nvar, 1000.*qrr(k,i), "qrr"//cty, 'qrr g/kg', '3d')
            call setGradsVar(jl, k, nvar, 1000.*qco(k,i), "qc"//cty, 'qc g/kg', '3d')
            call setGradsVar(jl, k, nvar, 1000.*qo_cup(k,i), "qcup"//cty, 'qc g/kg', '3d')
            call setGradsVar(jl, k, nvar, 1000.*qeso_cup(k,i), "qescup"//cty, 'qc g/kg', '3d')

            !~ call set_grads_var(jl,k,nvar,aa0(i),"a0"//cty,'aa0','2d')
            !~ call set_grads_var(jl,k,nvar,aa1_fa(i),"aa1fa"//cty,'aa1fa','2d')
            !~ call set_grads_var(jl,k,nvar,aa1_bl(i),"aa1bl"//cty,'aa1bl','2d')
            !~ call set_grads_var(jl,k,nvar,aa0_bl(i),"aa0bl"//cty,'aa0bl','2d')
            !~ call set_grads_var(jl,k,nvar,aa1(i),"a1"//cty,'aa1','2d')
            !~ call set_grads_var(jl,k,nvar,aa1(i)/(1.e-6+tau_ecmwf(i)),"mb13"//cty,'aa0','2d')
            !~ call set_grads_var(jl,k,nvar,xaa0(i),"xa0"//cty,'xaa0','2d')
            !~ call set_grads_var(jl,k,nvar,(XAA0(I)-AA1(I))/MBDT(I),"xk"//cty,'xk','2d')
         end do
         if (wrtgrads .and. .not. p_use_gate) then
            call wrtBinCtl(1, kte, po(1, 1:kte), cumulus)
         end if
      end do

   end subroutine gateSoundings
  ! ---------------------------------------------------------------------------------------------------
   subroutine setGradsVar(i_in, k_in, nvar, f_in, name1, name2, name3)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: procedureName = 'setGradsVar' ! Subroutine Name
   
      !Variables (input, output, inout)
      integer, intent(in)    :: i_in, k_in
      
      real, intent(in) :: f_in
      
      character(len=*), intent(in) :: name1
      character(len=*), intent(in) :: name2
      character(len=*), intent(in) :: name3

      integer, intent(inout) :: nvar

      cupout(nvar)%varp(i_in, k_in) = f_in
      cupout(nvar)%varn(1) = name1
      cupout(nvar)%varn(2) = name2
      cupout(nvar)%varn(3) = name3
      nvar = nvar + 1
      if (nvar > p_nvar_grads) stop 'nvar>nvar_grads'

   end subroutine setGradsVar
   ! ------------------------------------------------------------------------------------
   subroutine wrtBinCtl(n, mzp, p2d, cumulus)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: procedureName = 'wrtBinCtl' ! Subroutine Name

      real, parameter :: p_undef = -9.99e33
   
      !Variables (input, output, inout)
      integer, intent(in):: n, mzp

      real, intent(in) :: p2d(:)

      character(len=*), intent(in) :: cumulus
      
      !Local variables:
      integer:: nvartotal, klevgrads(200), jk, int_byte_size, nvar, maxklevgrads
      real   :: real_byte_size
      integer :: nrec, rec_size

      nrec = 0
      maxklevgrads = min(60, mzp)
      runname = '15geos5_'//cumulus
      runlabel = runname

      print *, "writing grads control file:',trim(runname)//'.ctl", ntimes
      call flush (6)
      !
      !number of variables to be written
      nvartotal = 0
      do nvar = 1, p_nvar_grads
         if (cupout(nvar)%varn(1) .ne. "xxxx") nvartotal = nvartotal + 1
         if (cupout(nvar)%varn(3) == "3d") klevgrads(nvar) = maxklevgrads
         if (cupout(nvar)%varn(3) == "2d") klevgrads(nvar) = 1
      end do

      !- binary file
      inquire (iolength=int_byte_size) real_byte_size  ! inquire by output list

      print *, 'opening grads file:', trim(runname)//'.gra'
      rec_size = size(cupout(nvar)%varp, 1)*real_byte_size
      if (ntimes == 1) then
         open (19, file=trim(runname)//'.gra', form='unformatted', &
               access='direct', status='replace', recl=rec_size)
      else
         open (19, file=trim(runname)//'.gra', form='unformatted', &
               access='direct', status='old', recl=rec_size)
      end if

      do nvar = 1, p_nvar_grads
         if (cupout(nvar)%varn(1) .ne. "xxxx") then
            do jk = 1, klevgrads(nvar)
               nrec = nrec + 1
               !write(19)          real((cupout(nvar)%varp(:,jk)),4)
               write (19, rec=nrec) real((cupout(nvar)%varp(:, jk)), 4)
            end do
         end if
      end do
      close (19)
      !-setting vertical dimension '0' for 2d var
      where (klevgrads == 1) klevgrads = 0
      !- ctl file
      open (20, file=trim(runname)//'.ctl', status='unknown')
      write (20, 2001) '^'//trim(runname)//'.gra'
      write (20, 2002) 'undef -9.99e33'
      write (20, 2002) 'options sequential byteswapped' ! zrev'
      write (20, 2002) 'title '//trim(runlabel)
      write (20, 2003) 1, 0., 1. ! units m/km
      write (20, 2004) n, 1., 1.
      write (20, 2005) maxklevgrads, (p2d(jk), jk=1, maxklevgrads)
      write (20, 2006) ntimes, '00:00Z24JAN1999', '10mn'
      write (20, 2007) nvartotal
      do nvar = 1, p_nvar_grads
         if (cupout(nvar)%varn(1) .ne. "xxxx") then
            !
            write (20, 2008) cupout(nvar)%varn(1) (1:len_trim(cupout(nvar)%varn(1))), klevgrads(nvar) &
               , cupout(nvar)%varn(2) (1:len_trim(cupout(nvar)%varn(2)))
         end if
      end do
      write (20, 2002) 'endvars'
      close (20)

2001  format('dset ', a)
2002  format(a)
2003  format('xdef ', i4, ' linear ', 2f15.3)
2004  format('ydef ', i4, ' linear ', 2f15.3)
2005  format('zdef ', i4, ' levels ', 60f8.3)
2006  format('tdef ', i4, ' linear ', 2a15)
2007  format('vars ', i4)
2008  format(a10, i4, ' 99 ', a40)!'[',a8,']')
!2055  format(60f7.0)
!133   format(1x, F7.0)

   end subroutine wrtBinCtl
   !---------------------------------------------------------------------------------------------------

   subroutine rain_evap_below_cloudbase(cumulus,itf,ktf, its,ite, kts,kte,ierr,kbcon,ktop     &
                                       ,xmb,psur,xland,qo_cup,t_cup,po_cup,qes_cup,pwavo,edto &
                                       ,pwevo,pwo,pwdo,pre,prec_flx,evap_flx,outt,outq,outbuoy,evap_bcb)

      implicit none
      real, parameter :: alpha1=5.44e-4 & !1/sec
                        ,alpha2=5.09e-3 & !unitless
                        ,alpha3=0.5777  & !unitless
                        ,c_conv=0.05      !conv fraction area, unitless

      character*(*)           ,intent(in)    :: cumulus
      integer                 ,intent(in)    :: itf,ktf, its,ite, kts,kte
      integer, dimension(:)   ,intent(in)    :: ierr,kbcon,ktop
      real,    dimension(:)   ,intent(in)    :: psur,xland,pwavo,edto,pwevo,xmb
      real,    dimension(:,:) ,intent(in)    :: po_cup,qo_cup,qes_cup,pwo,pwdo,t_cup
      real,    dimension(:)   ,intent(inout) :: pre
      real,    dimension(:,:) ,intent(inout) :: outt,outq,outbuoy,prec_flx,evap_flx
      real,    dimension(:,:) ,intent(out)   :: evap_bcb

      !-- locals
      integer :: i,k,vtp_index
      real    :: RH_cr , del_t,del_q,dp,q_deficit, pqsat, temp_pre
      real    :: RH_cr_OCEAN,RH_cr_LAND
      real,    dimension(its:ite) :: tot_evap_bcb,eff_c_conv

      if(cumulus == 'shallow') then
         RH_cr_OCEAN   = 1.
         RH_cr_LAND    = 1.
         eff_c_conv(:) = min(0.2,max(xmb(:),c_conv))
      else
         RH_cr_OCEAN   = 0.99 !test 0.90
         RH_cr_LAND    = 0.99 
         eff_c_conv(:) = c_conv
      endif

      prec_flx     = 0.0
      evap_flx     = 0.0
      tot_evap_bcb = 0.0
      if(c0 < 1.e-6 ) return

      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)

         !-- critical rel humidity  - check this, if the value is too small, not evapo will take place.
         RH_cr=RH_cr_OCEAN*xland(i)+RH_cr_LAND*(1.0-xland(i))

         !if(xland(i)  < 0.90 ) then !- over land
         !  RH_cr = RH_cr_LAND
         !else
         !  RH_cr = RH_cr_OCEAN
         !endif

         do k=ktop(i),kts,-1

            dp = 100.*(po_cup(k,i)-po_cup(k+1,i))

            !p_liq_ice(k,i) = FractLiqF(tempco(k,i))

            !---rainfall evaporation below cloud base
            if(k <= kbcon(i)) then
               q_deficit = max(0.,(RH_cr*qes_cup(k,i) -qo_cup(k,i)))
               !pqsat=satur_spec_hum(t_cup(k,i),po_cup(k,i))

               !--units here: kg[water]/kg[air}/sec
               evap_bcb(k,i) = eff_c_conv(i) * alpha1 * q_deficit * &
                  ( sqrt(po_cup(k,i)/psur(i))/alpha2 * prec_flx(k+1,i)/eff_c_conv(i) )**alpha3

               !--units here: kg[water]/kg[air}/sec * kg[air]/m3 * m = kg[water]/m2/sec
               evap_bcb(k,i)= evap_bcb(k,i)*dp/c_grav

            else

               evap_bcb(k,i)=0.0

            endif

            !-- before anything check if the evaporation already consumed all precipitation
            temp_pre = pre(i) - evap_bcb(k,i)
            if (temp_pre < 0.)  evap_bcb(k,i) =  pre(i)

            !-- get the net precitation flux after the local evaporation and downdraft
            prec_flx(k,i) = prec_flx(k+1,i) - evap_bcb(k,i) + xmb(i)*(pwo(k,i) + edto(i)*pwdo(k,i))
            prec_flx(k,i) = max(0.,prec_flx(k,i))

            evap_flx(k,i) = evap_flx(k+1,i) + evap_bcb(k,i) - xmb(i)*edto(i)*pwdo(k,i)
            evap_flx(k,i) = max(0.,evap_flx(k,i))

            tot_evap_bcb(i) = tot_evap_bcb(i)+evap_bcb(k,i)

            !-- feedback
            del_q =  evap_bcb(k,i)*c_grav/dp          ! > 0., units: kg[water]/kg[air}/sec
            del_t = -evap_bcb(k,i)*c_grav/dp*(c_alvl/c_cp) ! < 0., units: K/sec

            outq   (k,i) = outq   (k,i) + del_q
            outt   (k,i) = outt   (k,i) + del_t
            !--- comment out 17nov
            !outbuoy(k,i) = outbuoy(k,i) + c_cp*del_t+c_alvl*del_q

            pre(i) = pre(i) - evap_bcb(k,i)


            !--for future use (rain and snow precipitation fluxes)
            !prec_flx_rain(k) = prec_flx(k,i)*(1.-p_liq_ice(k))
            !prec_flx_snow(k) = prec_flx(k,i)*    p_liq_ice(k)
         enddo

         if(pre(i)<0.) then
            print*,"prec evap neg for cumulus=",pre(i),trim(cumulus)
            call flush(6)
            !stop '@subroutine rain_evap_below_cloudbase'
         endif

      enddo

   end subroutine rain_evap_below_cloudbase
   !---------------------------------------------------------------------------------------------------
   subroutine get_precip_fluxes(cumulus,klcl,kbcon,ktop,k22,ierr,xland,pre,xmb             &
                               ,pwo,pwavo,edto,pwevo,pwdo,t_cup,tempco,prec_flx,evap_flx   &
                               ,itf,ktf,its,ite, kts,kte)

      implicit none
      character *(*)         ,intent (in) :: cumulus
      integer                ,intent (in) :: itf,ktf,its,ite,kts,kte
      integer, dimension(:)  ,intent (in) :: kbcon,ktop,k22,klcl,ierr
      real,    dimension(:)  ,intent (in) :: xland,pwavo,pwevo,edto,pre,xmb
      real,    dimension(:,:),intent (in) :: pwo,pwdo,t_cup,tempco
      real,    dimension(:,:),intent (out):: prec_flx,evap_flx !-- units kg[water]/m2/s

      !-- locals
      integer :: i,k,vtp_index
      prec_flx = 0.0
      evap_flx = 0.0
      if(c0 < 1.e-6 ) return

      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)

         do k=ktop(i),kts,-1

            !--- precipitation flux (at 'cup' levels), units: kg[water]/m2/s
            prec_flx(k,i) = prec_flx(k+1,i) + xmb(i)*(pwo(k,i) + edto(i)*pwdo(k,i))
            prec_flx(k,i) = max(0.,prec_flx(k,i))

            !--- evaporation flux (at 'cup' levels), units: kg[water]/m2/s
            evap_flx(k,i) = evap_flx(k+1,i) - xmb(i)*edto(i)*pwdo(k,i)
            evap_flx(k,i) = max(0.,evap_flx(k,i))

            !
            !--for future use (rain and snow precipitation fluxes)
            !p_liq_ice(k,i) = FractLiqF(tempco(k,i))
            !prec_flx_rain(k) = prec_flx(k,i)*(1.-p_liq_ice(k))
            !prec_flx_snow(k) = prec_flx(k,i)*    p_liq_ice(k)
         enddo

            !if(prec_flx   (kts,i) .ne. pre(i)) then
            !print*,"error=",100.*(prec_flx   (kts,i) - pre(i))/(1.e-16+pre(i)),pre(i),prec_flx   (kts,i)
            !STOP 'problem with water balance'
            !endif
      enddo
   end   subroutine get_precip_fluxes
   !------------------------------------------------------------------------------------
   real function satur_spec_hum(pt,press) result(pqsat)
      implicit none
      real   , intent(in ) :: pt,press ! Kelvin, hPa
      !---locals
      real :: zew,zqs,zcor,foealfcu,foeewmcu
      real, parameter ::           &
          c_rd=287.06              &
         ,c_rv=461.52              &
         ,c_rtt=273.16             &
         ,c_retv=c_rv/c_rd-1.0     &
         ,c_r2es=611.21*c_rd/c_rv  &
         ,c_r3les=17.502           &
         ,c_r3ies=22.587           &
         ,c_r4les=32.19            &
         ,c_r4ies=-0.7             &
         ,c_rtwat=c_rtt            &
         ,c_rtice=c_rtt-23.        &
         ,c_rticecu=c_rtt-23.      &
         ,c_rtwat_rtice_r=1./(c_rtwat-c_rtice)   &
         ,c_rtwat_rticecu_r=1./(c_rtwat-c_rticecu)

      foealfcu = min(1.0,((max(c_rticecu,min(c_rtwat,pt))-c_rticecu)*c_rtwat_rticecu_r)**2)
      foeewmcu = c_r2es *(    foealfcu *exp(c_r3les*(pt-c_rtt)/(pt-c_r4les))+ &
                         (1.0-foealfcu)*exp(c_r3ies*(pt-c_rtt)/(pt-c_r4ies)))

      zew  = foeewmcu
      zqs  = zew/(100.*press)
      if(1.0-c_retv*zqs > 0. )then
         zcor = 1.0/(1.0-c_retv*zqs)  ! divide by zero
         pqsat= zqs*zcor
      else
         pqsat= p_max_qsat
      endif

   end function satur_spec_hum
   !---------------------------------------------------------------------------------------------------
   subroutine get_jmin(cumulus,itf,ktf,its,ite, kts,kte,ierr,kdet,ktop,kbcon,jmin,ierrc  &
                      ,beta,depth_min,heso_cup,zo_cup,melting_layer)

      implicit none
      character *(*)             ,intent (in)    :: cumulus
      real                       ,intent (in)    :: depth_min
      integer                    ,intent (in)    :: itf,ktf,its,ite,kts,kte
      integer, dimension(:)      ,intent (in)    :: ktop,kbcon
      real,    dimension(:,:)    ,intent (in)    :: heso_cup,zo_cup,melting_layer

      integer, dimension(:)      ,intent (inout) :: ierr,jmin,kdet
      real                       ,intent (out)   :: beta
      character*128,              intent (out)   :: ierrc(:)

      !-- local vars
      integer :: i,k,jmini,ki,vtp_index
      real    :: dh,dz
      real,    dimension(kts:kte,its:ite)  ::  hcdo
      logical :: keep_going

      if(cumulus == 'shallow'  ) then
         beta    = 0.02
         jmin(:) = 0
         return
      endif

      if(cumulus == 'deep') beta=0.05
      if(cumulus == 'mid' ) beta=0.02


      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)

         if(cumulus == 'deep' .and. melt_glac) jmin(i)=max(jmin(i),maxloc(melting_layer(:,i),1))

         !
         !--- check whether it would have buoyancy, if there where
         !--- no entrainment/detrainment
         !
         jmini = jmin(i)
         keep_going = .true.
         do while ( keep_going )
            keep_going = .false.
            if ( jmini - 1 .lt. kdet(i)   ) kdet(i) = jmini-1
            if ( jmini     .ge. ktop(i)-1 ) jmini = ktop(i) - 2
            ki = jmini
            hcdo(ki,i)=heso_cup(ki,i)
            dz=zo_cup(ki+1,i)-zo_cup(ki,i)
            dh=0.
            do k=ki-1,1,-1
               hcdo(k,i)=heso_cup(jmini,i)
               dz=zo_cup(k+1,i)-zo_cup(k,i)
               dh=dh+dz*(hcdo(k,i)-heso_cup(k,i))

               if(dh > 0.)then
                  jmini=jmini-1
                  if ( jmini .gt. 5 ) then
                     keep_going = .true.
                  else
                     ierr(i)=9 ; is_removed = remove(vec_ok, i)
                     ierrc(i) = "could not find jmini9"
                     exit
                  endif
               endif
            enddo
         enddo
         jmin(i) = jmini
         if ( jmini .le. 5 ) then
            ierr(i)=4 ; is_removed = remove(vec_ok, i)
            ierrc(i) = "could not find jmini4"
         endif
      enddo
      !
      ! - must have at least depth_min m between cloud convective base and cloud top.
      !
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         if ( jmin(i) - 1 < kdet(i)) kdet(i) = jmin(i)-1
         if (-zo_cup(kbcon(i),i)+zo_cup(ktop(i),i) < depth_min)then
            ierr(i)=6 ; is_removed = remove(vec_ok, i)
            ierrc(i)="cloud depth very shallow"
         endif
      enddo

   end subroutine get_jmin
   !------------------------------------------------------------------------------------
   subroutine precip_cwv_factor(itf,ktf,its,ite,kts,kte,ierr,t,po,qo,po_cup,cumulus,p_cwv_ave)
      implicit none
      !Parameters:
      character(len=*), parameter :: procedureName = 'precipCwvFactor' ! Subroutine Name

      character *(*), intent (in)             :: cumulus
      integer  ,intent (in )                  :: itf,ktf, its,ite, kts,kte
      integer  ,intent (in ), dimension(:)    :: ierr
      real     ,intent (in ), dimension(:,:)  :: t,po,qo,po_cup
      real     ,intent (out), dimension(:)    :: p_cwv_ave

      !--locals
      integer :: i,k, vtp_index
      real    :: dp, trash
      real, dimension(its:ite) :: w_col,w_ccrit,t_troposph
      real, parameter :: fpkup=0.8  !-- 90% of precip occurs above 80% of critical w

      p_cwv_ave = 0.0
      if(cumulus /= 'deep') return
      !
      !-- get the pickup of ensemble ave prec, following Neelin et al 2009.
      !
      w_col     (its:itf) = 0.       
      w_ccrit   (its:itf) = 0.    
      t_troposph(its:itf) = 0. 

      do vtp_index = 1, get_num_elements(vec_ok); i = get_data_value(vec_ok,vtp_index)
         trash=0.
         loopN:    do k=kts,ktf
            if(po(k,i) .lt. 200.) exit loopN

            dp=100.*(po_cup(k,i)-po_cup(k+1,i))
            trash=trash+dp/c_grav

            w_col     (i)= w_col     (i) + qo(k,i)*dp/c_grav ! unit mm
            t_troposph(i)= t_troposph(i) + t (k,i)*dp/c_grav
         enddo loopN
         !--average temperature
         t_troposph(i) =   t_troposph(i)/(1.e-8+trash)! unit K
         !
         !--- wcrit given by Neelin et al 2009.
         w_ccrit(i) = max(0.,56.2 + 2.1*(t_troposph(i)-268.)) ! unit mm
         !
         !--- pickup (normalized by the factor 'a')
         !-- <p>=a[(w-w_c)/w_c]**beta, a=0.15, beta=0.23
         !
         p_cwv_ave(i) = (max(0.,w_col(i)-fpkup*w_ccrit(i))/(1.e-8+fpkup*w_ccrit(i)))**0.23
         p_cwv_ave(i) = max(0., min(1., p_cwv_ave(i)))

           !print*,"NEE=",i,w_col(i),t_troposph(i),w_ccrit(i),p_cwv_ave    (i)
           !print*,"=================================================="
      enddo
   end subroutine precip_cwv_factor
   !------------------------------------------------------------------------------------
   subroutine get_wetbulb(jmin,qo_cup,t_cup,po_cup ,q_wetbulb,t_wetbulb)

      implicit none
      integer ,intent (in   ) :: jmin
      real    ,intent (in   ) :: qo_cup,t_cup,po_cup
      real    ,intent (inout) :: q_wetbulb,t_wetbulb

      !---locals
      real ::  zqp, zcond, zcond1, zcor, zqsat
      real :: psp, pt , pq
      real :: z3es,   z4es, z5alcp, zaldcp
      real :: ptare, evap
      real :: foedelta,foeewmcu,foealfcu,foedemcu,foeldcpmcu

      !-- for testing
      !              PSP                   TEMP        Q                     ZCOND1
      ! input   85090.0000000000        289.140030372766     1.105078557441815E-002
      ! output  85090.0000000000        287.230570412846     1.181792062536557E-002 -2.761256206705639E-005
      ! PT  = 289.140030372766
      ! PQ  = 1.105078557441815E-002
      ! PSP = 85090.
      !----------------------

      !-- environmental values
      PT  = t_cup       ! K
      PQ  = qo_cup      ! kg/kg
      psp = po_cup*100. ! hpa

      if (pt > c_rtt) then
         z3es=c_r3les
         z4es=c_r4les
         z5alcp=c_r5alvcp
         zaldcp=c_ralvdcp
      else
         z3es=c_r3ies
         z4es=c_r4ies
         z5alcp=c_r5alscp
         zaldcp=c_ralsdcp
      endif

      !--- get wet bulb thermo properties --------------------------
      ptare = pt
      zqp    =1.0/psp

      foealfcu = min(1.0,((max(c_rticecu,min(c_rtwat,ptare))-c_rticecu)*c_rtwat_rticecu_r)**2)
      foeewmcu = c_r2es *(foealfcu *exp(c_r3les*(ptare-c_rtt)/(ptare-c_r4les))+&
                     (1.0-foealfcu)*exp(c_r3ies*(ptare-c_rtt)/(ptare-c_r4ies)))
      zqsat=foeewmcu*zqp

      zqsat=min(p_max_qsat,zqsat)
      zcor=1.0/(1.0-c_retv  *zqsat)
      zqsat=zqsat*zcor

      foedemcu =  foealfcu *c_r5alvcp*(1.0/(ptare-c_r4les)**2)+&
             (1.0-foealfcu)*c_r5alscp*(1.0/(ptare-c_r4ies)**2)

      zcond=(pq-zqsat)/(1.0+zqsat*zcor*foedemcu)

      zcond=min(zcond,0.0)

      foeldcpmcu= foealfcu*c_ralvdcp+(1.0-foealfcu)*c_ralsdcp
      pt=pt+foeldcpmcu*zcond

      pq=pq-zcond

      !--update ptare
      ptare = pt

      foealfcu = min(1.0,((max(c_rticecu,min(c_rtwat,ptare))-c_rticecu)*c_rtwat_rticecu_r)**2)
      foeewmcu = c_r2es *(foealfcu *exp(c_r3les*(ptare-c_rtt)/(ptare-c_r4les))+&
                     (1.0-foealfcu)*exp(c_r3ies*(ptare-c_rtt)/(ptare-c_r4ies)))
      zqsat=foeewmcu*zqp

      zqsat=min(0.5,zqsat)
      zcor=1.0/(1.0-c_retv  *zqsat)
      zqsat=zqsat*zcor

      foedemcu =  foealfcu *c_r5alvcp*(1.0/(ptare-c_r4les)**2)+&
             (1.0-foealfcu)*c_r5alscp*(1.0/(ptare-c_r4ies)**2)
      zcond1=(pq-zqsat)/(1.0+zqsat*zcor*foedemcu)

      if(zcond == 0.0)zcond1=min(zcond1,0.0)
      foeldcpmcu= foealfcu*c_ralvdcp+(1.0-foealfcu)*c_ralsdcp
      pt=pt+foeldcpmcu*zcond1
      pq=pq-zcond1

      !-- set output --------------------------
      q_wetbulb =  pq
      t_wetbulb =  pt
      evap      = -ZCOND1 != q_wetbulb-qo_cup, source for water vapor
   end subroutine get_wetbulb
   !------------------------------------------------------------------------------------
   subroutine cup_forcing_ens_shal(itf,ktf,its,ite,kts,kte,dtime,ichoice,ierrc,ierr  &
                                  ,klcl,kpbl,kbcon,k22,ktop,xmb,tsur,cape,h_sfc_flux &
                                  ,le_sfc_flux,zws,po, hco, heo_cup,po_cup,t_cup,dhdt&
                                  ,rho,xff_shal2d,xf_dicycle,tke_pbl,wlpool          &
                                  ,xf_coldpool,ke_gustfront)

      implicit none
      integer                     ,intent (in)   :: itf,ktf,its,ite, kts,kte,ichoice
      integer ,dimension (:)      ,intent (in)   :: klcl,kpbl,kbcon,k22,ktop
      real                        ,intent (in)   :: dtime
      real    ,dimension (:)      ,intent (in)   :: tsur,cape,h_sfc_flux,le_sfc_flux &
                                                   ,zws,tke_pbl,wlpool,ke_gustfront
      real    ,dimension (:,:)    ,intent (in)   :: po,hco,heo_cup,po_cup,t_cup,dhdt,rho
 
      integer ,dimension (:)      ,intent (in)   :: ierr
      character*128,dimension (:) ,intent (in)   :: ierrc
      real    ,dimension (:)      ,intent (inout):: xmb,xf_dicycle,xf_coldpool
      real    ,dimension (:,:)    ,intent (out  ):: xff_shal2d

      !---local vars
      real   ,dimension (its:ite)    :: xmbmax
      integer :: i,k,kbase,vtp_index
      real    :: blqe,trash,tcold,fin,fsum,efic,thot,dp
      real    ,dimension (shall_closures)  :: xff_shal
      !-- tuning numbers for the TKE-based closure for shallow convection   
      real,parameter :: p_k1 = 1.2, p_cloud_area = 0.15
      !
      xmb       (:)     = 0.
      xf_dicycle(:)     = 0.

      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)

         xmbmax(i)=100.*(po(kbcon(i),i)-po(kbcon(i)+1,i))/(c_grav*dtime)

         !- limiting the mass flux at cloud base
         xmbmax(i)=min(p_xmbmaxshal,xmbmax(i))
         
         !- cloud base
         kbase=kbcon(i)
         !kbase=klcl(i)

         !--- closure from Grant (2001): ichoice = 1
         xff_shal(1)=.030*zws(i)*rho(kpbl(i),i)        
         xff_shal(2)=xff_shal(1)
         xff_shal(3)=xff_shal(1)
         
         !--- closure from the heat-engine principle : ichoice = 4
         !- Renno and Ingersoll(1996), Souza et al (1999)
         !- get the averaged environment temperature between cloud base
         !- and cloud top
         tcold=0.
         do k=kbase,ktop(i)
            dp   = po_cup(k,i)-po_cup(k+1,i)
            tcold= tcold + t_cup(k,i)*dp
         enddo
         tcold=tcold/(po_cup(kbase,i)-po_cup(ktop(i)+1,i))

         !-surface temperature
         thot=tsur(i)  ! + ztexec(i)
         !- thermodynamic eficiency
         !efic = max(0.05, (thot-tcold)/thot )
         efic = max(0.0, (thot-tcold)/thot )

         !- total heat flux from surface
         fin = max(0.0, h_sfc_flux(i)+le_sfc_flux(i))

         !--- mass flux at cloud base
         !if(cape(i) > 0.0 .and. h_sfc_flux(i) >0.0 ) then
         if(cape(i) > 0.0  ) then
            xff_shal(4) = efic * fin / cape(i)
         else
            xff_shal(4) = 0.0
         endif
         xff_shal(5)=xff_shal(4)
         xff_shal(6)=xff_shal(4)

         !--- closure from boundary layer QE (Raymond 1995): ichoice = 7
         blqe=0.
         trash=0.
         if(k22(i).lt.kpbl(i)+1)then
            do k=kts,kbase
               blqe=blqe+100.*dhdt(k,i)*(po_cup(k,i)-po_cup(k+1,i))/c_grav
            enddo
            trash = max((hco(kbase,i)-heo_cup(kbase,i)),1.e1)
            xff_shal(7)=max(0.,blqe/trash)
         else
            xff_shal(7)=0.0
         endif
         xff_shal(8)= xff_shal(7)
         xff_shal(9)= xff_shal(7)

         !--- new closure based on the PBL TKE mean (Zheng et al, 2020 GRL): ichoice = 10
         !-- shallow cumulus active area is for now keept by 0.15 (Zheng 2021 p. commun.)
         !-- k1 is 'slope' of the curve between Wb x (TKE_PBL)**0.5
         !--        and varies between 1.2 (from lidar) to 1.6 (from WRF and SAM models)
        !xff_shal(10) = p_cloud_area * rho(kbase,i) * p_k1 * sqrt(max(tke_pbl(i),ke_gustfront(i)))
         xff_shal(10) = p_cloud_area * rho(kbase,i) * p_k1 * sqrt(tke_pbl(i))
         xff_shal(11) = xff_shal(10)
         xff_shal(12) = xff_shal(10)

         !--- store all closures for later.
         xff_shal2d(i,:) = xff_shal(:)

      enddo
   end subroutine cup_forcing_ens_shal
   !------------------------------------------------------------------------------------
   subroutine cup_up_lightning(itf,ktf,its,ite, kts,kte, ierr, kbcon,ktop,xland,cape &
                              ,zo,zo_cup,t_cup,t,tempco,qrco,po_cup,rho,prec_flx     &
                              ,lightn_dens)

      !=====================================================================================
      !- Lightning parameterization based on:
      !- "A Lightning Parameterization for the ECMWF Integrated Forecasting System"
      !-  P. Lopez, 2016 MWR
      !
      !- Coded/adapted to the GF scheme by Saulo Freitas (10-Aug-2019)
      !=====================================================================================
      implicit none
      integer                 ,intent(in)  :: itf,ktf, its,ite, kts,kte
      integer, dimension(:)   ,intent(in)  :: ierr,kbcon,ktop
      real,    dimension(:)   ,intent(in)  :: cape,xland
      real,    dimension(:,:) ,intent(in)  :: po_cup,zo_cup,t_cup,t,tempco,zo &
                                             ,qrco,rho,prec_flx

      real,    dimension(:)   ,intent(out) :: lightn_dens ! lightning flash density
                                                          ! rate (units: 1/km2/day)

      !-- locals
      real, parameter :: p_V_graup     = 3.0  ! m/s
      real, parameter :: p_V_snow      = 0.5  ! m/s
      real, parameter :: p_beta_land   = 0.70 ! 1
      real, parameter :: p_beta_ocean  = 0.45 ! 1
      real, parameter :: p_alpha       = 37.5 ! 1
      real, parameter :: p_t_initial   =  0.0 + 273.15 ! K
      real, parameter :: p_t_final     = -25. + 273.15 ! K

      integer :: i, k, k_initial, k_final,vtp_index
      real    :: Q_R, z_base,beta,prec_flx_fr,dz
      real,    dimension(kts:kte) :: p_liq_ice, q_graup,q_snow
      
      lightn_dens(:) = 0.0
      
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)

         beta= xland(i)*p_beta_ocean + (1.-xland(i))*p_beta_land

         q_graup(:) = 0.
         q_snow (:) = 0.

         do k=kts,ktop(i)

            p_liq_ice(k) = FractLiqF(tempco(k,i))

            prec_flx_fr=   p_liq_ice(k)*prec_flx(k,i)/rho(k,i)

            q_graup(k) =      beta *prec_flx_fr/p_V_graup ! - graupel mixing ratio (kg/kg)
            q_snow (k) =  (1.-beta)*prec_flx_fr/p_V_snow  ! - snow    mixing ratio (kg/kg)

         enddo

         k_initial = minloc(abs(tempco(kbcon(i):ktop(i),i)-p_t_initial),1)+kbcon(i)-1
         k_final   = minloc(abs(tempco(kbcon(i):ktop(i),i)-p_t_final  ),1)+kbcon(i)-1

         Q_R = 0.0
         do k = k_initial, k_final
            dz  = zo(k,i)-zo(k-1,i)
            Q_R = Q_R + dz*rho(k,i)*(q_graup(k)*(qrco(k,i)+q_snow(k)))
           !print*,"qr=",q_r,tempco(k,i)-273.15,k,tempco(k,i)-t_initial
         enddo

         z_base = zo_cup(kbcon(i),i)/1000. ! km

         !---
         !--- lightning flash density (units: number of flashes/km2/day) - equation 5
         !--- (to compare with Lopez 2016's results, convert to per year: lightn_dens*365)
         !
         lightn_dens(i) = p_alpha * Q_R *sqrt (max(0.,cape(i))) * min(z_base,1.8)**2
        !
      enddo
   end subroutine cup_up_lightning
   !------------------------------------------------------------------------------------
   subroutine cup_up_rain(cumulus,klcl,kbcon,ktop,k22,ierr,xland,zo_cup,qco,qrco,pwo,pwavo &
                         ,po,p_cup,t_cup,tempco,zuo,up_massentr,up_massdetr,vvel2d,rho     &
                         ,qrs,itf,ktf,its,ite, kts,kte)

      implicit none
      character *(*)         ,intent (in) :: cumulus
      integer                ,intent (in) :: itf,ktf,its,ite,kts,kte
      integer, dimension(:)  ,intent (in) :: kbcon,ktop,k22,klcl,ierr
      real,    dimension(:)  ,intent (in) :: xland,pwavo
      real,    dimension(:,:),intent (in) :: zo_cup,qco,qrco,pwo,po,p_cup,t_cup,zuo       &
                                            ,up_massentr,up_massdetr,vvel2d,tempco,rho
      !--for future use (rain water mixing ratio)
      real,    dimension(:,:),intent (out) :: qrs      !-- units kg[water]/kg[air]

      !-- locals
      integer :: i,k,vtp_index
      real :: tmp
      integer, dimension(its:ite) :: start_level
      real :: dz,z1,zrnew,zc,zd,zint,z2,zrold,denom,fall_fact,wup, exp1,R_vr
      real,    dimension(kts:kte) :: prec_flx_rain,prec_flx_snow
      real,    dimension(kts:kte) :: pw,p_liq_ice ! - rainfall source
      real,    parameter :: rho1000mb = 1.2 , rhow = 1000., N_r = 0.1 ! cm^-3, rainfall drops number concen
      real,    parameter :: exp_KR    = 1./5. &! Kuo & Raymond 1983
         ,exp_SB    = 2./3.  ! Seifert & Beheng 2006 eq 27

      qrs = 0.
      if(c0 < 1.e-6 ) return

      !--- rain water mixing ratio
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)

         do k=ktop(i),kts,-1

            p_liq_ice(k) = FractLiqF(tempco(k,i))

            !--- transport + mixing
            denom = zuo(k+1,i)-.5*up_massdetr(k,i)+up_massentr(k,i)
            if(denom > 0.) then
               qrs(k,i) = (qrs(k+1,i)*zuo(k+1,i)-.5*up_massdetr(k,i)* qrs(k+1,i))/ denom
            else
               qrs(k,i) =  qrs(k+1,i)
            endif

            !--- rain source
            pw(k)= pwo(k,i)/(1.e-16+zuo(k,i))

            !-- rainfall sedimentation
            !-- Kuo & Raymond 1983
            !-- fallout of rain (21.18 * qrs^0.2 have m/s as units with qrs in kg/kg)
            !---                                half velocity for ice phase
            fall_fact = 21.18 *( p_liq_ice(k) + 0.5*(1.-p_liq_ice(k) ))
            exp1      = exp_KR
            !
            !-- Seifert & Beheng 2006 eq 27, units= m/s
            ! fall_fact = 159.*sqrt(rho1000mb/rho(k,i))*( p_liq_ice(k) + 0.5*(1.-p_liq_ice(k) ))
            ! exp1      = exp_SB

            !-- Kogan 2013
            R_vr = (4.*c_pi*rhow/(3.*rho(k,i)))**(-1./3.) * (qrs(k,i) + pw(k))**(1./3.) * &
               (N_r)**(-1./3.)! cm, N_r is the rainfall drops number concentration, and not CCN or N_c

            R_vr = max(40., R_vr * 1.e-2 * 1.e+6 )   ! micrometer
            fall_fact = 1.e-2*(2.4*R_vr-62.0)        ! m/s
            exp1      = 0.

            wup      = min(15.,max(2.,vvel2d(k,i)))
            z2       = fall_fact/wup

            !--exact solution
            if(qrs(k,i) + pw(k)>0.) then
               !-- this is the sedimentation speed divided by W_up
               zd= z2* (qrs(k,i) + pw(k))**exp1
               zint= EXP(-zd)
               zrold=qrs(k,i)
               zc=pw(k)
               zrnew=zrold*zint+zc/zd*(1.-zint)
               zrnew=max(0.0,min(qrs(k,i) + pw(k),zrnew))
            else
               zrnew=0.
            endif
            qrs(k,i)=zrnew

          !--solution with rk3
          !z1= qrs(k,i)
          !z1= qrs(k,i) +  1./3.*(pw(k) - z2*(z1)**exp1 * z1)
          !z1= qrs(k,i) +  0.500*(pw(k) - z2*(z1)**exp1 * z1)
          !z1= qrs(k,i) +  1.000*(pw(k) - z2*(z1)**exp1 * z1)
          !qrs(k,i)=z1
          !---
          !--- solution 2
          !qrs(k,i)= qrs(k,i) + pw(k) - (fall_fact*(qrs(k,i) + pw(k))**exp1) / wup * qrs(k,i)

          !--- solution 3
          !tmp = 0.5*(qrs(k,i) + pw(k) + qrs(k+1,i) + pw(k+1))
          !qrs(k,i)= qrs(k,i) + pw(k) - z2*(tmp)**exp1 * 0.5*(qrs(k,i)+qrs(k+1,i))

          !print*,"rr3=",k,zo_cup(k,i), pwo(k,i)*1000.,qrs(k,i)*1000.,fall_fact*(qrs(k,i) + pw(k))**exp1
         enddo
      enddo

   end subroutine cup_up_rain
   !------------------------------------------------------------------------------------
   subroutine get_interp(q_old,t_old,po_cup,q_new,t_new)
      implicit none
      real    ,intent (in   ) :: po_cup ! original
      real    ,intent (inout) :: q_old,t_old,q_new,t_new ! extrapolated

      !---locals
      real ::  zqp, zcond1, zcor, zqsat
      real ::  psp, pt , pq, ptare
      real ::  foealfcu, foeewmcu,foedemcu,foeldcpmcu
      
      !real, parameter :: &
      !  c_rd=287.06                             &
      ! ,c_rv=461.52                             &
      ! ,c_rcpd=1004.71                          &
      ! ,c_rtt=273.16                            &
      ! ,c_rhoh2o=1000.                          &
      ! ,c_rlvtt=2.5008e+6                       &
      ! ,c_rlstt=2.8345e+6                       &
      ! ,c_retv = c_rv/c_rd-1.0                  &
      ! ,c_rlmlt= c_rlstt-c_rlvtt                &
      ! ,c_rcpv=4.*c_rv                          &
      ! ,c_r2es=611.21*c_rd/c_rv                 &
      ! ,c_r3les=17.502                          &
      ! ,c_r3ies=22.587                          &
      ! ,c_r4les=32.19                           &
      ! ,c_r4ies=-0.7                            &
      ! ,c_r5les= c_r3les*(c_rtt-c_r4les)        &
      ! ,c_r5ies= c_r3ies*(c_rtt-c_r4ies)        &
      ! ,c_r5alvcp= c_r5les*c_rlvtt/c_rcpd       &
      ! ,c_r5alscp= c_r5ies*c_rlstt/c_rcpd       &
      ! ,c_ralvdcp= c_rlvtt/c_rcpd               &
      ! ,c_ralsdcp= c_rlstt/c_rcpd               &
      ! ,c_ralfdcp= c_rlmlt/c_rcpd               &
      ! ,c_rtwat= c_rtt                          &
      ! ,c_rtber= c_rtt-5.                       &
      ! ,c_rtbercu= c_rtt-5.0                    &
      ! ,c_rtice= c_rtt-23.                      &
      ! ,c_rticecu= c_rtt-23.                    &
      ! ,c_rtwat_rtice_r=1./(c_rtwat-c_rtice)    &
      ! ,c_rtwat_rticecu_r=1./(c_rtwat-c_rticecu)&
      ! ,c_rvtmp2= c_rcpv/c_rcpd-1.              &
      ! ,c_zqmax=0.5

      integer :: i

      pt  = t_old       ! k
      pq  = q_old       ! kg/kg
      psp = po_cup*100. ! hpa

      !-- for testing
      !              psp                   temp        q                     zcond1
      ! input    27940.0000000000        236.604976804749       3.220181796223121e-004
      ! output   27940.0000000000        236.361132108860       4.084506812610067e-004
      !  pt  = 236.604976804749      ! k
      !  pq  = 3.220181796223121e-004       ! kg/kg
      !  psp = 27940. ! hpa
      !----------------------
      !print*,"1",psp,pt,pq

      zqp   =1.0/psp
      do i=1,2
         ptare = pt

         foealfcu = min(1.0,((max(c_rticecu,min(c_rtwat,ptare))-c_rticecu)*c_rtwat_rticecu_r)**2)
         foeewmcu = c_r2es *(foealfcu *exp(c_r3les*(ptare-c_rtt)/(ptare-c_r4les))+&
                  (1.0-foealfcu)*exp(c_r3ies*(ptare-c_rtt)/(ptare-c_r4ies)))
         zqsat=foeewmcu*zqp

         !    if(1.0-retv  *zqsat == 0.) then
         !
         !      print*,"zqsat=",zqp,foeewmcu,q_old,t_old,po_cup,q_new,t_new
         !3.5491847e-02   46.36052      0.5000000       249.8219
         !  0.2817549      0.5000000       249.8219
         !      call flush(6)
         !      stop 3333
         !    endif

         zcor=1.0/(1.0-c_retv  *zqsat)
         zqsat=zqsat*zcor

         foedemcu =  foealfcu     *c_r5alvcp*(1.0/(ptare-c_r4les)**2)+&
                    (1.0-foealfcu)*c_r5alscp*(1.0/(ptare-c_r4ies)**2)


         zcond1=(pq-zqsat)/(1.0+zqsat*zcor*foedemcu)

         foeldcpmcu= foealfcu*c_ralvdcp+(1.0-foealfcu)*c_ralsdcp
         pt=pt+foeldcpmcu*zcond1
         pq=pq-zcond1
      enddo
      !-- FINAL --------------------------
      q_new =  PQ
      t_new =  PT
      !print*,"2",PSP,PT,PQ
      !print*,"E",100*(PT-236.361132108860)/236.361132108860,100*(PQ-4.084506812610067E-004)/4.084506812610067E-004
   end subroutine get_interp
   !------------------------------------------------------------------------------------
   real function FractLiqF(temp2) ! temp2 in Kelvin, fraction between 0 and 1.
      implicit none
      real,intent(in)  :: temp2 ! K
      real             :: temp,ptc
      real, parameter  :: max_temp = 46. !Celsius
      select case(FRAC_MODIS)

         case (1)
            temp = temp2-273.16 !Celsius
            temp = min(max_temp,max(-max_temp,temp))
            ptc  = 7.6725 + 1.0118 *temp    + 0.1422 *temp**2 + &
                            0.0106 *temp**3 + 3.39e-4*temp**4    + &
                            3.95e-6*temp**5
            FractLiqF = 1./(1.+exp(-ptc))

         !WMP skew ice fraction for deep convective clouds
         !       FractLiqF = FractLiqF**4
         !WMP
         case default
            FractLiqF =  min(1., (max(0.,(temp2-c_Tice))/(c_T00-c_Tice))**2)

         end select

   end function
   !------------------------------------------------------------------------------------
   subroutine sound(part, cumulus, int_time, dtime, ens4, itf, its, ite, kts, kte, xlats, xlons, jcol, whoami_all &
                    , z, qes, he, hes, t, q, po, z1, psur, zo, qeso, heo, heso, tn, qo, us, vs, omeg, xz, h_sfc_flux, le_sfc_flux &
                    , tsur, dx, stochastic_sig, zws, ztexec, zqexec, xland, kpbl, k22, klcl, kbcon, ktop, aa0, aa1, sig, xaa0, hkb &
                    , xmb, pre, edto, zo_cup, dhdt, rho, zuo, zdo, up_massentro, up_massdetro, outt, outq, outqc, outu, outv)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: procedureName = 'sound' ! Subroutine Name
      real, parameter :: p_latsnd = -10., p_lonsnd = 301., p_deltx = 0.2
      !      real, parameter :: LATSND= -8.72, LONSND= 186.6, DELTX=0.2
   
      !Variables (input, output, inout)
      integer, intent(in) ::ens4, itf, its, ite, kts, kte, jcol, whoami_all, part

      integer, intent(in) :: k22(:)
      integer, intent(in) :: klcl(:)
      integer, intent(in) :: kbcon(:)
      integer, intent(in) :: ktop(:)

      real, intent(in) :: zo_cup(:, :)
      real, intent(in) :: zuo(:, :)
      real, intent(in) :: zdo(:, :)
      real, intent(in) :: up_massentro(:, :)
      real, intent(in) :: up_massdetro(:, :)
      real, intent(in) :: outt(:, :)
      real, intent(in) :: outq(:, :)
      real, intent(in) :: outqc(:, :)
      real, intent(in) :: outu(:, :)
      real, intent(in) :: outv(:, :)
      real, intent(in) :: stochastic_sig(:)
      real, intent(in) :: xland(:)
      real, intent(in) :: aa0(:)
      real, intent(in) :: aa1(:)
      real, intent(in) :: xaa0(:)
      real, intent(in) :: hkb(:)
      real, intent(in) :: xmb(:)
      real, intent(in) :: pre(:)
      real, intent(in) :: edto(:)
      real, intent(in) :: sig(:)

      real, intent(in) :: int_time, dtime

      character(len=*), intent(in)    :: cumulus

      integer, intent(inout) :: kpbl(:)

      real, intent(inout) :: h_sfc_flux(:)
      real, intent(inout) :: le_sfc_flux(:)
      real, intent(inout) :: tsur(:)
      real, intent(inout) :: dx(:)
      real, intent(inout) :: zws(:)
      real, intent(inout) :: ztexec(:)
      real, intent(inout) :: zqexec(:)
      real, intent(inout) :: xlats(:)      
      real, intent(inout) :: xlons(:)      
      real, intent(inout) :: z1(:)
      real, intent(inout) :: psur(:)
      real, intent(inout) :: qes(:, :)
      real, intent(inout) :: he(:, :)
      real, intent(inout) :: hes(:, :)
      real, intent(inout) :: t(:, :)
      real, intent(inout) :: q(:, :)
      real, intent(inout) :: po(:, :)
      real, intent(inout) :: zo(:, :)
      real, intent(inout) :: heo(:, :)
      real, intent(inout) :: heso(:, :)
      real, intent(inout) :: tn(:, :)
      real, intent(inout) :: qo(:, :)
      real, intent(inout) :: us(:, :)
      real, intent(inout) :: vs(:, :)
      real, intent(inout) :: dhdt(:, :)
      real, intent(inout) :: omeg(:, :, :)

      real, intent(out) :: z(:, :)
      real, intent(out) :: xz(:, :)
      real, intent(out) :: qeso(:, :)
      real, intent(out) :: rho(:, :)

      !---locals
      integer :: i, k, x_kte, x_i, x_jcol, x_k
      real :: x_time
      real, dimension(its:ite) :: x_stochastic_sig, x_xland
      
      character(len=200) :: lixo


      if (trim(rundata) == "NONE") then
         if (mod(int_time, 3600.) < dtime) then
            open (15, file="dataLXXX.dat_"//trim(cumulus), status='unknown', position="APPEND")
            if (part == 1) then
               do i = its, itf
                  if (xlats(i) > p_latsnd - p_deltx .and. xlats(i) < p_latsnd + p_deltx) then
                     if (xlons(i) > p_lonsnd - p_deltx .and. xlons(i) < p_lonsnd + p_deltx) then

                        print *, "==============================================="
                        print *, "00>", i, jcol, xlats(i), xlons(i), whoami_all, int_time/3600.
                        call flush (6)

                        write (15, *) "====begin====="
                        write (15, *) "i,jcol,xlats(i),xlons(i),int_time/3600."
                        write (15, *) i, jcol, xlats(i), xlons(i), int_time/3600.

                        write (15, *) "kte,z1(i),psur(i),tsur(i),xland(i)"
                        write (15, *) kte, z1(i), psur(i), tsur(i), xland(i)

                        write (15, *) "h_sfc_flux(i),le_sfc_flux(i),ztexec(i),zqexec(i)"
                        write (15, *) h_sfc_flux(i), le_sfc_flux(i), ztexec(i), zqexec(i)

                        write (15, *) "stochastic_sig(i), dx(i),zws(i),kpbl(i)"
                        write (15, *) stochastic_sig(i), dx(i), zws(i), kpbl(i)

                        write (15, *) "=>k zo po t tn-t q qo-q us vs qes he hes qeso-qes heo-he heso-hes dhdt omeg"
                        do k = kts, kte
                           write (15, 100) k, zo(k,i), po(k,i), t(k,i) &
                                         , tn(k,i) - t(k,i), q(k,i), qo(k,i) - q(k,i) &
                                         , us(k,i), vs(k,i), qes(k,i), he(k,i) &
                                         , hes(k,i), qeso(k,i) - qes(k,i) &
                                         , heo(k,i) - he(k,i), heso(k,i) - hes(k,i) &
                                         , dhdt(k,i), omeg(k, i, 1:ens4)
                        end do

                     end if
                  end if
               end do
            else
               do i = its, itf
                  if (xlats(i) > p_latsnd - p_deltx .and. xlats(i) < p_latsnd + p_deltx) then
                     if (xlons(i) > p_lonsnd - p_deltx .and. xlons(i) < p_lonsnd + p_deltx) then

                        write (15, *) "====outputs======="
                        write (15, *) "L=", i, jcol, xlats(i), xlons(i), whoami_all
                        write (15, *) "A=", aa0(i), aa1(i), xaa0(i), sig(i)
                        write (15, *) "K=", k22(i), klcl(i), kpbl(i), kbcon(i), ktop(i)
                        write (15, *) "Z=", zo_cup(i, k22(i)) - z1(i), zo_cup(klcl(i),i) - z1(i), zo_cup(kpbl(i),i) - z1(i) &
                           , zo_cup(kbcon(i),i) - z1(i), zo_cup(ktop(i),i) - z1(i)
                        write (15, *) "H=", hkb(i)/real(c_cp), edto(i)
                        write (15, *) "T=", maxval(outt(1:ktop(i),i))*86400., maxval(outq(1:ktop(i),i))*86400.*1000., &
                           minval(outt(1:ktop(i),i))*86400., minval(outq(1:ktop(i),i))*86400.*1000.
                        write (15, *) "P=", xmb(i)*1000., 'g/m2/s', 3600*pre(i), 'mm/h'
                        if (xmb(i) > 0.0) then
                           write (15, *) "=> k zo po zuo,zdo,up_massentro,up_massdetro,outt, outq,outqc,outu,outv"
                           do k = kts, kte
                              write (15, 101) k, zo(k,i), po(k,i), zuo(k,i), zdo(k,i), up_massentro(k,i), up_massdetro(k,i) &
                                            , outt(k,i)*86400., outq(k,i)*86400.*1000., outqc(k,i)*86400.*1000., outu(k,i) &
                                            *86400., outv(k,i)*86400.

                           end do
                        end if
                        write (15, *) "=====end=========="
                     end if
                  end if
               end do
            end if
            close (15)
         end if
      else
         if (part == 1) then
            open (15, file=trim(rundata), status='old')
            i = 1
            read (15, *) lixo
            read (15, *) lixo
            read (15, *) x_i, x_jcol, xlats(i), xlons(i), x_time
            read (15, *) lixo
            read (15, *) x_kte, z1(i), psur(i), tsur(i), x_xland(i)
            !-- check
            if (x_kte .ne. kte) stop " X_kte .ne. kte "
            read (15, *) lixo
            read (15, *) h_sfc_flux(i), le_sfc_flux(i), ztexec(i), zqexec(i)
            read (15, *) lixo
            read (15, *) x_stochastic_sig(i), dx(i), zws(i), kpbl(i)
            read (15, *) lixo
            do k = kts, kte
               read (15, 100) x_k, zo(k,i), po(k,i), t(k,i), tn(k,i), q(k,i), qo(k,i), us(k,i), vs(k,i) , qes(k,i) &
                            , he(k,i), hes(k,i), qeso(k,i), heo(k,i), heso(k,i), dhdt(k,i), omeg( k, i,1:ens4)
            end do
            close (15)
            !---settings
            tn(:,i) = t(:,i) + tn(:,i) ! input is delta(T)
            qo(:,i) = q(:,i) + qo(:,i) ! input is delta(Q)
            qeso(:,i) = qes(:,i) + qeso(:,i) ! input is delta(Q)
            heo(:,i) = he(:,i) + heo(:,i) ! input is delta(H)
            heso(:,i) = hes(:,i) + heso(:,i) ! input is delta(HO)
            xz(:,i) = zo(:,i)
            z(:,i) = zo(:,i)
            rho(:,i) = 1.e2*po(:,i)/(c_rgas*t(:,i))
         else
            do i = its, itf
               if (xlats(i) > p_latsnd - p_deltx .and. xlats(i) < p_latsnd + p_deltx) then
                  if (xlons(i) > p_lonsnd - p_deltx .and. xlons(i) < p_lonsnd + p_deltx) then

                     print *, "====outputs======="
                     print *, "A=", aa0(i), aa1(i), xaa0(i), sig(i)
                     print *, "K=", k22(i), klcl(i), kpbl(i), kbcon(i), ktop(i)
                     print *, "Z=", zo_cup(k22(i),i) - z1(i), zo_cup(klcl(i),i) - z1(i), zo_cup(kpbl(i),i) - z1(i) &
                        , zo_cup(kbcon(i),i) - z1(i), zo_cup(ktop(i),i) - z1(i)
                     print *, "H=", hkb(i)/real(c_cp), edto(i)
                     print *, "T=", maxval(outt(1:ktop(i),i))*86400., maxval(outq(1:ktop(i),i))*86400.*1000. &
                                  , minval(outt(1:ktop(i),i))*86400., minval(outq(1:ktop(i),i))*86400.*1000.
                     print *, "P=", xmb(i)*1000., 'g/m2/s', 3600*pre(i), 'mm/h'
                     if (xmb(i) > 0.0) then
                        print *, "=> k zo po zuo,zdo,up_massentro,up_massdetro,outt, outq,outqc,outu,outv"
                        do k = kts, kte
                           write (*, 101) k, zo(k,i), po(k,i) &
                              , zuo(k,i), zdo(k,i), up_massentro(k,i), up_massdetro(k,i), outt(k,i)*86400. &
                              , outq(k,i)*86400.*1000., outqc(k,i)*86400.*1000., outu(k,i)*86400., outv(k,i)*86400.
                        end do
                     end if
                  end if
               end if
            end do
         end if
      end if
100   format(1x, i4, 16e16.8)
101   format(1x, i4, 11e16.8)

   end subroutine sound
   !------------------------------------------------------------------------------------
   subroutine cloud_dissipation(cumulus,itf,ktf, its,ite, kts,kte,ierr,kbcon,ktop,dtime,xmb,xland&
                               ,qo_cup,qeso_cup,po_cup,outt,outq,outqc,zuo,vvel2d       &
                               ,qrco,sig,tempco,qco,tn_cup,heso_cup,zo,zo_cup)

      implicit none
      character*(*)                ,intent(in)    :: cumulus
      integer                      ,intent(in)    :: itf,ktf, its,ite, kts,kte
      real                         ,intent(in)    :: dtime
      integer, dimension(:)        ,intent(in)    :: ierr,kbcon,ktop
      real,    dimension(:)        ,intent(in)    :: xmb,xland,sig
      real,    dimension(:,:)      ,intent(in)    :: po_cup,qo_cup,qeso_cup,zuo,vvel2d &
                                                    ,tempco,qco,tn_cup,heso_cup,zo,zo_cup
      real,    dimension(:,:)      ,intent(inout) :: outt,outq,outqc,qrco

      !-- locals
      integer :: i,k,vtp_index
      real    ::  del_t,del_q,dp,frh,rho_hydr
      real    :: qrc_diss,fractional_area,outqc_diss,outq_mix,outt_diss,outt_mix,tempx,qvx
      real, parameter :: cloud_lifetime= 1800.

      integer, parameter :: versionx = 2
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)

         do k=ktop(i),kbcon(i),-1

            !--- get air density at full layer (model levels) by hydrostatic balance (kg/m3)
            rho_hydr=100.*(po_cup(k,i)-po_cup(k+1,i))/(zo_cup(k+1,i)-zo_cup(k,i))/c_grav

            !--- cloud liq/ice remained in the convection plume
            qrc_diss = max(0., qrco(k,i) - outqc(k,i) * dtime)

            !dp  = 100.*(po_cup(k,i)-po_cup(k+1,i))

            !--- get relative humidity
            frh = 0. !min(qo_cup(k,i)/qeso_cup(k,i),1.)

            !--- estimation of the fractional area
            fractional_area = (xmb(i)/sig(i)) * zuo(k,i) / (rho_hydr*vvel2d(k,i))

            !--- source of enviroment moistening/cooling due to the 'remained' cloud dissipation into it.
            outqc_diss = ( qrc_diss * (1.-frh) ) / cloud_lifetime

            if(versionx==1 .or. COUPL_MPHYSICS .eqv. .false.) then

               outt_diss  = -outqc_diss*(c_alvl/c_cp) !--- cooling

               !--- source of enviroment moistening/warming due to the 'remained' in-cloud water vapor mixing into it.
               !  qvx   = qco   (k,i)
               !  tempx = tempco(k,i)
               qvx   = qeso_cup(k,i)
               tempx =(heso_cup(k,i)-c_grav*zo(k,i)-c_alvl*qeso_cup(k,i))/c_cp

               outq_mix = ( qvx   - qo_cup(k,i) ) / cloud_lifetime

               outt_mix = ( tempx - tn_cup(k,i) ) / cloud_lifetime

               !-- feedback
               del_q = (outqc_diss + outq_mix) * use_cloud_dissipation * fractional_area ! units: kg[water]/kg[air}/sec
               del_t = (outt_diss  + outt_mix) * use_cloud_dissipation * fractional_area ! units: K/sec

               outq (k,i) = outq (k,i) + del_q
               outt (k,i) = outt (k,i) + del_t

            else

               outqc(k,i) = outqc(k,i) + outqc_diss* fractional_area * use_cloud_dissipation

            endif

             !print*,"diss2=",k,real(outqc_diss*86400.*1000),real(sqrt(1.-sig(i)),4),real( fractional_area*100.,4)

            qrco    (k,i) = max(0., qrco(k,i) - outqc_diss * use_cloud_dissipation * fractional_area *dtime)
           !if(qrco (k,i) <0.) print*,"qrc<0",trim(cumulus),qrco(k,i)

         enddo
      enddo

   end subroutine cloud_dissipation
   !------------------------------------------------------------------------------------
   subroutine readGFConvParNML(mynum)
      implicit none
      integer,intent(in) :: mynum
      integer :: nlunit = 4
      character (len = 64) :: fn_nml = 'GF_ConvPar_nml'
      logical :: exists
      namelist /GF_NML/ icumulus_gf, closure_choice,use_scale_dep,dicycle     &
         ,use_tracer_transp, use_tracer_scaven,use_flux_form, use_tracer_evap &
         ,downdraft,use_fct  &
         ,use_rebcb, vert_discr, clev_grid,apply_sub_mp, alp1         &
         ,sgs_w_timescale, lightning_diag,autoconv, overshoot,use_wetbulb &
         ,c0_deep, qrc_crit,lambau_deep,lambau_shdn,c0_mid                     &
         ,cum_max_edt_land  ,cum_max_edt_ocean, cum_hei_down_land                 &
         ,cum_hei_down_ocean,cum_hei_updf_land, cum_hei_updf_ocean                &
         ,cum_entr_rate ,tau_deep,tau_mid                                         &
         ,use_momentum_transp ,moist_trigger,frac_modis                           &
         ,cum_use_excess,cum_ave_layer,adv_trigger,use_smooth_prof                &
         ,use_cloud_dissipation,use_smooth_tend,use_gustiness, use_random_num     &
         ,dcape_threshold,beta_sh,c0_shal,use_linear_subcl_mf,liq_ice_number_conc &
         ,cap_maxs,sig_factor,cum_fadj_massflx,lcl_trigger       &
         ,rh_dicycle,cum_t_star, convection_tracer, tau_ocea_cp, tau_land_cp      &
         ,use_memory, add_coldpool_prop ,mx_buoy1, mx_buoy2,max_tq_tend           &
         ,add_coldpool_clos,add_coldpool_trig,add_coldpool_diff,n_cldrop,output_sound&
         ,use_pass_cloudvol,use_lcl_ctrl_entr,use_rhu_ctrl_entr

      inquire (file = trim (fn_nml), exist = exists)
      if (.not. exists) then
         write (6, *) 'GF_convpar_nml :: namelist file: ', trim (fn_nml), ' does not exist'
         stop 31415
      else
         open (nlunit,file=fn_nml,status='old',form='formatted')
         read (nlunit,nml=GF_NML)
         close(nlunit)
      endif
      if(mynum==1) then
         !- print the namelist
         print*,"           "
         print*,"------------- GF ConvPar namelist -------------"
         print*,"!---- the main controls"
         print*, 'icumulus_gf        ' , icumulus_gf
         print*, 'cum_entr           ' , real(cum_entr_rate      ,4)
         print*, 'closure_choice     ' , closure_choice
         print*, 'use_scale_dep      ' , use_scale_dep
         print*, 'sig_factor         ' , real(sig_factor         ,4)
         print*, 'dicycle            ' , dicycle
         print*, 't_star             ' , real(cum_t_star         ,4)
         print*, 'rh_dicycle         ' , rh_dicycle
         print*, 'cap_maxs           ' , real(cap_maxs           ,4)
         print*, 'moist_trigger      ' , moist_trigger
         print*, 'adv_trigger        ' , adv_trigger
         print*, 'lcl_trigger        ' , lcl_trigger
         print*, 'dcape_threshold    ' , real(dcape_threshold    ,4)
         print*, 'tau_deep,tau_mid   ' , real(tau_deep,4),real(tau_mid,4)
         print*, 'sgs_w_timescale    ' , sgs_w_timescale
         print*, 'convection_tracer  ' , convection_tracer
         print*, 'add_coldpool_prop  ' , add_coldpool_prop
         print*, 'add_coldpool_clos  ' , add_coldpool_clos
         print*, 'add_coldpool_trig  ' , add_coldpool_trig
         print*, 'add_coldpool_diff  ' , add_coldpool_diff
         print*, 'tau_ocea_cp        ' , tau_ocea_cp
         print*, 'tau_land_cp        ' , tau_land_cp
         print*, 'mx_buoy1 - kJ/kg   ' , mx_buoy1*1.e-3
         print*, 'mx_buoy2 - kJ/kg   ' , mx_buoy2*1.e-3
         print*, 'use_memory         ' , use_memory
         print*, 'use_pass_cloudvol  ' , use_pass_cloudvol
         print*, 'use_lcl_ctrl_entr  ' , use_lcl_ctrl_entr
         print*, 'use_rhu_ctrl_entr  ' , use_rhu_ctrl_entr


         print*,'!--- controls rainfall evaporation'
         print*, 'use_rebcb          ' , use_rebcb
         print*, 'downdraft          ' , downdraft
         print*, 'max_edt_land       ' , real(cum_max_edt_land   ,4)
         print*, 'max_edt_ocean      ' , real(cum_max_edt_ocean  ,4)

         print*,'!---- boundary condition specification'
         print*, 'cum_use_excess     ' , cum_use_excess
         print*, 'cum_ave_layer      ' , real(cum_ave_layer      ,4)

         print*,'!---- for mass flux profiles - (deep ,shallow ,congestus)'
         print*, 'hei_down_land      ' , real(cum_hei_down_land  ,4)
         print*, 'hei_down_ocean     ' , real(cum_hei_down_ocean ,4)
         print*, 'hei_updf_land      ' , real(cum_hei_updf_land  ,4)
         print*, 'hei_updf_ocean     ' , real(cum_hei_updf_ocean ,4)
         print*, 'beta_sh            ' , real(beta_sh            ,4)
         print*, 'use_linear_subcl_mf' , use_linear_subcl_mf
         print*, 'use_smooth_prof    ' , use_smooth_prof
         print*, 'use_smooth_tend    ' , use_smooth_tend
         print*, 'use_random_num     ' , use_random_num

         print*,'!---- the cloud microphysics'
         print*, 'autoconv           ' , autoconv
         print*, 'c0_deep            ' , real(c0_deep            ,4)
         print*, 'c0_mid             ' , real(c0_mid             ,4)
         print*, 'c0_shal            ' , real(c0_shal            ,4)
         print*, 'qrc_crit           ' , real(qrc_crit           ,4)
         print*, 'n_cldrop           ' , real(n_cldrop           ,4)

         print*, '!--- for momentum transport'
         print*, 'use_momentum_trans ' , use_momentum_transp
         print*, 'lambau_deep        ' , real(lambau_deep        ,4)
         print*, 'lambau_shdn        ' , real(lambau_shdn        ,4)

         print*, '!--- for tracer transport'
         print*, 'use_tracer_transp  ' , use_tracer_transp
         print*, 'use_tracer_scaven  ' , use_tracer_scaven
         print*, 'use_flux_form      ' , use_flux_form
         print*, 'use_fct            ' , use_fct
         print*, 'use_tracer_evap    ' , use_tracer_evap
         print*, 'apply_sub_mp       ' , apply_sub_mp
         print*, 'alp1               ' , real(alp1               ,4)

         print*,'!---- couplings w/ other parameterizations'
         print*, 'lightning_diag     ' , lightning_diag
         print*, 'overshoot          ' , real(overshoot          ,4)
         print*, 'liq_ice_number_conc' , liq_ice_number_conc
         print*, 'use_gustiness      ' , use_gustiness

         print*, '!----misc controls'
         print*, 'output_sound       ' , output_sound
         print*, 'frac_modis         ' , frac_modis
         print*, 'use_cloud_dissipat ' , real(use_cloud_dissipation,4)
         print*, 'use_wetbulb        ' , use_wetbulb
         print*, 'clev_grid          ' , clev_grid
         print*, 'vert_discr         ' , vert_discr
         print*, 'max_tq_tend        ' , real(max_tq_tend,4)
         print*, 'cum_fadj_massflx   ' , real(cum_fadj_massflx     ,4)
         print*,"========================================================================"
         call flush(6)
      endif
   end subroutine readGFConvParNML
   !------------------------------------------------------------------------------------
   subroutine get_liq_ice_number_conc(itf,ktf,its,ite, kts,kte,ierr,ktop    &
                                     ,dtime,rho,outqc,tempco,outnliq,outnice)

      implicit none
      integer,   intent (in )  :: itf,ktf,its,ite,kts,kte
      real,      intent (in )  :: dtime

      integer, dimension (:)   ,intent (in )  :: ierr,ktop
      real,    dimension (:,:) ,intent (in )  :: outqc,tempco,rho
      real,    dimension (:,:) ,intent (out)  :: outnliq,outnice

      integer :: i,k,vtp_index
      real    :: fr,tqliq,tqice,dtinv

      real,   dimension (kts:kte,its:ite) :: nwfa   ! in the future set this as NCPL
      real,   dimension (kts:kte,its:ite) :: nifa   ! in the future set this as NCPI


      nwfa(:,:) =  99.e7  ! in the future set this as NCPL
      nifa(:,:) =  0.     ! in the future set this as NCPI
      dtinv    = 1./dtime
      
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)

         do k=kts,ktop(i)+1

            fr    = FractLiqF(tempco(k,i))
            tqliq = dtime * outqc(k,i)* rho(k,i) * fr
            tqice = dtime * outqc(k,i)* rho(k,i) * (1.-fr)

            outnice(k,i) = max(0.0,  makeIceNumber    (tqice, tempco(k,i))/rho(k,i))
            outnliq(k,i) = max(0.0,  makeDropletNumber(tqliq, nwfa  (k,i))/rho(k,i))

         enddo
         !-- convert in tendencies
         outnice = outnice * dtinv ! unit [1/s]
         outnliq = outnliq * dtinv ! unit [1/s]
      !--- for update
      ! nwfa =nwfa + outnliq*dtime
      ! nifa =nifa + outnice*dtime

      enddo

   end subroutine get_liq_ice_number_conc
   !+---+-----------------------------------------------------------------+
   !+---+-----------------------------------------------------------------+
   !----- module_mp_thompson_make_number_concentrations
   !- Developed by H. Barnes @ NOAA/OAR/ESRL/GSL Earth Prediction Advancement Division
   !-----------------------------------------------------------------------
   !      Q_ice              is cloud ice mixing ratio, units of kg/m3
   !      Q_cloud            is cloud water mixing ratio, units of kg/m3
   !      Q_rain             is rain mixing ratio, units of kg/m3
   !      temp               is air temperature in Kelvin
   !      makeIceNumber     is cloud droplet number mixing ratio, units of number per m3
   !      makeDropletNumber is rain number mixing ratio, units of number per kg of m3
   !      make_RainNumber    is rain number mixing ratio, units of number per kg of m3
   !      qnwfa              is number of water-friendly aerosols in number per kg

   !+---+-----------------------------------------------------------------+
   !+---+-----------------------------------------------------------------+
   elemental real function makeIceNumber (Q_ice, temp)

      implicit none
      real, parameter:: ice_density = 890.0
      real, parameter:: c_pi = 3.1415926536
      real, intent(in):: q_ice, temp
      integer :: idx_rei
      real :: corr, reice, deice
      double precision :: lambda

      !+---+-----------------------------------------------------------------+
      !..Table of lookup values of radiative effective radius of ice crystals
      !.. as a function of Temperature from -94C to 0C.  Taken from WRF RRTMG
      !.. radiation code where it is attributed to Jon Egill Kristjansson
      !.. and coauthors.
      !+---+-----------------------------------------------------------------+

      real, dimension(95), parameter:: retab = (/                       &
         5.92779, 6.26422, 6.61973, 6.99539, 7.39234,                   &
         7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930,          &
         10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319,          &
         15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955,          &
         20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125,          &
         27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943,          &
         31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601,          &
         34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078,          &
         38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635,          &
         42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221,          &
         50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898,          &
         65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833,          &
         93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424,          &
         124.954, 130.630, 136.457, 142.446, 148.608, 154.956,          &
         161.503, 168.262, 175.248, 182.473, 189.952, 197.699,          &
         205.728, 214.055, 222.694, 231.661, 240.971, 250.639 /)

      if (Q_ice == 0) then
         makeIceNumber = 0
         return
      end if

      !+---+-----------------------------------------------------------------+
      !..From the model 3D temperature field, subtract 179K for which
      !.. index value of retab as a start.  Value of corr is for
      !.. interpolating between neighboring values in the table.
      !+---+-----------------------------------------------------------------+

      idx_rei = int(temp-179.)
      idx_rei = min(max(idx_rei,1),94)
      corr = temp - int(temp)
      reice = retab(idx_rei)*(1.-corr) + retab(idx_rei+1)*corr
      deice = 2.*reice * 1.e-6

      !+---+-----------------------------------------------------------------+
      !..Now we have the final radiative effective size of ice (as function
      !.. of temperature only).  This size represents 3rd moment divided by
      !.. second moment of the ice size distribution, so we can compute a
      !.. number concentration from the mean size and mass mixing ratio.
      !.. The mean (radiative effective) diameter is 3./Slope for an inverse
      !.. exponential size distribution.  So, starting with slope, work
      !.. backwords to get number concentration.
      !+---+-----------------------------------------------------------------+

      lambda = 3.0 / deice
      makeIceNumber = Q_ice * lambda*lambda*lambda / (C_PI*Ice_density)

      !+---+-----------------------------------------------------------------+
      !..Example1: Common ice size coming from Thompson scheme is about 30 microns.
      !.. An example ice mixing ratio could be 0.001 g/kg for a temperature of -50C.
      !.. Remember to convert both into MKS units.  This gives N_ice=357652 per kg.
      !..Example2: Lower in atmosphere at T=-10C matching ~162 microns in retab,
      !.. and assuming we have 0.1 g/kg mixing ratio, then N_ice=28122 per kg,
      !.. which is 28 crystals per liter of air if the air density is 1.0.
      !+---+-----------------------------------------------------------------+

      return
   end function makeIceNumber

   !+---+-----------------------------------------------------------------+
   !+---+-----------------------------------------------------------------+

   elemental real function makeDropletNumber (Q_cloud, qnwfa)

      implicit none
      real, intent(in):: q_cloud, qnwfa
      real, parameter:: am_r = c_pi*1000./6.
      real, dimension(15), parameter:: g_ratio = (/24,60,120,210,336,   &
         504,720,990,1320,1716,2184,2730,3360,4080,4896/)
      double precision:: lambda, qnc
      real:: q_nwfa, x1, xDc
      integer:: nu_c

      if (Q_cloud == 0) then
         makeDropletNumber = 0
         return
      end if

      !+---+

      q_nwfa = MAX(99.e6, MIN(qnwfa,5.e10))
      nu_c = MAX(2, MIN(NINT(2.5e10/q_nwfa), 15))

      x1 = MAX(1., MIN(q_nwfa*1.e-9, 10.)) - 1.
      xDc = (30. - x1*20./9.) * 1.e-6

      lambda = (4.0D0 + nu_c) / xDc
      qnc = Q_cloud / g_ratio(nu_c) * lambda*lambda*lambda / am_r
      makeDropletNumber = SNGL(qnc)

      return
   end function makeDropletNumber

   !+---+-----------------------------------------------------------------+
   !+---+-----------------------------------------------------------------+

   elemental real function make_RainNumber (Q_rain, temp)

      implicit none

      real, intent(in):: q_rain, temp
      double precision:: lambda, n0, qnr
      real, parameter:: am_r = c_pi*1000./6.

      if (Q_rain == 0) then
         make_RainNumber = 0
         return
      end if

      !+---+-----------------------------------------------------------------+
      !.. Not thrilled with it, but set Y-intercept parameter to Marshal-Palmer value
      !.. that basically assumes melting snow becomes typical rain. However, for
      !.. -2C < T < 0C, make linear increase in exponent to attempt to keep
      !.. supercooled collision-coalescence (warm-rain) similar to drizzle rather
      !.. than bigger rain drops.  While this could also exist at T>0C, it is
      !.. more difficult to assume it directly from having mass and not number.
      !+---+-----------------------------------------------------------------+

      N0 = 8.e6

      if (temp .le. 271.15) then
         N0 = 8.e8
      elseif (temp .gt. 271.15 .and. temp.lt.273.15) then
         N0 = 8. * 10**(279.15-temp)
      endif

      lambda = SQRT(SQRT(N0*am_r*6.0/Q_rain))
      qnr = Q_rain / 6.0 * lambda*lambda*lambda / am_r
      make_RainNumber = SNGL(qnr)

      return
   end function make_RainNumber
   !----------------------------------------------------------------------
   subroutine rh_controls(mynum,itf,ktf,its,ite,kts,kte,ierr,t,po,qo,qeso,po_cup,cumulus  &
                       ,rh_entr_factor,rh_dicycle_fct,entr_rate_input, entr_rate,xlons,dt)
           
      implicit none
      character *(*), intent (in)             :: cumulus
      integer  ,intent (in )                  :: mynum,itf,ktf, its,ite, kts,kte
      integer  ,intent (in )  ,dimension(:)   :: ierr
      real     ,intent (in )                  :: entr_rate_input,dt
      real     ,intent (in )  ,dimension(:)   :: xlons
      real     ,intent (in )  ,dimension(:,:) :: t,po,qo,po_cup,qeso
      real     ,intent (inout),dimension(:)   :: entr_rate 
      real     ,intent (inout),dimension(:)   :: rh_entr_factor,rh_dicycle_fct

      !--locals
      integer :: i,k,vtp_index
      real*8  :: y,x
      real    :: dpg, trash, dayhr, p_start = 1000.
      real    ,dimension(its:ite) :: frh,dayhrr
      real    ,parameter :: ref_local_time = 8., ftun3=0.25
      logical ,parameter :: free_troposphere = .true. 
      
   if(moist_trigger > 1) then 
      
      !-- ave rh from 1000 -> 450 hPa, following Tian et al 2022 GRL.
      ! OR
      !-- ave rh from 800 -> 450 hPa accounts only for the ´free troposphere'
      if(free_troposphere) p_start = 800. 

      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         
         frh(i) = 0.
         trash  = 0.
         loopN:    do k=kts,ktf
            if( po(k,i) .gt. p_start .and. po(k,i) .lt. 450.) cycle loopN
            dpg=100.*(po_cup(k,i)-po_cup(k+1,i))/c_grav
            trash=trash+dpg
            frh(i)= frh(i) + (qo(k,i)/qeso(k,i))*dpg 
         enddo loopN
         
         !--average relative humidity
         frh(i) =   100.*frh(i)/(1.e-8+trash) ! no unit
         frh(i) = max(1., min(100., frh(i)))
         !
         !--- this is for the moist_trigger = 2
         x = dble(frh(i))          
         y = 9.192833D0 - 0.2529055D0*x + 0.002344832d0*x**2 &
           - 0.000007230408D0*x**3
         
         !--- if RH increases, rh_entr_factor decreases => entr decreases
         if(moist_trigger == 2) rh_entr_factor(i) = max(0.5,min(1.0,real(y,4)))
         
         !--- if RH increases, rh_entr_factor increases => entr increases
         if(moist_trigger == 3) rh_entr_factor(i) = max(0.1,min(1.0,frh(i)*0.01))

         entr_rate(i) = entr_rate(i) * rh_entr_factor(i)    
      enddo
   endif  
   if(rh_dicycle == 1) then 
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
          !--- local time
          dayhr  = (time_in / 3600. + float(itime1_in/100)+float(mod(itime1_in,100))/60.) 
          dayhrr(i) = mod(dayhr+xlons(i)/15.+24., 24.)

         if(abs ( dayhrr(i) - ref_local_time) < 1. .or. time_in < dt+1. ) &  
           !--- ftun3 controls the domain of the dicycle closure 
           !    ftun3 => 1 the closure is insensitive to the mean tropospheric RH
           !    ftun3 => 0 it depends on the RH, with mininum = ftun3
           rh_dicycle_fct(i) = ftun3 +(1. - ftun3)*&
                            (1.-(atan((frh(i)-55.)/10.)+atan(55.))/3.1016)
      enddo
   endif
      !-- code to test the atan function
      !  do i = 1 ,100 !relative humidity
      !      y = 0.25 +0.75*(1.-(atan((float(i)-60.)/10.)+atan(50.))/3.1016)/0.9523154
      !      print*,i,y
      !  enddo
      
      !print*,"FRH",maxval(frh),minval(frh),maxval(rh_dicycle_fct),minval(rh_dicycle_fct)
      !call flush(6)
   end subroutine rh_controls
   !------------------------------------------------------------------------------------
   real function coldPoolStart(CNV_TR)
      implicit none
      real,intent(in)  :: CNV_TR
      real             :: f1
      real, parameter  :: p_reduction_factor =  0.3 
      real, parameter  :: p_limiar = 1000.
   
      f1= min (mx_buoy2,CNV_TR)
      !--- f1 > 1000 => coldPoolStart ---> 1.0 
      !--- f1 < 1000 => coldPoolStart ---> 1.0 - reduction_factor
      coldPoolStart =  1.-min(p_reduction_factor,max(0.,(p_reduction_factor &
                        /(p_limiar-mx_buoy1))*(f1-mx_buoy1)))
   end function
   !------------------------------------------------------------------------------------
   subroutine sas_rainevap(cumulus,itf,ktf, its,ite, kts,kte,ierr,kbcon,ktop,dtime  &
                          ,sig,xmb,pre,xland,edto,pwo,pwdo,po_cup,qo,tn ,outt, outq,outbuoy &
                          ,rntot,delqev,delq2,qevap,rn,qcond, rainevap,qeso)
         implicit none

         character*(*)            ,intent(in)    :: cumulus
         integer                  ,intent(in)    :: itf,ktf, its,ite, kts,kte
         integer, dimension(:)    ,intent(in)    :: ierr,kbcon,ktop
         real,    dimension(:)    ,intent(inout) :: xland,edto,pre,xmb,sig

         real,    dimension(:,:)  ,intent(in)    :: po_cup,tn,qo,pwo,pwdo,qeso
         real,    dimension(:,:)  ,intent(inout) :: outt,outq,outbuoy
         real,                     intent(in)    :: dtime 
         real,    dimension(:)    ,intent(inout) :: rntot,delqev,delq2,qevap,rn&
                                                   ,qcond, rainevap

         logical,    dimension(its:ite)  :: flg
         integer :: i,k,vtp_index
         real :: dp, t1,q1, rain,evef
         real, parameter :: elocp=c_alvl/c_cp
         real, parameter :: el2orc=c_alvl*c_alvl/(c_rm*c_cp)
         real, parameter :: evfact=0.25 ! .4
         real, parameter :: evfactl=0.25 ! .2

         if(cumulus == 'shallow') return 
         rntot (:) = 0.
         delqev(:) = 0.
         delq2 (:) = 0.
         rn    (:) = 0.
         rntot (:) = 0.
         qevap (:) = 0.
         flg   (:) = .true.
         !
         do vtp_index = get_num_elements(vec_ok),1,-1
          i = get_data_value(vec_ok,vtp_index)
          rain  = 0.
            do k = ktop(i), kts, -1
                 rain     =  pwo(k,i) + edto(i) * pwdo(k,i)
                 rntot(i) = rntot(i) + rain * xmb(i)* .001 * dtime
            enddo
         enddo
         !
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            evef = edto(i) * evfact * sig(i)**2
            if(xland(i) > 0.5 .and. xland(i) < 1.5) evef = edto(i) * evfactl * sig(i)**2

            do k = ktop(i), 1, -1
               rain =  pwo(k,i) + edto(i) * pwdo(k,i)
               rn(i) = rn(i) + rain * xmb(i) * .001 * dtime

               if(flg(i))then
                     q1=qo(k,i)+(outq(k,i))*dtime
                     t1=tn(k,i)+(outt(k,i))*dtime
                     qcond(i) = evef * (q1 - qeso(k,i))            &
                              / (1. + el2orc * qeso(k,i) / t1**2)
                     dp = -100.*(po_cup(k+1,i)-po_cup(k,i))
                     if(rn(i) > 0. .and. qcond(i) < 0.) then
                       qevap(i) = -qcond(i) * (1.-exp(-.32*sqrt(dtime*rn(i))))
                       qevap(i) = min(qevap(i), rn(i)*1000.*c_grav/dp)
                       delq2(i) = delqev(i) + .001 * qevap(i) * dp / c_grav
                     endif
                     if(rn(i) > 0. .and. qcond(i) < 0. .and. delq2(i) > rntot(i)) then
                       qevap(i) = 1000.* c_grav * (rntot(i) - delqev(i)) / dp
                       flg(i) = .false.
                     endif
                     if(rn(i) > 0. .and. qevap(i) > 0.) then
                       outq(k,i) = outq(k,i) + qevap(i)/dtime
                       outt(k,i) = outt(k,i) - elocp * qevap(i)/dtime
                       rn  (i)   = max(0.,rn(i) - .001 * qevap(i) * dp / c_grav)
                       pre (i)   = pre(i) - qevap(i) * dp /c_grav/dtime
                       pre (i)   = max(pre(i),0.)
                       delqev(i) = delqev(i) + .001*dp*qevap(i)/c_grav
                     endif
               endif
            enddo
         enddo
         do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               rainevap(i)=delqev(i)
         enddo
   end subroutine sas_rainevap
   !------------------------------------------------------------------------------------
   subroutine coldPoolConvMem(cumulus,its, itf, kts, kte, ktf, aa2_, entr_rate_input, ztexec &
                            , zqexec, xland, po, buoy_exc, ierr, cap_max, wlpool, aa3_       &
                            , x_add_buoy, min_entr_rate, entr_rate)
      !! ## cold pool parameterization and convective memory
      !!
      !! Author: Saulo Freitas [SRF]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>
      !!
      !! Date: 20Fevereiro2023 13:42
      !!
      !! #####Version: version
      !!
      !! ---
      !! **Full description**:
      !!
      !! cold pool parameterization and convective memory
      !!
      !! ** History**:
      !!
      !! --- 
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!

      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_name = 'coldPoolConvMem' 
      !! subroutine name

      ! Variables (input, output, inout)
      integer, intent(in) :: its, itf, kts, kte, ktf
      character*(*), intent (in) :: cumulus
      real, intent(in)    :: aa2_(:)
      real, intent(in)    :: entr_rate_input
      real, intent(inout) :: ztexec(:)
      real, intent(inout) :: zqexec(:)
      real, intent(in)    :: xland(:)
      real, intent(in)    :: po(:,:)
      real, intent(in)    :: buoy_exc(:,:)

      integer, intent(inout) :: ierr(:)

      real, intent(inout) :: cap_max(:)
      real, intent(inout) :: wlpool(:)

      real, intent(out)   :: aa3_(:)
      real, intent(out)   :: x_add_buoy(:)
      real, intent(inout) :: min_entr_rate
      real, intent(inout) :: entr_rate(:)
      
      ! Local variables:
      integer :: i, vtp_index

      if(USE_MEMORY >= 0) then 
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)

            call get_cloud_bc(cumulus,ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i) &
                             ,buoy_exc (kts:kte,i),x_add_buoy (i),kts)
            !- for gate soundings
            !x_add_buoy (i) = float(JL*40)

            x_add_buoy (i) = min (x_add_buoy (i), 0.5*mx_buoy2)

            if(x_add_buoy (i) > mx_buoy1 ) then 
                zqexec(i) = 0.
                ztexec(i) = 0.
            endif
         enddo
      endif 

      !-- avoid extra-buoyancy where rained before
      if(USE_MEMORY == 1) then
               where(AA2_ > 10./3600.) 
                  x_add_buoy = 0.0
                  wlpool     = 0.0
               end where
      endif
      !-- change entrainment due improved cloud organization
      if(USE_MEMORY == 2 .or. USE_MEMORY == 222) then  
            !-- reduce entr rate, where cold pools exist
            do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               !entr_rate(i) = max(0.7, 1.-coldPoolStart_orig(x_add_buoy(i))) * entr_rate(i) 
                entr_rate(i) =             coldPoolStart     (x_add_buoy(i))  * entr_rate(i) 
            enddo
      endif

      !-- increase capmax 
      if(add_coldpool_trig >= 1)  then 
            do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
              !cap_max(i) = cap_max(i) + coldPoolStart_orig(x_add_buoy(i)) * 35.
               cap_max(i) = cap_max(i) + (1.-coldPoolStart(x_add_buoy(i))) * 35./0.3
            enddo
      endif
      !--- temporary for output
      AA3_(:)=x_add_buoy(:) 
            
      if(use_gustiness == 1 )  then 
                x_add_buoy(:) = 0.0
                AA3_(:)=c_cp*ztexec(:)+c_alvl*zqexec(:)
      endif
   end subroutine coldPoolConvMem
   !------------------------------------------------------------------------------------

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!-------- move the content below to another file named "Lib_ConvPar_GF_GEOS5.F90" --------

   !----------------------------------------------------------------------
   subroutine gen_random(its,ite,use_random_num,random)
      implicit none
      integer, intent(in)  :: its,ite
      real,    intent(in)  :: use_random_num
      real,    intent(out) :: random(:)

      !-local vars
      integer   :: i
      integer(8) :: iran, ranseed = 0

      call system_clock(ranseed)
      ranseed=mod(ranseed,2147483646)+1 !seed between 1 and 2^31-2
      iran = -ranseed

      !-- ran1 produces numbers between [ 0,1]
      !-- random        will be between [-1,1]
      !-- with use_random_num the interval will be [-use_random_num,+use_random_num]
      do i=its,ite
         random(i) = use_random_num * 2.0*(0.5-real(RAN1(IRAN),4))
        !print*,"ran=",i,random(i)
      enddo

      if(maxval(abs(random)) > use_random_num) stop "random > use_random_num"

   end subroutine 
   !----------------------------------------------------------------------
   real(8) function ran1(idum)

      ! This is contributed code standardized by Yong Wang
      ! Random number generator taken from Press et al.
      !
      ! Returns numbers in the range 0-->1
      !
      ! Their description...
      ! "Minimal" random number generator of Park and Miller with Bays-Durham
      ! shuffle and added safeguards. Returns a uniform deviate between 0.0 and 1.0
      ! (exclusive of the endpoint values). Call with idum a negative integer to
      ! initialize; thereafter, do not alter idum between successive calls in a
      ! sequence. RNMX should approximate the largest floating value that is less
      ! than 1.

      !use shr_kind_mod,    only: r8 => shr_kind_r8, i8 => shr_kind_i8
      implicit none
      integer(8), parameter:: ntab = 32,iq = 127773,ia = 16807,ir = 2836
      integer(8), parameter:: im = 2147483647,ndiv = 1+(im-1)/ntab
      real(8)   , parameter:: am = 1.0/im,eps = 1.2e-7,rnmx = 1.0-eps

      integer(8), intent(inout):: idum

      integer(8):: iy
      integer(8), dimension(ntab):: iv
      !save iv,iy
      data iv /ntab*0/, iy /0/
      integer(8):: j,k

      !
      if (idum.le.0.or.iy.eq.0) then
         ! initialize
         idum = max(-idum,1)
         do j = ntab+8,1,-1
            k = idum/iq
            idum = ia*(idum-k*iq)-ir*k
            if (idum.lt.0) idum = idum+im
            if (j.le.ntab) iv(j) = idum
         end do
         iy = iv(1)
      end if
      !
      k = idum/iq
      ! compute idum = mod(ia*idum,im) without overflows by schrage's method
      idum = ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum = idum+im
      ! j will be in the range 1-->ntab
      j = 1+iy/ndiv
      ! output previously stored value and refill the shuffle table
      iy = iv(j)
      iv(j) = idum
      ran1 = min(am*iy,rnmx)

   end function ran1
   !---------------------------------------------------------------------------------------------------
   subroutine get_kbmax_kdet_k22(cumulus,itf,ktf,its,ite,kts,kte,ierr,ierrc,z1 &
                                ,zo_cup,heo_cup,depth_min,z_detr,zkbmax&
                                ,kbmax,kdet,k22,kstabm)
      implicit none
      character*(*), intent (in)                   :: cumulus
      character*(*), intent (inout), dimension(:)  :: ierrc
      integer  ,intent (in   )                 :: itf,ktf, its,ite, kts,kte
      real     ,intent (inout)                 :: z_detr,zkbmax,depth_min
      integer  ,intent (inout), dimension(:)   :: ierr
      real     ,intent (in   ), dimension(:)   :: z1
      real     ,intent (in   ), dimension(:,:) :: zo_cup,heo_cup

      integer  ,intent (out  ), dimension(:)   :: kbmax,kdet,k22,kstabm
      
      !-local vars
      integer :: i,k,start_k22,vtp_index

       
      kbmax  (:) = kts
      kstabm (:) = ktf-1
      kdet   (:) = kts
      k22    (:) = kts
      !--- minimum depth (m), clouds must have
      !
      if(cumulus == 'deep'                         ) depth_min=1000.
      if(cumulus == 'mid' .or. cumulus == 'shallow') depth_min=500.
      !
      !--- max height(m) above ground where updraft air can originate
      !
      if(cumulus == 'deep'                         ) zkbmax=4000.
      if(cumulus == 'mid' .or. cumulus == 'shallow') zkbmax=3000.
      !
      !
      !--- depth(m) over which downdraft detrains all its mass
      !
      z_detr=1000.
      if(cumulus == 'deep'                         ) z_detr= 1000.
      if(cumulus == 'mid' .or. cumulus == 'shallow') z_detr= 300.

       do i=its,itf
         do k=kts,ktf
               if(zo_cup(k,i).gt.zkbmax+z1(i))then
                  kbmax(i)=k
                  exit
               endif
          end do
          !--- level where detrainment for downdraft starts
          do k=kts,ktf
               if(zo_cup(k,i).gt.z_detr+z1(i))then
                  kdet(i)=k
                  exit
               endif
          end do
      end do
      !
      !--- determine level with highest moist static energy content - k22
      if(cumulus == 'shallow') then
         start_k22 = 1
      else
         start_k22 = 2
      endif
      !
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         
         k22(i)=maxloc(heo_cup(start_k22:kbmax(i)+1,i),1)+start_k22-1
         k22(i)=max(k22(i),start_k22)
         if(cumulus == 'shallow') then
            k22(i)=min(2,k22(i))

            if(k22(i) > kbmax(i))then
               !is_removed  = remove(vec_ok, i)
               !is_inserted = insert_unique(vec_removed, i)
               ierr(i)=2
               ierrc(i)="could not find k22"
            endif

         else

            if(k22(i) > kbmax(i))then
               !- let's try k22=start_k22 for the cases k22>kbmax
               k22(i)= start_k22
               cycle
            endif

         endif
      enddo
   end subroutine get_kbmax_kdet_k22  
   !---------------------------------------------------------------------------------------------------
   subroutine get_edt(cumulus,itf,ktf,its,ite,kts,kte,xland,edtmin,edtmax &
                     ,MAX_EDT_OCEAN,MAX_EDT_LAND,c0_mid)
      implicit none
      character*(*), intent (in)             :: cumulus
      integer  ,intent (in )                 :: itf,ktf, its,ite, kts,kte
      real     ,intent (in )                 :: MAX_EDT_OCEAN,MAX_EDT_LAND,c0_mid
      real     ,intent (in ), dimension(:)   :: xland
      real     ,intent (out), dimension(:)   :: edtmin,edtmax
     
      !-local vars
      integer :: i,k
    !
      !--- max/min allowed value for epsilon (ratio downdraft base mass flux/updraft
      !    base mass flux
      !--  note : to make the evaporation stronger => increase "edtmin"
      if( cumulus == 'shallow') then
         edtmin(:)=0.0
         edtmax(:)=0.0
      endif
      if(cumulus == 'mid'     ) then
         do i=its,itf
            if(xland(i) > 0.99 ) then !- over water
               edtmin(i)=0.1
               edtmax(i)=MAX_EDT_OCEAN  !
            else!- over land
               edtmin(i)=0.1
               edtmax(i)=MAX_EDT_LAND  !
            endif
         enddo
         if(c0_mid < 1.e-8) edtmin(:)=0.0
      endif
      if(cumulus == 'deep'    ) then
         do i=its,itf
            if(xland(i) > 0.99 ) then !- over water
               edtmin(i)=0.2
               edtmax(i)=MAX_EDT_OCEAN  !
            else!- over land
               edtmin(i)=0.2
               edtmax(i)=MAX_EDT_LAND  !
            endif
         enddo
      endif
   end subroutine get_edt  
   !---------------------------------------------------------------------------------------------------
   subroutine get_lambdaU(cumulus,itf,ktf,its,ite,kts,kte,lambau_dp,lambau_dn &
                         ,lambau_deep,lambau_shdn,pgcon)

      implicit none
      character*(*), intent (in)             :: cumulus
      integer  ,intent (in )                 :: itf,ktf, its,ite, kts,kte
      real     ,intent (in )                 :: lambau_deep,lambau_shdn,pgcon
      real     ,intent (out), dimension(:)   :: lambau_dp, lambau_dn
      !-local vars
      integer :: i,k

      !--- lambda_U parameter for momentum transport
      !
      if(cumulus == 'deep'   ) then
         lambau_dp (:) = lambau_deep
         lambau_dn (:) = lambau_shdn
      endif
      if(cumulus == 'mid'    ) then
         lambau_dp (:) = lambau_shdn
         lambau_dn (:) = lambau_shdn
      endif
      if(cumulus == 'shallow') then
         lambau_dp (:) = lambau_shdn
         lambau_dn (:) = lambau_shdn
      endif

      if(pgcon .ne. 0.) then
         lambau_dp (:) = 0.
         lambau_dn (:) = 0.
      endif

   end subroutine get_lambdaU

   !---------------------------------------------------------------------------------------------------
   subroutine get_capmax(cumulus,itf,ktf,its,ite,kts,kte,cap_max_inc &
                        ,cap_max_increment,cap_max,cap_maxs,MOIST_TRIGGER)
      implicit none
      character*(*), intent (in)             :: cumulus
      integer  ,intent (in )                 :: itf,ktf, its,ite, kts,kte
      integer  ,intent (in )                 :: MOIST_TRIGGER
      real     ,intent (out)                 :: cap_max_inc,cap_maxs
      real     ,intent (out), dimension(:)   :: cap_max_increment,cap_max

      !-local vars
      integer :: i,k
      !
      !--- maximum depth (mb) of capping inversion (larger cap = no convection)
      if( MOIST_TRIGGER==0) then
         if(cumulus == 'deep'   ) then
            cap_max_inc=20.
         endif 
         if(cumulus == 'mid'    ) then
            cap_max_inc=10.
         endif 
         if(cumulus == 'shallow') then
            cap_max_inc=25.
         endif
      else
         if(cumulus == 'deep'   ) then
            cap_max_inc=90.
         endif
         if(cumulus == 'mid'    ) then
            cap_max_inc=90.
         endif 
         if(cumulus == 'shallow') then
            cap_max_inc=10.
         endif 
      endif
      cap_max_increment(:)= cap_max_inc
      cap_max          (:)= cap_maxs

   end subroutine get_capmax
   !---------------------------------------------------------------------------------------------------
   subroutine reset_1d(its,ite,ierrc,xland,xland1,aa0,aa1,aa2,aa3&
                   ,aa1_bl,aa1_fa,aa0_bl,q_adv,aa1_radpbl,aa1_adv,alpha_adv,cin1&
                   ,xk_x,edt,edto,tau_bl,q_wetbulb,t_wetbulb,tau_ecmwf,xf_dicycle&
                   ,x_add_buoy,xf_coldpool,wlpool_bcon,ke_gustfront,random,mbdt)
      implicit none       
      integer  ,intent (in )                      :: its,ite
      character*(*), intent (out), dimension(:)   :: ierrc
      real          ,intent (in ), dimension(:)   :: xland       
      real          ,intent (out), dimension(:)   ::   &       
         xland1,& 
         aa0   ,& 
         aa1   ,& 
         aa2   ,& 
         aa3   ,& 
         aa1_bl,& 
         aa1_fa,& 
         aa0_bl,& 
         q_adv ,& 
         aa1_radpbl,&
         aa1_adv   ,&
         alpha_adv ,&
         cin1      ,&
         xk_x      ,&
         edt       ,&
         edto      ,&
         tau_bl    ,&
         q_wetbulb ,&
         t_wetbulb ,&
         tau_ecmwf ,&
         xf_dicycle,&
         x_add_buoy  ,&
         xf_coldpool ,&
         wlpool_bcon ,&
         ke_gustfront,&
         random      ,&
         mbdt
         !---
         ierrc  (:) = "ierrtxt"
         
         xland1 (:) = xland(:) ! 1.
       
         aa0    (:) = 0.0
         aa1    (:) = 0.0
         aa2    (:) = 0.0
         aa3    (:) = 0.0
         aa1_bl (:) = 0.0
         aa1_fa (:) = 0.0
         aa0_bl (:) = 0.0
         q_adv  (:) = 0.0
         aa1_radpbl (:) = 0.0
         aa1_adv    (:) = 0.0
         alpha_adv  (:) = 0.0
         cin1       (:) = 0.0
         xk_x       (:) = 0.0
         edt        (:) = 0.0
         edto       (:) = 0.0
         tau_bl     (:) = 0.0
         q_wetbulb  (:) = 0.0
         t_wetbulb  (:) = 0.0
         tau_ecmwf  (:) = 0.0
         xf_dicycle (:) = 0.0
         x_add_buoy (:) = 0.0
         xf_coldpool(:) = 0.0
         wlpool_bcon(:) = 0.0
         ke_gustfront(:)= 0.0
         random      (:)= 0.0
         !
         !--- mbdt ~ xmb * timescale
         mbdt(:)= 0.1
        !mbdt(i)= 100.*(p_cup(kbcon(i),i)-p_cup(i,kbcon(i)+1))/(c_grav*dtime)
         !-- default flag for get_cloud_bc (only get_lcl uses it as 'true')
         ave_from_surface = .false. 

   end subroutine reset_1d
   !---------------------------------------------------------------------------------------------------
   subroutine reset_2d(its,ite,kts,kte,zo,z,xz,hcdo,cupclw,qrcdo&
                      ,hcot,xf_ens,pr_ens,evap_bcb,uc,vc,hc,hco,zuo,zdo,zenv)
      implicit none       
      integer  ,intent (in )                   :: kts,kte,its,ite
      real     ,intent (in ), dimension(:,:)   ::   zo   
      real     ,intent (out), dimension(:,:)   ::   & 
           z     ,&
           xz    ,&
           hcdo  ,&
           cupclw,&
           qrcdo ,&
           hcot  ,&
           xf_ens,&
           pr_ens,&
           evap_bcb, &
           uc  ,&
           vc  ,&
           hc  ,&
           hco ,&
           zuo ,&
           zdo ,&
           zenv
 
      z     (:,:) = zo(:,:)
      xz    (:,:) = zo(:,:)
      hcdo  (:,:) = 0.0
      uc    (:,:) = 0.0
      vc    (:,:) = 0.0
      hc    (:,:) = 0.0
      hco   (:,:) = 0.0
      zuo   (:,:) = 0.0
      zdo   (:,:) = 0.0
      zenv  (:,:) = 0.0
      cupclw(:,:) = 0.0
      qrcdo (:,:) = 0.0
      hcot  (:,:) = 0.0
      xf_ens(:,:) = 0.0
      pr_ens(:,:) = 0.0
      evap_bcb(:,:) = 0.0

   end subroutine reset_2d
   !---------------------------------------------------------------------------------------------------
   subroutine get_lcl(cumulus,convection_tracer,use_memory,ave_layer,its,ite,itf,kts,kte,ktf,ierr &
                     ,x_add_buoy,zqexec,ztexec,xland,po,t_cup,p_cup,z_cup,q_cup,k22,klcl)
      implicit none
      character*(*), intent (in)             :: cumulus
      integer  ,intent (in )                 :: itf,ktf, its,ite, kts,kte
      integer  ,intent (in )                 :: convection_tracer, use_memory
      integer  ,intent (in ), dimension(:)   :: k22,ierr
      real     ,intent (in )                 :: ave_layer
      real     ,intent (in ), dimension(:)   :: xland,x_add_buoy,zqexec,ztexec
      real     ,intent (in ), dimension(:,:) :: t_cup,p_cup,z_cup,po,q_cup

      integer  ,intent (out), dimension(:)   :: klcl

      real :: tlll,plll,rlll,tlcl,plcl,dzlcl,zlll,x_add
      integer :: i,k,vtp_index
      !--default value
      klcl(:) = k22(:)

      !-- tlll, rlll,plll - temp, water vapor and pressure of the source air parcel
      ! do i=its,itf  
      !   if(ierr(i) /= 0) cycle
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
           
         if( cumulus == 'deep' .and. convection_tracer == 1 .and. use_memory == 333) then
              x_add =  x_add_buoy(i)/c_alvl
         else
              x_add = max(0.,zqexec(i))
         endif

         call get_cloud_bc(cumulus,ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i),q_cup (kts:kte,i),rlll,k22(i),x_add)
           
         if( cumulus == 'deep' .and. convection_tracer == 1 .and. use_memory == 333) then
              x_add = 0.
         else
              x_add = max(0.,ztexec(i))
         endif

         call get_cloud_bc(cumulus,ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i),t_cup (kts:kte,i),tlll,k22(i),x_add)
         call get_cloud_bc(cumulus,ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i),p_cup (kts:kte,i),plll,k22(i))

         call calc_lcl(tlll,100.*plll,rlll,tlcl,plcl,dzlcl)

         !-get LCL
         if(dzlcl >= 0.) then ! LCL found (if dzlcl<0 => not found)
               call get_cloud_bc(cumulus,ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i),z_cup (kts:kte,i),zlll,k22(i))
               loop0:  do k=kts,ktf
                         if(z_cup(k,i).gt.zlll+dzlcl)then
                             klcl(i)=max(k,k22(i))
                             exit loop0
                         endif
               enddo loop0
               klcl(i)=min(klcl(i),ktf-4)
         endif
         
         !write(12,111)'MDlcl',tlcl,plcl,dzlcl,klcl(i),ierr(i)
         !111      format(1x,A5,3F10.2,2i4)
      enddo

   end subroutine get_lcl
   !---------------------------------------------------------------------------------------------------
   subroutine calc_lcl(t0,pp0,r0,tlcl,plcl,dzlcl)
      implicit none
      real,intent(in ) :: t0,pp0,r0
      real,intent(out) :: tlcl,plcl,dzlcl
      real :: ttd

      ttd=td(pp0,r0)
      tlcl=ttd-(0.001296*ttd+0.1963)*(t0-ttd)
      plcl=pp0*(tlcl/t0)**c_cpor
      dzlcl=127*(t0-ttd)
      if(dzlcl.le.0.)dzlcl=-999.
   
   end subroutine calc_lcl
   !---------------------------------------------------------------------------------------------------
   real function td(p,rs)
      implicit none
      real :: rr,rs,es,esln,p
      rr=rs+1e-8
      es=p*rr/(.622+rr)
      esln=log(es)
      td=(35.86*esln-4947.2325)/(esln-23.6837)
      return
   end function td
  !---------------------------------------------------------------------------------------------------
   subroutine get_cloud_bc(cumulus,ave_layer,kts,kte,ktf,xland,po,array,x_aver,k22,add,Tpert)
      implicit none
      character *(*)   ,intent (in) :: cumulus
      integer,intent(in)            :: kts,kte,ktf,k22
      real   ,intent(in)            :: array(:),po(:),xland,ave_layer
      real   ,optional ,intent(in)  :: add
      real   ,optional ,intent(in)  :: Tpert(:)
      real   ,intent(out)           :: x_aver
      !-- local vars
      integer                       :: i,local_order_aver,order_aver, i_beg,i_end,ic
      real,    parameter            :: frac_ave_layer_ocean= 0.3
      real                          :: count,dp,dp_layer,effec_frac,x_ave_layer

      if(ave_from_surface) then
            x_ave_layer = ave_layer
            i_beg = kts
            i_end = minloc(abs(po(kts:ktf)-(po(kts)-x_ave_layer)),1)
      else
            effec_frac  = (1.-xland) +xland*frac_ave_layer_ocean
            x_ave_layer = ave_layer*effec_frac
            i_beg = minloc(abs(po(kts:ktf)-(po(k22)+0.5*x_ave_layer)),1)
            i_end = minloc(abs(po(kts:ktf)-(po(k22)-0.5*x_ave_layer)),1)
      endif
      i_beg = min(ktf,max(i_beg,kts))
      i_end = min(ktf,max(i_end,kts))

      if(i_beg >= i_end) then
            x_aver   = array(k22)
            dp_layer = 0.
            ic       = i_beg

      else
            dp_layer = 1.e-06
            x_aver   = 0.
            ic       = 0
            do i = i_beg,ktf
               dp = -(po(i+1)-po(i))
               if(dp_layer + dp <= x_ave_layer)  then
                  dp_layer =  dp_layer  + dp
                  x_aver   =  x_aver    + array(i)*dp

               else
                  dp       =  x_ave_layer - dp_layer
                  dp_layer =  dp_layer    + dp
                  x_aver   =  x_aver      + array(i)*dp
                  exit
               endif
            enddo
            x_aver = x_aver/dp_layer
            ic  = max(i_beg,i)
      endif
      !-- this perturbation is included only for MSE
      if(present(Tpert)) x_aver = x_aver + c_cp*maxval(Tpert(i_beg:ic))  ! version 2 - maxval in the layer

      if(present(add)) x_aver = x_aver + add

   end subroutine get_cloud_bc
   !---------------------------------------------------------------------------------------------------
   
   subroutine get_lcl2(cumulus,ave_layer,its,ite,itf,kts,kte,ktf,ierr,zqexec,ztexec,xland &
                     ,po,t_cup,p_cup,z_cup,q_cup,k22,klcl,kpbl,psur,zlcl_sfc)
      implicit none
      character*(*), intent (in)             :: cumulus
      integer  ,intent (in )                 :: itf,ktf, its,ite, kts,kte
      integer  ,intent (in ), dimension(:)   :: k22,ierr,kpbl
      real     ,intent (in )                 :: ave_layer
      real     ,intent (in ), dimension(:)   :: xland,zqexec,ztexec,psur
      real     ,intent (in ), dimension(:,:) :: t_cup,p_cup,z_cup,po,q_cup

      integer  ,intent (out   ), dimension(:) :: klcl
      real     ,intent (inout ), dimension(:) :: zlcl_sfc


      real :: tlll,plll,rlll,tlcl,plcl,dzlcl,zlll,x_add,local_ave_layer,fr=0.
      integer :: i,k,vtp_index,k_beg

      ave_from_surface = .true. 
      !--default value
      klcl(:) = k22(:)
      x_add = 0.
      ! if(use_lcl_ctrl_entr == 1) fr = 0.1
      ! if(use_lcl_ctrl_entr == 2) fr = 0.5
      ! if(use_lcl_ctrl_entr == 3) fr = 0.2
      
      fr = 0.5
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         
         local_ave_layer = (psur(i) - po(kpbl(i),i)) * fr

         !x_add =  max(0.,zqexec(i))
         call get_cloud_bc(cumulus,local_ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i),q_cup (kts:kte,i),rlll,kts,x_add)
           
         !x_add = max(0.,ztexec(i))
         call get_cloud_bc(cumulus,local_ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i),t_cup (kts:kte,i),tlll,kts,x_add)
         
         call get_cloud_bc(cumulus,local_ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i),p_cup (kts:kte,i),plll,kts)

         call calc_lcl(tlll,100.*plll,rlll,tlcl,plcl,dzlcl)

         zlcl_sfc(i) = dzlcl
         
         !print*,"zlcl",dzlcl,zlcl_sfc(i),local_ave_layer
         

         !-get LCL
         ! if(dzlcl >= 0.) then ! LCL found (if dzlcl<0 => not found)
         !       call get_cloud_bc(cumulus,ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i),z_cup (kts:kte,i),zlll,k22(i))
         !       loop0:  do k=kts,ktf
         !                 if(z_cup(k,i).gt.zlll+dzlcl)then
         !                     klcl(i)=max(k,k22(i))
         !                     exit loop0
         !                 endif
         !       enddo loop0
         !       klcl(i)=min(klcl(i),ktf-4)
         ! endif
         
         !write(12,111)'MDlcl',tlcl,plcl,dzlcl,klcl(i),ierr(i)
         !111      format(1x,A5,3F10.2,2i4)
      enddo
      !
      !-- set back the default cloud_bc method
      ave_from_surface = .false. 

   end subroutine get_lcl2
   !---------------------------------------------------------------------------------------------------
   subroutine LCL_and_PBL_ctrl_on_entrainment(cumulus,its,ite,itf,kts,kte,ktf,min_entr_rate &
                                             ,entr_rate,zlcl_sfc,turb_len_scale,AA0_,AA1_,AA1_BL_)
      implicit none
      character *(*)   ,intent (in)          :: cumulus  
      integer  ,intent (in   )               :: its,ite,itf,kts,kte,ktf
      real     ,intent (in   )               :: min_entr_rate
      real     ,intent (in   ), dimension(:) :: zlcl_sfc
      real     ,intent (in   ), dimension(:,:) :: turb_len_scale
      real     ,intent (inout), dimension(:)   :: entr_rate,AA0_,AA1_,AA1_BL_
      !-- local vars
      integer         :: i
      real            :: k_lcl_entr,k_tls_entr,ax,bx,ax_lcl,ax_tls,tls_max
      bx = 1.0

      ax_tls = 0.5 ! 1 for NN, 0.5 MY
      ax_lcl = 2.0
      if(trim(cumulus) == 'deep') ax_lcl = 3.0
      do i=its,itf
         
         k_lcl_entr = min(bx, max(0.1, 1.e+3*0.2*ax_lcl/max(50., zlcl_sfc(i))))

         tls_max    = maxval(turb_len_scale(kts:ktf,i))
         k_tls_entr = min(bx, max(0.1, 1.e+3*0.2*ax_tls/max(1.0, tls_max)))
         
         if    (use_lcl_ctrl_entr == 1) then
            entr_rate(i) = k_tls_entr * entr_rate(i)
         
         elseif(use_lcl_ctrl_entr == 2) then
            entr_rate(i) = k_lcl_entr * entr_rate(i)
         
         elseif(use_lcl_ctrl_entr == 3) then
            entr_rate(i) = min(k_lcl_entr,k_tls_entr) * entr_rate(i)
         
         elseif(use_lcl_ctrl_entr == 4) then
            entr_rate(i) = 2./ ( 1.0/k_lcl_entr + 1.0/k_tls_entr ) * entr_rate(i)
         
         elseif(use_lcl_ctrl_entr == 5) then
            entr_rate(i) = max(k_lcl_entr,k_tls_entr) * entr_rate(i)
         
         endif

         if(trim(cumulus)=='deep') then 
          AA1_BL_(i)= entr_rate(i)/0.8 !1000.*2./ ( 1.0/k_lcl_entr + 1.0/k_tls_entr )
          AA0_(i)   = k_tls_entr*1000.
          AA1_(i)   = k_lcl_entr*1000.
         endif
      enddo
     
   end subroutine LCL_and_PBL_ctrl_on_entrainment
   
   !---------------------------------------------------------------------------------------------------
   subroutine set_entr_detr_rates(cumulus,its,ite,itf,kts,kte,ktf,ierr,klcl               &
                                 ,min_entr_rate,entr_rate,entr_rate_2d,cd,mentrd_rate,cdd &
                                 ,qo_cup, qeso_cup, cnvcf)
      implicit none
      character *(*)   ,intent (in) :: cumulus  
      integer  ,intent (in )                 :: its,ite,itf,kts,kte,ktf
      integer  ,intent (in ), dimension(:)   :: ierr, klcl
      real     ,intent (in )                 :: min_entr_rate
      real     ,intent (in ), dimension(:)   :: entr_rate
      real     ,intent (in ), dimension(:,:) :: qo_cup, qeso_cup,cnvcf

      real     ,intent (out), dimension(:)   :: mentrd_rate
      real     ,intent (out), dimension(:,:) :: entr_rate_2d ,cd, cdd

      !-- local vars
      integer :: i,k,vtp_index
      real    :: crh1 = 1.3, crh2 = 1.6
      real, dimension(kts:kte,its:ite) :: rh2d,entr_frh2d,detr_frh2d

      do i=its,itf
        entr_rate_2d(:,i) = entr_rate  (i)
        cd          (:,i) = entr_rate  (i)
      enddo
      
      if(use_rhu_ctrl_entr == 1) then

         !-- get the environmental relative humidity
         do vtp_index = get_num_elements(vec_ok),1,-1
           i = get_data_value(vec_ok,vtp_index)
           rh2d(:,i) = min(qo_cup(:,i)/qeso_cup(:,i),1.)
         enddo
      
         if(use_pass_cloudvol == 1 .or. use_pass_cloudvol == 12 .or. use_pass_cloudvol == 3) then
            crh1 = 1.6 ; crh2 = 1.9
            !-- The Role of Passive Cloud Volumes in the Transition 
            !-- From Shallow to Deep Atmospheric Convection - GRL 2023
            do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               !--effective RH =   environmental         + saturated 
               rh2d(:,i) = (1. - cnvcf(:,i)) * rh2d(:,i) + cnvcf(:,i)*1.0 
            enddo
         endif

         !--- get the relative humidity factor following Bechtold et al (2008)
         do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               entr_frh2d(:,i) = crh1 - rh2d(:,i)
               detr_frh2d(:,i) = crh2 - rh2d(:,i)
         enddo
      
      else ! no RH control on entr/detr rates

         entr_frh2d(:,:) = 1.0
         detr_frh2d(:,:) = 1.0
      endif

      !do i=its,itf
      !   if(ierr(i) /= 0) cycle
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)

         if(cumulus /= 'shallow') then

            do k=kts,ktf
               
               if(k >= klcl(i)) then
                     entr_rate_2d(k,i)=entr_rate(i)*entr_frh2d(k,i)*(qeso_cup(k,i)/qeso_cup(klcl(i),i))**1.25
               else
                     entr_rate_2d(k,i)=entr_rate(i)*entr_frh2d(k,i)
               endif
               cd(k,i)=0.75e-4*detr_frh2d(k,i)               
               entr_rate_2d(k,i) = max(entr_rate_2d(k,i),min_entr_rate)
            enddo

         else
            !-- for shallow 
            do k=kts,ktf
               entr_rate_2d(k,i)=entr_rate(i)*entr_frh2d(k,i)*max(min(1.,(qeso_cup(max(k,klcl(i)),i)&
                                                                         /qeso_cup(klcl(i),i))**1) ,0.1)
               cd(k,i)=0.75*entr_rate_2d(k,i)!+0.5e-3
            enddo
         endif

        !--- define entrainment/detrainment profiles for downdrafts
        mentrd_rate(i  ) =   entr_rate(i)*0.3
        cdd        (:,i) = mentrd_rate(i)
      enddo

   end subroutine set_entr_detr_rates
   !---------------------------------------------------------------------------------------------------
   subroutine set_scale_dep_factor(cumulus,use_scale_dep,downdraft,its,ite,itf,kts,kte,ktf&
                                  ,ierr,sig_factor,stochastic_sig,dx,sig,sigd)
      implicit none
      character*(*), intent (in)             :: cumulus
      integer  ,intent (in )                 :: itf,ktf, its,ite, kts,kte
      integer  ,intent (in )                 :: use_scale_dep,downdraft
      integer  ,intent (in ), dimension(:)   :: ierr
      real     ,intent (in )                 :: sig_factor
      real     ,intent (in ), dimension(:)   :: stochastic_sig,dx

      real     ,intent (out), dimension(:)   :: sig,sigd

      !-- local vars
      integer :: i,k,vtp_index
      
      !- scale dependence factor for updraft
      if( use_scale_dep == 0  .or.  cumulus == 'shallow') then
         sig(:) = 1.0
      else
         sig(:) = 0.0    
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)

            !-- for similar curve as in ifs/ec, use sig_factor = 0.22
            sig(i) = 1.0 - exp(-sig_factor*(dx(i)/1000.))

            if (stochastic_sig(i) /= 1.0) then
               sig(i) = sig(i)**(stochastic_sig(i)*max(0.9,0.9*sig(i)))
            end if
            sig(i)= max(0.001,min(sig(i),1.))
         end do
      endif
      !- scale dependence factor for downdraft
      if(downdraft == 0) then 
        sigd(:) = 0.0
      else
        sigd(:) = 1.0
      end if
   end subroutine set_scale_dep_factor
  !---------------------------------------------------------------------------------------------------
   subroutine get_cloud_top_by_inv_layers(cumulus,use_inv_layers,its,ite,itf,kts,kte,ktf &
                                         ,ierr,ierrc,psur,po_cup,tn_cup,zo_cup,mid,shal  &
                                         ,kbcon,ktop)
      implicit none
      character*(*), intent (in)             :: cumulus
      character*(*), intent (inout), dimension(:)  :: ierrc
      integer  ,intent (in )                 :: itf,ktf, its,ite, kts,kte,mid,shal
      logical  ,intent (in )                 :: use_inv_layers
      integer  ,intent (inout), dimension(:) :: ierr,kbcon
      integer  ,intent (inout), dimension(:) :: ktop
      real     ,intent (in   ), dimension(:) :: psur
      real     ,intent (in   ), dimension(:,:) :: tn_cup, zo_cup, po_cup


      !-- local vars
     integer :: i,k,vtp_index
     real :: min_deep_top, min_shall_top      
     integer, dimension (kts:kte,its:ite) ::  k_inv_layers
     real,    dimension (kts:kte,its:ite) ::  dtempdz

     if(cumulus == 'mid') then
         if(use_inv_layers) then
            !
            !--- get inversion layers
            call get_inversion_layers(cumulus,ierr,psur,po_cup,tn_cup,zo_cup,k_inv_layers,&
                                      dtempdz,itf,ktf,its,ite, kts,kte,mid,shal)
            
            do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               ktop(i) = min(ktop(i),k_inv_layers(mid,i))
             enddo
         endif
         !
         !-- check if ktop is above 450hpa layer for mid convection
         !
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            if(po_cup(ktop(i),i) < 450.) then
               ierr(i)=25 ; is_removed = remove(vec_ok, i)
               ierrc(i)='mid convection with cloud top above 450 hpa (~ 7km asl)'
            endif
         enddo
         !-- check if ktop is below 750hpa layer for mid convection
         !
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            if(po_cup(ktop(i),i) > 750.) then
               ierr(i)=55 ; is_removed = remove(vec_ok, i)
               ierrc(i)='ktop too low for mid'
            endif
         enddo
      endif

      if(cumulus == 'shallow') then
         if(use_inv_layers) then
            call get_inversion_layers(cumulus,ierr,psur,po_cup,tn_cup,zo_cup,k_inv_layers,&
                                      dtempdz,itf,ktf,its,ite, kts,kte,mid,shal)
            do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               ktop(i) = min(ktop(i),k_inv_layers(shal,i))
            enddo
         endif

         !--- check if ktop is above 700hpa layer for shallow convection
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            min_shall_top=600.
            !if(icumulus_gf(mid) == 0) min_shall_top=500.
            if(po_cup(ktop(i),i) < min_shall_top) then
               ierr(i)=26 ; is_removed = remove(vec_ok, i)
               ierrc(i)='shallow convection wit h cloud top above min_shall_top hpa' 
            endif
         enddo
      endif

      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         if(ktop(i) <= kbcon(i))then
            ierr(i)=5 ; is_removed = remove(vec_ok, i)
            ierrc(i)='ktop too small'
         endif
      enddo

      if(cumulus == 'deep') then
            min_deep_top=500.
            if(icumulus_gf(mid) == 0) min_deep_top=600.
            do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               if(po_cup(ktop(i),i) > min_deep_top) then
                  ierr(i)=55 ; is_removed = remove(vec_ok, i)
                  ierrc(i)='ktop too low for deep'
               endif
            enddo
      endif

   end subroutine get_cloud_top_by_inv_layers
   !------------------------------------------------------------------------------------
   subroutine get_inversion_layers(cumulus,ierr,psur,po_cup,to_cup,zo_cup,k_inv_layers &
                                  ,dtempdz,itf,ktf,its,ite, kts,kte,mid,shal)

      implicit none
      integer,                  intent (in )  :: itf,ktf,its,ite,kts,kte,mid,shal
      character *(*),           intent (in )  :: cumulus
      integer, dimension (:)  , intent (inout):: ierr
      real,    dimension (:)  , intent (in )  :: psur
      real,    dimension (:,:), intent (in )  :: po_cup,to_cup,zo_cup
      real,    dimension (:,:), intent (out)  :: dtempdz
      integer, dimension (:,:), intent (out)  :: k_inv_layers
      ! -- local vars
      integer, parameter :: ist = 3
      real :: dzm,delp
      integer:: i,k,ilev,kk,k1,ix,k800,k550,vtp_index
      integer, parameter :: extralayer = 0 !- makes plume top higher
      integer, dimension (kts:kte,its:ite)  :: local_k_inv_layers
      real   , dimension (kts:kte)          :: first_deriv,sec_deriv,distance
      
      !-initialize k_inv_layers as 1 (non-existent layer)_
      k_inv_layers= 1 !integer
      dtempdz     = 0.0
      first_deriv = 0.0
      sec_deriv   = 0.0
      distance    = 0.0
      local_k_inv_layers=1
      

      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)

         !- displacement from local surface pressure level
         delp=1000.-psur(i)

         !- 2nd method
         ! DO k = kts+1,ktf-2
         !dtempdz(k,i)=   ( deriv3(vtp_index,zo_cup(k,i), zo_cup(i,kts:ktf), to_cup(i,kts:ktf), ktf-kts+1, 1,ierr(i)))
         !!! sec_deriv(k)=abs( deriv3(vtp_index,zo_cup(k,i), zo_cup(i,kts:ktf), to_cup(i,kts:ktf), ktf-kts+1, 2))
           !print*,"2=",k,z_cup(k,i),dtempdz(k,i),
         ! ENDDO
         ! if(ierr(i) /= 0) cycle

         !-1st method
         !-  get the 1st derivative
         do k = kts+ist,ktf-ist
            first_deriv(k)  = (to_cup(k+1,i)-to_cup(k-1,i))/(zo_cup(k+1,i)-zo_cup(k-1,i))
         enddo
         first_deriv(kts      :kts+ist-1)  =first_deriv(kts+ist)
         first_deriv(ktf-ist+1:kte      )  =first_deriv(ktf-ist)

         dtempdz  (:,i)  = first_deriv(:)

         !-  get the abs of the 2nd derivative
         do k = kts+ist+1,ktf-ist-1
            sec_deriv(k)= abs((first_deriv(k+1)-first_deriv(k-1))/(zo_cup(k+1,i)-zo_cup(k-1,i)))
         enddo
         sec_deriv(kts    :kts+ist)=sec_deriv(kts+ist+1)
         sec_deriv(ktf-ist:kte    )=sec_deriv(ktf-ist-1)

         ix=1
         do kk=kts+ist+2,ktf-ist-2
            if(sec_deriv(kk) < sec_deriv(kk+1) .and. sec_deriv(kk) < sec_deriv(kk-1)) then
               local_k_inv_layers(ix,i)=kk
               ix  =ix+1
            endif
         enddo

         !- 2nd criteria
         do k=kts+ist+2,ktf-ist-2
            kk=local_k_inv_layers(k,i)
            if(kk == 1) cycle
            if( dtempdz(kk,i) < dtempdz(kk-1,i) .and. dtempdz(kk,i) < dtempdz(kk+1,i) ) then ! the layer is not a local maximum
               local_k_inv_layers(k,i) = 1
            endif
         enddo

      enddo


      !- find the locations of inversions around 800 and 550 hPa
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         !----------------
         !k_inv_layers(mid,i)=1
         !----------------
         !- displacement from local surface pressure level
         delp=1000.-psur(i)
         !----------------
         !k_inv_layers(mid,i)=21
         !cycle
         !----------------
         if( trim(cumulus)=='shallow') then
            !- now find the closest layers of 800 and 550 hPa.
            !- this is for shallow convection k800
            do k=kts,ktf
               distance(k)=abs(po_cup(local_k_inv_layers(k,i),i)-(750.-delp))
            enddo
            k800=minloc(abs(distance(kts:ktf)),1)

            if( k800 <= kts .or. k800 >= ktf - 4) then
               k_inv_layers(shal,i)= ktf
               !ierr(i)=8; is_removed = remove(vec_ok, i)
            else
               !-save k800 in the k_inv_layers array
               k_inv_layers(shal,i)=local_k_inv_layers(k800,i) +extralayer
            endif
 
         elseif( trim(cumulus)=='mid') then
            !- this is for mid/congestus convection k500
            do k=kts,ktf
               distance(k)=abs(po_cup(local_k_inv_layers(k,i),i)-(550.-delp))
            enddo
            k550=minloc(abs(distance(kts:ktf)),1)

            if( k550 <= kts .or. k550 >= ktf - 4) then
               k_inv_layers(mid,i) = 1
               ierr(i)=8 ; is_removed = remove(vec_ok, i)
            else
               !-save k550 in the k_inv_layers array
               k_inv_layers(mid,i)=local_k_inv_layers(k550,i) +extralayer
            endif
            if(  k_inv_layers(mid,i) <= kts .or. k_inv_layers(mid,i) >= ktf-4) then
               ierr(i)=12 ; is_removed = remove(vec_ok, i)
            endif
         else
            k_inv_layers(:,i)=1
            ierr(i)=88 ; is_removed = remove(vec_ok, i)
         endif

      enddo

   contains
      real function deriv3(vtp_index,xx, xi, yi, ni, m,ierr)
         !====================================================================
         ! Evaluate first- or second-order derivatives
         ! using three-point Lagrange interpolation
         ! written by: Alex Godunov (October 2009)
         !--------------------------------------------------------------------
         ! input ...
         ! xx    - the abscissa at which the interpolation is to be evaluated
         ! xi()  - the arrays of data abscissas
         ! yi()  - the arrays of data ordinates
         ! ni - size of the arrays xi() and yi()
         ! m  - order of a derivative (1 or 2)
         ! output ...
         ! deriv3  - interpolated value
         !============================================================================*/

         implicit none
         integer, parameter :: n=3
         real   , intent(in):: xx
         integer, intent(in):: ni, m,vtp_index
         real   , intent(in) :: xi(ni), yi(ni)
         real:: x(n), f(n)
         integer i, j, k, ix
         integer, intent(inout) :: ierr

         ! exit if too high-order derivative was needed,
         if (m > 2) then
            deriv3 = 0.0
            return
         end if

         ! if x is ouside the xi(1)-xi(ni) interval set deriv3=0.0
         if (xx < xi(1) .or. xx > xi(ni)) then
            deriv3 = 0.0
            ierr=8; is_removed = remove(vec_ok, vtp_index)
            !stop "problem with 2nd derivative-deriv3 routine"
            return
         endif

         ! a binary (bisectional) search to find i so that xi(i-1) < x < xi(i)
         i = 1
         j = ni
         do while (j > i+1)
            k = (i+j)/2
            if (xx < xi(k)) then
               j = k
            else
               i = k
            endif
         enddo

         ! shift i that will correspond to n-th order of interpolation
         ! the search point will be in the middle in x_i, x_i+1, x_i+2 ...
         i = i + 1 - n/2

         ! check boundaries: if i is ouside of the range [1, ... n] -> shift i
         if (i < 1) i=1
         if (i + n > ni) i=ni-n+1

         !  old output to test i
         !  write(*,100) xx, i
         !  100 format (f10.5, I5)

         ! just wanted to use index i
         ix = i

         ! initialization of f(n) and x(n)
         do i=1,n
            f(i) = yi(ix+i-1)
            x(i) = xi(ix+i-1)
         end do

         ! calculate the first-order derivative using Lagrange interpolation
         if (m == 1) then
            deriv3 =          (2.0*xx - (x(2)+x(3)))*f(1)/((x(1)-x(2))*(x(1)-x(3)))
            deriv3 = deriv3 + (2.0*xx - (x(1)+x(3)))*f(2)/((x(2)-x(1))*(x(2)-x(3)))
            deriv3 = deriv3 + (2.0*xx - (x(1)+x(2)))*f(3)/((x(3)-x(1))*(x(3)-x(2)))
         ! calculate the second-order derivative using Lagrange interpolation
         else
            deriv3 =          2.0*f(1)/((x(1)-x(2))*(x(1)-x(3)))
            deriv3 = deriv3 + 2.0*f(2)/((x(2)-x(1))*(x(2)-x(3)))
            deriv3 = deriv3 + 2.0*f(3)/((x(3)-x(1))*(x(3)-x(2)))
         endif
      end function deriv3
   end subroutine get_inversion_layers
   !------------------------------------------------------------------------------------
   subroutine get_1st_guess_MSE_profile(cumulus,ave_layer,its,ite,itf,kts,kte,ktf,ierr     &
                                       ,start_level,k22,ktop,hkbo,zqexec,ztexec,x_add_buoy &
                                       ,xland,zu,zuo,po,heo_cup,tpert,up_massentro,up_massdetro&
                                       ,heso_cup,heo,hco,cnvcf,heso)
      implicit none
      character*(*), intent (in)               :: cumulus
      integer  ,intent (in   )                 :: itf,ktf, its,ite, kts,kte
      integer  ,intent (in   ), dimension(:)   :: ierr,start_level,ktop,k22
      real     ,intent (in   )                 :: ave_layer
      real     ,intent (in   ), dimension(:)   :: zqexec,ztexec,x_add_buoy,xland
      real     ,intent (in   ), dimension(:,:) :: zu,zuo,up_massentro,up_massdetro,heso_cup,heo &
                                                 ,po,heo_cup,tpert,cnvcf,heso
      real     ,intent (inout), dimension(:)   :: hkbo
      real     ,intent (out  ), dimension(:,:) :: hco
      !-- local vars
      integer :: i,k,vtp_index
      real :: denom,x_add,h_env_eff

      !
      !--- 1st guess for moist static energy and dbyo (not including ice phase)
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         
         !-- update hkb/hkbo in case of k22 is redefined in 'cup_kbon'
         x_add = (c_alvl*zqexec(i)+c_cp*ztexec(i)) +  x_add_buoy(i)
         call get_cloud_bc(cumulus,ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i),heo_cup(kts:kte,i)&
                          ,hkbo(i),k22(i),x_add,Tpert(kts:kte,i))

         do k=kts,start_level(i)
               hco(k,i) = hkbo(i)
         enddo

         do k=start_level(i) + 1,ktop(i) + 1 
            denom=(zu(k-1,i)-.5*up_massdetro(k-1,i)+up_massentro(k-1,i)) + 1.e-12

            if(use_pass_cloudvol == 2 .or. use_pass_cloudvol == 12 .or. use_pass_cloudvol == 3) then
                  h_env_eff = (1.-cnvcf(k-1,i))*heo(k-1,i) + cnvcf(k-1,i)*heso(k-1,i)
                  
                  hco(k,i)=(hco(k-1,i)*zuo(k-1,i)-.5*up_massdetro(k-1,i)*hco(k-1,i) + &
                                                     up_massentro(k-1,i)*h_env_eff)/ denom
            
            else
                  hco(k,i)=(hco(k-1,i)*zuo(k-1,i)-.5*up_massdetro(k-1,i)*hco(k-1,i) + &
                                                     up_massentro(k-1,i)*heo(k-1,i))/ denom
            endif

            if(k==start_level(i)+1) then
                  !x_add = (c_alvl*zqexec(i)+c_cp*ztexec(i)) +  x_add_buoy(i)
                  hco(k,i)= hco(k,i) + x_add*up_massentro(k-1,i)/denom
            endif
         enddo
         do k=ktop(i)+2,ktf
            hco (k,i)=heso_cup(k,i)
         enddo
      enddo

   end subroutine get_1st_guess_MSE_profile
   !------------------------------------------------------------------------------------
   subroutine get_updraft_profile(cumulus,ave_layer,use_linear_subcl_mf,its,ite,itf,kts,kte,ktf     &
                                 ,ierr,start_level,k22,ktop,hkb,hkbo,zqexec,ztexec,x_add_buoy,xland &
                                 ,pgcon,po,u_cup,v_cup,he,heo,he_cup,heo_cup,heso_cup ,hes_cup      &
                                 ,p_liq_ice,qrco,Tpert ,us,vs,zuo,up_massdetru,up_massentru         &
                                 ,up_massdetro,up_massentro,hc,uc,vc,hco,cnvcf,hes,heso)
      implicit none
      character*(*), intent (in)             :: cumulus
      integer  ,intent (in )                 :: itf,ktf, its,ite, kts,kte,use_linear_subcl_mf
      integer  ,intent (in ), dimension(:)   :: ierr,start_level,ktop,k22
      real     ,intent (in )                 :: ave_layer, pgcon
      real     ,intent (in ), dimension(:)   :: zqexec,ztexec,x_add_buoy,xland
      real     ,intent (in ), dimension(:,:) :: up_massentro,up_massdetro,up_massdetru,up_massentru

      real     ,intent (in ), dimension(:,:) :: po,u_cup,v_cup,he,heo,he_cup,heo_cup,heso_cup,hes_cup&
                                               ,p_liq_ice,qrco,Tpert,us,vs,cnvcf,hes,heso

      real     ,intent (inout), dimension(:)   :: hkb,hkbo
      real     ,intent (inout), dimension(:,:) :: zuo,hc,uc,vc,hco
      !-- local vars
      integer :: i,k,vtp_index
      real    :: denom,denomU,x_add,h_env_eff

      !
      !-- updraft moist static energy + momentum budget
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         x_add = (c_alvl*zqexec(i)+c_cp*ztexec(i)) +  x_add_buoy(i)
         call get_cloud_bc(cumulus,ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i),he_cup (kts:kte,i),hkb (i)    ,k22(i),x_add,Tpert(kts:kte,i))
         call get_cloud_bc(cumulus,ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i),heo_cup(kts:kte,i),hkbo(i)    ,k22(i),x_add,Tpert(kts:kte,i))
         call get_cloud_bc(cumulus,ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i),u_cup  (kts:kte,i),uc  (kts,i),k22(i))
         call get_cloud_bc(cumulus,ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i),v_cup  (kts:kte,i),vc  (kts,i),k22(i))
      enddo

      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         do k=kts,start_level(i)
            hc (k,i) =hkb (i)
            hco(k,i) =hkbo(i)
            uc (k,i) =uc  (kts,i)
            vc (k,i) =vc  (kts,i)
         enddo
       enddo

      !--- option to produce linear fluxes in the sub-cloud layer.
      if(cumulus == 'shallow' .and. use_linear_subcl_mf == 1) then
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            call get_delmix(cumulus,kts,kte,ktf,xland(i),start_level(i),po(kts:kte,i) &
                           ,he_cup (kts:kte,i), hc (kts:kte,i))
            call get_delmix(cumulus,kts,kte,ktf,xland(i),start_level(i),po(kts:kte,i) &
                           ,heo_cup(kts:kte,i), hco(kts:kte,i))
         enddo
      endif

      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)

         do k=start_level(i) + 1 , ktop(i) + 1  
            denom =(zuo(k-1,i)-.5*up_massdetro(k-1,i)+up_massentro(k-1,i)) + 1.e-12
            denomU=(zuo(k-1,i)-.5*up_massdetru(k-1,i)+up_massentru(k-1,i)) + 1.e-12
           
            if(use_pass_cloudvol == 2 .or. use_pass_cloudvol == 12 .or. use_pass_cloudvol == 3) then
               h_env_eff = (1.-cnvcf(k-1,i))*he(k-1,i) + cnvcf(k-1,i)*hes(k-1,i)

               hc (k,i)=(hc (k-1,i)*zuo(k-1,i)-.5*up_massdetro(k-1,i)*hc (k-1,i) + &
                                                  up_massentro(k-1,i)*h_env_eff)/ denom
            
               h_env_eff = (1.-cnvcf(k-1,i))*heo(k-1,i) + cnvcf(k-1,i)*heso(k-1,i)
               
               hco(k,i)=(hco(k-1,i)*zuo(k-1,i)-.5*up_massdetro(k-1,i)*hco(k-1,i)+ &
                                                  up_massentro(k-1,i)*h_env_eff)/ denom
            else

               hc (k,i)=(hc (k-1,i)*zuo(k-1,i)-.5*up_massdetro(k-1,i)*hc (k-1,i) + &
                                                  up_massentro(k-1,i)*he (k-1,i))/ denom

               hco(k,i)=(hco(k-1,i)*zuo(k-1,i)-.5*up_massdetro(k-1,i)*hco(k-1,i)+ &
                                                  up_massentro(k-1,i)*heo(k-1,i))/ denom

            endif 


            if(k==start_level(i)+1) then
                  x_add = (c_alvl*zqexec(i)+c_cp*ztexec(i)) +  x_add_buoy(i)
                  hco(k,i)= hco(k,i) + x_add*up_massentro(k-1,i)/denom
                  hc (k,i)= hc (k,i) + x_add*up_massentro(k-1,i)/denom
            endif
           
            uc(k,i)=(uc(k-1,i)*zuo(k-1,i)-.5*up_massdetru(k-1,i)*uc(k-1,i)+     &
                                             up_massentru(k-1,i)*us(k-1,i)      &
                  -pgcon*.5*(zuo(k,i)+zuo(k-1,i))*(u_cup(k,i)-u_cup(k-1,i))) /  denomU

            vc(k,i)=(vc(k-1,i)*zuo(k-1,i)-.5*up_massdetru(k-1,i)*vc(k-1,i)+     &
                                             up_massentru(k-1,i)*vs(k-1,i)      &
                  -pgcon*.5*(zuo(k,i)+zuo(k-1,i))*(v_cup(k,i)-v_cup(k-1,i))) /  denomU

            !- includes glaciation effects on HC,HCO
            !                    ------ ice content --------
            hc (k,i)= hc (k,i)+(1.-p_liq_ice(k,i))*qrco(k,i)*c_xlf
            hco(k,i)= hco(k,i)+(1.-p_liq_ice(k,i))*qrco(k,i)*c_xlf

         enddo
         !
         do k=ktop(i)+2,ktf
            hc  (k,i)= hes_cup(k,i)
            uc  (k,i)=   u_cup(k,i)
            vc  (k,i)=   v_cup(k,i)
            hco (k,i)=heso_cup(k,i)
          enddo
      enddo

   end subroutine get_updraft_profile
   !----------------------------------------------------------------------
   subroutine get_delmix(cumulus,kts,kte,ktf,xland,subcl_level,po,ain,aout)
      implicit none
      character *(*)   ,intent (in) :: cumulus
      integer,intent(in)            :: kts,kte,ktf,subcl_level
      real   ,intent(in)            :: ain(kts:kte),po(kts:kte),xland
      real   ,intent(inout)         :: aout(kts:kte)

      !-- local var
      real :: x1,x2,dp,del,qc
      integer :: k

      !-
      qc = aout(kts)

      x2=0.
      x1=0.
      do k = kts,subcl_level
         dp = po(k+1)-po(k)
         x2 = x2 + dp
         x1 = x1 + dp*ain(k)
      enddo
      del = abs(qc-x1/(x2+1.e-12))
      aout(kts:subcl_level) =  ain(kts:subcl_level) + del

   end subroutine get_delmix
   !----------------------------------------------------------------------
   subroutine get_downdraft_profile(cumulus,ave_layer,pgcon,use_wetbulb &
                                   ,itf,ktf, its,ite, kts,kte           &
                                   ,ierrc,ierr,jmin,t_wetbulb,q_wetbulb &
                                   ,hc,heo,heso_cup,u_cup,v_cup,us,vs   &
                                   ,zdo,dd_massdetro,dd_massentro       &
                                   ,zo_cup,dd_massdetru,dd_massentru    &
                                   ,bud,hcdo,ucd,vcd,dbydo)

      implicit none
      character*(*), intent (in)             :: cumulus
      integer  ,intent (in )                 :: itf,ktf, its,ite, kts,kte,use_wetbulb
      integer  ,intent (in ), dimension(:)   :: jmin
      real     ,intent (in )                 :: ave_layer, pgcon
      real     ,intent (in ), dimension(:)   :: t_wetbulb,q_wetbulb
      real     ,intent (in ), dimension(:,:) :: hc,heo,heso_cup,u_cup,v_cup,us,vs,zo_cup,zdo&
                                               ,dd_massdetro,dd_massentro       &
                                               ,dd_massdetru,dd_massentru       
     
      character*(*), intent (inout), dimension(:)  :: ierrc
      integer  ,intent (inout), dimension(:)   :: ierr      
      real     ,intent (inout), dimension(:)   :: bud
      real     ,intent (inout), dimension(:,:) :: hcdo,ucd,vcd,dbydo
      !-- local vars
      integer :: i,ki,i_wb,vtp_index
      real :: denom,denomU,dzo
      !
      !--- downdraft moist static energy and momentum
      !
      do i=its,itf
         bud  (i)  = 0.
         hcdo (:,i)= heso_cup(:,i)
         ucd  (:,i)=    u_cup(:,i)
         vcd  (:,i)=    v_cup(:,i)
         dbydo(:,i)= 0.
      enddo
      if(cumulus == 'shallow') return

      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         i_wb=0
         if(use_wetbulb==1) then
            hcdo(jmin(i),i)=0.5*(c_cp*t_wetbulb(i)+c_alvl*q_wetbulb(i)+zo_cup(jmin(i),i)*c_grav + hc(jmin(i),i))
            i_wb=1
         endif

         dbydo(jmin(i),i)= hcdo (jmin(i),i)-heso_cup(jmin(i),i)
         bud  (i)        = dbydo(jmin(i),i)*(zo_cup(jmin(i)+1,i)-zo_cup(jmin(i),i))

         do ki=jmin(i) - i_wb ,kts,-1 
            denom = zdo(ki+1,i)-0.5*dd_massdetro(ki,i)+dd_massentro(ki,i)+1.e-12
            denomU= zdo(ki+1,i)-0.5*dd_massdetru(ki,i)+dd_massentru(ki,i)+1.e-12
            
            dzo=zo_cup(ki+1,i)-zo_cup(ki,i)

            ucd(ki,i)=(ucd(ki+1,i)*zdo(ki+1,i)-.5*dd_massdetru(ki,i)*ucd(ki+1,i)+ &
                                                  dd_massentru(ki,i)*us (ki,i)    &
                     -pgcon*zdo(ki+1,i)*(us(ki+1,i)-us(ki,i)))   /  denomU
            
            vcd(ki,i)=(vcd(ki+1,i)*zdo(ki+1,i)-.5*dd_massdetru(ki,i)*vcd(ki+1,i)+ &
                                                  dd_massentru(ki,i)*vs (ki,i)    &
                     -pgcon*zdo(ki+1,i)*(vs(ki+1,i)-vs(ki,i)))   /  denomU

            hcdo(ki,i)=(hcdo(ki+1,i)*zdo(ki+1,i)-.5*dd_massdetro(ki,i)*hcdo(ki+1,i)+     &
                                                    dd_massentro(ki,i)*heo(ki,i))  /denom

            dbydo(ki,i)=hcdo(ki,i)-heso_cup(ki,i)
            bud(i)=bud(i)+dbydo(ki,i)*dzo
            
         end do
         if(bud(i) > 0.)then
            ierr(i)=7 ; is_removed = remove(vec_ok, i)
            ierrc(i)='downdraft is not negatively buoyant '
         endif
      end do
   end subroutine get_downdraft_profile
   !----------------------------------------------------------------------
   subroutine check_mass_conserv(itf,ktf, its,ite, kts,kte,ierr,ktop,edto &
                                ,zo,zuo,zdo,dd_massdetro,dd_massentro     &
                                ,up_massentro,up_massdetro)
      implicit none
      integer  ,intent (in)                 :: itf,ktf, its,ite, kts,kte
      integer  ,intent (in), dimension(:)   :: ierr, ktop    
      real     ,intent (in), dimension(:)   :: edto
      real     ,intent (in), dimension(:,:) :: zo,zuo,zdo,dd_massdetro,dd_massentro &
                                              ,up_massentro,up_massdetro
      !-- local vars
      integer :: i,k,vtp_index
      real    :: detupk,subin,subdown,detdo,entdo,entup,detup,entupk,entdoj&
                ,totmas

      !-- only for debug
      !return
      !--- check mass conservation
      do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            do k=kts,ktop(i)
               ! these three are only used at or near mass detrainment and/or entrainment levels
               entupk=0.
               detupk=0.
               entdoj=0.
               ! detrainment and entrainment for downdrafts
               detdo=edto(i)*dd_massdetro(k,i)
               entdo=edto(i)*dd_massentro(k,i)
               ! entrainment/detrainment for updraft
               entup=up_massentro(k,i)
               detup=up_massdetro(k,i)
               ! subsidence by downdrafts only
               subin=-zdo(k+1,i)*edto(i)
               subdown=-zdo(k,i)*edto(i)
               if(k == ktop(i))then
                  detupk=zuo(ktop(i),i)  
                  subin=0.
                  subdown=0.
                  detdo=0.
                  entdo=0.
                  entup=0.
                  detup=0.
               endif
               totmas=subin-subdown+detup-entup-entdo+detdo-entupk-entdoj+detupk+zuo(k+1,i)-zuo(k,i)
               if(abs(totmas) > 1.e-6)then
                  write(6,*)'**mass cons: k,ktop,zo(ktop),totmas,subin,subdown,detup,entup,detdo,entdo,entupk,detupk'
                  write(6,123)'mass*1.e+6',k,ktop(i),zo(ktop(i),i),totmas*1.e+6,subin*1.e+6,subdown*1.e+6,detup*1.e+6,entup*1.e+6&
                              ,detdo*1.e+6,entdo*1.e+6,entupk*1.e+6,detupk*1.e+6
123               format(1X,A11,2i5,10e12.5)
               endif
            end do   ! k
      end do
   end subroutine check_mass_conserv
   !----------------------------------------------------------------------
   subroutine smooth_tend(use_smooth_tend,itf,ktf, its,ite, kts, kte, ierr, ktop &
                         ,po_cup ,dellu, dellv, dellah, dellat, dellaq, dellaqc)
      implicit none
      integer  ,intent (in )                   :: itf,ktf, its,ite, kts,kte
      integer  ,intent (in )                   :: use_smooth_tend
      integer  ,intent (in )  , dimension(:)   :: ierr,ktop
      real     ,intent (in )  , dimension(:,:) :: po_cup
      real     ,intent (inout), dimension(:,:) :: dellu, dellv, dellah, dellat &
                                                , dellaq, dellaqc
      !-- local vars
      integer :: i,k,kk,vtp_index
      real    :: dp, rcount
      real,    dimension (kts:kte,8)       ::  tend2d
      real,    dimension (8)               ::  tend1d

      do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               tend2d = 0.
 
               do k=kts,ktop(i)
                  rcount = 1.e-8
                  tend1d = 0.
                  do kk= max(kts,k-use_smooth_tend),min(ktop(i),k+use_smooth_tend)
                     dp         = (po_cup(kk,i)-po_cup(kk+1,i))
                     rcount     = rcount     +  dp
                     tend1d(1)  = tend1d(1)  +  dp* dellah  (kk,i)
                     tend1d(2)  = tend1d(2)  +  dp* dellaq  (kk,i)
                     tend1d(3)  = tend1d(3)  +  dp* dellaqc (kk,i)
                     tend1d(4)  = tend1d(4)  +  dp* dellu   (kk,i)
                     tend1d(5)  = tend1d(5)  +  dp* dellv   (kk,i)
                  end do
                  tend2d(k,1:5)  = tend1d(1:5) /rcount
               end do
               !--- get the final/smoother tendencies
               do k=kts,ktop(i)
                  dellah  (k,i) = tend2d(k,1)
                  dellaq  (k,i) = tend2d(k,2)
                  dellaqc (k,i) = tend2d(k,3)
                  dellu   (k,i) = tend2d(k,4)
                  dellv   (k,i) = tend2d(k,5)
               end do
      end do
   end  subroutine smooth_tend
   !----------------------------------------------------------------------
   subroutine get_env_change(coupl_mphysics,itf,ktf, its,ite, kts,kte,ierr &
                            ,heo,qo,tn,p_liq_ice,mbdt,dellu, dellv, dellah &
                            ,dellat, dellaq,dellaqc, subten_h, subten_q    &
                            ,subten_t ,xhe,xq,xt)
      implicit none
      logical  ,intent (in )                 :: coupl_mphysics 
      integer  ,intent (in )                 :: itf,ktf, its,ite, kts,kte
     
      integer  ,intent (in ), dimension(:)   :: ierr
      real     ,intent (in ), dimension(:)   :: mbdt
      real     ,intent (in ), dimension(:,:) :: heo,qo,tn,p_liq_ice,dellu   &
                                               ,dellv, dellah, subten_h, subten_q

      real     ,intent (inout), dimension(:,:) :: dellaqc, dellaq 

       real    ,intent (out  ), dimension(:,:) :: dellat,subten_t,xhe,xq,xt

      !-- local vars
      integer :: i,k,vtp_index

      dellat(:,:) = 0.0

     
      do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            do k=kts,ktf         !
               xhe(k,i)=(dellah(k,i)             )*mbdt(i)+heo(k,i)
               xq (k,i)=(dellaq(k,i)+dellaqc(k,i))*mbdt(i)+ qo(k,i)
               xq (k,i)=max(1.e-8,xq (k,i))
               
               !- do not feed dellat with dellaqc if the detrainment of liquid water
               !- will be used as a source for cloud microphysics
               if(coupl_mphysics) then
                  dellat(k,i)=(1./c_cp)*(dellah(k,i)-c_alvl*dellaq(k,i))

               else

                  dellat (k,i) = (1./c_cp)*( dellah(k,i) - c_alvl*(dellaq(k,i) + dellaqc(k,i))* & 
                                 (1.+(c_xlf/c_alvl)*(1.-p_liq_ice(k,i))))
                
                  !-adding dellaqc to dellaq:
                  dellaq (k,i)= dellaq(k,i)+dellaqc(k,i)
                  dellaqc(k,i)= 0.0
               endif

               xt(k,i)=((1./c_cp)*dellah(k,i)-(c_alvl/c_cp)*(dellaq(k,i)+dellaqc(k,i) * &
                        (1.+(c_xlf/c_alvl)*(1.-p_liq_ice(k,i)))))*mbdt(i) + tn(k,i)

               !--- temp tendency due to the environmental subsidence
               subten_t(k,i)=(1./c_cp)*(subten_h(k,i)-c_alvl*subten_q(k,i))

            enddo
      enddo
      !do vtp_index = get_num_elements(vec_ok),1,-1
      !      i = get_data_value(vec_ok,vtp_index)
      !      xhe(ktf,i)=heo(ktf,i)
      !      xq (ktf,i)=qo(ktf,i)
      !      xt (ktf,i)=tn(ktf,i)
      !      xq (ktf,i)=max(1.e-8,xq (ktf,i))
      !enddo
   end subroutine get_env_change
   !----------------------------------------------------------------------
   subroutine get_dellas(vert_discr,use_fct,alp1,dtime,itf,ktf, its,ite, kts,kte  &
                        ,ierr,ktop,edto                                           &
                        ,po_cup,hco,heo,heo_cup,hcdo,qo,qco,qrco,qrcdo,qo_cup     &
                        ,qcdo,pwo,pwdo,p_liq_ice,melting,dbydo                    &
                        ,us,vs,uc,vc,u_cup,v_cup,ucd,vcd                          &
                        ,zuo,zdo,zenv,dd_massdetro,dd_massentro                   &
                        ,up_massentro,up_massdetro                                &
                        ,dellu, dellv, dellah, dellat, dellaq                     &
                        ,dellaqc, dellabuoy, subten_H, subten_Q, subten_T)

      implicit none
      
      integer  ,intent (in )                 :: itf,ktf, its,ite, kts,kte
      integer  ,intent (in )                 :: vert_discr,use_fct
      integer  ,intent (in ), dimension(:)   :: ierr,ktop
   
      real     ,intent (in )                 :: alp1,dtime
      real     ,intent (in ), dimension(:)   :: edto 
      real     ,intent (in ), dimension(:,:) :: po_cup,hco,heo,heo_cup,hcdo,qo,qco,qrco,qrcdo &
                                               ,qo_cup,qcdo,pwo,pwdo,p_liq_ice,melting,dbydo  &
                                               ,us,vs,uc,vc,u_cup,v_cup,ucd,vcd               &
                                               ,zuo,zdo,zenv,dd_massdetro,dd_massentro        &
                                               ,up_massentro,up_massdetro   
      
      real     ,intent (out ), dimension(:,:) :: dellu, dellv, dellah, dellat, dellaq         &
                                               , dellaqc, dellabuoy, subten_H, subten_Q, subten_T

      !-- local vars
      integer :: i,k,vtp_index
      real :: trash,trash2,dp, C_up, E_dn,G_rain,dtime_max,detup,beta1,alp0
      real, dimension (kts:kte)             ::  aa,bb,cc,ddu,ddv,ddh,ddq,fp,fm
      real, dimension (kts:kte,its:ite)     ::  massflx
      real, dimension (kts:kte,1)           ::  trcflx_in,sub_tend
      
      !--- change per unit mass that a model cloud would modify the environment
      !
      !--- 1. in bottom layer
      !
      dellu     =0.
      dellv     =0.
      dellah    =0.
      dellat    =0.
      dellaq    =0.
      dellaqc   =0.
      dellabuoy =0.
      subten_H  =0.
      subten_Q  =0.
      subten_T  =0.

      !
      !----------------------------------------------  cloud level ktop
      !
      !- - - - - - - - - - - - - - - - - - - - - - - - model level ktop-1
      !      .               .                 .
      !      .               .                 .
      !      .               .                 .
      !      .               .                 .
      !      .               .                 .
      !      .               .                 .
      !
      !----------------------------------------------  cloud level k+2
      !
      !- - - - - - - - - - - - - - - - - - - - - - - - model level k+1
      !
      !----------------------------------------------  cloud level k+1
      !
      !- - - - - - - - - - - - - - - - - - - - - - - - model level k
      !
      !----------------------------------------------  cloud level k
      !
      !      .               .                 .
      !      .               .                 .
      !      .               .                 .
      !      .               .                 .
      !      .               .                 .
      !      .               .                 .
      !      .               .                 .
      !      .               .                 .
      !      .                .                 .
      !      .               .                 .
      !
      !----------------------------------------------  cloud level 3  _cup
      !
      !- - - - - - - - - - - - - - - - - - - - - - - - model level 2
      !
      !----------------------------------------------  cloud level 2  _cup
      !
      !- - - - - - - - - - - - - - - - - - - - - - - - model level 1

      !
      !------------------------------------------------------------------------------------------
      if(VERT_DISCR == 0) then

      do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            do k=kts,ktop(i)

               dp=100.*(po_cup(k,i)-po_cup(k+1,i))

               dellu(k,i) =-(zuo(k+1,i)*(uc (k+1,i)-u_cup(k+1,i) ) -             &
                             zuo(k,i  )*(uc (k,i  )-u_cup(k,i  ) ) )*c_grav/dp        &
                           +(zdo(k+1,i)*(ucd(k+1,i)-u_cup(k+1,i)) -              &
                             zdo(k,i  )*(ucd(k,i  )-u_cup(k,i  )) )*c_grav/dp*edto(i)

               dellv(k,i) =-(zuo(k+1,i)*(vc (k+1,i)-v_cup(k+1,i) ) -             &
                             zuo(k,i  )*(vc (k,i  )-v_cup(k,i  ) ) )*c_grav/dp        &
                           +(zdo(k+1,i)*(vcd(k+1,i)-v_cup(k+1,i) ) -             &
                             zdo(k,i  )*(vcd(k,i  )-v_cup(k,i  ) ) )*c_grav/dp*edto(i)

            enddo   ! k
         enddo

         do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               trash  = 0.0
               trash2 = 0.0
               do k=kts,ktop(i)

                  dp=100.*(po_cup(k,i)-po_cup(k+1,i))

                  dellah(k,i) =-(zuo(k+1,i)*(hco (k+1,i)-heo_cup(k+1,i) ) -                 &
                                 zuo(k,i  )*(hco (k,i  )-heo_cup(k,i  ) ) )*c_grav/dp            &
                               +(zdo(k+1,i)*(hcdo(k+1,i)-heo_cup(k+1,i) ) -                 &
                                 zdo(k,i  )*(hcdo(k,i  )-heo_cup(k,i  ) ) )*c_grav/dp*edto(i)

                  !---meltglac-------------------------------------------------
                  dellah(k,i) = dellah(k,i) + c_xlf*((1.-p_liq_ice(k,i))*0.5*(qrco(k+1,i)+qrco(k,i)) &
                              - melting(k,i))*c_grav/dp

                  !-- for output only
                  subten_H(k,i) = -(zuo(k+1,i)*(-heo_cup(k+1,i)) - zuo(k,i)*(-heo_cup(k,i)))*c_grav/dp       &
                                  +(zdo(k+1,i)*(-heo_cup(k+1,i)) - zdo(k,i)*(-heo_cup(k,i)))*c_grav/dp*edto(i)

                   !- check H conservation
                  trash2 = trash2+ (dellah(k,i))*dp/c_grav

                   !-- take out cloud liquid/ice water for detrainment
                  detup=up_massdetro(k,i)
                  dellaqc(k,i) = detup*0.5*(qrco(k+1,i)+qrco(k,i)) *c_grav/dp
                  !---
                  G_rain=  0.5*(pwo (k,i)+pwo (k+1,i))*c_grav/dp
                  E_dn  = -0.5*(pwdo(k,i)+pwdo(k+1,i))*c_grav/dp*edto(i) ! pwdo < 0 and E_dn must > 0
                  !
                  !print*,"eva=",k,pwdo(k,i),E_dn,zdo(k,i  ),G_rain
                  !
                  !-- condensation source term = detrained + flux divergence of
                  !-- cloud liquid/ice water (qrco) + converted to rain

                  C_up = dellaqc(k,i)+(zuo(k+1,i)* qrco(k+1,i) -       &
                     zuo(k,i  )* qrco(k,i  )  )*c_grav/dp + G_rain

                  !-- water vapor budget
                  !-- = flux divergence z*(Q_c - Q_env)_up_and_down &
                  !--   - condensation term + evaporation
                  dellaq(k,i) =-(zuo(k+1,i)*(qco (k+1,i)-qo_cup(k+1,i) ) -                 &
                                 zuo(k,i  )*(qco (k,i  )-qo_cup(k,i  ) ) )*c_grav/dp            &
                               +(zdo(k+1,i)*(qcdo(k+1,i)-qo_cup(k+1,i) ) -                 &
                                 zdo(k,i  )*(qcdo(k,i  )-qo_cup(k,i  ) ) )*c_grav/dp*edto(i)    &
                               - C_up + E_dn

                  !-- for output only
                  subten_Q(k,i) =-(zuo(k+1,i)*(-qo_cup(k+1,i)) - zuo(k,i)*(-qo_cup(k,i)))*c_grav/dp       &
                                 +(zdo(k+1,i)*(-qo_cup(k+1,i)) - zdo(k,i)*(-qo_cup(k,i)))*c_grav/dp*edto(i)

                  !- check water conservation liq+condensed (including rainfall)
                  trash= trash+ (dellaq(k,i)+dellaqc(k,i)+ G_rain-E_dn)*dp/c_grav

                  !---
                  dellabuoy(k,i) = edto(i)*dd_massdetro(k,i)*0.5*(dbydo(k+1,i)+dbydo(k,i))*c_grav/dp
                  !---
               enddo   ! k
         enddo

      elseif(VERT_DISCR == 1) then

         !---- convective transport of momentum
         if(alp1 == 0.) then !-- fully time explicit
            do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               do k=kts,ktop(i)
                  dp=100.*(po_cup(k,i)-po_cup(k+1,i))

                  dellu(k,i) =-(zuo(k+1,i)*(uc (k+1,i)-u_cup(k+1,i) ) -             &
                                zuo(k,i  )*(uc (k,i  )-u_cup(k,i  ) ) )*c_grav/dp        &
                              +(zdo(k+1,i)*(ucd(k+1,i)-u_cup(k+1,i)) -              &
                                zdo(k,i  )*(ucd(k,i  )-u_cup(k,i  )) )*c_grav/dp*edto(i)

                  dellv(k,i) =-(zuo(k+1,i)*(vc (k+1,i)-v_cup(k+1,i) ) -             &
                                zuo(k,i  )*(vc (k,i  )-v_cup(k,i  ) ) )*c_grav/dp        &
                              +(zdo(k+1,i)*(vcd(k+1,i)-v_cup(k+1,i) ) -             &
                                zdo(k,i  )*(vcd(k,i  )-v_cup(k,i  ) ) )*c_grav/dp*edto(i)
               enddo   ! k
            enddo

         elseif (alp1 > 0.) then              !-- time alp0*explict + alp1*implicit + upstream

            alp0=1.-alp1
            do vtp_index = get_num_elements(vec_ok),1,-1
               i = get_data_value(vec_ok,vtp_index)
               do k=kts,ktop(i)+1
                  fp(k) = 0.5*(zenv(k,i)+abs(zenv(k,i)))
                  fm(k) = 0.5*(zenv(k,i)-abs(zenv(k,i)))
               enddo

               do k=kts,ktop(i)
                  dp=100.*(po_cup(k,i)-po_cup(k+1,i))

                  beta1 = dtime*c_grav/dp
                  aa(k) =    alp1*beta1*fm(k)
                  bb(k) = 1.+alp1*beta1*(fp(k)-fm(k+1))
                  cc(k) =   -alp1*beta1*fp(k+1)

                  ddu(k) = us(k,i)-( zuo(k+1,i)*uc (k+1,i)-zuo(k,i)*uc (k,i) )*beta1 + &
                                   ( zdo(k+1,i)*ucd(k+1,i)-zdo(k,i)*ucd(k,i) )*beta1*edto(i)

                  ddu(k) = ddu(k) + alp0*beta1*(-fm(k)*us(max(kts,k-1),i) +(fm(k+1)-fp(k))*us(k,i) +fp(k+1)*us(k+1,i))


                  ddv(k) = vs(k,i)-( zuo(k+1,i)*vc (k+1,i)-zuo(k,i)*vc (k,i) )*beta1 + &
                                   ( zdo(k+1,i)*vcd(k+1,i)-zdo(k,i)*vcd(k,i) )*beta1*edto(i)

                  ddv(k) = ddv(k) + alp0*beta1*(-fm(k)*vs(max(kts,k-1),i) +(fm(k+1)-fp(k))*vs(k,i) +fp(k+1)*vs(k+1,i))

                 !print*,"trX4=",k,aa(k),bb(k),cc(k)!, 1.+alp1*beta1*zenv(k,i  ), -alp1*beta1*zenv(k+1,i)
               enddo
               call tridiag (ktop(i),aa (kts:ktop(i)),bb (kts:ktop(i)),cc (kts:ktop(i)),ddu (kts:ktop(i)))
               dellu(kts:ktop(i),i)=(ddu(kts:ktop(i))-us(kts:ktop(i),i))/dtime

               call tridiag (ktop(i),aa (kts:ktop(i)),bb (kts:ktop(i)),cc (kts:ktop(i)),ddv (kts:ktop(i)))
               dellv(kts:ktop(i),i)=(ddv(kts:ktop(i))-vs(kts:ktop(i),i))/dtime
            enddo
         endif


         !--- convective transport of MSE and Q/Qc
         !if(USE_FLUX_FORM == 1) then
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            !--- moist static energy : flux form + source/sink terms + time explicit
            !
            !   if(use_fct == 0 .or. adjustl(cumulus) == 'shallow') then
            if(use_fct == 0 ) then

               do k=kts,ktop(i)
                  dp=100.*(po_cup(k,i)-po_cup(k+1,i))
                  dellah(k,i) =-(zuo(k+1,i)*(hco (k+1,i)-heo_cup(k+1,i) ) -          &
                                 zuo(k,i  )*(hco (k,i  )-heo_cup(k,i  ) ) )*c_grav/dp     &
                               +(zdo(k+1,i)*(hcdo(k+1,i)-heo_cup(k+1,i) ) -          &
                                 zdo(k,i  )*(hcdo(k,i  )-heo_cup(k,i  ) ) )*c_grav/dp*edto(i)

                  dellah(k,i) = dellah(k,i) + c_xlf*((1.-p_liq_ice(k,i))* &
                              0.5*(qrco(k+1,i)+qrco(k,i)) - melting(k,i))*c_grav/dp

                  !--- for output only
                  subten_H(k,i) = -(zuo(k+1,i)*(-heo_cup(k+1,i)) - zuo(k,i)*(-heo_cup(k,i)))*c_grav/dp       &
                                  +(zdo(k+1,i)*(-heo_cup(k+1,i)) - zdo(k,i)*(-heo_cup(k,i)))*c_grav/dp*edto(i)
               enddo   ! k

            else

               !-- FCT scheme for the subsidence transport: d(M_env*S_env)/dz
               sub_tend (:,1) = 0. ! dummy array
               trcflx_in(:,1) = 0. ! dummy array
               massflx  (:,i) = 0.
               dtime_max      = dtime

               do k=kts,ktop(i)
                  dp=100.*(po_cup(k,i)-po_cup(k+1,i))
                  trcflx_in (k,1) =-(zuo(k,i)  -edto(i)*zdo(k,i))*heo_cup(k,i) !* xmb(i)
                  massflx   (k,i) =-(zuo(k,i)  -edto(i)*zdo(k,i))        !* xmb(i)
                  dtime_max=min(dtime_max,.5*dp)
               enddo
               call fct1d3 (ktop(i),kte,dtime_max,po_cup(:,i),heo(:,i),massflx(:,i),trcflx_in(:,1),sub_tend(:,1))

               do k=kts,ktop(i)
                  dp=100.*(po_cup(k,i)-po_cup(k+1,i))
                  dellah(k,i) =-( zuo(k+1,i)*hco (k+1,i) - zuo(k,i)*hco (k,i) )*c_grav/dp      &
                               +( zdo(k+1,i)*hcdo(k+1,i) - zdo(k,i)*hcdo(k,i) )*c_grav/dp*edto(i)

                  dellah(k,i) = dellah(k,i) + c_xlf*((1.-p_liq_ice(k,i))* &
                              0.5*(qrco(k+1,i)+qrco(k,i)) - melting(k,i))*c_grav/dp
                  !- update with subsidence term from the FCT scheme
                  dellah(k,i) = dellah(k,i) + sub_tend(k,1)
                  !--- for output only
                  subten_H(k,i) = sub_tend(k,1)
               enddo   ! k
            endif
         enddo

         !     elseif(USE_FLUX_FORM == 2) THEN
         !
         !        !- flux form + source/sink terms + time explicit + upstream with anti-diffusion step (Smolarkiewicz 1983)
         !        alp0=1.
         !        do i=its,itf
         !          if(ierr(i) /= 0) cycle
         !    do istep=1,-1, -2
         !
         !      if(istep == 1) then
         !         ddu(:) = heo(:,i)
         !     do k=kts,ktop(i)+1
         !       fp(k) = 0.5*(zenv(k,i)+abs(zenv(k,i)))
         !       fm(k) = 0.5*(zenv(k,i)-abs(zenv(k,i)))
         !     enddo
         !       else
         !         ddu(kts:ktop(i)+1) = heo(i,kts:ktop(i)+1) + dellah(i,kts:ktop(i)+1)*dtime
         !         zenv_diff(1,kts) = 0.
         !         do k=kts,ktop(i)+1
         !        dp = 100.*(po_cup(k,i)-po_cup(k+1,i))
         !        zenv_diff (1,k+1) = 1.06* ( dp*abs(zenv(k+1,i))/g - dtime*zenv(k+1,i)**2 )/dp/c_grav &
         !                * (ddu(k+1) - ddu(k)) /(ddu(k+1) + ddu(k) + 1.e-16)
         !         enddo
         !         do k=kts,ktop(i)+1
         !        fp(k) = 0.5*(zenv_diff(1,k)+abs(zenv_diff(1,k)))
         !        fm(k) = 0.5*(zenv_diff(1,k)-abs(zenv_diff(1,k)))
         !         enddo
         !       endif
         !       do k=kts,ktop(i)
         !         dp=100.*(po_cup(k,i)-po_cup(k+1,i))
         !         beta1 = dtime*c_grav/dp
         !         ddh(k) = ddu(k) + alp0*beta1*( -fm(k)*ddu(max(kts,k-1)) + (fm(k+1)-fp(k))*ddu(k) + fp(k+1)*ddu(k+1) )
         !       enddo
         !
         !       dellah(i,kts:ktop(i)+1)=(ddh(kts:ktop(i)+1)-heo(i,kts:ktop(i)+1))/dtime
         !
         !     enddo
         !
         !     do k=kts,ktop(i)
         !       dp=100.*(po_cup(k,i)-po_cup(k+1,i))
         !       beta1 = g/dp
         !
         !       ddh(k) =  -( zuo(k+1,i)*hco (k+1,i) - zuo(k,i)*hco (k,i) )*beta1       &
         !            +( zdo(k+1,i)*hcdo(k+1,i) - zdo(k,i)*hcdo(k,i) )*beta1*edto(i)
         !
         !       ddh(k) = ddh(k) + c_xlf*((1.-p_liq_ice(k,i))* &
         !              0.5*(qrco(k+1,i)+qrco(k,i)) - melting(k,i))*beta1
         !
         !       dellah(k,i) =  dellah(k,i) + ddh(k)
         !
         !     enddo
         !       enddo
         !     endif
         !-------------------------
         !--- water vapor + condensates : flux form + source/sink terms + time explicit
         do vtp_index = get_num_elements(vec_ok),1,-1
            i = get_data_value(vec_ok,vtp_index)
            !
            do k=kts,ktop(i)
               dp=100.*(po_cup(k,i)-po_cup(k+1,i))

                !-- take out cloud liquid/ice water for detrainment
               detup        = up_massdetro(k,i)
               dellaqc(k,i) = detup*0.5*(qrco(k+1,i)+qrco(k,i)) *c_grav/dp
               !---
               G_rain=  0.5*(pwo (k,i)+pwo (k+1,i))*c_grav/dp
               E_dn  = -0.5*(pwdo(k,i)+pwdo(k+1,i))*c_grav/dp*edto(i) ! pwdo < 0 and E_dn must > 0
               !
               !-- condensation source term = detrained + flux divergence of
               !-- cloud liquid/ice water (qrco) + converted to rain

               C_up = dellaqc(k,i)+(zuo(k+1,i)*qrco(k+1,i) - zuo(k,i)* qrco(k,i))*c_grav/dp + G_rain

               !-- water vapor budget
               !-- = flux divergence z*(Q_c - Q_env)_up_and_down  - condensation term + evaporation
               dellaq(k,i) =-(zuo(k+1,i)*qco (k+1,i) - zuo(k,i)*qco (k,i))*c_grav/dp     &
                            +(zdo(k+1,i)*qcdo(k+1,i) - zdo(k,i)*qcdo(k,i))*c_grav/dp*edto(i)    &
                            - C_up + E_dn

               !--- source of cold pools
               dellabuoy(k,i)=edto(i)*dd_massdetro(k,i)*0.5*(dbydo(k+1,i)+dbydo(k,i))*c_grav/dp

            enddo
            !        if(use_fct == 0 .or. adjustl(cumulus) == 'shallow') then
            if(use_fct == 0 ) then
               do k=kts,ktop(i)
                  dp=100.*(po_cup(k,i)-po_cup(k+1,i))
                  sub_tend(k,1) =-(zuo(k+1,i)*(-qo_cup(k+1,i)) - zuo(k,i)*(-qo_cup(k,i)))*c_grav/dp       &
                                 +(zdo(k+1,i)*(-qo_cup(k+1,i)) - zdo(k,i)*(-qo_cup(k,i)))*c_grav/dp*edto(i)
               enddo
            else
               !-- FCT scheme for the subsidence transport: d(M_env*S_env)/dz
               sub_tend (:,1) = 0. ! dummy array
               trcflx_in(:,1) = 0. ! dummy array
               massflx  (:,i) = 0.
               dtime_max      = dtime

               do k=kts,ktop(i)
                  dp=100.*(po_cup(k,i)-po_cup(k+1,i))
                  trcflx_in (k,1) =-(zuo(k,i)  -edto(i)*zdo(k,i))*qo_cup(k,i) !* xmb(i)
                  massflx   (k,i) =-(zuo(k,i)  -edto(i)*zdo(k,i))             !* xmb(i)
                  dtime_max=min(dtime_max,.5*dp)
               enddo
               call fct1d3 (ktop(i),kte,dtime_max,po_cup(:,i),qo(:,i),massflx(:,i),trcflx_in(:,1),sub_tend(:,1))
            endif

            !--- add the contribuition from the environ subsidence
            dellaq(kts:ktop(i),i) = dellaq(kts:ktop(i),i) + sub_tend(kts:ktop(i),1)

             !--- for output only
            subten_Q (kts:ktop(i),i) = sub_tend(kts:ktop(i),1)

            !     do k=kts,ktop(i)
            !     print*,"delq=",use_fct,k,dellaq(k,i) , sub_tend(1,k)
            !     enddo

            !- check H and water conservation liq+condensed (including rainfall)
            trash  = 0.
            trash2 = 0.0
            do k=kts,ktop(i)
               dp     = 100.*(po_cup(k,i)-po_cup(k+1,i))
               G_rain =  0.5*(pwo (k,i)+pwo (k+1,i))*c_grav/dp
               E_dn   = -0.5*(pwdo(k,i)+pwdo(k+1,i))*c_grav/dp*edto(i)
               trash  = trash + (dellaq(k,i) + dellaqc(k,i)+ G_rain-E_dn)*dp/c_grav
               trash2 = trash2+  dellah(k,i)*c_grav/dp + c_xlf*((1.-p_liq_ice(k,i))*0.5*(qrco(k+1,i)+qrco(k,i)) &
                  - melting(k,i))*c_grav/dp
            enddo   ! k
            !--- test only with double precision:
            !write(0,*)'=>H/W-FINAL= ',real(trash2,4),real(trash,4),k22(i),kbcon(i),ktop(i)
            !if(abs(trash)>1.e-6 .or. abs(trash2) > 1.e-6) then
            !    write(0,*)'=> not water mass or H cons for deep= ',i,trash,trash2
            !    !stop 33
            !endif

         enddo


      endif ! vertical discretization formulation

   end subroutine get_dellas

  !-----------------------------------------------------------------------------------------
   subroutine fct1d3 (ktop,n,dt,z,tracr,massflx,trflx_in,del_out)

      ! --- modify a 1-D array of tracer fluxes for the purpose of maintaining
      ! --- monotonicity (including positive-definiteness) in the tracer field
      ! --- during tracer transport.

      ! --- the underlying transport equation is   (d tracr/dt) = - (d trflx/dz)
      ! --- where  dz = |z(k+1)-z(k)| (k=1,...,n) and  trflx = massflx * tracr

      ! --- note: tracr is carried in grid cells while z and fluxes are carried on
      ! --- interfaces. interface variables at index k are at grid location k-1/2.
      ! --- sign convention: mass fluxes are considered positive in +k direction.

      ! --- massflx and trflx_in  must be provided independently to allow the
      ! --- algorithm to generate an auxiliary low-order (diffusive) tracer flux
      ! --- as a stepping stone toward the final product trflx_out.

      implicit none
      integer,intent(in) :: n,ktop                        ! number of grid cells
      real   ,intent(in) :: dt                            ! transport time step
      real   ,intent(in) :: z(n+0)                        ! location of cell interfaces
      real   ,intent(in) :: tracr(n)                      ! the transported variable
      real   ,intent(in) :: massflx  (n+0)                ! mass flux across interfaces
      real   ,intent(in) :: trflx_in (n+0)                ! original tracer flux
      real   ,intent(out):: del_out  (n+0)                ! modified tracr flux
      real               :: trflx_out(n+0)                ! modified tracr flux
      integer k,km1,kp1
      logical :: NaN, error=.false., vrbos=.false.
      real dtovdz(n),trmax(n),trmin(n),flx_lo(n+0),antifx(n+0),clipped(n+0),  &
         soln_hi(n),totlin(n),totlout(n),soln_lo(n),clipin(n),clipout(n),arg
      real,parameter :: epsil=1.e-22           ! prevent division by zero
      real,parameter :: damp=1.                ! damper of antidff flux (1=no damping)

      logical, parameter :: hi_order = .false.

      NaN(arg) = .not. (arg.ge.0. .or. arg.le.0.)        ! NaN detector
      soln_lo(:)=0.
      antifx (:)=0.
      clipout(:)=0.
      flx_lo (:)=0.

      do k=1,ktop
         dtovdz(k)=.01*dt/abs(z(k+1)-z(k))                ! time step / grid spacing
      !     if (z(k).eq.z(k+1)) error=.true.
      end do
      if (vrbos .or. error) print '(a/(8es10.3))','(fct1d) dtovdz =',dtovdz(1:ktop)

      do k=2,ktop
         if (massflx(k) > 0.) then
            flx_lo(k)=massflx(k)*tracr(k-1)              ! low-order flux, upstream
         else
            flx_lo(k)=massflx(k)*tracr(k)                ! low-order flux, upstream
         endif
         antifx(k)=trflx_in(k)-flx_lo(k)                ! antidiffusive flux
      end do
      flx_lo(  1)   =trflx_in(  1)
      flx_lo(ktop+1)=trflx_in(ktop+1)
      antifx(  1)   =0.
      antifx(ktop+1)=0.
      ! --- clip low-ord fluxes to make sure they don't violate positive-definiteness
      do k=1,ktop
         totlout(k)=max(0.,flx_lo(k+1))-min(0.,flx_lo(k  ))         ! total flux out
         clipout(k)=min(1.,tracr(k)/max(epsil,totlout(k))/ (1.0001*dtovdz(k)))
      end do

      do k=2,ktop
         if (massflx(k).ge.0.)  then
            flx_lo(k)=flx_lo(k)*clipout(k-1)
         else
            flx_lo(k)=flx_lo(k)*clipout(k)
         endif
      end do
      if (massflx(1)     .lt.0.) flx_lo(1)     =flx_lo(1)     *clipout(1)
      if (massflx(ktop+1).gt.0.) flx_lo(ktop+1)=flx_lo(ktop+1)*clipout(ktop)

      ! --- a positive-definite low-order (diffusive) solution can now be  constructed

      do k=1,ktop
         soln_lo  (k)=tracr(k)-(flx_lo(k+1)-flx_lo(k))*dtovdz(k)        ! low-ord solutn
         del_out  (k)=-c_grav*(flx_lo(k+1)-flx_lo(k))*dtovdz(k)/dt
      end do

      if(.not. hi_order) return

      soln_hi  (:)=0.
      clipin   (:)=0.
      trmin    (:)=0.
      trmax    (:)=0.
      clipped  (:)=0.
      trflx_out(:)=0.


      do k=1,ktop
         km1=max(1,k-1)
         kp1=min(n,k+1)
         trmax(k)=       max(soln_lo(km1),soln_lo(k),soln_lo(kp1),        &
            tracr  (km1),tracr  (k),tracr  (kp1))        ! upper bound
         trmin(k)=max(0.,min(soln_lo(km1),soln_lo(k),soln_lo(kp1),        &
            tracr  (km1),tracr  (k),tracr  (kp1)))       ! lower bound
      end do

      do k=1,ktop
         totlin (k)=max(0.,antifx(k  ))-min(0.,antifx(k+1))                ! total flux in
         totlout(k)=max(0.,antifx(k+1))-min(0.,antifx(k  ))                ! total flux out

         clipin (k)=min(damp,(trmax(k)-soln_lo(k))/max(epsil,totlin (k)) / (1.0001*dtovdz(k)))
         clipout(k)=min(damp,(soln_lo(k)-trmin(k))/max(epsil,totlout(k)) / (1.0001*dtovdz(k)))

         if (NaN(clipin (k))) print *,'(fct1d) error: clipin is NaN,  k=',k
         if (NaN(clipout(k))) print *,'(fct1d) error: clipout is NaN,  k=',k

         if (clipin(k).lt.0.) then
            print 100,'(fct1d) error: clipin < 0 at k =',k,                        &
               'clipin',clipin(k),'trmax',trmax(k),'soln_lo',soln_lo(k),        &
               'totlin',totlin(k),'dt/dz',dtovdz(k)
            error=.true.
         end if
         if (clipout(k).lt.0.) then
            print 100,'(fct1d) error: clipout < 0 at k =',k,                        &
               'clipout',clipout(k),'trmin',trmin(k),'soln_lo',soln_lo(k),        &
               'totlout',totlout(k),'dt/dz',dtovdz(k)
            error=.true.
         end if
100      format (a,i3/(4(a10,"=",es9.2)))
      end do

      do k=2,ktop
         if (antifx(k).gt.0.)  then
            clipped(k)=antifx(k)*min(clipout(k-1),clipin(k))
         else
            clipped(k)=antifx(k)*min(clipout(k),clipin(k-1))
         end if
         trflx_out(k)=flx_lo(k)+clipped(k)
         if (NaN(trflx_out(k)))  then
            print *,'(fct1d) error: trflx_out is NaN,  k=',k
            error=.true.
         end if
      end do

      trflx_out(     1)=trflx_in(     1)
      trflx_out(ktop+1)=trflx_in(ktop+1)
      do k=1,ktop
         soln_hi(k)=tracr(k)-(trflx_out(k+1)-trflx_out(k))*dtovdz(k)
         del_out(k) =     -c_grav*(trflx_out(k+1)-trflx_out(k))*dtovdz(k)/dt
        !write(32,*)'3',k,soln_lo(k),soln_hi(k)
      end do

      if (vrbos .or. error) then
         do k=2,ktop
            write(32,99)k,                   &
               'tracr(k)', tracr(k),            &
               'flx_in(k)', trflx_in(k),        &
               'flx_in(k+1)', trflx_in(k+1),    &
               'flx_lo(k)', flx_lo(k),          &
               'flx_lo(k+1)', flx_lo(k+1),      &
               'soln_lo(k)', soln_lo(k),        &
               'trmin(k)', trmin(k),            &
               'trmax(k)', trmax(k),            &
               'totlin(k)', totlin(k),          &
               'totlout(k)', totlout(k),        &
               'clipin(k-1)', clipin(k-1),      &
               'clipin(k)', clipin(k),          &
               'clipout(k-1)', clipout(k-1),    &
               'clipout(k)', clipout(k),        &
               'antifx(k)', antifx(k),          &
               'antifx(k+1)', antifx(k+1),      &
               'clipped(k)', clipped(k),        &
               'clipped(k+1)', clipped(k+1),    &
               'flx_out(k)', trflx_out(k),      &
               'flx_out(k+1)', trflx_out(k+1),  &
               'dt/dz(k)', dtovdz(k),           &
               'final', tracr(k)-(trflx_out(k+1)-trflx_out(k))*dtovdz(k)
99          format ('(trc1d)   k =',i4/(3(a13,'=',es13.6)))
         end do
         if (error) stop '(fct1d error)'
      end if
   end subroutine fct1d3
   !---------------------------------------------------------------------------------------------------
   subroutine tridiag (m,a,b,c,f)
      !-- this routine solves the problem: aa*f(k-1,t+1) + bb*f(k,t+1) + cc*f(k+1,t+1) = dd
      !-- an updated "f" at time t+1 is the output
      implicit none
      integer, intent(in) :: m
      real, dimension(m), intent(inout) :: a,b,c
      real, dimension(m), intent(inout) :: f
      !--locals
      real, dimension(m) :: q
      integer :: k
      real :: p

      c(m)=0.
      q(1)=-c(1)/b(1)
      f(1)= f(1)/b(1)
      do k=2,m
         p  = 1./( b(k)+a(k)*q(k-1) )
         q(k) = -c(k)*p
         f(k) = p*(f(k) - a(k)*f(k-1))
      enddo
      do k=m-1,1,-1
         f(k) = f(k) +q(k)*f(k+1)
      enddo
   end subroutine tridiag
   !---------------------------------------------------------------------------------------------------
   subroutine bidiag (m,b,c,f)
      !-- this routine solves the problem:  bb*f(k,t+1) + cc*f(k+1,t+1) = dd
      !-- an updated "f" at time t+1 is the output
      implicit none
      integer, intent(in) :: m
      real, dimension(m), intent(inout) :: b,c
      real, dimension(m), intent(inout) :: f
      !--locals
      real, dimension(m) :: q
      integer :: k
      real :: p

      c(m)=0.
      q(1)=-c(1)/b(1)
      f(1)= f(1)/b(1)
      do k=2,m
         p  = 1./b(k)
         q(k) = -c(k)*p
         f(k) =  f(k)*p
      enddo
      do k=m-1,1,-1
         f(k) = f(k) +q(k)*f(k+1)
      enddo
   end subroutine bidiag

   !------------------------------------------------------------------------------------
   subroutine BeckerEtAlClosure(cumulus,ave_layer,its, itf, ite, kts, ktf, kte, k22, start_level, ktop, klcl, kbcon &
                              , t_in, tn, tn_adv , q_in, qo, qo_adv, psur, tsur,us, vs, zqexec, ztexec, x_add_buoy  &
                              , xland, tpert, zu, up_massdetro, up_massentro, zuo, p_liq_ice, qrco, ierr, zo, po, z1 &
                              , aa1_radpbl, aa1_adv)
      !! ## Becker et Al closure
      !!
      !! Author: Saulo R Freitas [SRF]
      !!
      !! E-mail: <mailto:saulo.freitas@inpe.br>
      !!
      !! Date: 22Fevereiro2023 13:46
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Becker et Al closure
      !!
      !! ** History**:
      !!
      !! --- 
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!

      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_name = 'beckerEtAlClosure' 
      !! subroutine name
   
      ! Variables (input, output, inout)
      character *(*)   ,intent (in) :: cumulus
      real   , intent(in) :: ave_layer
      integer, intent(in) :: its, itf, ite, kts, ktf, kte
      integer, intent(in) :: k22(:)
      integer, intent(in) :: start_level(:)
      integer, intent(in) :: ktop(:)
      integer, intent(in) :: klcl(:)
      integer, intent(in) :: kbcon(:)
      real, intent(in) :: t_in(:,:)
      real, intent(in) :: tn(:,:)
      real, intent(in) :: tn_adv(:,:)
      real, intent(in) :: q_in(:,:)
      real, intent(in) :: qo(:,:)
      real, intent(in) :: qo_adv(:,:)
      real, intent(in) :: psur(:)      
      real, intent(in) :: tsur(:)
      real, intent(in) :: us(:, :)
      real, intent(in) :: vs(:, :)
      real, intent(in) :: zqexec(:)
      real, intent(in) :: ztexec(:)
      real, intent(in) :: x_add_buoy(:)
      real, intent(in) :: xland(:)
      real, intent(in) :: tpert(:,:)
      real, intent(in) :: zu(:,:)
      real, intent(in) :: up_massdetro(:,:)
      real, intent(in) :: up_massentro(:,:)
      real, intent(in) :: zuo(:,:)
      real, intent(in) :: p_liq_ice(:,:)
      real, intent(in) :: qrco(:,:)
      integer, intent(inout) :: ierr(:)
      
      real, intent(in ) :: zo(:,:)
      real, intent(in ) :: po(:, :)
      real, intent(in ) :: z1(:)
      real, intent(out) :: aa1_radpbl(:)
      real, intent(out) :: aa1_adv(:)
      
      ! Local variables:
      integer ::ierr_dummy(its:ite) 
      real :: hkbo_x      (its:ite)
      real :: qeso_x      (kts:kte,its:ite)
      real :: heo_x       (kts:kte,its:ite)
      real :: heso_x      (kts:kte,its:ite)
      real :: qeso_cup_x  (kts:kte,its:ite)
      real :: qo_cup_x    (kts:kte,its:ite)
      real :: heo_cup_x   (kts:kte,its:ite)
      real :: u_cup_x     (kts:kte,its:ite)
      real :: v_cup_x     (kts:kte,its:ite)
      real :: heso_cup_x  (kts:kte,its:ite)
      real :: zo_cup_x    (kts:kte,its:ite)
      real :: po_cup_x    (kts:kte,its:ite)
      real :: gammao_cup_x(kts:kte,its:ite)
      real :: tn_cup_x    (kts:kte,its:ite)
      real :: hco_x       (kts:kte,its:ite)
      real :: dbyo_x      (kts:kte,its:ite)
      real :: tn_x        (kts:kte,its:ite)
      real :: qo_x        (kts:kte,its:ite)
      integer :: ki, i, k
      real :: x_add, denom
      
      do ki = 1,2
         if(DICYCLE==2 .and. ki==2) cycle 

         if(ki==1) then
          !-- get the cloud work function for updrafts associated only with RAD + PBL
          tn_x = t_in + tn - tn_adv
          qo_x = q_in + qo - qo_adv
          !-- to check => aa1_radpbl=aa1
          !! tn_x = tn
          !! qo_x = qo
         endif

         if(ki==2) then
            !-- get the cloud work function for updrafts associated only with Qv-advection
           !tn_x = tn_adv  ! orig 
           tn_x = t_in        ! v2
           qo_x = qo_adv
         endif

         ierr_dummy=ierr

         call cup_env     (zo,qeso_x,heo_x,heso_x,tn_x,qo_x,po,z1,psur,ierr_dummy,-1,itf,ktf, its,ite, kts,kte)
         call cup_env_clev(tn_x,qeso_x,qo_x,heo_x,heso_x,zo,po,qeso_cup_x,qo_cup_x,heo_cup_x,us,vs   &
                          ,u_cup_x,v_cup_x,heso_cup_x,zo_cup_x,po_cup_x,gammao_cup_x,tn_cup_x,psur,tsur  &
                          ,ierr_dummy,z1,itf,ktf,its,ite, kts,kte)

         !--- get MSE
         do i=its,itf
            if(ierr_dummy(i) /= 0) cycle
            x_add = (c_alvl*zqexec(i)+c_cp*ztexec(i)) +  x_add_buoy(i)
            call get_cloud_bc(cumulus,ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i),heo_cup_x(kts:kte,i)&
                             ,hkbo_x(i),k22(i),x_add,Tpert(kts:kte,i))
            hco_x (kts:start_level(i),i) = hkbo_x(i)

            do k=start_level(i)+1 , ktop(i)+1  ! mass cons option
               denom =(zu(k-1,i)-.5*up_massdetro (k-1,i)+up_massentro (k-1,i)) + 1.e-12
               
               hco_x(k,i)=(hco_x(k-1,i)*zuo(k-1,i)-.5*up_massdetro(k-1,i)*hco_x(k-1,i)+ &
                                                      up_massentro(k-1,i)*heo_x(k-1,i))/ denom
               if(k==start_level(i)+1) then
                     !x_add = (xlv*zqexec(i)+cp*ztexec(i)) +  x_add_buoy(i)
                     hco_x(k,i)= hco_x(k,i) + x_add*up_massentro(k-1,i)/denom
               endif
               !
               !- includes glaciation effects on HCO_X
               hco_x(k,i)= hco_x(k,i)+(1.-p_liq_ice(k,i))*qrco(k,i)*c_xlf
            enddo
            hco_x (ktop(i)+2:ktf,i)=heso_cup_x(ktop(i)+2:ktf,i)
         enddo

         call get_buoyancy(itf,ktf, its,ite, kts,kte,ierr_dummy,klcl,kbcon,ktop &
                          ,hco_x,heo_cup_x,heso_cup_x,dbyo_x,zo_cup_x)

         if(ki==1) &  ! RAD+PBL only
            call cup_up_aa0(aa1_radpbl,zo_cup_x,zuo,dbyo_x,GAMMAo_CUP_x,tn_cup_x   &
                           ,k22,klcl,kbcon,ktop,ierr_dummy,itf,ktf,its,ite, kts,kte)
         !-- get AA1_ADV
         !aa1_adv = aa1 + aa0 - aa1_radpbl
         

         if(ki==2) & ! ADV of Qv only
         call cup_up_aa0(aa1_adv,zo_cup_x,zuo,dbyo_x,GAMMAo_CUP_x,tn_cup_x   &
                        ,k22,klcl,kbcon,ktop,ierr_dummy,itf,ktf,its,ite, kts,kte)
         !Observe that : 
         !aa1 ~ aa0 + (aa1_radpbl-aa0) + (aa1_adv-aa0)

      enddo ! ki   
   end subroutine BeckerEtAlClosure
   !------------------------------------------------------------------------------------
   subroutine get_Qadv(cumulus,itf,ktf,its,ite,kts,kte,ierr,dt,q,qo,qo_adv,po,po_cup &
                      ,qeso, Q_adv,col_sat_adv,alpha_adv,tau_bl,zo_cup,kbcon,ktop)

      implicit none
      real    :: alpha_adv_tuning      = 0.8  != tuning parameter for the Becker et al (2021) closure
      real    :: col_sat_adv_threshold = 0.94 != suppress Qadv closure for col_sat_adv > col_sat_adv_threshold

      character *(*), intent (in)                      :: cumulus
      integer ,intent (in)                             :: itf,ktf, its,ite, kts,kte
      real    ,intent (in)                             :: dt
      integer ,intent (in) ,dimension(:)   :: ierr,kbcon,ktop

      real    ,intent (in) ,dimension(:,:) :: q,qo,qo_adv,po_cup,qeso,po,zo_cup
      real    ,intent (in) ,dimension(:)   :: tau_bl

      real    ,intent (inout) ,dimension(:)   :: Q_adv,col_sat_adv,alpha_adv

      !--locals
      integer :: i,k
      real :: dp, layer,H_cloud,dz
      real ,parameter :: ptop = 60.
      !-- get the advective moisture tendency scaled with the relative humidity
      !--  Q_adv = integral( q/q*  DQv/Dt_adv dp), see Eq 1 Becker et al(2021 QJRMS)
      !-- units here are "J m^-3" _or_  "J kg^-1"

      do i=its,itf
         col_sat_adv(i) = 0.   !check if it needs be inout, perhavps only local var

         if(ierr(i) /= 0) cycle

         alpha_adv  (i) = alpha_adv_tuning
         layer = 0.

         loopN: do k=kts,ktf
            if(po(k,i) < ptop) exit loopN

            !dp=100.*(po_cup(k+1,i)-po_cup(k,i)) ! dp < 0.
            dz=      zo_cup(k+1,i)-zo_cup(k,i)  ! dz > 0

            !-- integral over dp
            !Q_adv(i) = Q_adv(i) + dp*(qo(k,i)/qeso(k,i))*(qo_adv(k,i)-q(k,i))/dt

            !-- integral over dz
            Q_adv(i) = Q_adv(i) + dz*(qo(k,i)/qeso(k,i))*(qo_adv(k,i)-q(k,i))/dt

            col_sat_adv(i) = col_sat_adv(i) + dz* qo(k,i)/qeso(k,i)

            layer = layer + dz

         enddo loopN
         !-- get the column-average saturation fraction
         col_sat_adv(i) = col_sat_adv(i)/(1.e-8+layer)

         !--check if the col-ave saturation fraction is over the threshold
         if(col_sat_adv(i) > col_sat_adv_threshold) then

            alpha_adv(i) = 0.0
            cycle

         endif

         !-- check if cloud top _OR_cloud layer   !<<<< check this
         H_cloud =  zo_cup(ktop(i),i)- zo_cup(kbcon(i),i)

         !-- convert Q_adv to units as in Eq (1) => J m^-3
         !Q_adv(i) = - Q_adv(i) * tau_bl(i) * c_alvl / (g * H_cloud)

         !-- convert Q_adv to units as in cloud work function => J kg-1
         Q_adv(i) =  Q_adv(i) * tau_bl(i) * c_alvl / (H_cloud)

          !if(abs(q_adv(i))>1.) print*,"Qadv=",i,q_adv(i),Q_adv_dz(i)call flush(6)
      enddo

   end subroutine get_Qadv

  !---------------------------------------------------------------------------------------------------
   subroutine get_pr_ens(cumulus,c0_mid,itf,ktf,its,ite,kts,kte,maxens3,ktop,ierr,ierrc &
                        ,pwo,pwdo,edto,pr_ens)
      implicit none
      character*(*), intent (in)                   :: cumulus
      character*(*), intent (inout), dimension(:)  :: ierrc

      integer  ,intent (in   )               :: itf,ktf, its,ite, kts,kte,maxens3
      integer  ,intent (in   ), dimension(:) :: ktop
      integer  ,intent (inout), dimension(:) :: ierr

      real     ,intent (in )                 :: c0_mid
      real     ,intent (in ), dimension(:)   :: edto
      real     ,intent (in ), dimension(:,:) :: pwo,pwdo
      real     ,intent (out), dimension(:,:) :: pr_ens
      
      !-- local vars
      integer :: i,k,nens3,vtp_index

      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         do k=kts,ktop(i)
               do nens3=1,maxens3
                  if(nens3.eq.7)then
                  !--- b=0
                      pr_ens(i,nens3)=pr_ens(i,nens3) +pwo(k,i)+edto(i)*pwdo(k,i)
                  !--- b=beta
                  else if(nens3.eq.8)then
                      pr_ens(i,nens3)=pr_ens(i,nens3) +pwo(k,i)+edto(i)*pwdo(k,i)
                  !--- b=beta/2
                  else if(nens3.eq.9)then
                      pr_ens(i,nens3)=pr_ens(i,nens3) +pwo(k,i)+edto(i)*pwdo(k,i)
                  else
                      pr_ens(i,nens3)=pr_ens(i,nens3) +pwo(k,i)+edto(i)*pwdo(k,i)
                  endif
               enddo
         enddo
         if(pr_ens(i,7).lt.1.e-6 .and. c0_mid > 0. .and.  cumulus /= 'shallow' )then
                  ierr(i)=18 ; is_removed = remove(vec_ok, i)
                  ierrc(i)="total normalized condensate too small"
                  do nens3=1,maxens3
                     pr_ens(i,nens3)=0.
                  enddo
         endif
         do nens3=1,maxens3
                  if(pr_ens(i,nens3) < 1.e-5) pr_ens(i,nens3)=0.
         enddo
      enddo

   end subroutine get_pr_ens

  !---------------------------------------------------------------------------------------------------
   subroutine apply_sub_microphys(cumulus,alp1,itf,ktf,its,ite,kts,kte,ktop,nmp,ierr,po_cup,zenv &
                                 ,mpql,mpqi,mpcf,dellampqi,dellampql,dellampcf)
                                   
      implicit none
      character*(*), intent (in)               :: cumulus

      integer  ,intent (in )                   :: itf,ktf, its,ite, kts,kte,nmp
      integer  ,intent (in ), dimension(:)     :: ktop
      integer  ,intent (in ), dimension(:)     :: ierr

      real     ,intent (in )                   :: alp1
      real     ,intent (in ), dimension(:,:)   :: po_cup,zenv
      real     ,intent (in ), dimension(:,:,:) :: mpql,mpqi,mpcf
      real     ,intent (out), dimension(:,:,:) :: dellampqi,dellampql,dellampcf
      
      !-- local vars
      integer :: i,k,kmp,vtp_index
      real    :: dp,env_mf,env_mf_m,env_mf_p,beta1,beta2,alp0
      real, dimension (kts:kte)     ::  aa,bb,cc,ddu,ddv,ddh,ddq,fp,fm
      real, dimension (nmp,kts:kte) ::  dd


      dellampqi =0.
      dellampql =0.
      dellampcf =0.

      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         do k=kts,ktop(i)
            dp=100.*(po_cup(k,i  )-po_cup(k+1,i))

            !--- apply environmental subsidence on grid-scale/anvil ice and liq water contents (Upwind scheme)
            !
            env_mf   = - 0.5* (zenv(k+1,i) + zenv(k,i))
            env_mf_m = min(env_mf,0.)*c_grav/dp
            env_mf_p = max(env_mf,0.)*c_grav/dp

            dellampqi(:,k,i) = - (  env_mf_m*(mpqi(:,i,k+1)-mpqi(:,k,i))  +        &
                                    env_mf_p*(mpqi(:,i,k  )-mpqi(:,i,max(k-1,kts))))
            dellampql(:,k,i) = - (  env_mf_m*(mpql(:,i,k+1)-mpql(:,k,i))  +         &
                                    env_mf_p*(mpql(:,i,k  )-mpql(:,i,max(k-1,kts))))

            !--- apply environmental subsidence on grid-scale/anvil cloud fraction
            dellampcf(:,k,i) = - (  env_mf_m*(mpcf(:,i,k+1)-mpcf(:,k,i))  +         &
                                    env_mf_p*(mpcf(:,i,k  )-mpcf(:,i,max(k-1,kts))))
         enddo

         !--- apply environmental subsidence on grid-scale and anvil cloud fraction using time implicit/explict method
         if(alp1 > 0.) then
            alp0=1.0-alp1
            do k=kts,ktop(i)
               dp=100.*(po_cup(k,i  )-po_cup(k+1,i))
               env_mf   = - 0.5* (zenv(k+1,i) + zenv(k,i))
               env_mf_m = min(env_mf,0.)*c_grav/dp
               env_mf_p = max(env_mf,0.)*c_grav/dp

               beta1 = -env_mf_m
               beta2 = -env_mf_p

               aa(k) =    alp1*beta2             ! coef of f(k-1,t+1),
               bb(k) = 1.+alp1*beta1-alp1*beta2  ! coef of f(k  ,t+1),
               cc(k) =   -alp1*beta1             ! coef of f(k+1,t+1),

               !-- this is the rhs of the discretization
               dd(:,k) = (1.-alp0*beta1+alp0*beta2)*mpcf(:,i,k  ) + &       ! coef of  f(k  ,t),
                             alp0*beta1 *mpcf(:,i,k+1)            - &       ! coef of  f(k+1,t),
                             alp0*beta2 *mpcf(:,i,max(kts,k-1))             ! coef of  f(k-1,t),
            enddo
            do kmp =1,nmp
                !-- this routine solves the problem: aa*f(k-1,t+1) + bb*f(k,t+1) + cc*f(k+1,t+1) = dd
               call tridiag (ktop(i),aa(kts:ktop(i)), bb(kts:ktop(i)), cc(kts:ktop(i)), dd(kmp,kts:ktop(i)))

               dellampcf(kmp,i,kts:ktop(i)) = dd(kmp,kts:ktop(i))-mpcf(kmp,i,kts:ktop(i))
            enddo
         endif
      enddo
   end subroutine apply_sub_microphys
  !---------------------------------------------------------------------------------------------------
   subroutine get_max_dd_height(cumulus,zcutdown,itf,ktf,its,ite,kts,kte,ktop,kzdown,ierr,z1,zo_cup)
                                   
      implicit none
      character*(*), intent (in)             :: cumulus

      integer  ,intent (in )                 :: itf,ktf, its,ite, kts,kte
      integer  ,intent (in ), dimension(:)   :: ktop
      integer  ,intent (in ), dimension(:)   :: ierr
      integer  ,intent (out), dimension(:)   :: kzdown

      real     ,intent (out)                 :: zcutdown
      real     ,intent (in ), dimension(:)   :: z1
      real     ,intent (in ), dimension(:,:) :: zo_cup
    
      !-- local vars
      integer :: i,k,vtp_index
      real :: zktop 

      !--- height(m) above which no downdrafts are allowed to originate
      !
      zcutdown  = 3000.
      kzdown(:) = kts
      if(trim(cumulus) == 'shallow') return
      !
      do vtp_index = get_num_elements(vec_ok),1,-1
         i = get_data_value(vec_ok,vtp_index)
         zktop=(zo_cup(ktop(i),i)-z1(i))*.6

         zktop=min(zktop+z1(i),zcutdown+z1(i))
         
         loopK: do k=kts,ktf
                   if(zo_cup(k,i) >= zktop)then
                     kzdown(i)=k
                     exit loopK 
                  endif
                enddo loopK
      enddo
   end subroutine get_max_dd_height
!---------------------------------------------------------------------------------------------------
   function initModConvParGF() result(is_init)
      !! # Initialize the module with values
      !!
      !! Author: Rodrigues, L.F. [LFR]
      !!
      !! E-mail: <mailto:luiz.rodrigues@inpe.br>
      !!
      !! Date: 26Janeiro2023 16:41
      !!
      !! #####Version: 0.1.0
      !!
      !! —
      !! **Full description**:
      !!
      !! Initialize all variables from module
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !!
      implicit none
      !Parameters:
      character(len=*), parameter :: procedureName = 'initModConvParGF' ! Function Name
   
      !Local variables:
      integer :: is_init
   
      !Code:
      ! Se o módulo já foi inicializado retorna -1
      if(modConvParGF_initialized) then
         is_init = -1
         return
      endif

      !-- inicializa as variáveis do módulo
      icumulus_gf      = (/1, 1, 1/)
      closure_choice   = (/10,10, 3/) 
      cum_entr_rate    = (/        &
                           6.3e-4 & !deep
                          ,1.0e-3 & !shallow
                          ,5.0e-4 & !mid
                         /)
      use_lcl_ctrl_entr = 0
      use_rhu_ctrl_entr = 1
      use_pass_cloudvol = 0
      use_tracer_transp = 1 
      use_tracer_scaven = 2 
      use_tracer_evap   = 1 
      use_flux_form     = 1 
      use_fct           = 1 
      clev_grid         = 1 
      alp1              = 1.0 
      vert_discr        = 1 

      convection_tracer = 1 
      use_memory        = 2 
      add_coldpool_prop = 3 
      add_coldpool_clos = 2
      add_coldpool_trig = 0
      add_coldpool_diff = 3 
      tau_ocea_cp       = 7200. 
      tau_land_cp       = 7200. 
      mx_buoy1          = (real(c_cp)*5.0 + real(c_alvl)*2.e-3)*0.025  
      mx_buoy2          = (real(c_cp)*10.+real(c_alvl)*4.e-3) 
      use_gustiness     = 0    
     
      use_scale_dep     = 1 
      sig_factor        = 0.22 
      
      dicycle           = 1 
      cum_t_star        = (/4., -99., -99./)
      rh_dicycle        = 0 
      
      downdraft         = 1   
      use_rebcb         = 1 
      irainevap         = 0
      cum_max_edt_land  = (/0.90, 0.00, 0.20/)
      cum_max_edt_ocean = (/0.90, 0.00, 0.20/)
      cum_use_excess    = (/2, 1, 1/)
      cum_ave_layer     = (/50., 25., 25./)
      sgs_w_timescale   = 1 
      tau_deep          = 3600.0  
      tau_mid           = 1200.0  
      lightning_diag    = 0 
      liq_ice_number_conc=0
      apply_sub_mp      = 0 
      use_wetbulb       = 0 
      overshoot         = 0.0

      autoconv          = 1     
      c0_deep           = 1.0e-3
      c0_mid            = 1.5e-3
      c0_shal           = 0.0    
      qrc_crit          = 6.e-4 
      n_cldrop          = 50.0  

      use_momentum_transp = 1   
      lambau_deep       = 0.0 
      lambau_shdn       = 2.0 
      
      max_tq_tend       = 500.0   
      use_smooth_prof   = 1      
      use_smooth_tend   = 1      
      cum_hei_down_land = (/0.40, 0.00, 0.35/)
      cum_hei_down_ocean= (/0.35, 0.00, 0.35/)
      cum_hei_updf_land = (/0.55, 0.10, 0.55/)
      cum_hei_updf_ocean= (/0.55, 0.10, 0.55/)
      cum_fadj_massflx  = (/1.00, 1.00, 1.00/)
      moist_trigger     = 0   
      frac_modis        = 1   
      adv_trigger       = 0   
      dcape_threshold   = 60.0
      lcl_trigger       = 0   
      use_random_num    = 0.0   
      cap_maxs          = 50.  
      beta_sh           = 2.2  
      use_linear_subcl_mf = 1    
      use_cloud_dissipation = 0.0 
 
      time_in        = 0.0
      int_time       = 0.0
      whoami_all     = 0
      jcol           = 0
      itime1_in      = 0
      ntimes         = 0
      output_sound   = 0   
      first_guess_w  = .false. 
      wrtgrads       = .false.

      chem_adj_autoc = 0.0
      ispc_co        = 0
      ind_chem       = 0
      chem_name_mask = 0
      chem_name_mask_evap = 0
      chem_name = ""
      if(allocated(hcts)) then
         hcts%hstar = 0.0
         hcts%dhr   = 0.0
         hcts%ak0   = 0.0
         hcts%dak   = 0.0
      endif
      ! informa que já inicializado
      modconvpargf_initialized = .true.
      ! retorna 0, foi inicializado dessa vez
      is_init = 0

   end function initModConvParGF
   
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
!--------------------------------------------------------------------------------------------!
!- section for atmospheric composition
!--------------------------------------------------------------------------------------------!
      subroutine tracerConvTransGF(its, itf, jl, jmin, ite, kts, kte, ktf, mtp, chem_name_mask, ave_layer, ierr, k22, ktop &
                           ,   start_level, dd_massdetro, dd_massentro, dtime, edto, fscav, po, po_cup, pw_up_chem &
                           ,   pwavo, pwevo, pwdo, pwo, qrco, sc_up_chem, tempco, tot_pw_up_chem &
                           ,   up_massdetro, up_massentro, vvel2d, xland, zdo, zo_cup, zuo, cumulus &
                           ,   se_chem, massf, out_chem, pw_dn_chem, sc_dn_chem, se_cup_chem, tot_pw_dn_chem, zenv) 
      !! ## section for atmospheric composition
      !!
      !! Author: Rodrigues, L. F. [LFR]
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 23Fevereiro2023 16:44
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! section for atmospheric composition
      !!
      !! ** History**:
      !!
      !! ---
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!

      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_name = 'atmosComposition'
      !! subroutine name

      ! Variables (input, output, inout)
      character(len=*), intent(in) :: cumulus
      integer, intent(in) :: its, itf, jl, ite, kts, kte, ktf, mtp
      integer, intent(in) :: chem_name_mask(:)
      integer, intent(in) :: ierr(:)
      integer, intent(in) :: jmin(:)
      integer, intent(in) :: k22(:)
      integer, intent(in) :: ktop(:)
      integer, intent(in) :: start_level(:)

      real, intent(in) :: ave_layer
      real, intent(in) :: dtime
      real, intent(in) :: edto(:)
      real, intent(in) :: fscav(:)
      real, intent(in) :: po(:,:)
      real, intent(in) :: po_cup(:,:)
      real, intent(in) :: pwavo(:)
      real, intent(in) :: pwevo(:)
      real, intent(in) :: pwdo(:,:)
      real, intent(in) :: pwo(:,:)
      real, intent(in) :: qrco(:,:)
      real, intent(in) :: tempco(:,:)
      real, intent(in) :: vvel2d(:,:)
      real, intent(in) :: xland(:)
      real, intent(in) :: zdo(:,:)
      real, intent(in) :: zo_cup(:,:)
      real, intent(in) :: zuo(:,:)
      real, intent(out) :: zenv(:,:)
      real, intent(in) :: dd_massdetro(:,:)
      real, intent(in) :: dd_massentro(:,:)
      real, intent(in) :: up_massdetro(:,:)
      real, intent(in) :: up_massentro(:,:)



      real, intent(inout) :: se_chem(:,:,:)

      real, intent(out) :: massf
      real, intent(out) :: out_chem(:,:,:)
      real, intent(out) :: pw_dn_chem(:,:,:)
      real, intent(out) :: pw_up_chem(:,:,:)
      real, intent(out) :: sc_dn_chem(:,:,:)
      real, intent(out) :: sc_up_chem(:,:,:)
      real, intent(out) :: se_cup_chem(:,:,:)
      real, intent(out) :: tot_pw_dn_chem(:,:)
      real, intent(out) :: tot_pw_up_chem(:,:)

      ! Local variables:
      real :: se_chem_update(3, its:ite, kts:kte)
      real :: massi, dp, beta1, evap,wetdep
      integer :: i,  k, vtp_index, ispc
      real, dimension(mtp) :: evap_, wetdep_, trash_ ,trash2_, residu_


      !--only for debug
      if (p_use_gate) then
         if (jl == 1) then
            se_chem_update(1,:,:) = se_chem(1,:,:)
         else
            se_chem(1,:,:) = se_chem_update(1,:,:)
         end if
         ! DE: manual if cycle remove
         do vtp_index = 1, get_num_elements(vec_ok) ; i=get_data_value(vec_ok, vtp_index) 
            massi = 0.
            do k = kts, ktop(i)
               dp = 100.*(po_cup(k,i) - po_cup(k+1,i))
               massi = massi + se_chem(1,k,i)*dp/c_grav
            end do
         end do
      end if
      !--only for debug

      !-1) get mass mixing ratios at the cloud levels

      call cupEnvClevChem(mtp, se_chem, se_cup_chem, ierr, itf, ktf, its, kts, kte)

      !-2) determine in-cloud tracer mixing ratios
      !
      ! a) chem - updraft
      !- note: here "sc_up_chem" stores the total in-cloud tracer mixing ratio (i.e., including the portion
      !        embedded in the condensates).
      call getInCloudScChemUP(cumulus, ave_layer, fscav, mtp, se_chem, se_cup_chem, sc_up_chem, pw_up_chem, tot_pw_up_chem, zo_cup &
                             , po, po_cup, qrco, tempco, pwo, zuo, up_massentro, up_massdetro, vvel2d, start_level, k22 &
                             , ktop, ierr, xland, itf, ktf, its, ite, kts, kte)

      ! b) chem - downdraft
      call getInCloudScChemDD(cumulus, mtp, se_chem, se_cup_chem, sc_dn_chem, pw_dn_chem &
                              , tot_pw_up_chem, tot_pw_dn_chem, po_cup, pwdo, pwevo, zdo, dd_massentro &
                              , dd_massdetro, pwavo, jmin, ierr, itf, its, kts)


      !-3) determine the vertical transport including mixing, scavenging and evaporation
      !
      !---a) change per unit mass that a model cloud would modify the environment
      ! DE: manual if cycle remove
      do vtp_index = 1, get_num_elements(vec_ok) ; i=get_data_value(vec_ok, vtp_index) 

         call fluxFormsSourceSink(kte, kts, mtp, ktop(i), chem_name_mask, dtime, edto(i), po_cup(:,i) &
                                 , pw_dn_chem(:,:,i), pw_up_chem(:,:,i), sc_dn_chem(:,:,i) &
                                 , sc_up_chem(:,:,i) ,se_chem(:,:,i) ,se_cup_chem(:,:,i) ,zo_cup(:,i) &
                                 , zdo(:,i), zuo(:,i) , cumulus, zenv(:,i) ,out_chem(:,:,i))
      end do ! loop 'i'


      if (p_use_gate) then
         !--only for debug
         do vtp_index = 1, get_num_elements(vec_ok) ; i=get_data_value(vec_ok, vtp_index) 
            massf = 0.
            do k = kts, ktop(i)
               se_chem_update(ispc_co,k,i) = se_chem_update(ispc_co,k,i) + out_chem(ispc_co,k,i)*dtime
               !se_chem_update(ispc_CO,i,k) = max(0.,se_chem_update(ispc_CO,i,k))
               dp = 100.*(po_cup(k,i) - po_cup(k+1,i))
               evap_(ispc_co) = -0.5*(zdo(k,i)*pw_dn_chem(ispc_co,k,i) + zdo(k+1,i) &
                              * pw_dn_chem(ispc_co,k+1,i)) *c_grav/dp*edto(i)
               wetdep_(ispc_co) = 0.5*(zuo(k,i)*pw_up_chem(ispc_co,k,i) + zuo(k+1,i)&
                                *pw_up_chem(ispc_co,k+1,i)) *c_grav/dp
               massf = massf + se_chem_update(1,k,i)*dp/c_grav + (-evap_(ispc_co) + wetdep_(ispc_co))*dp/c_grav
            end do
            if (abs((massf - massi)/(1.e-12 + massi)) > 1.e-6) print *, "mass con=>", (massf - massi)/(1.e-12 + massi)
         end do
!   19    format(1x, I3, 1x, 5e14.3)
!   18    format(1x, I3, 1x, 4e14.3)
!   20    format(1x, I3, 1x, 11e16.6)
         !--only for debug
      end if


      !--- check mass conservation for tracers
      do ispc = 1, mtp
               if(CHEM_NAME_MASK (ispc) == 0 ) cycle
               trash_ (:) = 0.
               trash2_(:) = 0.
               evap_  (:) = 0.
               wetdep_(:) = 0.
               residu_(:) = 0.
               do k=kts,ktop(i)
                  dp=100.*(po_cup(k,i)-po_cup(k+1,i))
                  evap   =  -0.5*(zdo(k,i)*pw_dn_chem(ispc,k,i)+zdo(k+1,i)*pw_dn_chem(ispc,k+1,i))*c_grav/dp*edto(i)
                  wetdep =   0.5*(zuo(k,i)*pw_up_chem(ispc,k,i)+zuo(k+1,i)*pw_up_chem(ispc,k+1,i))*c_grav/dp

                  evap_  (ispc) =   evap_  (ispc) + evap  *dp/c_grav
                  wetdep_(ispc) =   wetdep_(ispc) + wetdep*dp/c_grav
                  residu_(ispc) =   residu_(ispc) + (wetdep - evap)*dp/c_grav

                  !trash_ (ispc) =   trash_ (ispc) + (out_chem (ispc,k,i) - evap + wetdep)*dp/c_grav
                  trash_ (ispc) =   trash_ (ispc) + (out_chem (ispc,k,i)                )*dp/c_grav

                  trash2_(ispc) =   trash2_(ispc) + se_chem(ispc,k,i)*dp/c_grav
               end do
               if(residu_(ispc) < 0.) then
                  beta1 = c_grav/(po_cup(kts,i)-po_cup(i,ktop(i)+1))
                  do k=kts,ktop(i)
                     out_chem(ispc,k,i)=out_chem(ispc,k,i)+residu_(ispc)*beta1
                  end do
               endif

               !if(evap_  (ispc) > wetdep_(ispc)) then
                !print*,"budget=",ispc,evap_  (ispc), wetdep_(ispc),trash_ (ispc),trim(CHEM_NAME(ispc))!,trash_ (ispc),trash2_(ispc)
                !call flush(6)
                !endif
                !if(evap_  (ispc) > wetdep_(ispc)) stop " eva<wet "
                !if(abs(trash_(ispc)) >1.e-6 ) then
                !  if (MAPL_AM_I_ROOT())  write(6,*)'=> mass_cons=',trash_(ispc),spacing(trash2_(ispc)),trim(CHEM_NAME(ispc)),trim(cumulus)
                !endif
               !enddo

      end do

   end subroutine tracerConvTransGF     
   ! ---------------------------------------------------------------------------------------------------
   subroutine cupEnvClevChem(mtp, se_chem, se_cup_chem, ierr, itf, ktf, its, kts, kte)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: procedureName = 'cupEnvClevChem' ! Subroutine Name

      integer, parameter ::  p_clev_option = 2 
      !! use option 2
   
      !Variables (input, output, inout)
      integer, intent(in) :: itf, ktf, its, kts, kte, mtp
      
      integer, intent(in) :: ierr(:)

      real, intent(in) :: se_chem(:,:,:)

      real, intent(out) :: se_cup_chem(:,:,:)

      !Local variables:
      integer :: i, k
      integer :: vtp_index
      
      if (p_clev_option == 1) then
         !-- original version
         do vtp_index = 1, get_num_elements(vec_ok) ; i=get_data_value(vec_ok, vtp_index) !BD_n
            do k = kts + 1, ktf
               se_cup_chem(1:mtp,k,i) = 0.5*(se_chem(1:mtp,k-1,i) + se_chem(1:mtp,k,i))
            end do
            se_cup_chem(1:mtp,kts,i) = se_chem(1:mtp,kts,i)
            se_cup_chem(1:mtp,kte,i) = se_chem(1:mtp,ktf,i)
         end do
      else
         !-- version 2: se_cup (k+1/2) = se(k) => smoother profiles
         do vtp_index = 1, get_num_elements(vec_ok) ; i=get_data_value(vec_ok, vtp_index) !BD_n
            do k = kts, ktf
               se_cup_chem(1:mtp,k,i) = se_chem(1:mtp,k,i)
            end do
         end do
      end if

   end subroutine cupEnvClevChem
   ! ------------------------------------------------------------------------------------
   subroutine getInCloudScChemUp(cumulus, ave_layer, fscav, mtp, se, se_cup, sc_up, pw_up, tot_pw_up_chem &
                               , z_cup, po, po_cup, qrco, tempco, pwo, zuo, up_massentro, up_massdetro &
                               , vvel2d, start_level, k22, ktop, ierr, xland, itf, ktf &
                               , its, ite, kts, kte)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: procedureName = 'getInCloudScChemUp' ! Subroutine Name
   
      real, parameter :: p_scav_eff = 0.6  
      !! for smoke : Chuang et al. (1992) J. Atmos. Sci.
      real, parameter :: p_cte_w_upd = 10. 
      !! m/s
      !    real, parameter :: kc = 5.e-3  
      !! s-1
      real, parameter :: p_kc = 2.e-3
      !! autoconversion parameter in GF is lower than what is used in GOCART s-1

      !Variables (input, output, inout)
       integer, intent(in)  :: itf, ktf, its, ite, kts, kte, mtp

      integer, intent(in) :: ierr(:)
      integer, intent(in) :: ktop(:)
      integer, intent(in) :: k22(:)
      integer, intent(in) :: start_level(:)

      real, intent(in) :: fscav(:)
      real, intent(in) :: se(:,:,:)
      real, intent(in) :: se_cup(:,:,:)
      real, intent(in) :: z_cup(:,:)
      real, intent(in) :: po_cup(:,:)
      real, intent(in) :: qrco(:,:)
      real, intent(in) :: tempco(:,:)
      real, intent(in) :: pwo(:,:)
      real, intent(in) :: zuo(:,:)
      real, intent(in) :: up_massentro(:,:)
      real, intent(in) :: up_massdetro(:,:)
      real, intent(in) :: po(:,:)
      real, intent(in) :: vvel2d(:,:)
      real, intent(in) :: xland(:)
      real, intent(in) :: ave_layer

      character(len=*), intent(in) :: cumulus

      real, intent(out) :: sc_up(:,:,:)
      real, intent(out) :: pw_up(:,:,:)
      real, intent(out) :: tot_pw_up_chem(:,:)

      !Local variables:
      real, dimension(mtp, its:ite) ::  sc_b
      real, dimension(mtp) :: conc_mxr
      real :: dz, xzz, xzd, xze, denom, henry_coef, w_upd, fliq, dp
      integer :: i, k, ispc
      real, dimension(mtp, its:ite, kts:kte) ::  factor_temp
      integer :: vtp_index

      !--initialization
      sc_up = se_cup
      pw_up = 0.0
      tot_pw_up_chem = 0.0

      if (USE_TRACER_SCAVEN == 2 .and. cumulus /= 'shallow') then
         factor_temp = 1.
          do vtp_index = 1, get_num_elements(vec_ok) ; i=get_data_value(vec_ok, vtp_index) !BD_n
            do ispc = 1, mtp
               ! - if tracer is type "carbon" then set coefficient to 0 for hydrophobic
               if (trim(chem_name(ispc) (1:len_trim('OCphobic'))) == 'OCphobic') factor_temp(ispc, :, :) = 0.0

               ! - suppress scavenging most aerosols at cold T except BCn1 (hydrophobic), dust, and HNO3
               if (trim(chem_name(ispc) (1:len_trim('BCphobic'))) == 'BCphobic') then
                  where (tempco < 258.) factor_temp(ispc, :, :) = 0.0
               end if

               if (trim(chem_name(ispc)) == 'sulfur' .or. &
                   trim(chem_name(ispc) (1:len_trim('ss'))) == 'ss' .or. & ! 'seasalt'
                   trim(chem_name(ispc)) == 'SO2' .or. &
                   trim(chem_name(ispc)) == 'SO4' .or. &
                   trim(chem_name(ispc)) == 'nitrate' .or. &
                   trim(chem_name(ispc)) == 'bromine' .or. &
                   trim(chem_name(ispc)) == 'NH3' .or. &
                   trim(chem_name(ispc)) == 'NH4a') then

                  where (tempco < 258.) factor_temp(ispc, :, :) = 0.0
               end if

            end do
         end do
      end if

      do vtp_index = 1, get_num_elements(vec_ok) ; i=get_data_value(vec_ok, vtp_index) !BD_n
         !start_level(i) = klcl(i)
         !start_level(i) = k22(i)

         do ispc = 1, mtp
            call get_cloud_bc(cumulus,ave_layer,kts,kte,ktf,xland(i), po(kts:kte,i),  se_cup(ispc,i,  kts:kte),sc_b(ispc,i),k22(i))

         end do
         do k = kts, start_level(i)
            sc_up(:,k,i) = sc_b(:,i)
            !sc_up   (:,k,i) = se_cup(:,k,i)
         end do
      end do

      do vtp_index = 1, get_num_elements(vec_ok) ; i=get_data_value(vec_ok, vtp_index) !BD_n
         loopk: do k = start_level(i) + 1, ktop(i) + 1

            !-- entr,detr, mass flux ...
            xzz = zuo(i,k-1)
            xzd = 0.5*up_massdetro(i,k-1)
            xze = up_massentro(i,k-1)
            denom = (xzz - xzd + xze) + 1.e-12

            !-- transport + mixing
            sc_up(:,k,i) = (sc_up(:,k-1,i)*xzz - sc_up(:,k-1,i)*xzd + se(:,k-1,i)*xze)/denom

            !-- scavenging section
            if (USE_TRACER_SCAVEN == 0 .or. cumulus == 'shallow') cycle loopk
            dz = z_cup(k,i) - z_cup(i,k-1)

            !-- in-cloud vert velocity for scavenging formulation 2
            !           w_upd = cte_w_upd
            !           w_upd = vvel1d(i)
            w_upd = vvel2d(k,i)

            do ispc = 1, mtp
               if (fscav(ispc) > 1.e-6) then ! aerosol scavenging

                  !--formulation 1 as in GOCART with RAS conv_par
                  if (USE_TRACER_SCAVEN == 1) pw_up(ispc,k,i) = max(0., sc_up(ispc,k,i)*(1.-exp(-fscav(ispc)*(dz/1000.))))

                  !--formulation 2 as in GOCART
                  if (USE_TRACER_SCAVEN == 2) pw_up(ispc,k,i) = max(0., sc_up(ispc,k,i) * (1.-exp(-chem_adj_autoc(ispc)*p_kc &
                                            * (dz/w_upd)))*factor_temp(ispc,k,i))

                  !--formulation 3 - orignal GF conv_par
                  if (USE_TRACER_SCAVEN == 3) then
                     !--- cloud liquid water tracer concentration
                     conc_mxr(ispc) = p_scav_eff*sc_up(ispc,k,i) !unit [kg(aq)/kg(air)]  for aerosol/smoke
                     !---   aqueous-phase concentration in rain water
                     pw_up(ispc,k,i) = conc_mxr(ispc)*pwo(k,i)/(1.e-8 + qrco(k,i))
                  end if

                  !---(in cloud) total mixing ratio in gas and aqueous phases
                  sc_up(ispc,k,i) = sc_up(ispc,k,i) - pw_up(ispc,k,i)

                  !
               elseif (hcts(ispc)%hstar > 1.e-6) then ! tracer gas phase scavenging

                  !--- equilibrium tracer concentration - Henry's law
                  henry_coef = Henry(ispc, tempco(k,i))

                  if (USE_TRACER_SCAVEN == 3) then
                     !--- cloud liquid water tracer concentration
                     conc_mxr(ispc) = (henry_coef*qrco(k,i)/(1.+henry_coef*qrco(k,i)))*sc_up(ispc,k,i)
                     !
                     !---   aqueous-phase concentration in rain water
                     pw_up(ispc,k,i) = conc_mxr(ispc)*pwo(k,i)/(1.e-8 + qrco(k,i))

                  else

                     !-- this the 'alpha' parameter in Eq 8 of Mari et al (2000 JGR) = X_aq/X_total
                     fliq = henry_coef*qrco(k,i)/(1.+henry_coef*qrco(k,i))

                     !---   aqueous-phase concentration in rain water
                     pw_up(ispc,k,i) = max(0., sc_up(ispc,k,i) &
                                       *(1.-exp(-fliq*chem_adj_autoc(ispc)*p_kc*dz/w_upd)))!*factor_temp(ispc,k,i))

                  end if

                  !---(in cloud) total mixing ratio in gas and aqueous phases
                  sc_up(ispc,k,i) = sc_up(ispc,k,i) - pw_up(ispc,k,i)

                  !
                  !---(in cloud)  mixing ratio in aqueous phase
                  !sc_up_aq(ispc,k,i) = conc_mxr(ispc) !if using set to zero at the begin.
               end if
            end do
            !
            !-- total aerosol/gas in the rain water
            dp = 100.*(po_cup(k,i) - po_cup(k+1,i))

            tot_pw_up_chem(:,i) = tot_pw_up_chem(:,i) + pw_up(:,k,i)*dp/c_grav
         end do loopk
         !
         !----- get back the in-cloud updraft gas-phase mixing ratio : sc_up(ispc,k)
         !          do k=start_level(i)+1,ktop(i)+1
         !            do ispc = 1,mtp
         !             sc_up(ispc,k,i) = sc_up(ispc,k,i) - sc_up_aq(ispc,k,i)
         !            enddo
         !          enddo
      end do
   end subroutine getInCloudScChemUp
    ! ---------------------------------------------------------------------------------------------------
   subroutine getInCloudScChemDd(cumulus, mtp, se, se_cup, sc_dn, pw_dn, tot_pw_up_chem &
                               , tot_pw_dn_chem, po_cup, pwdo, pwevo, zdo, dd_massentro &
                               , dd_massdetro, pwavo, jmin, ierr, itf, its, kts)
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
      implicit none
      !Parameters:
      character(len=*), parameter :: procedureName = 'getInCloudScChemDd' ! Subroutine Name
   
      !Variables (input, output, inout)
      integer, intent(in)  :: itf, its, kts, mtp

      integer, intent(in) :: ierr(:)
      integer, intent(in) :: jmin(:)

      real, intent(in) :: se(:,:,:)
      real, intent(in) :: se_cup(:,:,:)
      real, intent(in) :: pwavo(:)
      real, intent(in) :: pwevo(:)
      real, intent(in) :: po_cup(:,:)
      real, intent(in) :: pwdo(:,:)
      real, intent(in) :: zdo(:,:)
      real, intent(in) :: dd_massentro(:,:)
      real, intent(in) :: dd_massdetro(:,:)
      real, intent(in) :: tot_pw_up_chem(:,:)

      character(len=*), intent(in)  :: cumulus

      real, intent(out) :: sc_dn(:,:,:)
      real, intent(out) :: pw_dn(:,:,:)
      real, intent(out) :: tot_pw_dn_chem(:,:)

      !Local variables:
      real :: xzz, xzd, xze, denom, pwdper, frac_evap, dp
      integer :: i, k, ispc
      integer :: vtp_index

      sc_dn = 0.0
      pw_dn = 0.0
      tot_pw_dn_chem = 0.0
      if (cumulus == 'shallow') return

      do vtp_index = 1, get_num_elements(vec_ok) ; i=get_data_value(vec_ok, vtp_index) !BD_n

         !--- fration of the total rain that was evaporated
         frac_evap = -pwevo(i)/(1.e-16 + pwavo(i))

         !--- scalar concentration in-cloud - downdraft

         !--- at k=jmim
         k = jmin(i)
         pwdper = pwdo(k,i)/(1.e-16 + pwevo(i))*frac_evap  ! > 0
         if (USE_TRACER_EVAP == 0) pwdper = 0.0

         dp = 100.*(po_cup(k,i) - po_cup(k+1,i))

         do ispc = 1, mtp
            !--downdrafts will be initiate with a mixture of 50% environmental and in-cloud concentrations
            sc_dn(ispc,k,i) = se_cup(ispc,k,i)
            !sc_dn(ispc,k,i) = 0.9*se_cup(ispc,k,i)+0.1*sc_up(ispc,k,i)

            pw_dn(ispc,k,i) = -pwdper*tot_pw_up_chem(ispc, i)*c_grav/dp
            sc_dn(ispc,k,i) = sc_dn(ispc,k,i) - pw_dn(ispc,k,i)
            tot_pw_dn_chem(ispc, i) = tot_pw_dn_chem(ispc, i) + pw_dn(ispc,k,i)*dp/c_grav
         end do
         !
         !--- calculate downdraft mass terms
         do k = jmin(i) - 1, kts, -1
            xzz = zdo(k+1,i)
            xzd = 0.5*dd_massdetro(k,i)
            xze = dd_massentro(k,i)
            denom = (xzz - xzd + xze) + 1.e-12

            !-- transport + mixing
            sc_dn(:,k,i) = (sc_dn(:,k+1,i)*xzz - sc_dn(:,k+1,i)*xzd + se(:,k,i)*xze)/denom

            !
            !-- evaporation term
            if (USE_TRACER_EVAP == 0) cycle

            dp = 100.*(po_cup(k,i) - po_cup(k+1,i))

            !-- fraction of evaporated precip per layer
            pwdper = pwdo(k,i)/(1.e-16 + pwevo(i))! > 0

            !-- fraction of the total precip that was actually evaporated at layer k
            pwdper = pwdper*frac_evap

            !-- sanity check
            pwdper = min(1., max(pwdper, 0.))

            do ispc = 1, mtp
               !-- amount evaporated by the downdraft from the precipitation
               pw_dn(ispc,k,i) = -pwdper*tot_pw_up_chem(ispc,i)*c_grav/dp ! < 0. => source term for the downdraft tracer concentration

               !if(ispc==1) print*,"pw=",pwdper,tot_pw_up_chem (ispc,i),pwevo(i)/pwavo(i),pwdo(k,i)/(1.e-16+pwo(k,i))

               !-- final tracer in the downdraft
               sc_dn(ispc,k,i) = sc_dn(ispc,k,i) - pw_dn(ispc,k,i) ! observe that -pw_dn is > 0.

               !-- total evaporated tracer
               tot_pw_dn_chem(ispc,i) = tot_pw_dn_chem(ispc,i) + pw_dn(ispc,k,i)*dp/c_grav

               !print*,"to=",k,tot_pw_dn_chem(ispc,i),pwdo(k,i)/(1.e-16+pwevo(i)),frac_evap,tot_pw_dn_chem(ispc,i)/tot_pw_up_chem (ispc,i)

            end do
         end do
         !
      end do
   end subroutine getInCloudScChemDd
   !---------------------------------------------------------------------------------------------------
   function henry(ispc,temp) result(henry_coef)
      !--- calculate Henry's constant for solubility of gases into cloud water
      !--- inputs : ak0(ispc), dak(ispc),  hstar(ispc), dhr(ispc)
      implicit none
      integer, intent(in) :: ispc
      real   , intent(in) :: temp
      real :: henry_coef
      real :: fct ,tcorr, corrh

      !--- define some constants!
      ! real, parameter:: c_rgas_atm = 8.205e-2 ! atm M^-1 K^-1 ! 8.314 gas constant [J/(mol*K)]
      ! real, parameter:: c_avogad= 6.022e23! Avogadro constant [1/mol]
      ! real, parameter:: c_rhoH2O= 999.9668! density of water [kg/m3]
      ! real, parameter:: c_temp0 = 298.15! standard temperature [K]
      ! real, parameter:: c_temp0i= 1./298.15! inverse of standard temperature [K]
      ! real, parameter:: c_MWH2O = 18.02! molecular mass of water [kg/kmol]
      ! real, parameter:: c_MWAIR = 28.97! effective molecular mass of air [kg/kmol]
      ! real, parameter:: c_conv3 = c_avogad / 1.0e6!  [mol(g)/m3(air)]  to [molec(g)/cm3(air)]
      ! real, parameter:: c_conv4 = 100.            !  [m]   to [cm]
      ! real, parameter:: c_conv5 = 1000.           !  [m^3]      to [l]
      ! real, parameter:: c_conv7 = 1/c_conv5       !  [l]    to [m^3]
      ! real, parameter:: c_conv6 = 1. / 101325.    !  [Pa]      to [atm]
      ! real, parameter:: c_hplus = 1.175e-4        ! for cloud water. pH is assumed to be 3.93: pH=3.93 
      !                                             !  =>hplus=10**(-pH)

      ! aqueous-phase concentrations XXXa [mol/m3(air)]!
      ! gas-phase concentrations XXXg [mol/m3(air)]!
      ! Henry constants XXXh for scavenging [mol/(l*atm)]!
      ! converted to [(mol(aq)/m3(aq))/(mol(g)/m3(air))], i.e. dimensionless!
      ! in equilibrium XXXa = XXXh * LWC * XXXg!
      tcorr = 1./temp - c_temp0i

      !-P. Colarco corrected the expression below
      !fct  = conv7 * c_rgas_atm * temp ! - for henry_coef in units 1/m3
      fct   =         c_rgas_atm * temp ! - for henry_coef dimensioless


      !-taking into account the acid dissociation constant
      ! ak=ak0*exp(dak*(1/t-1/298))
      corrh=1.+Hcts(ispc)%ak0 * exp(Hcts(ispc)%dak * tcorr)/c_hplus

      !-- for concentration in mol[specie]/mol[air] - Eq 5 in 'Compilation of Henry's law constants (version 4.0) for
      !-- water as solvent, R. Sander, ACP 2015'.
      henry_coef =  Hcts(ispc)%hstar* exp(Hcts(ispc)%dhr*tcorr) * fct * corrh

   end function henry
   !---------------------------------------------------------------------------------------------------
   subroutine fluxFormsSourceSink(kte, kts, mtp, ktop, chem_name_mask, dtime, edto, po_cup, pw_dn_chem, pw_up_chem, sc_dn_chem &
                              ,   sc_up_chem ,se_chem ,se_cup_chem ,zo_cup ,zdo ,zuo, cumulus, zenv ,out_chem)
      !! ## brief
      !!
      !! Author: Rodrigues, L. F. [LFR]
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 24Fevereiro2023 14:47
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! brief
      !!
      !! ** History**:
      !!
      !! ---
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!

      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_name = 'fluxFormsSourceSink'
      !! subroutine name

      ! Variables (input, output, inout)
      integer, intent(in) ::  kte, kts, mtp, ktop
      integer, intent(in) :: chem_name_mask(:)

      real, intent(in) :: dtime
      real, intent(in) :: edto
      real, intent(in) :: po_cup(:)
      real, intent(in) :: pw_dn_chem(:,:)
      real, intent(in) :: pw_up_chem(:,:)
      real, intent(in) :: sc_dn_chem(:,:)
      real, intent(in) :: sc_up_chem(:,:)
      real, intent(in) :: se_chem(:,:)
      real, intent(in) :: se_cup_chem(:,:)
      real, intent(in) :: zo_cup(:)
      real, intent(in) :: zdo(:)
      real, intent(in) :: zuo(:)

      character(len=*), intent(in) :: cumulus

      real, intent(inout) :: zenv(:)

      real, intent(out) :: out_chem(:,:)

      ! Local variables:
      integer :: k, ispc, istep, lstep
      real :: dp, alp0, dtime_max, beta1, wetdep, evap, dz
      real, dimension(mtp, kts:kte) :: trcflx_in, sub_tend, zenv_diff
      real, dimension(kts:kte) :: massflx
      real, dimension(kts:kte) :: aa, bb, cc, fp, fm
      real, dimension(mtp, kts:kte) :: ddtr, ddtr_upd, fp_mtp, fm_mtp
      real, dimension(mtp) :: trash_, trash2_, evap_, residu_, wetdep_

      !- flux form + source/sink terms + time explicit + FCT
      if (USE_FLUX_FORM == 1 .and. ALP1 == 0.) then

         if (USE_FCT == 0) then
            do k = kts, ktop
               dp=100.*(po_cup(k)-po_cup(k+1))

               out_chem(:,k) =-(zuo(k+1)*(sc_up_chem(:,k+1)-se_cup_chem(:,k+1) ) -                 &
                                zuo(k  )*(sc_up_chem(:,k  )-se_cup_chem(:,k  ) ))*c_grav/dp        &
                              +(zdo(k+1)*(sc_dn_chem(:,k+1)-se_cup_chem(:,k+1) ) -                 &
                                zdo(k  )*(sc_dn_chem(:,k  )-se_cup_chem(:,k  ) ))*c_grav/dp*edto
            end do

         else

            !-- FCT scheme for the subsidence transport: d(M_env*S_env)/dz
            sub_tend  = 0.
            trcflx_in = 0.
            dtime_max = dtime
            massflx(:)= 0.

            do k=kts+1,ktop+1
               dp              = 100.*(po_cup(k)-po_cup(k+1))
               trcflx_in (:,k) =-(zuo(k)  -edto*zdo(k))*se_cup_chem(:,k) !* xmb
               massflx     (k) =-(zuo(k)  -edto*zdo(k))           !* xmb
               dtime_max       = min(dtime_max,.5*dp)
            end do
            !- if dtime_max<dtime => needs a loop to update from t to t+dtime (check this!)
            !if( dtime_max < dtime ) stop "dtime_max < dtime in GF scheme"

            do ispc = 1, mtp
              ! sub_tend(ispc, :) = Fct1d3(ktop, kte, dtime_max, po_cup( :), se_chem(ispc,  :) &
              !                         , massflx( :), trcflx_in(ispc, :))
              call fct1d3 (ktop,kte,dtime_max,po_cup(:),se_chem(ispc,:),massflx(:),trcflx_in(ispc,:),sub_tend(ispc,:))

            end do

            do k = kts, ktop
               dp = 100.*(po_cup(k) - po_cup(k+1))
               out_chem(:,k) = -(zuo(k+1)*(sc_up_chem(:,k+1)) - zuo(k)*(sc_up_chem(:,k)))*c_grav/dp     &
                               +(zdo(k+1)*(sc_dn_chem(:,k+1)) - zdo(k)*(sc_dn_chem(:,k)))*c_grav/dp*edto

               !- update with the subsidence term from FCT scheme
               out_chem(:,k) = out_chem(:,k) + sub_tend(:,k)

            end do
         end if

         !- include evaporation (this term must not be applied to the tracer 'QW')
         if (USE_TRACER_EVAP == 1 .and. trim(cumulus) /= 'shallow') then
            do k = kts, ktop
               dp = 100.*(po_cup(k) - po_cup(k+1))
               out_chem(:,k) = out_chem(:,k) - 0.5*edto*(zdo(k)*pw_dn_chem(:,k) &
                             + zdo(k+1) * pw_dn_chem(:,k+1)) * c_grav/dp   ! evaporated ( pw_dn < 0 => E_dn > 0)
               !*chem_name_mask_evap(:) !-- to avoid the "Dry Mass Violation"
            end do
         end if

         !- include scavenging
         if (USE_TRACER_SCAVEN > 0 .and. trim(cumulus) /= 'shallow') then
            do k = kts, ktop
               dp = 100.*(po_cup(k) - po_cup(k+1))
               out_chem(:,k) = out_chem(:,k) - 0.5*(zuo(k)*pw_up_chem(:,k) &
                             + zuo(k+1)*pw_up_chem(:,k+1))*c_grav/dp  ! incorporated in rainfall (<0)
            end do
         end if
      end if ! IF(USE_FLUX_FORM == 1 .and. ALP1 == 0. )

      !- flux form + source/sink terms + time explicit/implicit + upstream
      if (USE_FLUX_FORM == 1 .and. ALP1 > 0.) then

         alp0 = 1.-ALP1
         do k = kts, ktop + 1
            fp(k) = 0.5*(zenv(k) + abs(zenv(k)))
            fm(k) = 0.5*(zenv(k) - abs(zenv(k)))
         end do

         do k = kts, ktop
            dp = 100.*(po_cup(k) - po_cup(k+1))
            beta1 = dtime*c_grav/dp
            aa(k) = ALP1*beta1*fm(k)
            bb(k) = 1.+ALP1*beta1*(fp(k) - fm(k+1))
            cc(k) = -ALP1*beta1*fp(k+1)

            ddtr(:,k) = se_chem(:,k) - (zuo(k+1)*sc_up_chem(:,k+1) - zuo(k)*sc_up_chem(:,k))*beta1      &
                                     + (zdo(k+1)*sc_dn_chem(:,k+1) - zdo(k)*sc_dn_chem(:,k))*beta1*edto


            !- include evaporation (this term must not be applied to the tracer 'QW')
            if (USE_TRACER_EVAP == 1 .and. trim(cumulus) /= 'shallow') then
               out_chem(:,k) = out_chem(:,k) - 0.5*edto*(zdo(k)*pw_dn_chem(:,k) &
                             + zdo(k+1) * pw_dn_chem(:,k+1)) * beta1 !&  ! evaporated ( pw_dn < 0 => E_dn > 0)
               !*chem_name_mask_evap(:) !-- to avoid the "Dry Mass Violation"
            end if

            !- include scavenging
            if (USE_TRACER_SCAVEN > 0 .and. trim(cumulus) /= 'shallow') then
               out_chem(:,k) = out_chem(:,k) - 0.5*(zuo(k)*pw_up_chem(:,k) &
                             + zuo(k+1)*pw_up_chem(:,k+1))*beta1  ! incorporated in rainfall (<0)
            end if

            ddtr(:, k) = ddtr(:, k) + out_chem(:,k) + alp0*beta1*(-fm(k) &
                       * se_chem(:,  max(kts, k-1)) + (fm(k+1) - fp(k))*se_chem(:,k) + fp(k+1)*se_chem(:,k+1))

         end do

         do ispc = 1, mtp
            if (chem_name_mask(ispc) == 0) cycle
            call tridiag(ktop, aa(kts:ktop), bb(kts:ktop), cc(kts:ktop), ddtr(ispc, kts:ktop))
            out_chem(ispc,  kts:ktop) = (ddtr(ispc, kts:ktop) - se_chem(ispc,  kts:ktop))/dtime
         end do
      end if !USE_FLUX_FORM == 1 .and. ALP1 > 0.


      !- flux form + source/sink terms + time explicit + upstream with anti-diffusion step (Smolarkiewicz 1983)
      if (USE_FLUX_FORM == 2 .or. USE_FLUX_FORM == 3) then
         if (USE_FLUX_FORM == 2) lstep = -1 ! upstream + anti-diffusion step
         if (USE_FLUX_FORM == 3) lstep =  1 ! only upstream
         alp0 = 1.

         !--- Zenv here have the following reference:  < 0 => downward motion
         zenv(:) = -(zuo(:) - edto*zdo(:))

         do istep = 1, lstep, -2

            if (istep == 1) then
               ddtr_upd(:,:) = se_chem(:,:)
               do k = kts, ktop + 1
                  fp_mtp(:,k) = 0.5*(zenv(k) + abs(zenv(k)))
                  fm_mtp(:,k) = 0.5*(zenv(k) - abs(zenv(k)))
               end do

            else

               ddtr_upd(:,kts:ktop + 1) = se_chem(:,kts:ktop + 1) + out_chem(:,kts:ktop + 1)*dtime
               zenv_diff(:,kts) = 0.
               do k = kts, ktop + 1
                  dz = zo_cup(k+1) - zo_cup(k)
                  zenv_diff(:, k+1) = 1.08*(dz*abs(zenv(k+1)) - dtime*zenv(k+1)**2) &
                                            *(ddtr_upd(:,k+1) - ddtr_upd(:,k)) &
                                           /((ddtr_upd(:,k+1) + ddtr_upd(:,k) + 1.e-16)*dz)
               end do
               do k = kts, ktop + 1
                  fp_mtp(:,k) = 0.5*(zenv_diff(:,k) + abs(zenv_diff(:,k)))
                  fm_mtp(:,k) = 0.5*(zenv_diff(:,k) - abs(zenv_diff(:,k)))
               end do

            end if

            do k = kts, ktop
               dp = -100.*(po_cup(k) - po_cup(k+1))
               beta1 = dtime*c_grav/dp
               
               ddtr(:,k) = ddtr_upd(:,k) + alp0*beta1*( &
                          (fp_mtp(:,k+1)*ddtr_upd(:,k)           +fm_mtp(:,k+1)*ddtr_upd(:,k+1)) &
                         -(fp_mtp(:,k  )*ddtr_upd(:,max(kts,k-1))+fm_mtp(:,k  )*ddtr_upd(:,k  )) )
            end do
            do ispc = 1, mtp
               if (chem_name_mask(ispc) == 0) cycle
               out_chem(ispc,kts:ktop) = (ddtr(ispc,kts:ktop) - se_chem(ispc,kts:ktop))/dtime
            end do

         end do ! anti-diff steps

         do k = kts, ktop
            dp = 100.*(po_cup(k) - po_cup(k+1))
            beta1 = c_grav/dp

            out_chem(:,k) = out_chem(:,k) - (zuo(k+1)*sc_up_chem(:,k+1) - zuo(k)*sc_up_chem(:,k))*beta1 &
                                          + (zdo(k+1)*sc_dn_chem(:,k+1) - zdo(k)*sc_dn_chem(:,k))*beta1*edto

            !- include evaporation (this term must not be applied to the tracer 'QW')
            if (USE_TRACER_EVAP == 1 .and. trim(cumulus) /= 'shallow') then
               out_chem(:,k) = out_chem(:,k) - 0.5*edto*(zdo(k)*pw_dn_chem(:,k) &
                                             +zdo(k+1) * pw_dn_chem(:,k+1)) * beta1 !&  ! evaporated ( pw_dn < 0 => E_dn > 0)
               !*chem_name_mask_evap(:) !-- to avoid the "Dry Mass Violation"
            end if

            !- include scavenging
            if (USE_TRACER_SCAVEN > 0 .and. trim(cumulus) /= 'shallow') then
               out_chem(:,k) = out_chem(:,k) - 0.5*(zuo(k)*pw_up_chem(:,k) &
                                        + zuo(k+1)*pw_up_chem(:,k+1))*beta1  ! incorporated in rainfall (<0)
            end if
         end do
      end if ! USE_FLUX_FORM == 2 .or. USE_FLUX_FORM == 3
      !
   end subroutine fluxFormsSourceSink
!--------------------------------------------------------------------------------------------!
!- end of section for atmospheric composition
!---------------------------------------------------------------------------------------------!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
! -- for GEOS-5   - save for later use
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

   !---------------------------------------------------------------------------------------------------
   subroutine GF_GEOS5_INTERFACE(mxp,myp,mzp,mtp,ITRCR,LONS,LATS,DT_MOIST          &
      ,T, PLE, PLO, ZLE, ZLO, PK,  U, V, OMEGA , KH      &
      ,TH1, Q1, U1, V1 ,QLCN ,QICN,QLLS,QILS, CNPCPRATE  &
      ,CNV_MF0, CNV_PRC3, CNV_MFD, CNV_DQLDT ,ENTLAM     &
      ,CNV_MFC, CNV_UPDF, CNV_CVW, CNV_QC, CLCN,CLLS     &
      ,QV_DYN_IN,PLE_DYN_IN,U_DYN_IN,V_DYN_IN,T_DYN_IN   &
      ,RADSW   ,RADLW ,DQDT_BL  ,DTDT_BL                 &
      ,FRLAND  ,AREA  ,USTAR ,TSTAR ,QSTAR ,T2M ,Q2M     &
      ,TA      ,QA    ,SH    ,EVAP  ,PHIS                &
      ,KPBLIN                                            &
      ,STOCHASTIC_SIG, SIGMA_DEEP, SIGMA_MID             &
      ,DQDT_GF,DTDT_GF,MUPDP,MUPSH,MUPMD,MDNDP           &
      ,MFDP,MFSH,MFMD,ERRDP,ERRSH,ERRMD                  &
      ,AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC      &
      ,DTDTDYN,DQVDTDYN                                  &
      ,NCPL, NCPI, CNV_NICE, CNV_NDROP,CNV_FICE,CLDMICRO &
      ,TRACER,FSCAV,CNAMES,QNAMES,DTRDT_GF               &
      ,RSU_CN,REV_CN, PFI_CN, PFL_CN                     &
      ,TPWI,TPWI_star,LIGHTN_DENS                        &
      ,VAR3d_a,VAR3d_b,VAR3d_c,VAR3d_d                   &
       !,CNV_TR &!,VAR2d,ZKBCON
      )

      implicit none
      !INCLUDE "mpif.h"
      real, parameter :: MAPL_GRAV   = 9.80665                ! m^2/s
      real, parameter :: MAPL_PI     = 3.14159265358979323846
      logical, parameter :: FEED_3DMODEL   = .true.  != set "false" to not feedback the AGCM with the
                                                     != heating/drying/transport conv tendencies

      character(len=*),intent(in) :: CLDMICRO !set two-moment microphysics

      integer ,intent(in) :: mxp,myp,mzp,mtp,ITRCR

      real                             ,intent(in)   :: DT_moist

      real   ,dimension(mxp,myp)       ,intent(in)   :: FRLAND ,AREA ,USTAR ,TSTAR ,QSTAR &
         ,T2M ,Q2M ,TA ,QA ,SH ,EVAP ,PHIS  &
         ,KPBLIN,LONS,LATS &
         ,STOCHASTIC_SIG
      real   ,dimension(mxp,myp)       ,intent(out)  :: SIGMA_DEEP, SIGMA_MID

      real   ,dimension(mxp,myp,0:mzp) ,intent(in)   :: PLE,ZLE,PLE_DYN_IN
      real   ,dimension(mxp,myp,mzp)   ,intent(in)   :: T ,U ,V ,ZLO ,PLO ,PK ,OMEGA, KH       &
         ,RADSW  ,RADLW  ,DQDT_BL  ,DTDT_BL      &
         ,QV_DYN_IN,U_DYN_IN,V_DYN_IN,T_DYN_IN   &
         ,DTDTDYN,DQVDTDYN

      real   ,dimension(mxp,myp,mzp)   ,intent(inout):: TH1,Q1,U1,V1,QLCN ,QICN, NCPL, NCPI    &
         ,QLLS,QILS,CLLS
      real   ,dimension(mxp,myp,mzp)   ,intent(inout):: RSU_CN,REV_CN,VAR3d_a,VAR3d_b,VAR3d_c,VAR3d_d

      real   ,dimension(mxp,myp,0:mzp) ,intent(inout):: PFI_CN, PFL_CN

      !   REAL   ,DIMENSION(mxp,myp,mzp)   ,INTENT(INOUT):: CNV_TR
      real   ,dimension(mxp,myp,mzp)                 :: CNV_TR

      real   ,dimension(mxp,myp,0:mzp) ,intent(out)  :: CNV_MFC
      real   ,dimension(mxp,myp,mzp)   ,intent(out)  :: CNV_MF0 ,CNV_PRC3,CNV_MFD,CNV_DQLDT  &
         ,CNV_UPDF, CNV_CVW, CNV_QC,CLCN,ENTLAM&
         ,CNV_NICE, CNV_NDROP, CNV_FICE
      !-for debug purposes
      real   ,dimension(mxp,myp,mzp)   ,intent(inout):: DQDT_GF,DTDT_GF,MUPDP,MDNDP,MUPSH,MUPMD ,DTRDT_GF
      real   ,dimension(mxp,myp)       ,intent(inout):: MFDP,MFSH,MFMD,ERRDP,ERRSH,ERRMD
      real   ,dimension(mxp,myp)       ,intent(inout):: AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC &
         ,TPWI, TPWI_star

      real   ,dimension(MXP,MYP)       ,intent(out)  :: CNPCPRATE ,LIGHTN_DENS
      real   ,dimension(MXP,MYP)                     :: VAR2d,ZKBCON

      real   ,dimension(mxp,myp,mzp,itrcr) ,intent(inout)   :: TRACER  !=XHO in grid_moist_comp.f90
      real   ,dimension(itrcr)             ,intent(in   )   :: FSCAV
      character(len=*)  ,dimension(mtp)    ,intent(in   )   :: CNAMES,QNAMES

      real   ,dimension(mxp,myp,mzp)                        :: frct_liq,AA1_ADV,AA1_RADPBL

      integer      :: ims,ime, jms,jme, kms,kme,    &
         its,ite, jts,jte, kts,kte,    &
         mynum

      real,  dimension(mzp , mxp, myp ) :: up       &
         ,vp       &
         ,wp       &
         ,rvap     &
         ,temp     &
         ,press    &
         ,zm3d     &
         ,zt3d     &
         ,dm3d     &
         ,curr_rvap&
         ,buoy_exc &
         ,khloc



      real,  dimension(mzp , mxp, myp ) ::          &
           gsf_t   & ! grid-scale forcing for temp
         , gsf_q   & ! advection forcing for rv
         ,advf_t   & ! advection forcing for temp
         ,sgsf_t   & ! sub-grid scale forcing for temp
         ,sgsf_q   & ! sub-grid scale forcing for rv
         ,SRC_T    & ! temp tendency      from convection
         ,SRC_Q    & ! rv tendency        from convection
         ,SRC_CI   & ! cloud/ice tendency from convection
         ,SRC_U    & ! U tendency         from convection
         ,SRC_V    & ! V tendency         from convection
         ,SRC_NI   & ! Ice     number tendency from convection
         ,SRC_NL   & ! Droplet number tendency from convection
         ,SRC_BUOY & ! buoyancy tendency from downdrafts
         ,REVSU_GF & ! evaporation_or_sublimation of_convective_precipitation kg kg-1 s-1
         ,PRFIL_GF & ! ice_or_liq convective_precipitation flux: kg m2 s-1 (deep only)
         ,VAR3d_aGF& ! dummy 3-d var for output
         ,VAR3d_bGF& ! dummy 3-d var for output
         ,VAR3d_cGF& ! dummy 3-d var for output
         ,VAR3d_dGF& ! dummy 3-d var for output
         ,qexcp    & ! placeholder for Q   ex from cold pool param
         ,hexcp    & ! placeholder for MSE ex from cold pool param
         ,cnvcf    & ! placeholder
         ,turb_len_scale  ! placeholder


      real,  dimension(nmp, mzp , mxp, myp ) ::     &
          mp_ice   &
         ,mp_liq   &
         ,mp_cf

      real,  dimension(nmp, mzp , mxp, myp ) ::     &
          SUB_MPQI & ! subsidence transport applied to ice mix ratio
         ,SUB_MPQL & ! subsidence transport applied to cloud mix ratio
         ,SUB_MPCF   ! subsidence transport applied to cloud fraction


      real,  allocatable, dimension(:,:,:,:) :: SRC_CHEM ! tracer mixing ratio tendencies from the parameterized convection

      real,  dimension(mtp)       :: FSCAV_INT
      character(len=100)          :: AER_CHEM_MECH

      real,    dimension(mxp,myp) :: CONPRR
      real                        :: time=0.
      real,    dimension(mxp,myp) ::  aot500  ,temp2m  ,sfc_press &
         ,sflux_r ,sflux_t ,topt      &
         ,xland   ,dx2d    ,water_bud &
         ,col_sat, tke_pbl,rh_dicy_fct&
         ,wlpool
      integer, dimension(mxp,myp) :: kpbl,do_this_column
      integer, dimension(mzp)     :: flip
      integer :: k,i,j,iens,ispc, itime1=0
      !- for convective transport
      integer, dimension(mxp,myp,maxiens) ::  &
          ierr4d                       &
         ,jmin4d                       &
         ,klcl4d                       &
         ,k224d                        &
         ,kbcon4d                      &
         ,ktop4d                       &
         ,kstabi4d                     &
         ,kstabm4d

      real,dimension(mxp,myp,maxiens)     ::  &
          cprr4d                       &
         ,xmb4d                        &
         ,edt4d                        &
         ,pwav4d                       &
         ,sigma4d
      real,dimension(mxp,myp,mzp,maxiens) ::  &
          pcup5d                       &
         ,up_massentr5d                &
         ,up_massdetr5d                &
         ,dd_massentr5d                &
         ,dd_massdetr5d                &
         ,zup5d                        &
         ,zdn5d                        &
         ,prup5d                       &
         ,prdn5d                       &
         ,clwup5d                      &
         ,tup5d                        &
         ,conv_cld_fr5d

      !-----------local var in GEOS-5 data structure
      real,  dimension(mxp, myp, mzp) :: T1, ZLO_N,PLO_N,PK_N,TH_N,Q_N,T_N,U_N,V_N,DM
      real,  dimension(mxp,myp,0:mzp) :: ZLE_N,PLE_N
      integer :: status,alloc_stat ,wantrank=-99999
      real    :: tem1,dz,air_dens, src_cnvtr,snk_cnvtr,dz_int,tau_cp
      character(len=10) :: ENV_SETTING='DEFAULT'! ! 'CURRENT'/'DEFAULT'
      integer, parameter :: itest=1!3
      real :: RL, RI, disp_factor,x1,x2

      !-- some initialization
      do_this_column =0
      ierr4d         =0
      jmin4d         =0
      klcl4d         =0
      k224d          =0
      kbcon4d        =0
      ktop4d         =0
      kstabi4d       =0
      kstabm4d       =0
      xmb4d          =0.
      cprr4d         =0.
      edt4d          =0.
      pwav4d         =0.
      sigma4d        =0.
      pcup5d         =0.
      up_massentr5d  =0.
      up_massdetr5d  =0.
      dd_massentr5d  =0.
      dd_massdetr5d  =0.
      zup5d          =0.
      zdn5d          =0.
      prup5d         =0.
      prdn5d         =0.
      clwup5d        =0.
      tup5d          =0.
      conv_cld_fr5d  =0.
      CNV_NDROP      =0.
      CNV_NICE       =0.
      CNV_FICE       =0.
      SRC_NI         =0.
      SRC_NL         =0.
      SRC_T          =0.
      SRC_Q          =0.
      SRC_CI         =0.
      SRC_U          =0.
      SRC_V          =0.
      CNPCPRATE      =0.
      SUB_MPQI       =0.
      SUB_MPQL       =0.
      SUB_MPCF       =0.
      LIGHTN_DENS    =0.
      rh_dicy_fct    =0.
      SRC_BUOY       =0.
      REVSU_GF       =0.
      PRFIL_GF       =0.
      VAR3d_aGF      =0.
      VAR3d_bGF      =0.
      VAR3d_cGF      =0.
      VAR3d_dGF      =0.
      VAR2d          =0.
      turb_len_scale =0.
      !-
      !---temporary settings for debugging purposes
      !- special setting for SCM runs
      if(mxp==1 .and. myp==1 .and. maxval(T2m) < 1.e-6) return

      !- special setting for SCM runs
      if(mxp>1 .and. myp>1) wrtgrads = .false.
      !call mpi_comm_rank(MPI_COMM_WORLD,WHOAMI_ALL,status)
      !----
      if(p_use_gate) stop "p_use_gate must be false for GEOS5 runs"
      !if( .not. p_use_gate .and. wrtgrads) call alloc_grads_arr(1,mzp,1,jl)
      !--------------------------------------------------------

      !- time counter
      ntimes = ntimes + 1
      !if(ntimes==1 .and. WHOAMI_ALL == 0) print *,'==> Applying GF convection scheme '
      mynum = -999
      call set_index_loops( ims,ime, jms,jme, kms,kme,    &
         its,ite, jts,jte, kts,kte,    &
         mxp,myp,mzp                   )

      !- define the vector "flip" to invert the z-axis orientation
      call flipz(flip,mzp)
      !
      if(.not.allocated(SRC_CHEM)) then
         allocate(SRC_CHEM(mtp, mzp, mxp, myp),stat=alloc_stat) !- tendency from convection
         if(alloc_stat==0) SRC_CHEM=0.0
      endif
      if(USE_TRACER_TRANSP==1) then
         AER_CHEM_MECH='GOCART' !in the future set as intent(in) from MoistGridComp
         call interface_aerchem(mtp,itrcr,aer_chem_mech, cnames,qnames, fscav, fscav_int)
      endif

      !
      !- 2-d input data
      aot500  (:,:) = 0.1  ! #
      !- as moist is called before surface, at the 1st time step all arrays
      !- from surface are zero
      if(maxval(T2m) < 1.e-6) then
         temp2m   (:,:) = T  (:,:,mzp) ! Kelvin
      else
         temp2m   (:,:) = T2M(:,:) ! or TA(:,:) ! Kelvin
      endif
      !- moisture flux from sfc
      sflux_r  (:,:) = EVAP(:,:) ! kg m-2 s-1
      !- sensible heat flux (sh) comes in W m-2, below it is converted to K m s-1
      !-(air_dens_sfc = ple(:,:,mzp)/( 287.04*TA(:,:)*(1.+0.608*QA(:,:)))))
      sflux_t  (:,:) = SH  (:,:) /(1004. * ple(:,:,mzp)/(287.04*T(:,:,mzp)*(1.+0.608*Q1(:,:,mzp)))) ! K m s-1
      !- topography height  (m)
      topt     (:,:) = PHIS(:,:)/MAPL_GRAV
      !- land/ocean fraction: land if < 1 ,ocean if = 1
      xland    (:,:) = 1.0-FRLAND(:,:)
      !
      !- grid length for the scale awareness
      dx2d(:,:) = sqrt(AREA(:,:)) ! meters
      !- special setting for SCM runs
      if(mxp==1 .and. myp==1) dx2d = 100000.

      !-pbl heigth index
      do j=1,myp
         do i=1,mxp
            if (nint(KPBLIN(i,j)) /= 0) then
               kpbl(i,j) = max(1, flip(min( nint(KPBLIN(i,j)), mzp)))
            else
               kpbl(i,j) = 1
            endif
            tke_pbl(i,j) = p_tkmin ! dummy 
         enddo
      enddo
      !
      !- 3-d input data
      !- current temperature T1 (after dyn+every physics process called before moist)
      T1 = PK * TH1
      !- any var with index "1" (and omega and pk) are already updated with dynamics
      !  tendencies and everything else (from physics) that was called before moist
      !
      if(trim(env_setting)=='CURRENT') then
         PLO_N = PLO
         T_N   = T1
         DM = ( PLE(:,:,1:mzp)-PLE(:,:,0:mzp-1) )*(1./MAPL_GRAV)
         !- 1st setting: enviromental state is the one already modified by dyn + physics
         do j=1,myp
            do i=1,mxp
               do k=1,mzp
                  temp  (k,i,j) = T1    (i,j,flip(k))
                  press (k,i,j) = PLO   (i,j,flip(k))*100. !Pa
                  rvap  (k,i,j) = Q1    (i,j,flip(k))!
                  up    (k,i,j) = U1    (i,j,flip(k))! already @ A-grid (m/s)
                  vp    (k,i,j) = V1    (i,j,flip(k))! already @ A-grid (m/s)
                  wp    (k,i,j) = OMEGA (i,j,flip(k))! Pa/s
                  zt3d  (k,i,j) = ZLO   (i,j,flip(k))! mid -layer level
                  zm3d  (k,i,j) = ZLE   (i,j,flip(k))! edge-layer level
                  dm3d  (k,i,j) = DM    (i,j,flip(k))
                  khloc (k,i,j) = KH    (i,j,flip(k))
                  curr_rvap(k,i,j) = Q1 (i,j,flip(k))!current rvap

                  mp_ice(lsmp,k,i,j) = QILS  (i,j,flip(k))
                  mp_liq(lsmp,k,i,j) = QLLS  (i,j,flip(k))
                  mp_cf (lsmp,k,i,j) = CLLS  (i,j,flip(k))
                  mp_ice(cnmp,k,i,j) = QICN  (i,j,flip(k))
                  mp_liq(cnmp,k,i,j) = QLCN  (i,j,flip(k))
                  mp_cf (cnmp,k,i,j) = CLCN  (i,j,flip(k))

               enddo
            enddo
         enddo
         !- sfc pressure (Pa)
         sfc_press(:,:) = PLE(:,:,mzp)
         !- Grid and sub-grid scale forcings for convection
         do j=1,myp
            do i=1,mxp
               do k=1,mzp
                  gsf_t (k,i,j) = 0.
                  gsf_q (k,i,j) = 0.
                  sgsf_t (k,i,j) = RADSW  (i,j,flip(k))+ RADLW(i,j,flip(k)) + DTDT_BL(i,j,flip(k))
                  sgsf_q (k,i,j) = DQDT_BL(i,j,flip(k))
                  advf_t (k,i,j) = 0.
               enddo
            enddo
         enddo
         !.. !- Tracer transport/scavenging section
         if(USE_TRACER_TRANSP==1) then
            where(TRACER<=0.0) TRACER = p_mintracer
         endif
      elseif(trim(env_setting)=='DEFAULT') then
         !-2nd setting: environmental state is that one before any tendency
         !- is applied (i.e, at begin of each time step).
         !- Get back the model state, heights and others variables at time N
         !- (or at the beggining of current time step)
         !- In physics, the state vars (T,U,V,PLE) are untouched and represent the
         !- model state after dynamics phase 1. But, "Q" is modified by physics, so
         !- depending on what was called before this subroutine, "Q" may be already
         !- changed from what it was just after dynamics phase 1. To solve this issue,
         !- "Q" just after dynamics is saved in the var named "QV_DYN_IN" in "GEOS_AgcmGridComp.F90".
         Q_N   =  QV_DYN_IN
         T_N   =   T_DYN_IN
         U_N   =   U_DYN_IN
         V_N   =   V_DYN_IN
         PLE_N = PLE_DYN_IN

         DM = ( PLE_N(:,:,1:mzp)-PLE_N(:,:,0:mzp-1) )*(1./MAPL_GRAV)
         !DM = ( PLE  (:,:,1:mzp)-PLE  (:,:,0:mzp-1) )*(1./MAPL_GRAV)
         call get_vars(mzp,mxp,myp,Q_N,T_N,PLE_N,ZLE_N,ZLO_N,PLO_N,PK_N,TH_N)
         !
         do j=1,myp
            do i=1,mxp
               do k=1,mzp
                  !
                  temp       (k,i,j) = T_N   (i,j,flip(k))! (K)
                  press      (k,i,j) = PLO_N (i,j,flip(k))*100.! (Pa) @ mid-layer level
                  rvap       (k,i,j) = Q_N   (i,j,flip(k))! water vapor mix ratio
                  up         (k,i,j) = U_N   (i,j,flip(k))! already @ A-grid (m/s)
                  vp         (k,i,j) = V_N   (i,j,flip(k))! already @ A-grid (m/s)
                  wp         (k,i,j) = OMEGA (i,j,flip(k))! (Pa/s)
                  zt3d       (k,i,j) = ZLO_N (i,j,flip(k))! mid -layer level (m)
                  zm3d       (k,i,j) = ZLE_N (i,j,flip(k))! edge-layer level (m)
                  dm3d       (k,i,j) = DM    (i,j,flip(k))
                  khloc      (k,i,j) = KH    (i,j,flip(k))
                  curr_rvap  (k,i,j) = Q1    (i,j,flip(k)) ! current rvap (dyn+phys)
                  mp_ice(lsmp,k,i,j) = QILS  (i,j,flip(k))
                  mp_liq(lsmp,k,i,j) = QLLS  (i,j,flip(k))
                  mp_cf (lsmp,k,i,j) = CLLS  (i,j,flip(k))
                  mp_ice(cnmp,k,i,j) = QICN  (i,j,flip(k))
                  mp_liq(cnmp,k,i,j) = QLCN  (i,j,flip(k))
                  mp_cf (cnmp,k,i,j) = CLCN  (i,j,flip(k))
               enddo
            enddo
         enddo
         !- sfc pressure (Pa)
         sfc_press(:,:) = PLE_N(:,:,mzp)
         !- Grid and sub-grid scale forcings for convection
         do j=1,myp
            do i=1,mxp
               do k=1,mzp
                  advf_t (k,i,j) = DTDTDYN (i,j,flip(k))
                  gsf_t  (k,i,j) = DTDTDYN (i,j,flip(k)) + RADLW(i,j,flip(k)) + RADSW  (i,j,flip(k))
                  gsf_q  (k,i,j) = DQVDTDYN(i,j,flip(k))

                  sgsf_t (k,i,j) = DTDT_BL (i,j,flip(k))
                  sgsf_q (k,i,j) = DQDT_BL (i,j,flip(k))
               enddo
            enddo
         enddo
         !
         !- Tracer transport/scavenging section
         if(USE_TRACER_TRANSP==1) then
            where(TRACER<=0.0) TRACER = p_mintracer
         endif
      else
         stop 'unknown env_setting at convpar_gf_geos5.F90'
      endif

      if(CONVECTION_TRACER==1) then
         do j=1,myp
            do i=1,mxp
               !--- saturation CWV
               col_sat(i,j) = TPWI(i,j)/(1.e-6+TPWI_star(i,j))
               col_sat(i,j) = min(1.,max(0.,col_sat(i,j)))
               !--- temporary to hold only CWV in mm
               ! col_sat(i,j) = TPWI(i,j)
               !
               do k=1,mzp
                  buoy_exc(k,i,j)=CNV_TR (i,j,flip(k))
               enddo
            enddo
         enddo
         qexcp = 0. 
         hexcp = 0. 
         wlpool= 0.
        !print*,"buoy_exc1=",maxval(buoy_exc),minval(buoy_exc)
      else
         buoy_exc = 0.0
      endif

      !- call the driver routine to apply the parameterization
      call convParGFDriver(mxp,myp,mzp,mtp,nmp, time, itime1   &
         ,ims,ime, jms,jme, kms,kme   &
         ,its,ite, jts,jte, kts,kte   &
         ,flip        &
         ,fscav_int   &
         ,mynum       &
         ,dt_moist    &
         ,dx2d        &
         ,stochastic_sig &
         ,zm3d        &
         ,zt3d        &
         ,dm3d        &
         !--- sfc inputs
         ,lons        &
         ,lats        &
         ,aot500      &
         ,temp2m      &
         ,sflux_r     &
         ,sflux_t     &
         ,qexcp       &
         ,hexcp       &
         ,wlpool      &
         ,topt        &
         ,xland       &
         ,sfc_press   &
         ,kpbl        &
         ,tke_pbl     &
         ,turb_len_scale &
         !--- atmos state
         ,col_sat     &
         ,up          &
         ,vp          &
         ,wp          &
         ,temp        &
         ,press       &
         ,rvap        &
         ,mp_ice      &
         ,mp_liq      &
         ,mp_cf       &
         ,curr_rvap   &
         !--- atmos composition state
         ,TRACER      & !- note: uses GEOS-5 data structure
         ,cnvcf       & ! conv cloud fraction
         !---- forcings---
         ,buoy_exc    &
         ,gsf_t       &
         ,gsf_q       &
         ,advf_t      &
         ,sgsf_t      &
         ,sgsf_q      &
         !---- output ----
         ,conprr      &
         ,lightn_dens &
         ,rh_dicy_fct &
         ,src_t       &
         ,src_q       &
         ,src_ci      &
         ,src_nl      &
         ,src_ni      &
         ,src_u       &
         ,src_v       &
         ,sub_mpqi    &
         ,sub_mpql    &
         ,sub_mpcf    &
         ,src_buoy    &
         ,src_chem    &
         ,revsu_gf    &
         ,prfil_gf    &
         !
         !
         ,do_this_column&
         ,ierr4d       &
         ,jmin4d       &
         ,klcl4d       &
         ,k224d        &
         ,kbcon4d      &
         ,ktop4d       &
         ,kstabi4d     &
         ,kstabm4d     &
         ,cprr4d       &
         ,xmb4d        &
         ,edt4d        &
         ,pwav4d       &
         ,sigma4d      &
         ,pcup5d       &
         ,up_massentr5d&
         ,up_massdetr5d&
         ,dd_massentr5d&
         ,dd_massdetr5d&
         ,zup5d        &
         ,zdn5d        &
         ,prup5d       &
         ,prdn5d       &
         ,clwup5d      &
         ,tup5d        &
         ,conv_cld_fr5d&
         !-- for debug/diagnostic
         ,AA0,AA1,AA1_ADV,AA1_RADPBL,AA1_BL,AA2,AA3,AA1_CIN,TAU_BL,TAU_EC &
         ,VAR2d,VAR3d_aGF,VAR3d_bGF,VAR3d_cGF,VAR3d_dGF)


      if(FEED_3DMODEL)then
         !--- vertical fraction of liq/ice water phases
         if(FRAC_MODIS==1 .and. icumulus_gf(DEEP)==ON) then
            do j=1,myp
               do i=1,mxp
                  do k=1,mzp
                     frct_liq(i,j,k) = FractLiqF(tup5d(i,j,flip(k),DEEP))
                  enddo
               enddo
            enddo
         else
            do j=1,myp
               do i=1,mxp
                  do k=1,mzp
                     frct_liq(i,j,k) = FractLiqF(T(i,j,k))
                  enddo
               enddo
            enddo
         endif
         !-- update GEOS-5 model state with the feedback from cumulus convection
         !- to include the tendencies from the convection,  update the vars th1,q1,v1 and u1
         do j=1,myp
            do i=1,mxp
               if(do_this_column(i,j) == 0) cycle
               !- conv precip rate: mm/s = kg m-2 s-1
               CNPCPRATE (i,j) =  CONPRR (i,j)

               if(ITEST==0) CNPCPRATE(i,j) =  0.

               !
               do k=1,mzp ! in the future, limit the vertical loop to ktop (DO k=mzp,flip(ktop),-1)?
                  !- convert from d temp/dt to d theta/dt using PK => d theta/dt = (1/pk)*d temp/dt
                  !- (think if PK must be the current one _or_ at the begin of the time step
                  TH1(i,j,k) = TH1(i,j,k) + DT_moist * SRC_T(flip(k),i,j) / PK(i,j,k)

                  Q1 (i,j,k) = Q1 (i,j,k) + DT_moist * SRC_Q(flip(k),i,j)
                  Q1 (i,j,k) = max(p_smaller_qv, Q1 (i,j,k))
                  !
                  !- simple splitting of condensate tendency into liq/ice phases
                  !- these are 'anvil' mixing ratio and not 'grid-scale' mix ratio
                  !- (the convective source will be applied in progno_cloud, routine "consrc")
                  ! QLCN (i,j,k) = QLCN (i,j,k) + DT_moist * SRC_CI(flip(k),i,j) * frct_liq(i,j,k)
                  ! QICN (i,j,k) = QICN (i,j,k) + DT_moist * SRC_CI(flip(k),i,j) * (1.0-frct_liq(i,j,k))

                  if(ITEST==3) then
                     !- simple splitting of condensate tendency into liq/ice phases
                     !- these are 'grid-scale' mix ratio
                     !- (the convective source will be set to zero, see below)
                     QLLS (i,j,k) = QLLS (i,j,k) + DT_moist * SRC_CI(flip(k),i,j) *  frct_liq(i,j,k)
                     QILS (i,j,k) = QILS (i,j,k) + DT_moist * SRC_CI(flip(k),i,j) *  (1.-frct_liq(i,j,k))
                  endif
               enddo
               !--- sublimation/evaporation tendencies (kg/kg/s)
               do k=1,mzp
                  !--- sublimation/evaporation tendencies (kg/kg/s)
                  RSU_CN (i,j,k) = REVSU_GF(flip(k),i,j)* (1.0-frct_liq(i,j,k))
                  REV_CN (i,j,k) = REVSU_GF(flip(k),i,j)*  frct_liq(i,j,k)
                  !--- preciptation fluxes (kg/kg/s)
                  PFI_CN (i,j,k) = PRFIL_GF(flip(k),i,j)* (1.0-frct_liq(i,j,k))
                  PFL_CN (i,j,k) = PRFIL_GF(flip(k),i,j)*  frct_liq(i,j,k)
               enddo

            enddo
         enddo
         !-----
         if(USE_MOMENTUM_TRANSP > 0) then
            do j=1,myp
               do i=1,mxp
                  if(do_this_column(i,j) == 0) cycle
                  do k=1,mzp
                     U1 (i,j,k) = U1 (i,j,k) + DT_moist * SRC_U(flip(k),i,j)
                     V1 (i,j,k) = V1 (i,j,k) + DT_moist * SRC_V(flip(k),i,j)
                  enddo
               enddo
            enddo
         endif


         if(APPLY_SUB_MP == 1) then
            do j=1,myp
               do i=1,mxp
                  if(do_this_column(i,j) == 0) cycle
                  do k=1,mzp
                     QLLS (i,j,k) = QLLS (i,j,k) + DT_moist * SUB_MPQL(LSMP,flip(k),i,j)
                     QILS (i,j,k) = QILS (i,j,k) + DT_moist * SUB_MPQI(LSMP,flip(k),i,j)
                     CLLS (i,j,k) = CLLS (i,j,k) + DT_moist * SUB_MPCF(LSMP,flip(k),i,j)
                     QLCN (i,j,k) = QLCN (i,j,k) + DT_moist * SUB_MPQL(CNMP,flip(k),i,j)
                     QICN (i,j,k) = QICN (i,j,k) + DT_moist * SUB_MPQI(CNMP,flip(k),i,j)
                     CLCN (i,j,k) = CLCN (i,j,k) + DT_moist * SUB_MPCF(CNMP,flip(k),i,j)
                  enddo
               enddo
            enddo
         endif

         if(USE_TRACER_TRANSP==1) then

            do j=1,myp
               do i=1,mxp
                  if(do_this_column(i,j) == 0) cycle
                  do k=1,mzp
                     !-special array for output of CO tendency
                     DTRDT_GF(i,j,k)=SRC_CHEM(ispc_CO,flip(k),i,j)

                     !- update tracer mass mixing ratios
                     do ispc=1,mtp

                        TRACER(i,j,k,ispc)=TRACER(i,j,k,ispc)+ DT_moist * SRC_CHEM(ispc,flip(k),i,j) * CHEM_NAME_MASK(ispc)!

                        !-- final check for negative tracer mass mixing ratio
                        TRACER(i,j,k,ispc)=max(p_mintracer, TRACER(i,j,k,ispc))
                     enddo
                  enddo
               enddo
            enddo
           !-- final check for negative tracer mass mixing ratio
           !where (TRACER < p_mintracer) TRACER = p_mintracer

         endif

         !--- for lightning flash rate
         !--- ice/liq precip fluxes
         if(icumulus_gf(deep) == ON) then
            do j=1,myp
               do i=1,mxp
                  if(ierr4d(i,j,deep) .ne. 0) cycle
                  ZKBCON(i,j) = ZLE(i,j,flip(kbcon4d(i,j,DEEP)))
                  do k=mzp,1,-1
                     PFI_CN (i,j,k) = prfil_gf(flip(k),i,j)* (1.0-frct_liq(i,j,k))
                     PFL_CN (i,j,k) = prfil_gf(flip(k),i,j)*  frct_liq(i,j,k)
                  enddo
               enddo
            enddo
         endif

         !--- for dummy output
         do j=1,myp
            do i=1,mxp
               do k=1,mzp
                  VAR3d_a  (i,j,k) = VAR3d_aGF(flip(k),i,j)
                  VAR3d_b  (i,j,k) = VAR3d_bGF(flip(k),i,j)
                  VAR3d_c  (i,j,k) = VAR3d_cGF(flip(k),i,j)
                  VAR3d_d  (i,j,k) = VAR3d_dGF(flip(k),i,j)
               enddo
            enddo
         enddo


         do IENS=1, maxiens
            if(icumulus_gf(IENS) == ON) then
               do j=1,myp
                  do i=1,mxp
                     if(ierr4d(i,j,IENS) .ne. 0) cycle
                     do k=mzp,flip(ktop4d(i,j,IENS))-1,-1

                        DZ       = -( ZLE(i,j,k) - ZLE(i,j,k-1) )
                        air_dens = 100.*plo_n(i,j,k)/(287.04*T_n(i,j,k)*(1.+0.608*Q_n(i,j,k)))

                        !- special treatment for CNV_DQLDT: 'convective_condensate_source',  UNITS     = 'kg m-2 s-1',
                        !- SRC_CI contains contributions from deep, shallow,... . So, do not need to be accumulated over  CNV_DQLDT
                        !- note the SRC_CI has different array structure (k,i,j) _not_ (i,j,k)
                        if(ITEST /= 3) &
                           CNV_DQLDT(i,j,k)=  SRC_CI(flip(k),i,j)*DZ * air_dens !units: kg[w]/(kg[air] s) * m * kg[air]/m^3 = kg[w]/(m^2 s)

                        !CNV_DQLDT(i,j,k)=  SRC_CI(flip(k),i,j)*DM(i,j,k)     !units: kg[w]/(kg[air] s) * m * kg[air]/m^3 = kg[w]/(m^2 s)
                        !
                        !-'detraining_mass_flux',UNITS     = 'kg m-2 s-1',
                        CNV_MFD (i,j,k) = CNV_MFD (i,j,k) + ( up_massdetr5d(i,j,flip(k),IENS) )


                        !-'cloud_base_mass_flux',    units = 'kg m-2 s-1',                                  &
                        CNV_MF0 (i,j,k) = CNV_MF0 (i,j,k) + zup5d(i,j,flip(k),IENS)

                        !-convective mass flux [kg m^{-2} s^{-1}] - with downdrafts
                        ! CNV_MFC (i,j,k) = CNV_MFC (i,j,k) +  ( zup5d(i,j,flip(k),IENS) + edt4d(i,j,IENS)*  zdn5d(i,j,flip(k),IENS))

                        !---only updraft
                        CNV_MFC (i,j,k) = CNV_MFC (i,j,k) + zup5d(i,j,flip(k),IENS)

                        if(zup5d(i,j,flip(k),IENS) > 1.0e-6) then
                           !-'entrainment parameter',  UNITS     ='m-1',
                           ENTLAM  (i,j,k) =  ENTLAM   (i,j,k) + (up_massentr5d(i,j,flip(k),IENS)/(DZ*zup5d(i,j,flip(k),IENS)))

                           !-'updraft_vertical_velocity',            UNITS     = 'hPa s-1',
                           CNV_CVW (i,j,k) = -0.2 ! hPa/s =>  4 m/s
                        endif

                        !-'grid_mean_convective_condensate', UNITS     ='kg kg-1'
                        CNV_QC  (i,j,k) = CNV_QC  (i,j,k) + clwup5d(i,j,flip(k),IENS)
                        !
                        !
                        !~ !----------------------------------------------------------------------------------------------------
                        !- not using progno-cloud to calculate the precip from the convective column
                        !- if CNV_PRC3 will be send to progno-cloud, set CNPCPRATE = zero
                        !-'convective_precipitation_from_GF',UNITS     = 'kg m-2 s-1',
                        !- JAN/17/2017 : the units above are wrong. The correct are kg[precip water]/kg[air]
                        CNV_PRC3(i,j,k) = CNV_PRC3(i,j,k) +                 (prup5d(i,j,flip(k),IENS) + &
                           edt4d(i,j,IENS)* prdn5d(i,j,flip(k),IENS) ) &
                           * DT_moist/(DZ*air_dens)

                        !-'updraft_areal_fraction',
                        if(zup5d(i,j,flip(k),IENS) > 1.0e-6) CNV_UPDF(i,j,k) = 0.033
                        !----------------------------------------------------------------------------------------------------
                        if(ITEST==2) then
                           !-'updraft_areal_fraction',
                           if(zup5d(i,j,flip(k),IENS) > 1.0e-6) then
                              CNV_UPDF(i,j,k) = 0.033
                           else
                              CNV_UPDF(i,j,k) = 0.0
                           endif

                           !-'convective_cloud_area_fraction', adimensional
                           !- Tiedtke formulation
                           CLCN(i,j,k) = CLCN(i,j,k) + (1.0-CLCN(i,j,k))*(up_massdetr5d(i,j,flip(k),IENS) &
                              * DT_moist/(DZ*air_dens)+CNV_UPDF(i,j,k))

                           !- Chab&Bechtold 2002/2005 formulation
                           ! CLCN(i,j,k) = CLCN(i,j,k) + (1.0-CLCN(i,j,k))*(conv_cld_fr5d(i,j,flip(k),IENS)+CNV_UPDF(i,j,k))
                           !
                           CLCN(i,j,k) = max(0.,min(CLCN(i,j,k),0.99))
                        endif
                     !----------------------------------------------------------------------------------------------------
                     enddo
                   !print*,"iens=",iens,maxval(conv_cld_fr5d(i,j,:,IENS)),minval(conv_cld_fr5d(i,j,:,IENS));call flush(6)
                  enddo
               enddo
            endif
         enddo
      endif

      if(adjustl(CLDMICRO) =="2MOMENT") then
         !- we adjust convective cloud condensate and number here
         do j=1,myp
            do i=1,mxp
               do k=1,mzp

                  !---obsolete
                  !tem1 = T(i,j,k)
                  !RL =   10.0  + (12.0*(283.0- tem1)/40.0)
                  !RL =   min(max(RL, 10.0), 30.0)*1.e-6
                  !RI =   100.0 + (80.0*(tem1- 253.0)/40.0)
                  !RI =   min(max(RI, 20.0), 250.0)*1.e-6
                  !tem1 = 1.- (tem1 - 235.0) /38.0
                  !tem1 =  min(max(0.0, tem1), 1.0)
                  ! make up some "number" sources. In the future this should depend explicitly on the convective mphysics
                  !disp_factor =  10.0 ! used to account somehow for the size dist
                  !SRC_NL(flip(k),i,j) = SRC_CI(flip(k),i,j)* (1.0-tem1) /(1.333 * MAPL_PI*RL*RL*RL*997.0*disp_factor)
                  !SRC_NI(flip(k),i,j)= SRC_CI(flip(k),i,j) * tem1 /(1.333 * MAPL_PI *RI*RI*RI*500.0*disp_factor)
                  !CNV_FICE (i, j, k)   =   tem1
                  !---obsolete

                  tem1 = frct_liq(i,j,k)

                  QLCN (i,j,k) = QLCN (i,j,k) + DT_moist * SRC_CI(flip(k),i,j) * (1.0-tem1)
                  QICN (i,j,k) = QICN (i,j,k) + DT_moist * SRC_CI(flip(k),i,j) * tem1

                  NCPL (i,j,k) = NCPL (i,j,k) + DT_moist * SRC_NL(flip(k),i,j)
                  NCPI (i,j,k) = NCPI (i,j,k) + DT_moist * SRC_NI(flip(k),i,j)

                  DZ       = -( ZLE(i,j,k) - ZLE(i,j,k-1) )
                  air_dens = 100.*PLO_n(i,j,k)/(287.04*T_n(i,j,k)*(1.+0.608*Q_n(i,j,k)))

                  if (CNV_MFD (i,j,k)  .gt. 0.) then
                     CNV_NICE (i,j,k)=   SRC_NI(flip(k),i,j)*DZ * air_dens/CNV_MFD (i,j,k)
                     CNV_NDROP(i,j,k)=   SRC_NL(flip(k),i,j)*DZ * air_dens/CNV_MFD (i,j,k)
                  endif
               enddo
            enddo
         enddo
         !--- special section for convective cloud fraction
         do iens=1,maxiens
            if(icumulus_gf(iens) .ne. ON) cycle
            do j=1,myp
               do i=1,mxp
                  if(ierr4d(i,j,IENS) .ne. 0) cycle
                  do k=mzp,flip(ktop4d(i,j,IENS))-1,-1

                     DZ       = -( ZLE(i,j,k) - ZLE(i,j,k-1) )
                     air_dens = 100.*plo_n(i,j,k)/(287.04*T_n(i,j,k)*(1.+0.608*Q_n(i,j,k)))
                     CLCN(i,j,k) = CLCN(i,j,k) + (1.0-CLCN(i,j,k))*(up_massdetr5d(i,j,flip(k),iens) &
                        * dt_moist/(dz*air_dens)+CNV_UPDF(i,j,k))
                     CLCN(i,j,k) = max(0.,min(CLCN(i,j,k),0.99))
                  enddo
               enddo
            enddo
         enddo
      endif !2 moment
      !
      !--- cold pool/"convection tracer"
      if(CONVECTION_TRACER==1) then
         DTRDT_GF=0. !temporary   for output only
         do j=1,myp
            do i=1,mxp

               tau_cp=FRLAND(i,j)*tau_land_cp + (1.0-FRLAND(i,j))*tau_ocea_cp

               do k=1,mzp

                  !- sink term (exp decay 1h)
                  snk_cnvtr =  dt_moist * abs(CNV_TR (i,j,k))/tau_cp

                  !DZ          = -( ZLE(i,j,k) - ZLE(i,j,k-1) )
                  !air_dens = 100.*PLO_N(i,j,k)/(287.04*T_n(i,j,k)*(1.+0.608*Q_n(i,j,k)))
                  !
                  !- downdraft convective mass flux [kg m^{-2} s^{-1}]
                  ! iens =?
                  !src_cnvtr =  edt4d(i,j,iens)*zdn5d(i,j,flip(k),iens)

                  !- downdraft detrainment mass flux [kg m^{-2} s^{-1}]
                  ! iens =?
                  !src_cnvtr =  edt4d(i,j,iens)*dd_massdetr5d(i,j,flip(k),iens)


                  !- source term
                  !- downdraft detrainment of buoyancy [ J/kg s^{-1}]
                  !- negative sign => source for updraft lifting
                  src_cnvtr = - dt_moist * min(0.,SRC_BUOY(flip(k),i,j))

                  !- 'continuity' equation = ADV + SRC - SINK
                  CNV_TR (i,j,k) = CNV_TR (i,j,k)  + src_cnvtr  - snk_cnvtr

                  !temporary for output only
                  DTRDT_GF(i,j,k)=SRC_BUOY(flip(k),i,j)

               enddo
            enddo
         enddo
        !print*,"buoy_exc2=",maxval(SRC_BUOY),minval(SRC_BUOY)
      endif

      !
      if(maxval(icumulus_gf)>0) then
         do IENS=1, maxiens
            if(icumulus_gf(IENS) == ON .and. IENS== DEEP) then
               do j=1,myp
                  do i=1,mxp
                     if(ierr4d(i,j,DEEP) /= 0) cycle
                     MFDP      (i,j)      =   xmb4d(i,j,DEEP)
                     SIGMA_DEEP(i,j)      = sigma4d(i,j,deep)
                     MUPDP     (i,j,1:mzp)=zup5d(i,j,flip(1):flip(mzp):-1,DEEP)
                     MDNDP     (i,j,1:mzp)=zdn5d(i,j,flip(1):flip(mzp):-1,DEEP) * edt4d(i,j,IENS)
                  enddo
               enddo
            elseif(icumulus_gf(IENS) == ON .and. IENS== SHAL) then
               do j=1,myp
                  do i=1,mxp
                     if(ierr4d(i,j,SHAL) /= 0) cycle
                     MFSH (i,j)      =xmb4d(i,j,SHAL)
                     MUPSH(i,j,1:mzp)=zup5d(i,j,flip(1):flip(mzp):-1,SHAL)
                  enddo
               enddo
            elseif(icumulus_gf(IENS) == ON .and. IENS== MID) then
               do j=1,myp
                  do i=1,mxp
                     if(ierr4d(i,j,MID) /= 0) cycle
                     MFMD      (i,j)      = cprr4d (i,j,MID) ! xmb4d(i,j,MID)temporary saving for mid precip
                     SIGMA_MID (i,j)      = sigma4d(i,j,MID )
                     MUPMD     (i,j,1:mzp)=zup5d(i,j,flip(1):flip(mzp):-1,MID)
                  enddo
               enddo
            endif
         enddo
         !- for output purposes, ierr=0 (convection is ON) will be changed to 1
         where(ierr4d==0)ierr4d=1
         where(ierr4d>1)ierr4d=0
         do j=1,myp
            do i=1,mxp
               DQDT_GF(i,j,1:mzp)=SRC_Q(flip(1):flip(mzp):-1,i,j)!note SQR_Q is (k,i,j)
               DTDT_GF(i,j,1:mzp)=SRC_T(flip(1):flip(mzp):-1,i,j)!note SQR_T is (k,i,j)

               ERRDP(i,j)=float(ierr4d(i,j,DEEP))
               ERRSH(i,j)=float(ierr4d(i,j,SHAL))
               ERRMD(i,j)=float(ierr4d(i,j,MID ))
            enddo
         enddo
      endif
      !
      if( allocated(src_chem))  deallocate(src_chem,stat=alloc_stat) !tendency   from convection

      !- for debugging purposes only
      !if(.not. p_use_gate .and. wrtgrads) call alloc_grads_arr(1,mzp,2,jl)

   end subroutine GF_GEOS5_INTERFACE
   !------------------------------------------------------------------------------------
   subroutine get_vars(LM,mxp,myp,Q,T,PLE,ZLE,ZLO,PLO,PK,TH)
      implicit none
      integer, intent(in) :: LM,mxp,myp
      real, intent(in) , dimension(mxp,myp,0:LM) :: PLE
      real, intent(in) , dimension(mxp,myp,LM)   :: T,Q

      real, intent(out), dimension(mxp,myp,0:LM) :: ZLE
      real, intent(out), dimension(mxp,myp,LM)   :: ZLO,PLO,PK,TH

      !-local var
      real, parameter:: MAPL_GRAV   = 9.80665                ! m^2/s
      real, parameter:: MAPL_AIRMW  = 28.965                 ! kg/Kmole
      real, parameter:: MAPL_H2OMW  = 18.015                 ! kg/Kmole
      real, parameter:: MAPL_RUNIV  = 8314.47                ! J/(Kmole K)
      real, parameter:: MAPL_RDRY   = MAPL_RUNIV/MAPL_AIRMW  ! J/(kg K)
      real, parameter:: MAPL_CPDRY  = 3.5*MAPL_RDRY          ! J/(kg K)
      real, parameter:: MAPL_KAPPA  = MAPL_RDRY/MAPL_CPDRY   ! (2.0/7.0)
      real, parameter:: MAPL_EPSILON= MAPL_H2OMW/MAPL_AIRMW  ! --
      real, parameter:: MAPL_RGAS   = MAPL_RDRY              ! J/(kg K) (DEPRECATED)
      real, parameter:: MAPL_CP     = MAPL_RGAS/MAPL_KAPPA   ! J/(kg K) (DEPRECATED)
      real, parameter:: MAPL_VIREPS = 1.0/MAPL_EPSILON-1.0   !          (DEPRECATED)
      integer :: L
      real, dimension(mxp,myp,0:LM) :: CNV_PLE,PKE

      CNV_PLE  = PLE*.01
      PLO      = 0.5*(CNV_PLE(:,:,0:LM-1) +  CNV_PLE(:,:,1:LM  ) )
      PKE      = (CNV_PLE/1000.)**(MAPL_RGAS/MAPL_CP)
      PK       = (PLO/1000.)**(MAPL_RGAS/MAPL_CP)
      TH       = T/PK

      ZLE(:,:,LM) = 0.
      do L=LM,1,-1
         ZLE(:,:,L-1) = TH (:,:,L) * (1.+MAPL_VIREPS*Q(:,:,L))
         ZLO(:,:,L  ) = ZLE(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PKE(:,:,L)-PK (:,:,L  ) ) * ZLE(:,:,L-1)
         ZLE(:,:,L-1) = ZLO(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PK (:,:,L)-PKE(:,:,L-1) ) * ZLE(:,:,L-1)
      end do

        !.. IF( MAPL_AM_I_ROOT()) then
             !.. print*,"1get-vars =============================================================="
             !.. do L=LM,1,-1
                 !.. print*,"PLE/PLO",L,PLO(1,1,L),PLE(1,1,L)
             !.. end do
             !.. print*,"PLE/PLO",0,PLO(1,1,1),PLE(1,1,0)
             !.. print*,"2get-vars =============================================================="
             !.. call flush(6)
        !.. ENDIF

   end subroutine get_vars
   !---------------------------------------------------------------------------------------------------
   subroutine interface_aerchem(mtp,itrcr,aer_chem_mech, cnames,qnames, fscav, fscav_int)
      implicit none
      integer, intent(in) :: mtp,ITRCR
      character(len=*),                intent(in)   :: AER_CHEM_MECH
      character(len=*),dimension(mtp) ,intent(in)   :: CNAMES,QNAMES
      real            ,dimension(ITRCR) ,intent(in) :: FSCAV

      real,    dimension(mtp) , intent(out)   :: FSCAV_INT
      !-local vars
      integer :: ispc, len_ACM, len_spc, irun=0
      character(len=100) :: TMP_AER_NAME
      ispc_CO             = 1
      CHEM_NAME_MASK      = 1
      CHEM_NAME_MASK_EVAP = 1
      CHEM_ADJ_AUTOC      = 1.0

      !- GOCART + PCHEM section
      if(AER_CHEM_MECH=='GOCART') then
         len_ACM=len(TRIM(AER_CHEM_MECH))
         do ispc=1,mtp
            FSCAV_INT(ispc)   = fscav(ispc)
            len_spc           = len(TRIM(cnames(ispc)))
            CHEM_name (ispc)  = TRIM(qnames(ispc))!(len_ACM+1:len_spc))
            if(TRIM(CHEM_name (ispc)) == 'CO') ispc_CO=ispc

            !-the tracers below are not being transported :
            if(trim(chem_name(ispc)) == "NCPI"    ) CHEM_NAME_MASK     (ispc) = 0
            if(trim(chem_name(ispc)) == "NCPL"    ) CHEM_NAME_MASK     (ispc) = 0
            if(trim(chem_name(ispc)) == "QGRAUPEL") CHEM_NAME_MASK     (ispc) = 0
            if(trim(chem_name(ispc)) == "QRAIN"   ) CHEM_NAME_MASK     (ispc) = 0
            if(trim(chem_name(ispc)) == "QSNOW"   ) CHEM_NAME_MASK     (ispc) = 0
            !if(trim(chem_name(ispc)) == "QW"      ) CHEM_NAME_MASK     (ispc) = 0
            !if(trim(chem_name(ispc)) == "QW"      ) CHEM_NAME_MASK_EVAP(ispc) = 0

            if(trim(chem_name(ispc)(1:3)) == "du0") CHEM_ADJ_AUTOC     (ispc) = 1.0
            if(trim(chem_name(ispc)(1:3)) == "ss0") CHEM_ADJ_AUTOC     (ispc) = 1.0
            if(trim(chem_name(ispc)(1:3)) == "OCp") CHEM_ADJ_AUTOC     (ispc) = 1.0
            if(trim(chem_name(ispc)(1:3)) == "BCp") CHEM_ADJ_AUTOC     (ispc) = 1.0
            if(trim(chem_name(ispc)(1:3)) == "SO4") CHEM_ADJ_AUTOC     (ispc) = 1.0
         enddo
         !-----------------------------------------------------------------------------------------
         !---temporary section to fill Henrys cts for N2O and CH4 of PCHEM chemical mechanism
         do ispc=1,mtp
            if (TRIM(CHEM_name (ispc)) == 'N2O' .or. &
               TRIM(CHEM_name (ispc)) == 'CH4'      &
               )then

               call getHenrysLawCts(TRIM(CHEM_name (ispc)), &
                  Hcts(ispc)%hstar,Hcts(ispc)%dhr,Hcts(ispc)%ak0,Hcts(ispc)%dak)
            endif
         enddo
        !***   all tracers not having wet deposition _must_ have all Hcts reset to zero here.
        !.. WHERE( Hcts(:)%hstar == -1.) Hcts(:)%hstar = 0.
        !.. WHERE( Hcts(:)%dhr   == -1.) Hcts(:)%dhr   = 0.
        !.. WHERE( Hcts(:)%ak0   == -1.) Hcts(:)%ak0   = 0.
        !.. WHERE( Hcts(:)%dak   == -1.) Hcts(:)%dak   = 0.
      !-----------------------------------------------------------------------------------------
      endif
      return
      !IF( MAPL_AM_I_ROOT() .and. irun == 0)THEN
      irun = 1
      print*,"=========================================================================";call flush(6)
      print*," the following tracers will be transport by GF scheme                    ";call flush(6)

      write(10,*)"================= table of tracers in the GF conv transport ==================="
      write(10,*)" SPC,  CHEM_NAME,  FSCAV          - the four Henrys law cts  -   Transport Flag - kc adjust"
      do ispc=1,mtp
         write(10,121) ispc,trim(chem_name(ispc)), FSCAV_INT(ispc),Hcts(ispc)%hstar &
            ,Hcts(ispc)%dhr,Hcts(ispc)%ak0,Hcts(ispc)%dak, CHEM_NAME_MASK(ispc),CHEM_ADJ_AUTOC (ispc)
         if( CHEM_NAME_MASK     (ispc) == 1) then
            print*,"GF is doing transp and wet removal of: ",trim(chem_name(ispc))
            call flush(6)
         endif
      enddo
      print*,"=========================================================================";call flush(6)
121   format(1x,I4,A10,5F14.5,I4,F14.5)
      !ENDIF
   end subroutine interface_aerchem
   !------------------------------------------------------------------------------------
   subroutine writetxt(mzp,t,ple,th1,pk,q1,u1,zle,zlo      &
      ,DYNF_Q ,DYNF_T ,DYNF_PLE, DYNF_UA  &
      !
      ,theta,pp,rv,up,zm3d,zt3d,vp,omega&
      ,sflux_r,sflux_t,topt,xland,sfc_press,dx,kpbl,temp2m,dt_moist&
      )

      implicit none
      integer, intent(in)    :: mzp
      real, dimension(mzp)   :: t ,th1 ,pk ,q1 ,u1 ,zlo &
         ,DYNF_Q ,DYNF_T , DYNF_UA

      real, dimension(mzp)   :: theta ,pp ,rv ,up ,zm3d ,zt3d,vp,omega
      real, dimension(0:mzp) :: ple ,zle ,DYNF_PLE
      real :: sflux_r,sflux_t,topt,xland,sfc_press,dx,temp2m,dt_moist

      integer :: k,kpbl

      write(8,*) "================================================"
      write(7,*) "================================================"
      write(7,*) kpbl,sflux_r,sflux_t,topt,xland,sfc_press,dx,temp2m,dt_moist
      do k=1,mzp
         write(8,10) k,ple(k),t(k),th1(k),pk(k),1000.*q1(k),u1(k),zle(k),zlo(k),86400.*DYNF_Q(k)
         write(7,11) k,theta(k),pp(k),1000.*rv(k),up(k),zm3d(k),zt3d(k),vp(k),omega(k)
      enddo
      call flush(7)
      call flush(8)
10    format(1x,i4,9F11.3)
11    format(1x,i4,8F11.3)
   end subroutine writetxt
  !------------------------------------------------------------------------------------

   subroutine FLIPZ(flip,mzp)
      implicit none
      integer, intent(in) :: mzp
      integer, dimension(mzp), intent(inout) :: flip
      integer :: m,k
      m=mzp
      do k=1,mzp
         flip(k)=m
         m=m-1
      enddo
   end subroutine FLIPZ
   !------------------------------------------------------------------------------------

   subroutine set_index_loops( ims,ime, jms,jme, kms,kme,    &
      its,ite, jts,jte, kts,kte,    &
      mxp,myp,mzp                   )

      implicit none
      integer, intent(in)         :: mxp,myp,mzp
      integer, intent(inout)      :: ims,ime, jms,jme, kms,kme,&
         its,ite, jts,jte, kts,kte


      ims=1
      ime=mxp
      jms=1
      jme=myp
      kms=1
      kme=mzp
      its=1
      ite=mxp
      jts=1
      jte=myp
      kts=1
      kte=mzp

   end subroutine set_index_loops
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
! -- END for GEOS-5   - save for later use
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

   !------------------------------------------------------------------------------------
   subroutine get_condensation(q_old,t_old,po_cup,q_new,t_new)

      !-- calculate condensation and adjust t and q accordingly
      implicit none
      real    ,intent (in   ) :: po_cup,q_old,t_old ! before condensation
      real    ,intent (inout) ::        q_new,t_new ! after  condensation

      !---locals
      real ::  zqp, zcond, zcond1, zcor, zqsat,zi,zl,zf
      real :: psp, pt , pq
      real :: z3es,   z4es, z5alcp, zaldcp
      real :: ptare, cond
      real :: foeewmcu,foealfcu,foedemcu,foeldcpmcu

      ! real, parameter :: &
      !     RD=287.06                             &
      !    ,RV=461.52                             &
      !    ,RCPD=1004.71                          &
      !    ,RTT=273.16                            &
      !    ,RHOH2O=1000.                          &
      !    ,RLVTT=2.5008e+6                       &
      !    ,RLSTT=2.8345e+6                       &
      !    ,RETV=RV/RD-1.0                        &
      !    ,RLMLT=RLSTT-RLVTT                     &
      !    ,RCPV=4.*RV                            &
      !    ,R2ES=611.21*RD/RV                     &
      !    ,R3LES=17.502                          &
      !    ,R3IES=22.587                          &
      !    ,R4LES=32.19                           &
      !    ,R4IES=-0.7                            &
      !    ,R5LES=R3LES*(RTT-R4LES)               &
      !    ,R5IES=R3IES*(RTT-R4IES)               &
      !    ,R5ALVCP=R5LES*RLVTT/RCPD              &
      !    ,R5ALSCP=R5IES*RLSTT/RCPD              &
      !    ,RALVDCP=RLVTT/RCPD                    &
      !    ,RALSDCP=RLSTT/RCPD                    &
      !    ,RALFDCP=RLMLT/RCPD                    &
      !    ,RTWAT=RTT                             &
      !    ,RTBER=RTT-5.                          &
      !    ,RTBERCU=RTT-5.0                       &
      !    ,RTICE=RTT-23.                         &
      !    ,RTICECU=RTT-23.                       &
      !    ,RTWAT_RTICE_R=1./(RTWAT-RTICE)        &
      !    ,RTWAT_RTICECU_R=1./(RTWAT-RTICECU)    &
      !    ,RVTMP2=RCPV/RCPD-1.                   &
      !    ,ZQMAX=0.5

      !----------------------
      !-- for testing
      !----------------------
      !                          PSP (hPA)            TEMP(K)              Q(kg/kg)                ZCOND1(kg/kg)
      ! input          1   98020.0000000000        295.163188640152     1.745679200989956E-002
      ! output         1   98020.0000000000        295.513490789916     1.731605637618801E-002 -9.779916453243843E-007
      !----------------------
      ! input       157   85090.0000000000        288.089188935407     1.399404805052166E-002
      ! output      157   85090.0000000000        289.294751760460     1.350970717999820E-002 -1.146268822756454E-005
      !----------------------
      ! PT  = 288.089188935407
      ! PQ  = 1.399404805052166E-002
      ! PSP = 85090
      !----------------------

      !-- initial values
      PT  = t_old       ! K
      PQ  = q_old       ! kg/kg
      PSP = po_cup*100. ! hPa


      !--- get condensation in moist ascent --------------------------
      PTARE = PT
      ZQP    =1.0/PSP

      ZL=1.0/(PT-c_R4LES)
      ZI=1.0/(PT-c_R4IES)

      FOEALFCU = MIN(1.0,((MAX(c_RTICECU,MIN(c_RTWAT,PTARE))-c_RTICECU)*c_RTWAT_RTICECU_R)**2)
      ZQSAT=c_R2ES *(FOEALFCU *EXP(c_R3LES*(PTARE-c_RTT)*ZL)+&
                (1.0-FOEALFCU)*EXP(c_R3IES*(PTARE-c_RTT)*ZI))

      ZQSAT=ZQSAT*ZQP
      ZQSAT=MIN(0.5,ZQSAT)
      ZCOR=1.0-c_RETV*ZQSAT

      ZF=FOEALFCU*c_R5ALVCP*ZL**2 + (1.0-FOEALFCU)*c_R5ALSCP*ZI**2
      ZCOND=(PQ*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*ZF)

      if(ZCOND > 0.0)then
         FOELDCPMCU= FOEALFCU*c_RALVDCP+(1.0-FOEALFCU)*c_RALSDCP
         PT=PT+FOELDCPMCU*ZCOND
         PTARE = PT
         PQ=PQ-ZCOND

         ZL=1.0/(PT-c_R4LES)
         ZI=1.0/(PT-c_R4IES)


         FOEALFCU = MIN(1.0,((MAX(c_RTICECU,MIN(c_RTWAT,PTARE))-c_RTICECU)*c_RTWAT_RTICECU_R)**2)
         ZQSAT=c_R2ES *(FOEALFCU* EXP(c_R3LES*(PT-c_RTT)*ZL)+&
                   (1.0-FOEALFCU)*EXP(c_R3IES*(PT-c_RTT)*ZI))

         ZQSAT=ZQSAT*ZQP
         ZQSAT=ZQSAT-0.5*(ABS(0.5-ZQSAT)-(0.5-ZQSAT))


         ZCOR=1.0-c_RETV*ZQSAT
         ZF=FOEALFCU*c_R5ALVCP*ZL**2 + (1.0-FOEALFCU)*c_R5ALSCP*ZI**2

         ZCOND1=(PQ*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*ZF)
         if(ZCOND ==  0.0)ZCOND1=0.0
         FOELDCPMCU= FOEALFCU*c_RALVDCP+(1.0-FOEALFCU)*c_RALSDCP
         PT=PT+FOELDCPMCU*ZCOND1
         PQ=PQ-ZCOND1
      endif

      !-- FINAL --------------------------
      q_new =  PQ
      t_new =  PT
      cond  = -ZCOND1 != q_old-qnew, source for the liquid water
   end subroutine get_condensation


!********************
!********************
!********************
!********************
!******************** OBSOLETOS *********************************




   !--- get cloud fraction
   !
   ! do i=its,itf
   !    clfrac(i,:)=0.
   !    if(ierr(i) /= 0) cycle
   !    dummy1(kts:ktf) = xmb(i)* zuo(i,kts:ktf)
   !    dummy2(kts:ktf) = 100.*po_cup(i,kts:ktf)
   !    call get_cloud_fraction(ktf,kts,ktf                                                   &
   !     ,dummy2(kts:ktf),zo_cup(i,kts:ktf),tn_cup(i,kts:ktf),qo_cup(i,kts:ktf) &
   !     ,qco (i,kts:ktf),  qrco(i,kts:ktf),  dummy1(kts:ktf),clfrac(i,kts:ktf) )
   ! enddo
   !--------------------------------------------------------------------------------------------!
 
   subroutine get_cloud_fraction(  mzp, kts, ktf, &
      PPABS, PZZ, PT, PRV, QCO, QRCO, PMFLX, PCLDFR ,PRC, PRI )
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
      real, parameter :: C_ALVLTT  = 2.5008e6    ! latent heat of vaporisation at XTT
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
      !*       0.1   Declarations of dummy arguments :
      !
      !
      integer,              intent(in)   :: mzp     ! vertical dimension
      integer,              intent(in)   :: kts     ! vertical  computations start at
      !                                             ! KTS that is at least 1
      integer,              intent(in)   :: ktf     ! vertical computations can be
                                                    ! limited to MZP + 1 - KTF
                                                    ! default=1
      real, dimension(mzp), intent(in)    :: PPABS  ! pressure (Pa)
      real, dimension(mzp), intent(in)    :: PZZ    ! height of model levels (m)
      real, dimension(mzp), intent(in)    :: PT     ! grid scale T  (K)
      real, dimension(mzp), intent(in)    :: PRV    ! grid scale water vapor mixing ratio (kg/kg)
      real, dimension(mzp), intent(in)    :: PMFLX  ! convective mass flux (kg/(s m^2))
      real, dimension(mzp), intent(in)    :: QRCO   ! sub-grid scale liq water mixing ratio (kg/kg)
      real, dimension(mzp), intent(in)    :: QCO    ! in-cloud water mixing ratio (kg/kg)
      real, dimension(mzp), intent(inout),optional :: PRC    ! grid scale r_c mixing ratio (kg/kg)
      real, dimension(mzp), intent(inout),optional :: PRI    ! grid scale r_i (kg/kg)
      real, dimension(mzp), intent(out)   :: PCLDFR ! fractional cloudiness (between 0 and 1)
      !
      !
      !*       0.2   Declarations of local variables :
      !
      integer  ::  JKT, JKP, JKM,K     ! loop index
      real, dimension(mzp) :: ZTLK, ZRT       ! work arrays for T_l, r_t
      real, dimension(mzp) :: ZL              ! length-scale
      integer   :: ITPL    ! top levels of tropopause/highest inversion
      real      :: ZTMIN   ! min Temp. related to ITPL
      real, dimension(mzp) :: LOC_PRC,LOC_PRI
      !
      real :: ZTEMP, ZLV, ZLS, ZTL, ZPV, ZQSL, ZPIV, ZQSI, ZFRAC, ZCOND, ZCPD ! thermodynamics
      real :: ZLL, DZZ, ZZZ ! length scales
      real :: ZAH, ZA, ZB, ZSBAR, ZQ1, ZSIGMA, ZDRW, ZDTL ! related to computation of Sig_s
      real :: ZSIG_CONV,  ZSIGMA_NOCONV,  ZQ1_NOCONV      ! convective part of Sig_s
      !
      !*       0.3  Definition of constants :
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
      logical :: ONLY_CONVECTIVE_CLOUD_FRACTION=.true. ! set .false. for the total cloud fraction
      !-------------------------------------------------------------------------------
      !RETURN
      !
      if(PRESENT(PRC)) then
         LOC_PRC(:)=PRC(:)
      else
         LOC_PRC(:)=0.0
      endif
      if(PRESENT(PRI)) then
         LOC_PRI(:)=PRI(:)
      else
         LOC_PRI(:)=0.0
      endif

      PCLDFR(:) = 0. ! Initialize values
      !

      JKT = MZP+1-KTS
      !-will limit the model vertical column to 60 hPa
      do K=KTF,KTS,-1
         if(PPABS(k) > 60.*100.) then
            JKT = k
            !PRINT*,"JKT=",K,MZP+1-KTS ;CALL FLUSH(6)
            exit
         endif
      enddo

      do K=KTS,JKT
         ZTEMP  = PT(k)
          !latent heat of vaporisation/sublimation
         ZLV    = C_ALVLTT + ( XCPV - XCL ) * ( ZTEMP - XTT )
         ZLS    = XLSTT + ( XCPV - XCI ) * ( ZTEMP - XTT )

         !store temperature at saturation and total water mixing ratio
         ZRT(k)   = PRV(k) + LOC_PRC(k) + LOC_PRI(k)
         ZCPD     = XCPD  + XCPV*PRV(k) + XCL*LOC_PRC(k) + XCI*LOC_PRI(k)
         ZTLK(k)  = ZTEMP - ZLV*LOC_PRC(k)/ZCPD - ZLS*LOC_PRI(k)/ZCPD
      end do

      !-------------------------------------------------------------------------------
      ! Determine tropopause/inversion  height from minimum temperature

      ITPL  = KTS+1
      ZTMIN = 400.
      do k = KTS+1,JKT-1
         if ( PT(k) < ZTMIN ) then
            ZTMIN = PT(k)
            ITPL  = K
         endif
      end do

      ! Set the mixing length scale - used for computing the "turbulent part" of Sigma_s

      ZL(:) = 20.
      do k = KTS+1,JKT

         ! free troposphere
         ZL(k) = ZL0
         JKP   = ITPL
         ZZZ   = PZZ(k) -  PZZ(KTS)
            ! approximate length for boundary-layer : linear increase
         if ( ZL0 > ZZZ )  ZL(k) = ZZZ
            ! gradual decrease of length-scale near and above tropopause/top inversion
         if ( ZZZ > 0.9*(PZZ(JKP)-PZZ(KTS)) ) &
            ZL(k) = .6 * ZL(K-1)
      end do
      !-------------------------------------------------------------------------------

      do k=KTS+1,JKT-1
         JKP=k+1
         JKM=k-1
         ZTEMP  = PT(k)

         !latent heat of vaporisation/sublimation
         ZLV    = C_ALVLTT + ( XCPV - XCL ) * ( ZTEMP - XTT )
         ZLS    = XLSTT + ( XCPV - XCI ) * ( ZTEMP - XTT )

         ZCPD   = XCPD + XCPV*PRV(k) + XCL*LOC_PRC(k) + XCI*LOC_PRI(k)
         !temperature at saturation
         ZTL    = ZTEMP - ZLV*LOC_PRC(k)/ZCPD - ZLS*LOC_PRI(k)/ZCPD

         !saturated water vapor mixing ratio over liquid water
         ZPV    = MIN(EXP( XALPW - XBETAW / ZTL - XGAMW * LOG( ZTL ) ),0.99*PPABS(k))
         ZQSL   = XRD / XRV * ZPV / ( PPABS(k) - ZPV )

          !saturated water vapor mixing ratio over ice
         ZPIV   = MIN(EXP( XALPI - XBETAI / ZTL - XGAMI * LOG( ZTL ) ),0.99*PPABS(k))
         ZQSI   = XRD / XRV * ZPIV / ( PPABS(k) - ZPIV )

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

         ZA   = 1. / ( 1. + ZLV/ZCPD * ZAH )
         ZB   = ZAH * ZA

         !- parameterize Sigma_s with first_order closure
         DZZ    =  PZZ (JKP)  - PZZ(JKM)
         ZDRW   =  ZRT (JKP)  - ZRT(JKM)
         ZDTL   =  ZTLK(JKP) - ZTLK(JKM) + XG/ZCPD * DZZ
         ZLL    =  ZL(k)

         !- standard deviation due to convection
         ZSIG_CONV = ZCSIG_CONV * PMFLX(k) / ZA

         !- turb + conv
         ZSIGMA = SQRT( MAX( 1.e-25, ZCSIGMA*ZCSIGMA* ZLL*ZLL/(DZZ*DZZ) * ( &
            ZA*ZA*ZDRW*ZDRW - 2.*ZA*ZB*ZDRW*ZDTL   &
            + ZB*ZB*ZDTL*ZDTL                      ) &
            + ZSIG_CONV * ZSIG_CONV ) )

         !- zsigma should be of order 4.e-4 in lowest 5 km of atmosphere
         ZSIGMA = MAX( ZSIGMA, 1.e-10 )

         !- normalized saturation deficit
         ZSBAR = ZA * ( ZRT (k) - ZQSL )
         !- "Q1" parameter
         ZQ1   = ZSBAR / ZSIGMA

         !- total cloud fraction
         PCLDFR(k) = MAX( 0., MIN(1.,0.5+0.36*ATAN(1.55*ZQ1)) )

         if(ONLY_CONVECTIVE_CLOUD_FRACTION) then
            !- get cloud fraction associated with ONLY the sub-grid scale convective part
            !- this sigma does not include the sub-grid scale convective part
            ZSIGMA_NOCONV = SQRT( MAX( 1.e-25, ZCSIGMA*ZCSIGMA* ZLL*ZLL/(DZZ*DZZ) * ( &
               ZA*ZA*ZDRW*ZDRW - 2.*ZA*ZB*ZDRW*ZDTL   &
               + ZB*ZB*ZDTL*ZDTL  )))
            !- zsigma should be of order 4.e-4 in lowest 5 km of atmosphere
            ZSIGMA_NOCONV = MAX( ZSIGMA_NOCONV, 1.e-10 )
            ZQ1_NOCONV = ZSBAR / ZSIGMA_NOCONV

            !- cloud fraction associated with ONLY convective part ("total-turb")
            PCLDFR(k) = 0.36*(ATAN(1.55*ZQ1)-ATAN(1.55*ZQ1_NOCONV))

            PCLDFR(k) = MAX( 0., MIN(1.,PCLDFR(k)) )

         endif
         !- newer formulation, see GMD 2015
         !PCLDFR(k) = MAX( 0., MIN(1.,0.5+0.34*ATAN(1.85*ZQ1+2.33)) )
         !- this is area fraction of cloud cores
         !PCLDFR(k) = MAX( 0., MIN(1.,0.292/ZQ1**2) )

         cycle
         !total condensate diagnostic (not being used)
         if (ZQ1 > 0. .and. ZQ1 <= 2. ) then
            !orig   ZCOND =     EXP(-1.)+.66*ZQ1+.086*ZQ1*ZQ1
            ZCOND = MIN(EXP(-1.)+.66*ZQ1+.086*ZQ1**2, 2.) ! We use the MIN function for continuity
         else if (ZQ1 > 2.) then
            ZCOND = ZQ1
         else
            ZCOND = EXP( 1.2*ZQ1-1. )
         end if
         ZCOND = ZCOND * ZSIGMA

         if ( zcond < 1.e-12) then
            zcond = 0.
            pcldfr(k) = 0.
         end if
         if ( pcldfr(k) == 0.) then
            zcond = 0.
         end if

         LOC_PRC(k) = ZFRAC * ZCOND ! liquid condensate
         if (LUSERI) then
            LOC_PRI(k) = (1.-ZFRAC) * ZCOND   ! solid condensate
         end if

      !---
      ! compute s'rl'/Sigs^2
      ! used in w'rl'= w's' * s'rl'/Sigs^2
      !  PSIGRC(k) = PCLDFR(k)   ! Gaussian relation
      !
      ! s r_c/ sig_s^2
      !    PSIGRC(JI,JJ,JK) = PCLDFR(JI,JJ,JK)  ! use simple Gaussian relation
      !
      !    multiply PSRCS by the lambda3 coefficient
      !
      !      PSIGRC(JI,JJ,JK) = 2.*PCLDFR(JI,JJ,JK) * MIN( 3. , MAX(1.,1.-ZQ1) )
      ! in the 3D case lambda_3 = 1.
      !      INQ1 = MIN( MAX(-22,FLOOR(2*ZQ1) ), 10)
      !      ZINC = 2.*ZQ1 - INQ1
      !
      !      PSIGRC(k) =  MIN(1.,(1.-ZINC)*ZSRC_1D(INQ1)+ZINC*ZSRC_1D(INQ1+1))
      !
      !      PSIGRC(k) = PSIGRC(k)* MIN( 3. , MAX(1.,1.-ZQ1) )
      !---
      end do
     !
   end subroutine get_cloud_fraction

 !------------------------------------------------------------------------------------
   function auto_rk(n,step,aux,xexp,qrc1) result(PW)
      integer, intent(in) :: n
      real   , intent(in) :: step,aux,qrc1,xexp
      real                :: PW

      PW=step*qrc1*(1.0-exp(-aux**xexp))/float(n)

   end function auto_rk
!+---+-----------------------------------------------------------------+
   !DSM {
   pure function intfuncgamma(x, y) result(z)
      real :: z
      real, intent(in) :: x, y

      z = x**(y-1.0) * exp(-x)
   end function intfuncgamma

   function gammaBrams(a) result(g)
      real :: g
      real, intent(in) :: a

      real, parameter :: small = 1.0e-4
      integer, parameter :: points = 100000

      real :: infty, dx, p, sp(2, points), x
      integer :: i
      logical :: correction

      x = a

      correction = .false.
      ! value with x<1 gives \infty, so we use
      ! \Gamma(x+1) = x\Gamma(x)
      ! to avoid the problem
      if ( x < 1.0 ) then
         correction = .true.
         x = x + 1
      end if

      ! find a "reasonable" infinity...
      ! we compute this integral indeed
      ! \int_0^M dt t^{x-1} e^{-t}
      ! where M is such that M^{x-1} e^{-M} ≤ \epsilon
      infty = 1.0e4
      do while ( intfuncgamma(infty, x) > small )
         infty = infty * 10.0
      end do

      ! using simpson
      dx = infty/real(points)
      sp = 0.0
      forall(i=1:points/2-1) sp(1, 2*i) = intfuncgamma(2.0*(i)*dx, x)
      forall(i=1:points/2) sp(2, 2*i - 1) = intfuncgamma((2.0*(i)-1.0)*dx, x)
      g = (intfuncgamma(0.0, x) + 2.0*sum(sp(1,:)) + 4.0*sum(sp(2,:)) + &
         intfuncgamma(infty, x))*dx/3.0

      if ( correction ) g = g/a

   end function gammaBrams
   !DSM}   !------------------------------------------------------------------------------------
   real function coldPoolStart_orig(CNV_TR)
      implicit none
      real,intent(in)  :: CNV_TR
      real             :: f1
      real,parameter   :: width = 100. !orig 100.
      real    :: mx_buoy1       = (c_cp*5.0 + c_alvl*2.e-3)*0.025  !=   250.5 J/kg
      real    :: mx_buoy2       = (c_cp*10. + c_alvl*4.e-3)        != 20004.0 J/kg: temp exc=10 K, q deficit=4 g/kg (=> mx_buoy ~ 20 kJ/kg)

      f1= min (mx_buoy2,CNV_TR)
      !--- f1 > mx_buoy1 => coldPoolStart ---> 1 
      coldPoolStart_orig =  (1.35+atan( (f1-mx_buoy1)/width))/2.8
      coldPoolStart_orig =  max(0.00,min(coldPoolStart_orig,1.00))
   end function
   !------------------------------------------------------------------------------------
   subroutine fix_scale_dep(cumulus,kts,kte,ktf,zo_cup,w,rho,v_ratio)
      implicit none
      character *(*),intent (in) :: cumulus
      integer       ,intent (in) :: kts,kte,ktf
      real          ,intent (in) :: zo_cup(kts:kte),w(kts:kte),rho(kts:kte)
      
      real          ,intent (out) :: v_ratio
      !-- local vars
      integer :: k
      real :: vert_int_w,total_dz,dz
      real :: vert_int_w_threshold  = 1.0  != grid-scale vertical velocity threshold (m/s)
      real :: vert_height_threshold = 3000.!= max height to calculate the mean grid-scale vertical velocity (m)

      total_dz    = 0.
      vert_int_w  = 0.

      loopK: do k = kts,ktf 
      
         dz         = zo_cup(k+1) - zo_cup(k)   
         total_dz   = total_dz    + dz
         vert_int_w = vert_int_w  + dz * w(k) ! vert veloc m/s
         if( zo_cup(k+1) >= vert_height_threshold) exit loopK 
         
      end do loopK
      vert_int_w =  vert_int_w / (1.e-12+total_dz)
      
      v_ratio    = (vert_int_w / vert_int_w_threshold)**4.
 
   end subroutine fix_scale_dep

!------------------------------------------------------------------------------------
   subroutine get_zu_zd_pdf_orig(draft,ierr,kb,kt,zs,zuf,ztop,zu,kts,kte,ktf)

      implicit none
      integer, intent(in) ::kb,kt,kts,kte,ktf
      real, intent(in) :: Zs,Zuf,Ztop
      real, intent(inout) :: zu(kts:kte)
      integer, intent(inout) :: ierr
      character*(*), intent(in) ::draft

      !- local var
      integer :: add,i,nrec=0,k,kb_adj
      real ::zumax,ztop_adj
      real ::beta, alpha,kratio,tunning

      !- kb cannot be at 1st level
      kb_adj=max(kb,2)

      !-- fill zu with zeros
      zu=0.0

      if(draft == "UP" .or. draft == "up" ) then
         if(kt<=kb_adj) then
            !stop "ktop must be larger than kbcon"
            ierr=99
            return
         endif
         !beta=4.  !=> must larger than 1
                   !=> higher makes the profile sharper
                   !=> around the maximum zu
         add=0     !=> additional levels above kbcon, where
                   !=> the maximum zu will resides
         kb_adj=kb_adj+add
         kb_adj=max(10,kb_adj)

         !- this alpha constrains the location of the maximun ZU to be at
         !- "kb_adj" vertical level
         !alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))

         !- 2nd approach for beta and alpha parameters
         !- the tunning parameter must be between 0.5 (low  level max zu)
         !-                                   and 1.5 (high level max zu)
         tunning = 0.6
         beta    = 2.0/tunning
         alpha   = tunning*beta

          !- Beta PDF
         do k=kts,min(kte,kt+1)
            kratio= float(k)/float(kt+1)

            zu(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
         enddo

      elseif(draft == "DOWN" .or. draft == "DOWNM") then
         add=0    !=> additional levels above kbcon, where
                  !=> the maximum zu will resides
         beta=4.  !=> must larger than 1
                  !=> higher makes the profile sharper
                  !=> around the maximum zu
         alpha= 0.25*beta

         !- for downdrafts kb = jmin(i)-levadj
         kb_adj=kb_adj+add

         !- 2nd approach for beta and alpha parameters
         !- the tunning parameter must be between 0.5 (low  level max zu)
         !-                                   and 1.5 (high level max zu)
         tunning = 1.
         beta    = 2.0/tunning
         alpha   = tunning*beta

         !- Beta PDF
         do k=kts,min(kte,kt)
            kratio= float(k)/float(kt)

            zu(k+1) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
         enddo

      elseif(draft == "shallow" .or. draft == "SHALLOW") then

         alpha= 3.
         beta = 2.*alpha
         kb_adj=1 ! level where mass flux starts

         !- Beta PDF
         do k=kts+kb_adj-1,min(kte,kt+1)
            kratio=float(k+1-kb_adj)/float(kt+1)  !-kb_adj+1)

            zu(k)=kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
         enddo

      else
         print*, "unknown type of flow" ,draft
         stop "routine get_zu_zd"

      endif

      !- normalize ZU
      zu(kts:min(kte,kt+1))= zu(kts:min(kte,kt+1))/ maxval(zu(kts:min(kte,kt+1)))

      !--- Sanity checks
      if(beta <= 1) stop "beta must be larger than 1"

      if(minval(zu(:)) < 0.0 ) then
         print*," zu has negative values for ", draft
         stop   " zu < zero"
      endif
      if(maxval(zu(:)) > 1.0 ) then
         print*," zu has values greater than 1 for ", draft
         stop   " zu  >  one"
      endif

      return

   !OPEN(19,FILE= 'zu.gra', FORM='unformatted',ACCESS='direct'&
   !       ,STATUS='unknown',RECL=4)
   ! DO k = kts,kte
   !    nrec=nrec+1
   !     WRITE(19,REC=nrec) zu(k)
   ! END DO
   !close (19)

   end subroutine get_zu_zd_pdf_orig

  !------------------------------------------------------------------------------------

   real function satvap(temp2)
      implicit none
      real :: temp2, temp, toot, toto, eilog, tsot,  &
         ewlog, ewlog2, ewlog3, ewlog4
      temp = temp2-273.155
      if (temp.lt.-20.) then   !!!! ice saturation
         toot = 273.16 / temp2
         toto = 1 / toot
         eilog = -9.09718 * (toot - 1) - 3.56654 * (log(toot) / &
            log(10.)) + .876793 * (1 - toto) + (log(6.1071) / log(10.))
         satvap = 10 ** eilog
      else
         tsot = 373.16 / temp2
         ewlog = -7.90298 * (tsot - 1) + 5.02808 * &
            (log(tsot) / log(10.))
         ewlog2 = ewlog - 1.3816e-07 * &
            (10 ** (11.344 * (1 - (1 / tsot))) - 1)
         ewlog3 = ewlog2 + .0081328 * &
            (10 ** (-3.49149 * (tsot - 1)) - 1)
         ewlog4 = ewlog3 + (log(1013.246) / log(10.))
         satvap = 10 ** ewlog4
      endif

   end function

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  !---------------------------------------------------------------------------------------------------
end  module modConvParGF



