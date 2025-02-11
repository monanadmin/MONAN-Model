! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver_convection
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_timer, only : mpas_timer_start, mpas_timer_stop

 use mpas_atmphys_constants
 use mpas_atmphys_utilities
 use mpas_atmphys_vars

!wrf physics:
 use module_cu_gf
 use module_cu_kfeta
 use module_cu_tiedtke
 use module_cu_ntiedtke,only: cu_ntiedtke_driver
 use module_cu_gf_monan,only: cu_gf_monan_driver                  &
                            , readGFConvParNML                    &
                            , modConvParGF_initialized            &
                            , initModConvParGF                    &
                            , tau_ocea_cp, tau_land_cp            &
                            , mx_buoy1, mx_buoy2                  &
                            , satur_spec_hum
 use mpas_atm_dimensions, only : nVertLevels
 use mpas_dmpar, only : IO_NODE
 implicit none
 !--srf
 ! Abstract interface for routine used to communicate halos of fields
 ! in a named group
 !
 abstract interface
      subroutine halo_exchange_routine(domain, halo_group, ierr)

         use mpas_derived_types, only : domain_type

         type (domain_type), intent(inout) :: domain
         character(len=*), intent(in) :: halo_group
         integer, intent(out), optional :: ierr

      end subroutine halo_exchange_routine
 end interface
 !--srf
 private
 public:: allocate_convection,     &
          deallocate_convection,   &
          init_convection,         &
          driver_convection,       &
          update_convection_step1, &
          update_convection_step2


!MPAS driver for parameterization of convection.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in mpas_atmphys_driver_convection:
! ----------------------------------------------
! allocate_convection     : allocate local arrays for parameterization of convection.
! deallocate_convection   : deallocate local arrays for parameterization of convection.
! init_convection         : initialization of individual convection scheme.
! driver_convection       : main driver (called from subroutine physics_driver).
! convection_from_MPAS    : initialize local arrays.
! convection_to_MPAS      : copy local arrays to MPAS arrays.
! update_convection_step1 : updates lifetime of deep convective clouds in Kain-Fritsch scheme.
! update_convection_step2 : updates accumulated precipitation output from convection schemes.
!
! WRF physics called from driver_convection:
! ------------------------------------------
! * module_cu_kfeta   : Kain-Fritsch convection scheme.
! * module_cu_tiedtke : Tiedtke convection scheme.

! add-ons and modifications to sourcecode:
! ----------------------------------------
! * removed the pre-processor option "do_hydrostatic_pressure" before call to the subroutines kf_eta_cps
!   and tiedtke. 
! * removed call to the updated Kain-Fritsch convection scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2013-05-29.
! * added the mean distance between cell centers in the call to the Kain-Fritsch parameterization of convection.
!   Laura D. Fowler (laura@ucar.edu) / 2013-08-22. 
! * in call to subroutine kf_eta_cps, replaced the variable g (that originally pointed to gravity) with gravity,
!   for simplicity.
!   Laura D. Fowler (laura@ucar.edu) / 2014-03-21.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * Modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * in ./physics_wrf, updated the Tiedtke convection scheme to that in WRF version 3.6.1. The call to subroutine 
!   cu_tiedtke has been updated accordingly to include the sensible heat flux.
!   Laura D. Fowler (laura@ucar.edu) / 2014-09-11.
! * renamed "kain_fritsch" with "cu_kain_fritsch" and "tiedtke" with "cu_tiedtke".
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-22.
! * added the implementation of the Grell-Freitas convection scheme (option cu_grell_freitas).
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-30.
!   code cleaner, now that we are using physics packages.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-01.
! * in the call to subroutine cu_grell_freitas, added the array pratec for the calculation of rainc.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-14.
! * in the call to subroutine cu_grell_freitas, removed the argument ktop_deep.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-20.
! * updated the call to the Tiedtke parameterization of convection to WRF version 3.8.1. This option is
!   triggered using the option cu_tiedtke.
!   Laura D. Fowler (laura@ucar.edu) / 2016-08-18.
! * added the call to the "new" Tiedtke parameterization of convection from WRF version 3.8.1. This option is
!   triggered using the option cu_ntiedtke.
!   Laura D. Fowler (laura@ucar.edu) / 2016-09-20.
! * for the kain_fritsch parameterization of convection, change the definition of dx_p to match that used in the
!   Grell-Freitas and "new Tiedtke" parameterization.
!   Laura D. Fowler (laura@ucar.edu) / 2016-10-18.
! * since we removed the local variable convection_scheme from mpas_atmphys_vars.F, now defines convection_scheme
!   as a pointer to config_convection_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2917-02-16.
! * removed f_qv,f_qr,and f_qs in the calls to cu_tiedtke and cu_ntiedtke. removed f_qv and f_qc in the call to
!   kf_eta_cps.
!   Laura D. Fowler (laura@ucar.edu) / 2024-02-13.
!---------  MONAN implementations -------------------------------------------------
! * implemented the Grell-Freitas convection parameterization (MONAN version). 
!   This option is triggered using the option 'cu_gf_monan'.
!   Saulo Freitas INPE/Brazil (contact: saulo.freitas@inpe.br) / 2024-04-10


 contains


!=================================================================================================================
 subroutine allocate_convection(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local variables and pointers:
 character(len=StrKIND),pointer:: convection_scheme
 integer,pointer:: config_gf_sub3d

 integer:: i,k,j

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)
 !print*,"---allocating for---------------------------------------------",trim(convection_scheme)

 if(.not.allocated(cu_act_flag)) allocate(cu_act_flag(ims:ime,jms:jme)       )
 if(.not.allocated(rthcuten_p) ) allocate(rthcuten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqvcuten_p) ) allocate(rqvcuten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqccuten_p) ) allocate(rqccuten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(rqicuten_p) ) allocate(rqicuten_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(pratec_p)   ) allocate(pratec_p(ims:ime,jms:jme)          )
 if(.not.allocated(raincv_p)   ) allocate(raincv_p(ims:ime,jms:jme)          )

 do i = its,ite
 do j = jts,jte
    pratec_p(i,j) = 0._RKIND
    raincv_p(i,j) = 0._RKIND
 enddo
 enddo

 do i = its,ite
 do k = kts,kte
 do j = jts,jte
    rthcuten_p(i,k,j) = 0._RKIND
    rqvcuten_p(i,k,j) = 0._RKIND
    rqccuten_p(i,k,j) = 0._RKIND
    rqicuten_p(i,k,j) = 0._RKIND
 enddo
 enddo
 enddo
 !print*,'in alloc',maxval(rthcuten_p(its:ite,kts:kte,jts:jte))*86400._RKIND

 convection_select: select case(convection_scheme)

   case ("cu_gf_monan")
       if(.not.allocated(ht_p)           ) allocate(ht_p(ims:ime,jms:jme)               )
       if(.not.allocated(kpbl_p)         ) allocate(kpbl_p(ims:ime,jms:jme)             )
 
       if(.not.allocated(dx_p)           ) allocate(dx_p(ims:ime,jms:jme)               )
       if(.not.allocated(area_p)         ) allocate(area_p(ims:ime,jms:jme)             )
       if(.not.allocated(hfx_p)          ) allocate(hfx_p(ims:ime,jms:jme)              )
       if(.not.allocated(qfx_p)          ) allocate(qfx_p(ims:ime,jms:jme)              )
       if(.not.allocated(xland_p)        ) allocate(xland_p(ims:ime,jms:jme)            )
       if(.not.allocated(t2m_p)          ) allocate(t2m_p(ims:ime,jms:jme)              )
       if(.not.allocated(cape_p)         ) allocate(cape_p(ims:ime,jms:jme)             )
       if(.not.allocated(cin_p)          ) allocate(cin_p(ims:ime,jms:jme)              )
       if(.not.allocated(wlpool_p)       ) allocate(wlpool_p(ims:ime,jms:jme)           )
       if(.not.allocated(xlat_p)         ) allocate(xlat_p(ims:ime,jms:jme)             )
       if(.not.allocated(xlon_p)         ) allocate(xlon_p(ims:ime,jms:jme)             )
       if(.not.allocated(rmfxdpcu_p)     ) allocate(rmfxdpcu_p(ims:ime,jms:jme)         )
       if(.not.allocated(rmfxdncu_p)     ) allocate(rmfxdncu_p(ims:ime,jms:jme)         )
       if(.not.allocated(rmfxmdcu_p)     ) allocate(rmfxmdcu_p(ims:ime,jms:jme)         )
       if(.not.allocated(rmfxshcu_p)     ) allocate(rmfxshcu_p(ims:ime,jms:jme)         )
       if(.not.allocated(rtopdpcu_p)     ) allocate(rtopdpcu_p(ims:ime,jms:jme)         )
       if(.not.allocated(rtopmdcu_p)     ) allocate(rtopmdcu_p(ims:ime,jms:jme)         )
       if(.not.allocated(rtopshcu_p)     ) allocate(rtopshcu_p(ims:ime,jms:jme)         )
       if(.not.allocated(rbotdpcu_p)     ) allocate(rbotdpcu_p(ims:ime,jms:jme)         )
       if(.not.allocated(v2d1_p)         ) allocate(v2d1_p(ims:ime,jms:jme)             )
       if(.not.allocated(v2d2_p)         ) allocate(v2d2_p(ims:ime,jms:jme)             )
       if(.not.allocated(lightn_dens_p)  ) allocate(lightn_dens_p(ims:ime,jms:jme)      )
       if(.not.allocated(sigma_deep_p)   ) allocate(sigma_deep_p(ims:ime,jms:jme)       )

       if(.not.allocated(rthblten_p)     ) allocate(rthblten_p(ims:ime,kms:kme,jms:jme) )
       if(.not.allocated(rthdynten_p)    ) allocate(rthdynten_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(rthratensw_p)   ) allocate(rthratensw_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(rthratenlw_p)   ) allocate(rthratenlw_p(ims:ime,kms:kme,jms:jme))

       if(.not.allocated(rqvblten_p)     ) allocate(rqvblten_p(ims:ime,kms:kme,jms:jme) )
       if(.not.allocated(rqvdynten_p)    ) allocate(rqvdynten_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(rucuten_p)      ) allocate(rucuten_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(rvcuten_p)      ) allocate(rvcuten_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(rcnvcfcuten_p)  ) allocate(rcnvcfcuten_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(rbuoyxcuten_p)  ) allocate(rbuoyxcuten_p(ims:ime,kms:kme,jms:jme))

       if(.not.allocated(tkepbl_p)       ) allocate(tkepbl_p  (ims:ime,kms:kme,jms:jme) )
       if(.not.allocated(elpbl_p)        ) allocate(elpbl_p   (ims:ime,kms:kme,jms:jme) )
       if(.not.allocated(rupmfxcu_p)     ) allocate(rupmfxcu_p(ims:ime,kms:kme,jms:jme) )
       if(.not.allocated(rdnmfxcu_p)     ) allocate(rdnmfxcu_p(ims:ime,kms:kme,jms:jme) )
       if(.not.allocated(v3d1_p)         ) allocate(v3d1_p    (ims:ime,kms:kme,jms:jme) )
 

       do i = its,ite
       do k = kts,kte
       do j = jts,jte
          rucuten_p(i,k,j)     = 0._RKIND
          rvcuten_p(i,k,j)     = 0._RKIND
          rcnvcfcuten_p(i,k,j) = 0._RKIND
          rbuoyxcuten_p(i,k,j) = 0._RKIND
          tkepbl_p(i,k,j)      = 0._RKIND
          elpbl_p(i,k,j)       = 0._RKIND
          rupmfxcu_p(i,k,j)    = 0._RKIND
          rdnmfxcu_p(i,k,j)    = 0._RKIND
       enddo
       enddo
       enddo

       do i = its,ite
       do j = jts,jte
          rmfxdpcu_p(i,j)   = 0._RKIND
          rmfxdncu_p(i,j)   = 0._RKIND
          rmfxmdcu_p(i,j)   = 0._RKIND
          rmfxshcu_p(i,j)   = 0._RKIND
          rtopdpcu_p(i,j)   = 0._RKIND
          rtopmdcu_p(i,j)   = 0._RKIND
          rtopshcu_p(i,j)   = 0._RKIND
          rbotdpcu_p(i,j)   = 0._RKIND
          v2d2_p(i,j)       = 0._RKIND
          v2d1_p(i,j)       = 0._RKIND
          sigma_deep_p(i,j) = 0._RKIND
       enddo
       enddo
       do i = its,ite
       do k = kts,kte
       do j = jts,jte
          v3d1_p(i,k,j) = 0._RKIND
       enddo
       enddo
       enddo

       call mpas_pool_get_config(configs,'config_gf_sub3d',config_gf_sub3d)
       if(config_gf_sub3d > 0) then
            if(.not.allocated(sub3d_rthcuten_p)) allocate(sub3d_rthcuten_p(ims:ime,kms:kme,jms:jme)  )
            if(.not.allocated(sub3d_rqvcuten_p)) allocate(sub3d_rqvcuten_p(ims:ime,kms:kme,jms:jme)  )
            if(.not.allocated(sub3d_rucuten_p))  allocate(sub3d_rucuten_p (ims:ime,kms:kme,jms:jme)  )
            if(.not.allocated(sub3d_rvcuten_p))  allocate(sub3d_rvcuten_p (ims:ime,kms:kme,jms:jme)  )
            do i = its,ite
            do k = kts,kte
            do j = jts,jte
               sub3d_rthcuten_p(i,k,j)  = 0._RKIND
               sub3d_rqvcuten_p(i,k,j)  = 0._RKIND
               sub3d_rucuten_p(i,k,j)   = 0._RKIND
               sub3d_rvcuten_p(i,k,j)   = 0._RKIND
            enddo
            enddo
            enddo
       else
            if(.not.allocated(sub3d_rthcuten_p)) allocate(sub3d_rthcuten_p(0,0,0)  )
            if(.not.allocated(sub3d_rqvcuten_p)) allocate(sub3d_rqvcuten_p(0,0,0)  )
            if(.not.allocated(sub3d_rucuten_p))  allocate(sub3d_rucuten_p (0,0,0)  )
            if(.not.allocated(sub3d_rvcuten_p))  allocate(sub3d_rvcuten_p (0,0,0)  )
       endif

    case ("cu_grell_freitas")
       if(.not.allocated(ht_p)           ) allocate(ht_p(ims:ime,jms:jme)               )
       if(.not.allocated(cubot_p)        ) allocate(cubot_p(ims:ime,jms:jme)            )
       if(.not.allocated(cutop_p)        ) allocate(cutop_p(ims:ime,jms:jme)            )
       if(.not.allocated(kpbl_p)         ) allocate(kpbl_p(ims:ime,jms:jme)             )
       if(.not.allocated(k22_shallow_p)  ) allocate(k22_shallow_p(ims:ime,jms:jme)      )
       if(.not.allocated(kbcon_shallow_p)) allocate(kbcon_shallow_p(ims:ime,jms:jme)    )
       if(.not.allocated(ktop_shallow_p) ) allocate(ktop_shallow_p(ims:ime,jms:jme)     )
       if(.not.allocated(kbot_shallow_p) ) allocate(kbot_shallow_p(ims:ime,jms:jme)     )
       if(.not.allocated(ktop_deep_p)    ) allocate(ktop_deep_p(ims:ime,jms:jme)        )

       if(.not.allocated(dx_p)           ) allocate(dx_p(ims:ime,jms:jme)               )
       if(.not.allocated(area_p)         ) allocate(area_p(ims:ime,jms:jme)             )
       if(.not.allocated(gsw_p)          ) allocate(gsw_p(ims:ime,jms:jme)              )
       if(.not.allocated(hfx_p)          ) allocate(hfx_p(ims:ime,jms:jme)              )
       if(.not.allocated(qfx_p)          ) allocate(qfx_p(ims:ime,jms:jme)              )
       if(.not.allocated(xland_p)        ) allocate(xland_p(ims:ime,jms:jme)            )

       if(.not.allocated(rthblten_p)     ) allocate(rthblten_p(ims:ime,kms:kme,jms:jme) )
       if(.not.allocated(rthdynten_p)    ) allocate(rthdynten_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(rthraten_p)     ) allocate(rthraten_p(ims:ime,kms:kme,jms:jme) )
       if(.not.allocated(rqvblten_p)     ) allocate(rqvblten_p(ims:ime,kms:kme,jms:jme) )
       if(.not.allocated(rqvdynten_p)    ) allocate(rqvdynten_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(rucuten_p)      ) allocate(rucuten_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(rvcuten_p)      ) allocate(rvcuten_p(ims:ime,kms:kme,jms:jme)  )

       if(.not.allocated(xmb_total_p)    ) allocate(xmb_total_p(ims:ime,jms:jme)        )
       if(.not.allocated(xmb_shallow_p)  ) allocate(xmb_shallow_p(ims:ime,jms:jme)      )

       if(.not.allocated(qccu_p)         ) allocate(qccu_p(ims:ime,kms:kme,jms:jme)     )
       if(.not.allocated(qicu_p)         ) allocate(qicu_p(ims:ime,kms:kme,jms:jme)     )
       do i = its,ite
       do k = kts,kte
       do j = jts,jte
          rucuten_p(i,k,j) = 0._RKIND
          rvcuten_p(i,k,j) = 0._RKIND
          qccu_p(i,k,j)    = 0._RKIND
          qicu_p(i,k,j)    = 0._RKIND
       enddo
       enddo
       enddo

    case ("cu_kain_fritsch")
       if(.not.allocated(dx_p)       ) allocate(dx_p(ims:ime,jms:jme)              )
       if(.not.allocated(area_p)     ) allocate(area_p(ims:ime,jms:jme)            )
       if(.not.allocated(nca_p)      ) allocate(nca_p(ims:ime,jms:jme)             )
       if(.not.allocated(cubot_p)    ) allocate(cubot_p(ims:ime,jms:jme)           )
       if(.not.allocated(cutop_p)    ) allocate(cutop_p(ims:ime,jms:jme)           )
       if(.not.allocated(w0avg_p)    ) allocate(w0avg_p(ims:ime,kms:kme,jms:jme)   )
       if(.not.allocated(rqrcuten_p) ) allocate(rqrcuten_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(rqscuten_p) ) allocate(rqscuten_p(ims:ime,kms:kme,jms:jme))

       do i = its,ite
       do j = jts,jte
          cubot_p(i,j) = DBLE(kte+1)
          cutop_p(i,j) = DBLE(kts)
       enddo
       enddo

       do i = its,ite
       do k = kts,kte
       do j = jts,jte
          rqrcuten_p(i,k,j) = 0._RKIND
          rqscuten_p(i,k,j) = 0._RKIND
       enddo
       enddo
       enddo

    case ("cu_tiedtke","cu_ntiedtke")
       if(.not.allocated(hfx_p)        ) allocate(hfx_p(ims:ime,jms:jme)                )
       if(.not.allocated(qfx_p)        ) allocate(qfx_p(ims:ime,jms:jme)                )
       if(.not.allocated(xland_p)      ) allocate(xland_p(ims:ime,jms:jme)              )
       if(.not.allocated(rucuten_p)    ) allocate(rucuten_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(rvcuten_p)    ) allocate(rvcuten_p(ims:ime,kms:kme,jms:jme)    )

       do i = its,ite
       do j = jts,jte
          hfx_p(i,j)   = 0._RKIND
          qfx_p(i,j)   = 0._RKIND
          xland_p(i,j) = 0._RKIND
       enddo
       enddo

       do i = its,ite
       do k = kts,kte
       do j = jts,jte
          rucuten_p(i,k,j)     = 0._RKIND
          rvcuten_p(i,k,j)     = 0._RKIND
       enddo
       enddo
       enddo

       cu_tiedtke_select: select case(convection_scheme)

          case ("cu_tiedtke")
             if(.not.allocated(rqvdynten_p)  ) allocate(rqvdynten_p(ims:ime,kms:kme,jms:jme)  )
             if(.not.allocated(rqvdynblten_p)) allocate(rqvdynblten_p(ims:ime,kms:kme,jms:jme))

             do i = its,ite
             do k = kts,kte
             do j = jts,jte
                rqvdynten_p(i,k,j)   = 0._RKIND
                rqvdynblten_p(i,k,j) = 0._RKIND
             enddo
             enddo
             enddo

          case("cu_ntiedtke")
             if(.not.allocated(dx_p)     ) allocate(dx_p(ims:ime,jms:jme)             )
             if(.not.allocated(rqvften_p)) allocate(rqvften_p(ims:ime,kms:kme,jms:jme))
             if(.not.allocated(rthften_p)) allocate(rthften_p(ims:ime,kms:kme,jms:jme))

             do i = its,ite
             do k = kts,kte
             do j = jts,jte
                rqvften_p(i,k,j) = 0._RKIND
                rthften_p(i,k,j) = 0._RKIND
             enddo
             enddo
             enddo

           case default

       end select cu_tiedtke_select

    case default

 end select convection_select

 end subroutine allocate_convection

!=================================================================================================================
 subroutine deallocate_convection(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: convection_scheme

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 if(allocated(cu_act_flag)) deallocate(cu_act_flag)
 if(allocated(rthcuten_p) ) deallocate(rthcuten_p )
 if(allocated(rqvcuten_p) ) deallocate(rqvcuten_p )
 if(allocated(rqccuten_p) ) deallocate(rqccuten_p )
 if(allocated(rqicuten_p) ) deallocate(rqicuten_p )
 if(allocated(pratec_p)   ) deallocate(pratec_p   )
 if(allocated(raincv_p)   ) deallocate(raincv_p   )

 convection_select: select case(convection_scheme)

    case ("cu_grell_freitas")
       if(allocated(ht_p)           ) deallocate(ht_p           )
       if(allocated(cubot_p)        ) deallocate(cubot_p        )
       if(allocated(cutop_p)        ) deallocate(cutop_p        )
       if(allocated(kpbl_p)         ) deallocate(kpbl_p         )
       if(allocated(k22_shallow_p)  ) deallocate(k22_shallow_p  )
       if(allocated(kbcon_shallow_p)) deallocate(kbcon_shallow_p)
       if(allocated(ktop_shallow_p) ) deallocate(ktop_shallow_p )
       if(allocated(kbot_shallow_p) ) deallocate(kbot_shallow_p )
       if(allocated(ktop_deep_p)    ) deallocate(ktop_deep_p    )

       if(allocated(dx_p)           ) deallocate(dx_p           )
       if(allocated(area_p)         ) deallocate(area_p         )
       if(allocated(gsw_p)          ) deallocate(gsw_p          )
       if(allocated(hfx_p)          ) deallocate(hfx_p          )
       if(allocated(qfx_p)          ) deallocate(qfx_p          )
       if(allocated(xland_p)        ) deallocate(xland_p        )
       if(allocated(rthblten_p)     ) deallocate(rthblten_p     )
       if(allocated(rthdynten_p)    ) deallocate(rthdynten_p    )
       if(allocated(rthraten_p)     ) deallocate(rthraten_p     )
       if(allocated(rqvblten_p)     ) deallocate(rqvblten_p     )
       if(allocated(rqvdynten_p)    ) deallocate(rqvdynten_p    )
       if(allocated(rucuten_p)      ) deallocate(rucuten_p      )
       if(allocated(rvcuten_p)      ) deallocate(rvcuten_p      )

       if(allocated(xmb_total_p)    ) deallocate(xmb_total_p    )
       if(allocated(xmb_shallow_p)  ) deallocate(xmb_shallow_p  )

       if(allocated(qccu_p)         ) deallocate(qccu_p         )
       if(allocated(qicu_p)         ) deallocate(qicu_p         )

   case ("cu_gf_monan")
       if(allocated(ht_p)           ) deallocate(ht_p           )
       if(allocated(kpbl_p)         ) deallocate(kpbl_p         )
 
       if(allocated(dx_p)           ) deallocate(dx_p           )
       if(allocated(area_p)         ) deallocate(area_p         )
       if(allocated(hfx_p)          ) deallocate(hfx_p          )
       if(allocated(qfx_p)          ) deallocate(qfx_p          )
       if(allocated(xland_p)        ) deallocate(xland_p        )
       if(allocated(t2m_p)          ) deallocate(t2m_p          )
       if(allocated(cape_p)         ) deallocate(cape_p         )
       if(allocated(cin_p)          ) deallocate(cin_p          )
       if(allocated(wlpool_p)       ) deallocate(wlpool_p       )
       if(allocated(xlat_p)         ) deallocate(xlat_p         )
       if(allocated(xlon_p)         ) deallocate(xlon_p         )
       if(allocated(rmfxmdcu_p)     ) deallocate(rmfxmdcu_p     )
       if(allocated(rmfxshcu_p)     ) deallocate(rmfxshcu_p     )
       if(allocated(rtopdpcu_p)     ) deallocate(rtopdpcu_p     )
       if(allocated(rtopmdcu_p)     ) deallocate(rtopmdcu_p     )
       if(allocated(rtopshcu_p)     ) deallocate(rtopshcu_p     )
       if(allocated(rbotdpcu_p)     ) deallocate(rbotdpcu_p     )
       if(allocated(v2d1_p)         ) deallocate(v2d1_p         )
       if(allocated(v2d2_p)         ) deallocate(v2d2_p         )
       if(allocated(lightn_dens_p)  ) deallocate(lightn_dens_p  )
       if(allocated(sigma_deep_p)   ) deallocate(sigma_deep_p   )

       if(allocated(rthblten_p)     ) deallocate(rthblten_p     )
       if(allocated(rthdynten_p)    ) deallocate(rthdynten_p    )
       if(allocated(rthraten_p)     ) deallocate(rthraten_p     )
       if(allocated(rqvblten_p)     ) deallocate(rqvblten_p     )
       if(allocated(rqvdynten_p)    ) deallocate(rqvdynten_p    )
       if(allocated(rucuten_p)      ) deallocate(rucuten_p      )
       if(allocated(rvcuten_p)      ) deallocate(rvcuten_p      )
       if(allocated(rcnvcfcuten_p)  ) deallocate(rcnvcfcuten_p  )
       if(allocated(rbuoyxcuten_p)  ) deallocate(rbuoyxcuten_p  )

       if(allocated(tkepbl_p)       ) deallocate(tkepbl_p       )
       if(allocated(elpbl_p)        ) deallocate(elpbl_p        )
       if(allocated(rupmfxcu_p)     ) deallocate(rupmfxcu_p     )
       if(allocated(rdnmfxcu_p)     ) deallocate(rdnmfxcu_p     )
       if(allocated(rmfxdpcu_p)     ) deallocate(rmfxdpcu_p     )
       if(allocated(rmfxdncu_p)     ) deallocate(rmfxdncu_p     )
       if(allocated(v3d1_p)         ) deallocate(v3d1_p         )
       if(allocated(sub3d_rthcuten_p)) deallocate(sub3d_rthcuten_p)
       if(allocated(sub3d_rqvcuten_p)) deallocate(sub3d_rqvcuten_p)
       if(allocated(sub3d_rucuten_p) ) deallocate(sub3d_rucuten_p )
       if(allocated(sub3d_rvcuten_p) ) deallocate(sub3d_rvcuten_p )

    case ("cu_kain_fritsch")
       if(allocated(dx_p)         ) deallocate(dx_p         )
       if(allocated(area_p)       ) deallocate(area_p       )
       if(allocated(nca_p)        ) deallocate(nca_p        )
       if(allocated(cubot_p)      ) deallocate(cubot_p      )
       if(allocated(cutop_p)      ) deallocate(cutop_p      )
       if(allocated(w0avg_p)      ) deallocate(w0avg_p      )
       if(allocated(rqrcuten_p)   ) deallocate(rqrcuten_p   )
       if(allocated(rqscuten_p)   ) deallocate(rqscuten_p   )

    case ("cu_tiedtke","cu_ntiedtke")
       if(allocated(hfx_p)        ) deallocate(hfx_p        )
       if(allocated(qfx_p)        ) deallocate(qfx_p        )
       if(allocated(xland_p)      ) deallocate(xland_p      )
       if(allocated(rucuten_p)    ) deallocate(rucuten_p    )
       if(allocated(rvcuten_p)    ) deallocate(rvcuten_p    )

       cu_tiedtke_select: select case(convection_scheme)

          case ("cu_tiedtke")
             if(allocated(rqvdynten_p)  ) deallocate(rqvdynten_p  )
             if(allocated(rqvdynblten_p)) deallocate(rqvdynblten_p)

          case ("cu_ntiedtke")
             if(allocated(dx_p)     ) deallocate(dx_p)
             if(allocated(rqvften_p)) deallocate(rqvften_p)
             if(allocated(rthften_p)) deallocate(rthften_p)

          case default

       end select cu_tiedtke_select

    case default

 end select convection_select

 end subroutine deallocate_convection

!=================================================================================================================
 subroutine init_convection(dminfo,mesh,configs,diag_physics)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs
 type(dm_info),intent(in):: dminfo

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 logical,pointer:: config_do_restart
 character(len=StrKIND),pointer:: convection_scheme

 integer,pointer:: nCells
 real(kind=RKIND),dimension(:),pointer:: nca

!local variables:
 integer:: iCell

!local var for gf_monan
 integer, pointer :: config_gf_cporg,config_gf_pcvol,config_gf_gustf,config_gf_sub3d
 integer :: mynum, init_stat
!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_config(configs,'config_do_restart'       ,config_do_restart)
 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 convection_select: select case(convection_scheme)

    case ("cu_kain_fritsch")
       call mpas_pool_get_array(diag_physics,'nca',nca)
       if(.not. config_do_restart) then
          do iCell = 1, nCells
             nca(iCell) = -100._RKIND
          enddo
       endif
       call kf_lutab(svp1,svp2,svp3,svpt0)

    case ("cu_gf_monan")
       call mpas_pool_get_config(configs,'config_gf_gustf',config_gf_gustf)
       call mpas_pool_get_config(configs,'config_gf_cporg',config_gf_cporg)
       call mpas_pool_get_config(configs,'config_gf_pcvol',config_gf_pcvol)
       call mpas_pool_get_config(configs,'config_gf_sub3d',config_gf_sub3d)
       modConvParGF_initialized = .false.
       init_stat = initModConvParGF()
       mynum = dminfo % my_proc_id
       !print*," mynum", dminfo % my_proc_id
       call readGFConvParNML(mynum,IO_NODE,config_gf_cporg,config_gf_pcvol,config_gf_gustf &
                            ,config_gf_sub3d)
    case default

 end select convection_select

 end subroutine init_convection

!=================================================================================================================
 subroutine driver_convection(itimestep,configs,mesh,sfc_input,diag_physics,tend_physics,diag,dminfo,its,ite &
                             ,exchange_halo_group,block)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(in):: diag
 type(dm_info),intent(in):: dminfo
 type(block_type), intent(in) :: block
 procedure (halo_exchange_routine) :: exchange_halo_group

 integer,intent(in):: its,ite
 integer,intent(in):: itimestep

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!local variables and arrays:
 logical:: log_convection

 integer:: i,j
 integer:: icount,initflag

 real(kind=RKIND):: dx

!local pointers:
 logical,pointer:: config_do_restart
 integer,pointer:: gfconv_closure_deep,gfconv_closure_shallow
 character(len=StrKIND),pointer:: convection_scheme
 real(kind=RKIND),pointer:: len_disp

!variables specific to Kain_Fritsch parameterization:
 logical:: warm_rain,adapt_step_flag
 integer:: ktau
 real(kind=RKIND):: curr_secs
 real(kind=RKIND):: cudt
 real(kind=RKIND):: cudtacttime

!CCPP-compliant flags:
 character(len=StrKIND):: errmsg
 integer:: errflg

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine driver_convection:')

!initialization of CCPP-compliant flags:
 errmsg = ' '
 errflg = 0

 call mpas_pool_get_config(configs,'config_gfconv_closure_deep',gfconv_closure_deep)
 call mpas_pool_get_config(configs,'config_gfconv_closure_shallow',gfconv_closure_shallow)
 call mpas_pool_get_config(configs,'config_len_disp'         ,len_disp         )
 call mpas_pool_get_config(configs,'config_do_restart'       ,config_do_restart)
 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

!initialize instantaneous precipitation, and copy convective tendencies from the dynamics to
!the physics grid:
 call convection_from_MPAS(dt_dyn,configs,mesh,sfc_input,diag_physics,tend_physics,diag,its,ite)

!... convert the convection time-step to minutes:
 cudt = dt_cu/60.

!... call to convection schemes:
 curr_secs   = -1
 cudtacttime = -1
 adapt_step_flag = .false.
 do j = jts, jte
 do i = its, ite
    cu_act_flag(i,j) = .false.
 enddo
 enddo

!... initialization of initflag needed in the calls cu_tiedtke and cu_ntiedtke:
 initflag = 1
 if(config_do_restart .or. itimestep > 1) initflag = 0

 convection_select: select case(convection_scheme)

    case ("cu_gf_monan")
       call mpas_timer_start('Grell-Freitas-MONAN')
       call cu_gf_monan_driver(                                                           &
             dt              = dt_dyn                 , confrq        = dt_cu           , &
             dx_p            = dx_p                   , areaCell      = area_p          , &
             lats            = xlat_p                 , lons          = xlon_p          , &
             u               = u_p                    , v             = v_p             , &
             w               = w_p                    , temp          = t_p             , &
             rvap            = qv_p                   , rho           = rho_p           , &
             press           = pres_hyd_p             , pi            = pi_p            , &
             p8w             = pres2_hyd_p            , dz8w          = dz_p            , &
             topt            = ht_p                   , xland         = xland_p         , &
             sflux_t         = hfx_p                  , sflux_r       = qfx_p           , &
             temp2m          = t2m_p                  , wlpool        = wlpool_p        , &
             mpas_cape       = cape_p                 , mpas_cin      = cin_p           , &
             kpbl            = kpbl_p                 , tke_pbl       = tkepbl_p        , &
             turb_len_scale  = elpbl_p                ,                                   &
             buoyx           = buoyx_p                , cnvcf         = cnvcf_p         , &
             rthblten        = rthblten_p             , rqvblten      = rqvblten_p      , &
             rthratenlw      = rthratenlw_p           , rthratensw    = rthratensw_p    , &
             rthdyten        = rthdynten_p            , rqvdyten      = rqvdynten_p     , &
             
             raincv          = raincv_p               , conprr        = pratec_p        , &
             lightn_dens     = lightn_dens_p          , sigma_deep    = sigma_deep_p    , &
             rthcuten        = rthcuten_p             , rqvcuten      = rqvcuten_p      , &
             rqccuten        = rqccuten_p             , rqicuten      = rqicuten_p      , &
             rucuten         = rucuten_p              , rvcuten       = rvcuten_p       , &
             rbuoyxcuten     = rbuoyxcuten_p          , rcnvcfcuten   = rcnvcfcuten_p   , &

             sub3d_rthcuten  = sub3d_rthcuten_p       , sub3d_rqvcuten= sub3d_rqvcuten_p, &
             sub3d_rucuten   = sub3d_rucuten_p        , sub3d_rvcuten = sub3d_rvcuten_p , &

             rupmfxcu        = rupmfxcu_p             , rdnmfxcu      = rdnmfxcu_p      , &
             
             rmfxdpcu        = rmfxdpcu_p             , rmfxdncu      = rmfxdncu_p      , &
             rmfxmdcu        = rmfxmdcu_p                                               , &
             rmfxshcu        = rmfxshcu_p             , rtopdpcu      = rtopdpcu_p      , &
             rtopmdcu        = rtopmdcu_p             , rtopshcu      = rtopshcu_p      , &
             rbotdpcu        = rbotdpcu_p                                               , &
             var2d1          = v2d1_p                 , var2d2        = v2d2_p          , &
             var3d1          = v3d1_p                                                   , &
             ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde      , &
             ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme      , &
             its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte      , &
             itimestep       = itimestep              , mynum         = dminfo%my_proc_id )

       call mpas_timer_stop('Grell-Freitas-MONAN')

    case ("cu_grell_freitas")
       call mpas_timer_start('Grell-Freitas')
       call cu_grell_freitas( &
             itimestep       = itimestep              , dt            = dt_dyn          , &
             dxCell          = dx_p                   , areaCell      = area_p          , &
             u               = u_p                    , v             = v_p             , &
             w               = w_p                    , t             = t_p             , &
             q               = qv_p                   , rho           = rho_p           , &
             p               = pres_hyd_p             , pi            = pi_p            , &
             p8w             = pres2_hyd_p            , dz8w          = dz_p            , &
             ht              = ht_p                   , xland         = xland_p         , &
             gsw             = gsw_p                  , xlv           = xlv             , &
             cp              = cp                     , g             = gravity         , &
             r_v             = R_v                    , hfx           = hfx_p           , &
             qfx             = qfx_p                  , rthblten      = rthblten_p      , &
             rqvblten        = rqvblten_p             , rthften       = rthdynten_p     , &
             rqvften         = rqvdynten_p            , rthraten      = rthraten_p      , &
             kpbl            = kpbl_p                 , raincv        = raincv_p        , &
             pratec          = pratec_p               , htop          = cutop_p         , &
             hbot            = cubot_p                , k22_shallow   = k22_shallow_p   , & 
             kbcon_shallow   = kbcon_shallow_p        , ktop_shallow  = ktop_shallow_p  , & 
             xmb_total       = xmb_total_p            , xmb_shallow   = xmb_shallow_p   , &
             gdc             = qccu_p                 , gdc2          = qicu_p          , &
             rthcuten        = rthcuten_p             , rqvcuten      = rqvcuten_p      , & 
             rqccuten        = rqccuten_p             , rqicuten      = rqicuten_p      , &
             ichoice_deep    = gfconv_closure_deep                                      , &
             ichoice_shallow = gfconv_closure_shallow                                   , &
             ishallow_g3     = ishallow                                                 , &
             ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde      , &
             ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme      , &
             its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte        &
                            )
       call mpas_timer_stop('Grell-Freitas')

    case ("cu_kain_fritsch")
       if(itimestep == 1) then
          ktau = itimestep
       else
          ktau = itimestep + 1
       endif

       call mpas_timer_start('Kain-Fritsch')
       call kf_eta_cps ( &
             pcps            = pres_hyd_p      , t         = t_p        ,            &
             dt              = dt_dyn          , ktau      = ktau       ,            &
             dxCell          = dx_p            , areaCell  = area_p     ,            &
             cudt            = cudt            , curr_secs = curr_secs  ,            & 
             adapt_step_flag = adapt_step_flag , rho       = rho_p      ,            &
             raincv          = raincv_p        , pratec    = pratec_p   ,            &
             nca             = nca_p           , u         = u_p        ,            &
             v               = v_p             , th        = th_p       ,            &
             pi              = pi_p            , w         = w_p        ,            &
             dz8w            = dz_p            , w0avg     = w0avg_p    ,            &
             xlv0            = xlv0            , xlv1      = xlv1       ,            &
             xls0            = xls0            , xls1      = xls1       ,            &
             cp              = cp              , r         = r_d        ,            &
             g               = gravity         , ep1       = ep_1       ,            &
             ep2             = ep_2            , svp1      = svp1       ,            &
             svp2            = svp2            , svp3      = svp3       ,            &
             svpt0           = svpt0           , stepcu    = n_cu       ,            &
             cu_act_flag     = cu_act_flag     , warm_rain = warm_rain  ,            &
             cutop           = cutop_p         , cubot     = cubot_p    ,            &
             qv              = qv_p            , f_qr      = f_qr       ,            &
             f_qi            = f_qi            , f_qs      = f_qs       ,            &
             rthcuten        = rthcuten_p      , rqvcuten  = rqvcuten_p ,            &
             rqccuten        = rqccuten_p      , rqrcuten  = rqrcuten_p ,            &
             rqicuten        = rqicuten_p      , rqscuten  = rqscuten_p ,            &
             ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde , &
             ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme , &
             its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte   &
                        )
       call mpas_timer_stop('Kain-Fritsch')

    case("cu_tiedtke")
       call mpas_timer_start('Tiedtke')
       call cu_tiedtke( &
             pcps        = pres_hyd_p  , p8w       = pres2_hyd_p   ,                 &
             znu         = znu_hyd_p   , t3d       = t_p           ,                 &
             dt          = dt_dyn      , itimestep = initflag      ,                 &
             stepcu      = n_cu        , raincv    = raincv_p      ,                 &
             pratec      = pratec_p    , qfx       = qfx_p         ,                 &
             u3d         = u_p         , v3d       = v_p           ,                 &
             w           = w_p         , qv3d      = qv_p          ,                 &
             qc3d        = qc_p        , qi3d      = qi_p          ,                 &
             pi3d        = pi_p        , rho3d     = rho_p         ,                 &
             qvften      = rqvdynten_p , qvpblten  = rqvdynblten_p ,                 &
             dz8w        = dz_p        , xland     = xland_p       ,                 &
             cu_act_flag = cu_act_flag , f_qc      = f_qc          ,                 &
             f_qi        = f_qi        , rthcuten  = rthcuten_p    ,                 &
             rqvcuten    = rqvcuten_p  , rqccuten  = rqccuten_p    ,                 &
             rqicuten    = rqicuten_p  , rucuten   = rucuten_p     ,                 &
             rvcuten     = rvcuten_p   ,                                             &
             ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde , &
             ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme , &
             its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte   &
                      )
       call mpas_timer_stop('Tiedtke')

    case("cu_ntiedtke")
       call mpas_timer_start('cu_ntiedtke')
       call cu_ntiedtke_driver( &
             pcps     = pres_hyd_p , p8w         = pres2_hyd_p ,                     &
             t3d      = t_p        , dz8w        = dz_p        ,                     &
             dt       = dt_dyn     , itimestep   = initflag    ,                     &
             stepcu   = n_cu       , raincv      = raincv_p    ,                     &
             pratec   = pratec_p   , qfx         = qfx_p       ,                     &
             hfx      = hfx_p      , xland       = xland_p     ,                     &
             dx       = dx_p       , u3d         = u_p         ,                     &
             v3d      = v_p        , w           = w_p         ,                     &
             qv3d     = qv_p       , qc3d        = qc_p        ,                     &
             qi3d     = qi_p       , pi3d        = pi_p        ,                     &
             rho3d    = rho_p      , qvften      = rqvften_p   ,                     &
             thften   = rthften_p  , cu_act_flag = cu_act_flag ,                     &
             f_qc     = f_qc       , f_qi        = f_qi        ,                     &
             rthcuten = rthcuten_p , rqvcuten    = rqvcuten_p  ,                     &
             rqccuten = rqccuten_p , rqicuten    = rqicuten_p  ,                     &
             rucuten  = rucuten_p  , rvcuten     = rvcuten_p   ,                     &
             grav     = gravity    , xlf         = xlf         ,                     &
             xls      = xls        , xlv         = xlv         ,                     &
             rd       = R_d        , rv          = R_v         ,                     &
             cp       = cp         ,                                                 &
             errmsg   = errmsg     , errflg      = errflg      ,                     &
             ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde , &
             ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme , &
             its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte   &
                              )
       call mpas_timer_stop('cu_ntiedtke')

    case default

 end select convection_select

!copy instantaneous and accumulated precipitation, convective tendencies, and "other" arrays
!specific to convection parameterization back to the dynamics grid:
 call convection_to_MPAS(configs,mesh,diag_physics,tend_physics,its,ite,exchange_halo_group,block)

!call mpas_log_write('--- end subroutine driver_convection.')
!call mpas_log_write('')

 end subroutine driver_convection

!=================================================================================================================
 subroutine convection_from_MPAS(dt_dyn,configs,mesh,sfc_input,diag_physics,tend_physics,diag,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(in):: diag_physics
 type(mpas_pool_type),intent(in):: tend_physics
 type(mpas_pool_type),intent(in):: diag

 integer,intent(in):: its,ite
 real(kind=RKIND),intent(in):: dt_dyn

!local pointers:
 character(len=StrKIND),pointer:: convection_scheme,config_pbl_scheme
 integer,dimension(:),pointer:: kpbl,k22_shallow,kbcon_shallow,ktop_shallow,ktop_deep
 real(kind=RKIND),dimension(:),pointer  :: areaCell,meshDensity
 real(kind=RKIND),dimension(:),pointer  :: nca,cubot,cutop,cuprec,raincv,wlpool,lightn_dens,sigma_deep
 real(kind=RKIND),dimension(:),pointer  :: gsw,hfx,qfx,xland,t2m,cape,cin
 real(kind=RKIND),dimension(:),pointer  :: xmb_total,xmb_shallow
 real(kind=RKIND),dimension(:),pointer  :: latCell,lonCell
 real(kind=RKIND),dimension(:,:),pointer:: zgrid
 real(kind=RKIND),dimension(:,:),pointer:: qc_cu,qi_cu,tke_pbl,el_pbl,qke
 real(kind=RKIND),dimension(:,:),pointer:: w0avg
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten,rqvcuten,rqccuten,rqicuten,rqrcuten,rqscuten 
 real(kind=RKIND),dimension(:,:),pointer:: rthblten,rthdynten,rthratenlw,rthratensw
 real(kind=RKIND),dimension(:,:),pointer:: rqvblten,rqvdynten,rucuten,rvcuten
 real(kind=RKIND),dimension(:,:),pointer:: rcnvcfcuten, rbuoyxcuten

!local variables:
 integer:: i,j,k
 integer:: iEdge
 real(kind=RKIND),pointer:: len_disp

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 if(trim(convection_scheme) /= "cu_gf_monan") then 
 call mpas_pool_get_array(diag_physics,'cuprec',cuprec)
 call mpas_pool_get_array(diag_physics,'raincv',raincv)

 call mpas_pool_get_array(tend_physics,'rthcuten',rthcuten)
 call mpas_pool_get_array(tend_physics,'rqvcuten',rqvcuten)
 call mpas_pool_get_array(tend_physics,'rqccuten',rqccuten)
 call mpas_pool_get_array(tend_physics,'rqicuten',rqicuten)

 do j = jts,jte
 do i = its,ite
    raincv_p(i,j) = raincv(i)
    pratec_p(i,j) = cuprec(i)
    do k = kts,kte
       rthcuten_p(i,k,j) = rthcuten(k,i)
       rqvcuten_p(i,k,j) = rqvcuten(k,i)
       rqccuten_p(i,k,j) = rqccuten(k,i)
       rqicuten_p(i,k,j) = rqicuten(k,i)
    enddo
 enddo
 enddo
 end if
 convection_select: select case(convection_scheme)

     case ("cu_gf_monan")

         call mpas_pool_get_config(configs,'config_len_disp',len_disp)

         call mpas_pool_get_array(mesh,'areaCell'   ,areaCell   )
         call mpas_pool_get_array(mesh,'meshDensity',meshDensity)
         call mpas_pool_get_array(mesh,'zgrid'      ,zgrid      )
         call mpas_pool_get_array(mesh,'latCell'    ,latCell    )
         call mpas_pool_get_array(mesh,'lonCell'    ,lonCell    )

         call mpas_pool_get_array(sfc_input,'xland',xland)

         call mpas_pool_get_array(diag,'cape'               ,cape    )
         call mpas_pool_get_array(diag,'cin'                ,cin     )
         call mpas_pool_get_array(diag_physics,'hfx'        ,hfx     )
         call mpas_pool_get_array(diag_physics,'qfx'        ,qfx     )
         call mpas_pool_get_array(diag_physics,'t2m'        ,t2m     )
         call mpas_pool_get_array(diag_physics,'kpbl'       ,kpbl    )
         call mpas_pool_get_array(diag_physics,'wlpool'     ,wlpool  )
         call mpas_pool_get_array(diag_physics,'lightn_dens',lightn_dens)
         call mpas_pool_get_array(diag_physics,'sigma_deep' ,sigma_deep)

         call mpas_pool_get_array(tend_physics,'rqvblten'   ,rqvblten  )
         call mpas_pool_get_array(tend_physics,'rqvdynten'  ,rqvdynten )
         call mpas_pool_get_array(tend_physics,'rthblten'   ,rthblten  )
         call mpas_pool_get_array(tend_physics,'rthdynten'  ,rthdynten )
         call mpas_pool_get_array(tend_physics,'rthratenlw' ,rthratenlw)
         call mpas_pool_get_array(tend_physics,'rthratensw' ,rthratensw)
         call mpas_pool_get_array(tend_physics,'rcnvcfcuten',rcnvcfcuten)
         call mpas_pool_get_array(tend_physics,'rbuoyxcuten',rbuoyxcuten)

         do j = jts,jte
         do i = its,ite
            dx_p(i,j) = len_disp / meshDensity(i)**0.25
            area_p(i,j)          = areaCell(i)
            ht_p(i,j)            = zgrid(1,i)
            xland_p(i,j)         = xland(i)
            hfx_p(i,j)           = hfx(i)
            qfx_p(i,j)           = qfx(i)
            t2m_p(i,j)           = t2m(i)
            kpbl_p(i,j)          = kpbl(i)
            cape_p(i,j)          = cape(i)
            cin_p(i,j)           = cin(i)
            wlpool_p(i,j)        = wlpool(i)
            lightn_dens_p(i,j)   = lightn_dens(i)
            xlat_p(i,j)          = latCell(i) / degrad
            xlon_p(i,j)          = lonCell(i) / degrad
            sigma_deep_p(i,j)    = sigma_deep(i)

            do k = kts,kte
 
               rthblten_p(i,k,j)    = rthblten(k,i)
               rthdynten_p(i,k,j)   = rthdynten(k,i)
               rthratenlw_p(i,k,j)  = rthratenlw(k,i) 
               rthratensw_p(i,k,j)  = rthratensw(k,i)
               rqvblten_p(i,k,j)    = rqvblten(k,i)
               rqvdynten_p(i,k,j)   = rqvdynten(k,i)
               rbuoyxcuten_p(i,k,j) = 0.0
               rcnvcfcuten_p(i,k,j) = 0.0
              
            enddo
         enddo
         enddo
         
         call mpas_pool_get_config(configs, 'config_pbl_scheme', config_pbl_scheme)
         if (trim(config_pbl_scheme) == "bl_mynn") then
               !call mpas_pool_get_array(diag_physics,'tke_pbl'   ,tke_pbl )
               call mpas_pool_get_array(diag_physics,'qke'       ,qke )
               call mpas_pool_get_array(diag_physics,'el_pbl'    ,el_pbl  )
               do j = jts,jte
               do i = its,ite               
                  do k = kts,kte
                     !--- srf: for some reason, in this MPAS version "tke_pbl" is always zero
                     !--- using 'QKE' which is defined as 2 * TKE 
                     !tkepbl_p(i,k,j) = tke_pbl(k,i)
                     tkepbl_p(i,k,j) = 0.5*qke(k,i)
                     elpbl_p(i,k,j)  = el_pbl(k,i)                    
                  enddo
               enddo
               enddo
         end if

    case ("cu_grell_freitas")
       call mpas_pool_get_config(configs,'config_len_disp',len_disp)

       call mpas_pool_get_array(mesh,'areaCell'   ,areaCell   )
       call mpas_pool_get_array(mesh,'meshDensity',meshDensity)
       call mpas_pool_get_array(mesh,'zgrid'      ,zgrid      )

       call mpas_pool_get_array(sfc_input,'xland',xland)

       call mpas_pool_get_array(diag_physics,'gsw'          ,gsw          )
       call mpas_pool_get_array(diag_physics,'hfx'          ,hfx          )
       call mpas_pool_get_array(diag_physics,'qfx'          ,qfx          )
       call mpas_pool_get_array(diag_physics,'kpbl'         ,kpbl         )
       call mpas_pool_get_array(diag_physics,'cubot'        ,cubot        )
       call mpas_pool_get_array(diag_physics,'cutop'        ,cutop        )
       call mpas_pool_get_array(diag_physics,'xmb_total'    ,xmb_total    )
       call mpas_pool_get_array(diag_physics,'xmb_shallow'  ,xmb_shallow  )
       call mpas_pool_get_array(diag_physics,'k22_shallow'  ,k22_shallow  )
       call mpas_pool_get_array(diag_physics,'kbcon_shallow',kbcon_shallow)
       call mpas_pool_get_array(diag_physics,'ktop_shallow' ,ktop_shallow )
       call mpas_pool_get_array(diag_physics,'ktop_deep'    ,ktop_deep    )
       call mpas_pool_get_array(diag_physics,'qc_cu'        ,qc_cu        )
       call mpas_pool_get_array(diag_physics,'qi_cu'        ,qi_cu        )

       call mpas_pool_get_array(tend_physics,'rqvblten'  ,rqvblten  )
       call mpas_pool_get_array(tend_physics,'rqvdynten' ,rqvdynten )
       call mpas_pool_get_array(tend_physics,'rthblten'  ,rthblten  )
       call mpas_pool_get_array(tend_physics,'rthdynten' ,rthdynten )
       call mpas_pool_get_array(tend_physics,'rthratenlw',rthratenlw)
       call mpas_pool_get_array(tend_physics,'rthratensw',rthratensw)
       call mpas_pool_get_array(tend_physics,'rucuten'   ,rucuten   )
       call mpas_pool_get_array(tend_physics,'rvcuten'   ,rvcuten   )

       do j = jts,jte
       do i = its,ite
          dx_p(i,j) = len_disp / meshDensity(i)**0.25
          area_p(i,j)          = areaCell(i)
          ht_p(i,j)            = zgrid(1,i)
          xland_p(i,j)         = xland(i)
          gsw_p(i,j)           = gsw(i)
          hfx_p(i,j)           = hfx(i)
          qfx_p(i,j)           = qfx(i)
          kpbl_p(i,j)          = kpbl(i)
          cubot_p(i,j)         = cubot(i)
          cutop_p(i,j)         = cutop(i)

          xmb_total_p(i,j)     = xmb_total(i)
          xmb_shallow_p(i,j)   = xmb_shallow(i)

          k22_shallow_p(i,j)   = k22_shallow(i)
          kbcon_shallow_p(i,j) = kbcon_shallow(i)
          ktop_shallow_p(i,j)  = ktop_shallow(i)
          ktop_deep_p(i,j)     = ktop_deep(i)

          do k = kts,kte
             qccu_p(i,k,j)      = qc_cu(k,i)
             qicu_p(i,k,j)      = qi_cu(k,i)

             rthblten_p(i,k,j)  = rthblten(k,i)
             rthdynten_p(i,k,j) = rthdynten(k,i)
             rthraten_p(i,k,j)  = rthratenlw(k,i) + rthratensw(k,i)
             rqvblten_p(i,k,j)  = rqvblten(k,i)
             rqvdynten_p(i,k,j) = rqvdynten(k,i)
             rucuten_p(i,k,j)   = rucuten(k,i)
             rvcuten_p(i,k,j)   = rvcuten(k,i)
          enddo
       enddo
       enddo

    case ("cu_kain_fritsch")
       call mpas_pool_get_config(configs,'config_len_disp',len_disp)

       call mpas_pool_get_array(mesh,'areaCell',areaCell)
       call mpas_pool_get_array(mesh,'meshDensity',meshDensity)

       call mpas_pool_get_array(diag_physics,'cubot'   ,cubot   )
       call mpas_pool_get_array(diag_physics,'cutop'   ,cutop   )
       call mpas_pool_get_array(diag_physics,'nca'     ,nca     )
       call mpas_pool_get_array(diag_physics,'w0avg'   ,w0avg   )

       call mpas_pool_get_array(tend_physics,'rqrcuten',rqrcuten)
       call mpas_pool_get_array(tend_physics,'rqscuten',rqscuten)

       do j = jts,jte
       do i = its,ite
          !area of grid-cell:       
          dx_p(i,j)    = len_disp / meshDensity(i)**0.25
          area_p(i,j)  = areaCell(i)
          cubot_p(i,j) = cubot(i)
          cutop_p(i,j) = cutop(i)

          do k = kts,kte
             rqrcuten_p(i,k,j) = rqrcuten(k,i)
             rqscuten_p(i,k,j) = rqscuten(k,i)
          enddo

          !decreases the characteristic time period that convection remains active. When nca_p
          !becomes less than the convective timestep, convective tendencies and precipitation
          !are reset to zero (note that this is also done in subroutine kf_eta_cps).
          nca_p(i,j) = nca(i)

          do k = kts,kte
             w0avg_p(i,k,j) = w0avg(k,i)
          enddo
       enddo
       enddo

    case ("cu_tiedtke","cu_ntiedtke")
       call mpas_pool_get_array(sfc_input,'xland',xland)
       call mpas_pool_get_array(diag_physics,'hfx',hfx)
       call mpas_pool_get_array(diag_physics,'qfx',qfx)

       call mpas_pool_get_array(tend_physics,'rqvblten' ,rqvblten )
       call mpas_pool_get_array(tend_physics,'rqvdynten',rqvdynten)
       call mpas_pool_get_array(tend_physics,'rucuten'  ,rucuten  )
       call mpas_pool_get_array(tend_physics,'rvcuten'  ,rvcuten  )

       do j = jts,jte
       do i = its,ite
          xland_p(i,j) = xland(i)
          hfx_p(i,j)   = hfx(i)
          qfx_p(i,j)   = qfx(i)
       enddo

       do k = kts,kte
       do i = its,ite
          rucuten_p(i,k,j) = rucuten(k,i)
          rvcuten_p(i,k,j) = rvcuten(k,i)
       enddo
       enddo
       enddo

       cu_tiedtke_select: select case(convection_scheme)

          case ("cu_tiedtke")
             do j = jts,jte
             do k = kts,kte
             do i = its,ite
                rqvdynblten_p(i,k,j) = rqvblten(k,i) 
                rqvdynten_p(i,k,j)   = rqvdynten(k,i)
             enddo
             enddo
             enddo
 
          case("cu_ntiedtke")
             call mpas_pool_get_config(configs,'config_len_disp',len_disp)
             call mpas_pool_get_array(mesh,'meshDensity',meshDensity)

             call mpas_pool_get_array(tend_physics,'rqvdynten' ,rqvdynten )
             call mpas_pool_get_array(tend_physics,'rqvblten'  ,rqvblten  )
             call mpas_pool_get_array(tend_physics,'rthdynten' ,rthdynten )
             call mpas_pool_get_array(tend_physics,'rthblten'  ,rthblten  )
             call mpas_pool_get_array(tend_physics,'rthratenlw',rthratenlw)
             call mpas_pool_get_array(tend_physics,'rthratensw',rthratensw)
             
             do j = jts,jte
             do i = its,ite
                dx_p(i,j) = len_disp / meshDensity(i)**0.25
             enddo

             do k = kts,kte
             do i = its,ite
                rqvften_p(i,k,j) = rqvdynten(k,i) + rqvblten(k,i)
                rthften_p(i,k,j) = (rthdynten(k,i) + rthblten(k,i) &
                                 + rthratenlw(k,i) + rthratensw(k,i)) * pi_p(i,k,j) 
             enddo
             enddo
             enddo

          case default

       end select cu_tiedtke_select

    case default

 end select convection_select

 end subroutine convection_from_MPAS

!=================================================================================================================
 subroutine convection_to_MPAS(configs,mesh,diag_physics,tend_physics,its,ite,exchange_halo_group,block)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(block_type), intent(in) :: block
 procedure (halo_exchange_routine) :: exchange_halo_group

 integer,intent(in):: its,ite

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!local variables:
 integer:: i,k,j,iCell,ineighbour

!local pointers:
 character(len=StrKIND),pointer:: convection_scheme
 integer,dimension(:),pointer:: k22_shallow,kbcon_shallow,ktop_shallow,ktop_deep
 real(kind=RKIND),dimension(:),pointer  :: nca,cubot,cutop,cuprec,raincv,lightn_dens,sigma_deep
 real(kind=RKIND),dimension(:),pointer  :: xmb_total,xmb_shallow
 real(kind=RKIND),dimension(:,:),pointer:: qc_cu,qi_cu
 real(kind=RKIND),dimension(:,:),pointer:: w0avg
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten,rqvcuten,rqccuten,rqicuten,rqrcuten,rqscuten
 real(kind=RKIND),dimension(:,:),pointer:: rucuten,rvcuten
 !local pointers specific for gf_monan convection scheme
 integer,pointer:: config_gf_sub3d,nCellsSolve
 integer, dimension(:), pointer:: nEdgesOnCell
 integer, dimension(:,:), pointer:: cellsOnCell
 real(kind=RKIND),dimension(:),pointer  :: rmfxdpcu,rmfxdncu,rmfxmdcu,rmfxshcu
 real(kind=RKIND),dimension(:),pointer  :: rtopdpcu,rtopmdcu,rtopshcu
 real(kind=RKIND),dimension(:),pointer  :: rbotdpcu,v2d1,v2d2
 real(kind=RKIND),dimension(:,:),pointer:: v3d1
 real(kind=RKIND),dimension(:,:),pointer:: rcnvcfcuten, rbuoyxcuten
 real(kind=RKIND),dimension(:,:),pointer:: rupmfxcu, rdnmfxcu
 real(kind=RKIND),dimension(:,:),pointer:: sub3d_rucuten,sub3d_rvcuten,sub3d_rthcuten,sub3d_rqvcuten

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 call mpas_pool_get_array(diag_physics,'cuprec',cuprec)
 call mpas_pool_get_array(diag_physics,'raincv',raincv)

 call mpas_pool_get_array(tend_physics,'rthcuten',rthcuten)
 call mpas_pool_get_array(tend_physics,'rqvcuten',rqvcuten)
 call mpas_pool_get_array(tend_physics,'rqccuten',rqccuten)
 call mpas_pool_get_array(tend_physics,'rqicuten',rqicuten)

 do j = jts,jte
 do i = its,ite
    raincv(i) = raincv_p(i,j)
    cuprec(i) = pratec_p(i,j)
    do k = kts, kte
       rthcuten(k,i) = rthcuten_p(i,k,j)
       rqvcuten(k,i) = rqvcuten_p(i,k,j)
       rqccuten(k,i) = rqccuten_p(i,k,j)
       rqicuten(k,i) = rqicuten_p(i,k,j)
    enddo
 enddo
 enddo

 convection_select: select case(convection_scheme)
    case ("cu_gf_monan")
       call mpas_pool_get_array(diag_physics,'rmfxdpcu'    ,rmfxdpcu   )
       call mpas_pool_get_array(diag_physics,'rmfxdncu'    ,rmfxdncu   )
       call mpas_pool_get_array(diag_physics,'rmfxmdcu'    ,rmfxmdcu   )
       call mpas_pool_get_array(diag_physics,'rmfxshcu'    ,rmfxshcu   )
       call mpas_pool_get_array(diag_physics,'rtopdpcu'    ,rtopdpcu   )
       call mpas_pool_get_array(diag_physics,'rtopmdcu'    ,rtopmdcu   )
       call mpas_pool_get_array(diag_physics,'rtopshcu'    ,rtopshcu   )
       call mpas_pool_get_array(diag_physics,'rbotdpcu'    ,rbotdpcu   )
       call mpas_pool_get_array(diag_physics,'v2d1'        ,v2d1       )
       call mpas_pool_get_array(diag_physics,'v2d2'        ,v2d2       )
       call mpas_pool_get_array(diag_physics,'lightn_dens' ,lightn_dens)
       call mpas_pool_get_array(diag_physics,'sigma_deep'  ,sigma_deep )

       call mpas_pool_get_array(diag_physics,'rupmfxcu'    ,rupmfxcu   )
       call mpas_pool_get_array(diag_physics,'rdnmfxcu'    ,rdnmfxcu   )
       call mpas_pool_get_array(diag_physics,'v3d1'        ,v3d1       )

       call mpas_pool_get_array(tend_physics,'rucuten'     ,rucuten    )
       call mpas_pool_get_array(tend_physics,'rvcuten'     ,rvcuten    )
       call mpas_pool_get_array(tend_physics,'rbuoyxcuten' ,rbuoyxcuten)
       call mpas_pool_get_array(tend_physics,'rcnvcfcuten' ,rcnvcfcuten)


       do j = jts,jte
       do i = its,ite
          do k = kts,kte
             rucuten(k,i)      = rucuten_p(i,k,j)
             rvcuten(k,i)      = rvcuten_p(i,k,j)
             rbuoyxcuten(k,i)  = rbuoyxcuten_p(i,k,j)
             rcnvcfcuten(k,i)  = rcnvcfcuten_p(i,k,j)
             rupmfxcu(k,i)     = rupmfxcu_p(i,k,j)
             rdnmfxcu(k,i)     = rdnmfxcu_p(i,k,j)
         enddo
       enddo
       enddo
       !- 2-d
       do j = jts,jte
       do i = its,ite
          lightn_dens(i) = lightn_dens_p(i,j)
          sigma_deep(i) = sigma_deep_p(i,j)
          rmfxdpcu(i) = rmfxdpcu_p(i,j)  
          rmfxdncu(i) = rmfxdncu_p(i,j)  
          rmfxmdcu(i) = rmfxmdcu_p(i,j)  
          rmfxshcu(i) = rmfxshcu_p(i,j)  
          rtopdpcu(i) = rtopdpcu_p(i,j)  
          rtopmdcu(i) = rtopmdcu_p(i,j)  
          rtopshcu(i) = rtopshcu_p(i,j)  
          rbotdpcu(i) = rbotdpcu_p(i,j)
       enddo
       enddo

       !- dummy vars for debug/output 
       do j = jts,jte
       do i = its,ite
          v2d1(i)  = v2d1_p(i,j)
         !--- used in sfclayer
         !v2d2(i)  = v2d2_p(i,j)
          do k = kts,kte
             v3d1(k,i) = v3d1_p(i,k,j)
          enddo
       enddo
       enddo

       !-- including option of the env subsidence lateral spreading              
       call mpas_pool_get_config(configs,'config_gf_sub3d',config_gf_sub3d)
       
       if(config_gf_sub3d > 0) then
            call mpas_pool_get_array(tend_physics,'sub3d_rucuten'  ,sub3d_rucuten )
            call mpas_pool_get_array(tend_physics,'sub3d_rvcuten'  ,sub3d_rvcuten )
            call mpas_pool_get_array(tend_physics,'sub3d_rqvcuten' ,sub3d_rqvcuten)
            call mpas_pool_get_array(tend_physics,'sub3d_rthcuten' ,sub3d_rthcuten)

            do j = jts,jte
            do i = its,ite
               do k = kts,kte
                  sub3d_rthcuten(k,i) = sub3d_rthcuten_p(i,k,j)
                  sub3d_rqvcuten(k,i) = sub3d_rqvcuten_p(i,k,j)
                  sub3d_rucuten(k,i)  = sub3d_rucuten_p(i,k,j)
                  sub3d_rvcuten(k,i)  = sub3d_rvcuten_p(i,k,j)
               enddo
            enddo
            enddo
            !
            !
            !-- part 1: remove the subsidence tendencies from the total tendencies
            !
            call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
            do i = 1, nCellsSolve
            do k = 1, nVertLevels
                 rucuten(k,i)  = rucuten(k,i)  - sub3d_rucuten(k,i)
                 rvcuten(k,i)  = rvcuten(k,i)  - sub3d_rvcuten(k,i)
                 rthcuten(k,i) = rthcuten(k,i) - sub3d_rthcuten(k,i)
                 rqvcuten(k,i) = rqvcuten(k,i) - sub3d_rqvcuten(k,i)
            enddo
            enddo

            !-- part 2: broadcast the subsidence tendencies
            call exchange_halo_group(block % domain, 'physics:sub3d_cuten')

            !-- part 3: spread the subsidence tendencies to the 1st neighboors
            call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
            call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)

            do iCell = 1, nCellsSolve
              do k = 1, nVertLevels
                 rthcuten(k,iCell) = rthcuten(k,iCell) + sigma_deep(iCell) * sub3d_rthcuten(k,iCell)
                 rqvcuten(k,iCell) = rqvcuten(k,iCell) + sigma_deep(iCell) * sub3d_rqvcuten(k,iCell)
                 rucuten(k,iCell)  = rucuten(k,iCell)  + sigma_deep(iCell) * sub3d_rucuten(k,iCell)
                 rvcuten(k,iCell)  = rvcuten(k,iCell)  + sigma_deep(iCell) * sub3d_rvcuten(k,iCell)
                 
                 !-- loop over the neighbours of the grid cell 'iCell'
                 do i = 1, nEdgesOnCell(iCell)
                   ineighbour = cellsOnCell(i,iCell)
                   rthcuten(k,iCell) = rthcuten(k,iCell) + &
                                       (1._RKIND - sigma_deep(ineighbour))*sub3d_rthcuten(k,ineighbour)/ &
                                       real(nEdgesOnCell(ineighbour),RKIND)  
                   rqvcuten(k,iCell) = rqvcuten(k,iCell) + &
                                       (1._RKIND - sigma_deep(ineighbour))*sub3d_rqvcuten(k,ineighbour)/ &
                                       real(nEdgesOnCell(ineighbour),RKIND)  
                   rucuten(k,iCell) = rucuten(k,iCell) + &
                                       (1._RKIND - sigma_deep(ineighbour))*sub3d_rucuten(k,ineighbour)/ &
                                       real(nEdgesOnCell(ineighbour),RKIND)  
                   rvcuten(k,iCell) = rvcuten(k,iCell) + &
                                       (1._RKIND - sigma_deep(ineighbour))*sub3d_rvcuten(k,ineighbour)/ &
                                       real(nEdgesOnCell(ineighbour),RKIND)
                 enddo
              enddo
            enddo

       endif

    case ("cu_grell_freitas")
       call mpas_pool_get_array(diag_physics,'cubot'        ,cubot        )
       call mpas_pool_get_array(diag_physics,'cutop'        ,cutop        )
       call mpas_pool_get_array(diag_physics,'xmb_total'    ,xmb_total    )
       call mpas_pool_get_array(diag_physics,'xmb_shallow'  ,xmb_shallow  )
       call mpas_pool_get_array(diag_physics,'k22_shallow'  ,k22_shallow  )
       call mpas_pool_get_array(diag_physics,'kbcon_shallow',kbcon_shallow)
       call mpas_pool_get_array(diag_physics,'ktop_shallow' ,ktop_shallow )
       call mpas_pool_get_array(diag_physics,'ktop_deep'    ,ktop_deep    )
       call mpas_pool_get_array(diag_physics,'qc_cu'        ,qc_cu        )
       call mpas_pool_get_array(diag_physics,'qi_cu'        ,qi_cu        )

       call mpas_pool_get_array(tend_physics,'rucuten',rucuten)
       call mpas_pool_get_array(tend_physics,'rvcuten',rvcuten)

       do j = jts,jte
       do i = its,ite
          cubot(i)         = cubot_p(i,j)
          cutop(i)         = cutop_p(i,j)
          xmb_total(i)     = xmb_total_p(i,j)
          xmb_shallow(i)   = xmb_shallow_p(i,j)
          k22_shallow(i)   = k22_shallow_p(i,j)
          kbcon_shallow(i) = kbcon_shallow_p(i,j)
          ktop_shallow(i)  = ktop_shallow_p(i,j)
          ktop_deep(i)     = ktop_deep_p(i,j)

          do k = kts,kte
             qc_cu(k,i)    = qccu_p(i,k,j)
             qi_cu(k,i)    = qicu_p(i,k,j)
             rucuten(k,i)  = rucuten_p(i,k,j)
             rvcuten(k,i)  = rvcuten_p(i,k,j)
          enddo
       enddo
       enddo

    case ("cu_kain_fritsch")
       call mpas_pool_get_array(diag_physics,'cubot',cubot)
       call mpas_pool_get_array(diag_physics,'cutop',cutop)
       call mpas_pool_get_array(diag_physics,'nca'  ,nca  )
       call mpas_pool_get_array(diag_physics,'w0avg',w0avg)

       call mpas_pool_get_array(tend_physics,'rqrcuten',rqrcuten)
       call mpas_pool_get_array(tend_physics,'rqscuten',rqscuten)

       do j = jts,jte
       do i = its,ite
          cubot(i) = cubot_p(i,j)
          cutop(i) = cutop_p(i,j)
          nca(i)   = nca_p(i,j)
          do k = kts, kte
             w0avg(k,i)    = w0avg_p(i,k,j)
             rqrcuten(k,i) = rqrcuten_p(i,k,j)
             rqscuten(k,i) = rqscuten_p(i,k,j)
          enddo                          
       enddo
       enddo

    case ("cu_tiedtke","cu_ntiedtke")
       call mpas_pool_get_array(tend_physics,'rucuten',rucuten)
       call mpas_pool_get_array(tend_physics,'rvcuten',rvcuten)

       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          rucuten(k,i) = rucuten_p(i,k,j)
          rvcuten(k,i) = rvcuten_p(i,k,j)
       enddo
       enddo
       enddo

    case default

 end select convection_select

 end subroutine convection_to_MPAS

!=================================================================================================================
 subroutine update_convection_step1(configs,diag_physics,tend_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 integer,intent(in):: its,ite

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!local pointers:
 character(len=StrKIND),pointer:: convection_scheme
 real(kind=RKIND),dimension(:),pointer  :: nca,cubot,cutop,cuprec,raincv
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten,rqvcuten,rqccuten,rqicuten,rqrcuten,rqscuten

!local variables and arrays:
 integer:: i,k

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 convection_select: select case(convection_scheme)

    case ("cu_kain_fritsch")
      call mpas_pool_get_array(diag_physics,'nca'   ,nca   )
      call mpas_pool_get_array(diag_physics,'cubot' ,cubot )
      call mpas_pool_get_array(diag_physics,'cutop' ,cutop )
      call mpas_pool_get_array(diag_physics,'cuprec',cuprec)
      call mpas_pool_get_array(diag_physics,'raincv',raincv)

      call mpas_pool_get_array(tend_physics,'rthcuten',rthcuten)
      call mpas_pool_get_array(tend_physics,'rqvcuten',rqvcuten)
      call mpas_pool_get_array(tend_physics,'rqccuten',rqccuten)
      call mpas_pool_get_array(tend_physics,'rqicuten',rqicuten)
      call mpas_pool_get_array(tend_physics,'rqrcuten',rqrcuten)
      call mpas_pool_get_array(tend_physics,'rqscuten',rqscuten)

       do i = its, ite
          !decreases the characteristic time period that convection remains active. When nca_p
          !becomes less than the convective timestep, convective tendencies and precipitation
          !are reset to zero (note that this is also done in subroutine kf_eta_cps).
          if(nca(i) .gt. 0.) then
             nca(i) = nca(i) - dt_dyn
             
             if(nca(i) .lt. 0.5*dt_dyn) then
                do k = kts,kte
                   rthcuten(k,i) = 0._RKIND
                   rqvcuten(k,i) = 0._RKIND
                   rqccuten(k,i) = 0._RKIND
                   rqrcuten(k,i) = 0._RKIND
                   rqicuten(k,i) = 0._RKIND
                   rqscuten(k,i) = 0._RKIND
                enddo
                raincv(i) = 0._RKIND
                cuprec(i) = 0._RKIND
                cubot(i)  = kte+1
                cutop(i)  = kts
             endif
          endif
       enddo

    case default

 end select convection_select

 end subroutine update_convection_step1

!=================================================================================================================
 subroutine update_convection_step2(configs,mesh,sfc_input,diag_physics,tend_physics,state,diag,time_lev,its,ite,  &
                                    exchange_halo_group,block)
!=================================================================================================================

!input arguments:
 type(block_type), intent(in) :: block
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: diag
 procedure (halo_exchange_routine) :: exchange_halo_group
 integer,intent(in):: its,ite
 integer,intent(in):: time_lev

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!local pointers:
 integer,pointer:: index_cnvcf, index_buoyx
 integer,pointer:: config_gf_cporg, config_gf_pcvol
 integer,dimension(:),pointer:: i_rainc
 integer,pointer:: nEdgesSolve,nCells,nCellsSolve,nEdges

!local pointers:
 character(len=StrKIND),pointer:: convection_scheme
 real(kind=RKIND),pointer:: bucket_rainc
 real(kind=RKIND),dimension(:),pointer:: cuprec,rainc,lightn_dens,acc_lightn_dens
 real(kind=RKIND),dimension(:),pointer:: xland
 real(kind=RKIND),dimension(:),pointer:: wlpool,vcpool,u_gustfront,v_gustfront

 real(kind=RKIND),dimension(:,:),pointer:: umcl,vmcl
 real(kind=RKIND),dimension(:,:),pointer:: rumcl
 real(kind=RKIND),dimension(:,:),pointer:: zgrid,rh
 real(kind=RKIND),dimension(:,:),pointer:: cnvcf,buoyx
 real(kind=RKIND),dimension(:,:),pointer:: rcnvcfcuten,rbuoyxcuten
 real(kind=RKIND),dimension(:,:),pointer:: rho_edge    ! diag rho_edge
 real(kind=RKIND),dimension(:,:),pointer:: rho_zz      ! state rho at center

 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

!local variables and arrays for the 'cu_gf_monan':
 integer:: i,k,j,k_amb_wind
 real(kind=RKIND),parameter:: alp = 0.5                      ! semi-implicit, 1=explicit
 real(kind=RKIND),parameter:: tau_cloud_diss0 = 1.0e-3       ! seconds^(-1)
 real(kind=RKIND),parameter:: width  = 1000.                 ! meters
 real(kind=RKIND),parameter:: turncf = 2000.                 ! meters above local terrain
 real(kind=RKIND),parameter:: Kfr = 0.9                      ! internal Froude number 
 real(kind=RKIND),parameter:: epsx = 100.                    ! threshold 
 real(kind=RKIND),parameter:: slope_pool = 45.*(3.1416/180.) ! using mean value of Reif et al 2020
 real(kind=RKIND),parameter:: startlev = 1000.               ! meters above local surface
 real(kind=RKIND),parameter:: endlev   = 4000.               ! meters above local surface

 real(kind=RKIND):: vert_scale_diss,tau_cloud_diss,htopx,snk,src,tau_cp,denom 
 real(kind=RKIND):: total_dz,H_env,gama,temp, rvap, press, qes, zlevel, dzrho, aux, MCL_speed
!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_bucket_rainc',bucket_rainc)

 call mpas_pool_get_array(diag_physics,'i_rainc',i_rainc)
 call mpas_pool_get_array(diag_physics,'cuprec' ,cuprec )
 call mpas_pool_get_array(diag_physics,'rainc'  ,rainc  )

!update the accumulated precipitation at the end of each dynamic time-step:
 do i = its, ite
    rainc(i) = rainc(i) + cuprec(i) * dt_dyn

    if(l_acrain .and. bucket_rainc.gt.0._RKIND .and. &
       rainc(i).gt.bucket_rainc) then
       i_rainc(i) = i_rainc(i) + 1
       rainc(i)   = rainc(i) - bucket_rainc
    endif
 enddo
 !CR: negative rainc values correction:
 !where (rainc < 0._RKIND)
 !  rainc = 0._RKIND
 !endwhere

 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 convection_select: select case(convection_scheme)
    case ("cu_gf_monan")      
      !--- update the accumulated lightning density at the end of each dynamic time-step:
      !--- if(lightning_diag == 1) then -- for future
      call mpas_pool_get_array(diag_physics,'lightn_dens'      ,lightn_dens)
      call mpas_pool_get_array(diag_physics,'acc_lightn_dens'  ,acc_lightn_dens)
      do i = its, ite
        acc_lightn_dens(i) = acc_lightn_dens(i) + lightn_dens(i) * dt_dyn
      enddo

      call mpas_pool_get_config(configs,'config_gf_pcvol',config_gf_pcvol)
      call mpas_pool_get_config(configs,'config_gf_cporg',config_gf_cporg)

      if(config_gf_pcvol > 0) then 
          !-- get cnvcf scalar
          call mpas_pool_get_dimension(state,'index_cnvcf',index_cnvcf)
          call mpas_pool_get_array(state,'scalars',scalars,time_lev)
          cnvcf => scalars(index_cnvcf,:,:)
          !-- get the cnvcf tendency from convection 
          call mpas_pool_get_array(tend_physics,'rcnvcfcuten',rcnvcfcuten)

          call mpas_pool_get_array(sfc_input,'xland' ,xland )
          call mpas_pool_get_array(mesh,     'zgrid' ,zgrid )
          call mpas_pool_get_array(diag,     'relhum',rh    )! pct 

          do i = its,ite
              htopx = zgrid(1,i) !- local terrain
              do k = kts,kte
                !-- this goes from 1 to 0.1 as the height goes from the local surface height to the upper levels
                vert_scale_diss = 1._RKIND/(1.5_RKIND+atan(2._RKIND*(zgrid(k,i) - htopx - turncf)/width))/3.      
                vert_scale_diss = max(min (vert_scale_diss,1.) , 0.1_RKIND) ! = 1. to turn off vertical variation
                tau_cloud_diss = 1._RKIND/(vert_scale_diss*tau_cloud_diss0)
                
                !---- sink term for cloud fraction
                snk = dt_dyn * (1._RKIND/tau_cloud_diss) * (1.0_RKIND - 0.01_RKIND*rh(k,i)) * abs(cnvcf(k,i))
 
                !---- source term for cloud fraction
                src = dt_dyn * rcnvcfcuten(k,i)
 
                !- using semi-impl formulation;
                !- f(t+dt)  = ( f(t) - (dt/tau)* (alp* f(t)) ] / (1+(dt/tau)* (1-alp))
                denom = 1._RKIND + (1._RKIND-alp)*src + (1.-alp)*(1._RKIND-0.01_RKIND*rh(k,i))*dt_dyn/tau_cloud_diss

                !---- update cnvcf
                cnvcf(k,i) = ( cnvcf(k,i) - alp*snk + (1._RKIND-alp*cnvcf(k,i))*src )/denom
                
                !---- fuly explicit (Tiedtke 1993)
                !cnvcf(k,i) = cnvcf(k,i) + (1.-cnvcf(k,i))*src - snk

                !---- fully explicit (alp=1) (ECMWF and Morcrette (2012) formulation)
                !snk = snk * (1.-cnvcf(k,i))
                !cnvcf(k,i) = cnvcf(k,i) + (1.-cnvcf(k,i)*src - snk

               !---- sanity check for the cloud fraction
                cnvcf(k,i) = min(1.0_RKIND, max(0._RKIND, cnvcf(k,i)))
             enddo
          enddo

      endif 
      
      if(config_gf_cporg == 1) then
          call mpas_pool_get_dimension(state,'index_buoyx',index_buoyx)
          call mpas_pool_get_array(state,'scalars',scalars,time_lev)
          buoyx => scalars(index_buoyx,:,:)
          call mpas_pool_get_array(tend_physics,'rbuoyxcuten',rbuoyxcuten)
          call mpas_pool_get_array(sfc_input,'xland' ,xland )
          call mpas_pool_get_array(mesh,     'zgrid' ,zgrid )

          do i = its,ite
             !<var name="xland" type="real" dimensions="nCells Time" units="unitless"
             !description="land-ocean mask (1=land including sea-ice ; 2=ocean)"/>
             if(xland(i) < 1.5_RKIND ) tau_cp = tau_land_cp
             if(xland(i) >= 1.5_RKIND) tau_cp = tau_ocea_cp

              do k = kts,kte                
                !---- sink term for the buoyancy excess
                snk = dt_dyn*(1./tau_cp) * abs(buoyx(k,i))

                !---- source term for the buoyancy excess
                src = dt_dyn*rbuoyxcuten(k,i)
                
                !- using semi-impl formulation;
                !- f(t+dt)  = ( f(t) - (dt/tau)* (alp* f(t)) ] / (1+(dt/tau)* (1-alp))
                buoyx(k,i) = (buoyx(k,i) + src - snk) / (1._RKIND + (dt_dyn/tau_cp)*(1._RKIND-alp))
               
                !---- sanity check for the buoyancy excess
                buoyx(k,i) = max(0._RKIND, buoyx(k,i))
             enddo
          enddo
          !-- Max vertical velocity at the leading of the cold pool (based on Reif et al 2020)
          call mpas_pool_get_array(diag_physics,'wlpool'  ,wlpool  )
          
          !-- Gust front horizontal velocity relative to the storm 
          !-- (eq 8.5 page 346 Storm and Cloud Dynamics, Cotton et al. 2nd edition.)
          call mpas_pool_get_array(diag_physics,'vcpool'  ,vcpool  )
          
          do j = jts,jte
          do i = its,ite
             wlpool(i) = 0.0_RKIND
             vcpool(i) = 0.0_RKIND
             total_dz  = 0.0_RKIND
             
             loopK: do k = kts, kte
               
                if(buoyx(k,i) < epsx .or. total_dz .gt. 2000._RKIND ) exit loopK 
                total_dz  = total_dz + dz_p(i,k,j)
                temp  = t_p(i,k,j)
                press = pres_p(i,k,j)
                rvap  = qv_p(i,k,j)
                
                qes=satur_spec_hum(temp,press)
                gama=(xlv/cp)*(xlv/(R_v*temp*temp))*qes
                gama= 1._RKIND/(1._RKIND+gama)

                vcpool(i) = vcpool(i) + (gravity* dz_p(i,k,j)* min(mx_buoy2, buoyx(k,i))*gama/(cp*temp))
                wlpool(i) = wlpool(i) + (gravity* dz_p(i,k,j)* min(mx_buoy2, buoyx(k,i))*gama/(cp*temp))

            enddo loopK
            vcpool(i) = min(20._RKIND, Kfr *sqrt(vcpool(i)))

            !-- estimating the slope of head cold pool, option using RKW theory
            ! rkw = exp(-abs(vcpool(i)-wshear(i))/(vcpool(i)+wshear(i)+1.e-6))
            ! slope_pool(i) = max(45.,abs(asind(rkw)))

            wlpool(i) = min(10._RKIND, Kfr *sin(slope_pool)*sqrt(wlpool(i)))
          enddo
          enddo
          !
          !--- mean cloud layer horizontal wind (MCL)
          !
          call mpas_pool_get_array(diag_physics,'umcl'  ,umcl )
          call mpas_pool_get_array(diag_physics,'vmcl'  ,vmcl )
          !
          !--- gust front speed at surface level (1st model layer)
          !
          call mpas_pool_get_array(diag_physics,'u_gustfront'  ,u_gustfront )
          call mpas_pool_get_array(diag_physics,'v_gustfront'  ,v_gustfront )

          do j = jts,jte
          do i = its,ite
            umcl(1,i)  = 0.0_RKIND
            vmcl(1,i)  = 0.0_RKIND
            total_dz = 0.0_RKIND
            
            loopK_2: do k = kts, kte
                     !-- height level above local surface 
                     zlevel = zgrid(k,i) - zgrid(1,i)
                     if(zlevel < startlev .or. zlevel > endlev) cycle loopK_2  
                     dzrho    = dz_p(i,k,j)*rho_p(i,k,j)
                     umcl(1,i)  = umcl(1,i)  + u_p(i,k,j)*dzrho
                     vmcl(1,i)  = vmcl(1,i)  + v_p(i,k,j)*dzrho
                     total_dz = total_dz + dzrho
            enddo loopK_2
            umcl(1,i) = umcl(1,i)/(total_dz+1.e-12_RKIND)
            vmcl(1,i) = vmcl(1,i)/(total_dz+1.e-12_RKIND)        
          enddo
          enddo
          !-- add cold pool gust velocity to the MCL wind, keeping the original MCL wind direction
          !-- add the influence of ambient flow to the MCL wind
          !-- 'factor' = 0.62 is suggested by Simpson and Britter (1980) to add the cold pool velocity 
          !-- the to ambient wind (see page 346 Storm and Cloud Dynamics, Cotton et al. 2nd edition)
          k_amb_wind = 2
          do j = jts,jte
          do i = its,ite
            MCL_speed = sqrt (umcl(1,i)**2 + vmcl(1,i)**2)
            aux = (MCL_speed + vcpool(i))/(MCL_speed+1.e-12_RKIND)
            umcl(1,i) = aux * umcl(1,i)
            vmcl(1,i) = aux * vmcl(1,i)
           
            !-- including the influence of the ambient flow to the MCL wind
            umcl(1,i) = umcl(1,i) + 0.62_RKIND * u_p(i,k_amb_wind,j)
            vmcl(1,i) = vmcl(1,i) + 0.62_RKIND * v_p(i,k_amb_wind,j)
          enddo
          enddo
          !-- save the gust front velocities
          u_gustfront(:) = umcl(1,:) 
          v_gustfront(:) = vmcl(1,:)                           
          !
          !-- get dry air density divided by d(zeta)/dz"
          call mpas_pool_get_array(state, 'rho_zz', rho_zz,time_lev)
          !-- copying the value @ level one to levels 2:nvertlevels
          !-- (the 1-d (ncells) arrays umcl,vmcl are converted to 2-d (nvert,ncells) because 
          !-- the current MPAS paralelization does not support exchange of halos of 1-d arrays)
          !-- convert to mass flux u*rho/zz (momentum)
          do j = jts,jte
          do i = its,ite
            do k = 2, nVertLevels  
               umcl(k,i) = umcl(1,i)*rho_zz(k,i)
               vmcl(k,i) = vmcl(1,i)*rho_zz(k,i)
            enddo
            !-- do not change this
            umcl(1,i) = umcl(1,i)*rho_zz(1,i)
            vmcl(1,i) = vmcl(1,i)*rho_zz(1,i)
          enddo
          enddo
          !
          call mpas_pool_get_dimension(mesh, 'nEdgesSolve', nEdgesSolve)
          !-
          !-- calculating the horizontal mean cloud layer momentum at cell edge (rho*u/zz)
          !-- used in the cold pool gust front 2-d advection for the convection 
          !---organization parameterization.
          !-- exchange halos of rho*umcl and rho*vmcl (see src/core_atmosphere/mpas_atm_core.F)
          call exchange_halo_group(block % domain, 'physics:coldpool')
          !-- array for the horizontal mean cloud layer momentum at cell edge
          call mpas_pool_get_array(diag, 'rumcl' ,rumcl)
          !-- calculate the horizontal mean cloud layer momentum at cell edge
          call UVMCL_toEdges(mesh,UMCL,VMCL,RUMCL)
          
          !----------
          !-- the below does not work; there is no reproducibility using a different
          !-- number of cores/threads.
          !
          !-- get rho_zz averaged from cell centers to the cell edge
          !call mpas_pool_get_array(diag , 'rho_edge', rho_edge)
          !do i = 1, nEdgesSolve
          !   do k = 1, nVertLevels  
          !      rumcl(k,i)=rumcl(k,i)*rho_edge(k,i)
          !   enddo
          !enddo
          !---------
      endif  !- (config_gf_cporg == 1)

    case default

 end select convection_select
 end subroutine update_convection_step2
!=================================================================================================================
 !-- alternative routine for "tend_toEdges" --!
 subroutine UVMCL_toEdges(mesh,UMCL,VMCL,RUMCL)
  !-- calculating the horizontal mean cloud layer momentum at cell edge (rho*u/zz)
  !-- used in the cold pool gust front advection for the convection
  !-- organization parameterization.
  
  use mpas_atm_dimensions

 !input arguments:
 ! type(block_type),intent(in),target:: block
  type(mpas_pool_type),intent(in):: mesh
  real(kind=RKIND),intent(in),dimension(:,:),target:: umcl,vmcl 

 !output arguments:
  real(kind=RKIND),intent(out),dimension(:,:):: rumcl

 !local variables:
  integer:: iCell,iEdge,k,j
  integer:: cell1, cell2
  integer,pointer:: nCells,nCellsSolve,nEdges
  integer,dimension(:,:),pointer:: cellsOnEdge

  real(kind=RKIND), dimension(:,:), pointer :: east, north, edgeNormalVectors

  call mpas_pool_get_dimension(mesh, 'nCells', nCells)
  call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
  call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)

  call mpas_pool_get_array(mesh, 'east', east)
  call mpas_pool_get_array(mesh, 'north', north)
  call mpas_pool_get_array(mesh, 'edgeNormalVectors', edgeNormalVectors)

  call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)

  !--srf: rumcl does not vary on vertical, see Freitas et al. (JAMES 2024) for explanation.
  !--     It is used to transport the scalar 'buoyx' representing the MSE anomalies in the 
  !--     leading edge of the cold pool.
  !print*,'ncells=',nCells,nCellsSolve, nEdges

  do iEdge = 1, nEdges
     cell1 = cellsOnEdge(1,iEdge)
     cell2 = cellsOnEdge(2,iEdge)
     rumcl(:,iEdge) =   umcl(:,cell1) * 0.5 * (edgeNormalVectors(1,iEdge) * east(1,cell1)   &
                                          +    edgeNormalVectors(2,iEdge) * east(2,cell1)   &
                                          +    edgeNormalVectors(3,iEdge) * east(3,cell1))  &
                      + vmcl(:,cell1) * 0.5 * (edgeNormalVectors(1,iEdge) * north(1,cell1)  &
                                          +    edgeNormalVectors(2,iEdge) * north(2,cell1)  &
                                          +    edgeNormalVectors(3,iEdge) * north(3,cell1)) &
                      + umcl(:,cell2) * 0.5 * (edgeNormalVectors(1,iEdge) * east(1,cell2)   &
                                          +    edgeNormalVectors(2,iEdge) * east(2,cell2)   &
                                          +    edgeNormalVectors(3,iEdge) * east(3,cell2))  &
                      + vmcl(:,cell2) * 0.5 * (edgeNormalVectors(1,iEdge) * north(1,cell2)  &
                                          +    edgeNormalVectors(2,iEdge) * north(2,cell2)  &
                                          +    edgeNormalVectors(3,iEdge) * north(3,cell2))
   
  end do
  
  end subroutine UVMCL_toEdges

!=================================================================================================================
 end module mpas_atmphys_driver_convection
!=================================================================================================================
