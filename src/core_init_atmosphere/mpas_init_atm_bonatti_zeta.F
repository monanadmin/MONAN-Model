! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_init_atm_bonatti_zeta

 use mpas_kind_types
 use mpas_derived_types, only : MPAS_LOG_WARN
 use mpas_log, only : mpas_log_write

 implicit none
 private

 public :: bonatti_zeta3_levels

 contains


!==================================================================================================
 subroutine bonatti_zeta3_levels(config_ztop, config_bonatti_zeta0, config_bonatti_fm, nz1, zw, dzw)
!==================================================================================================
!
! Cubic zeta3 vertical coordinate proposed by J. P. Bonatti (paulo.bonatti@inpe.br,
! CPTEC/INPE) for the MONAN 2.0 model, presented at the "Reuniao GCC MONAN" of
! 09/09/2026 ("Proposta para Coordenadas Verticais Zeta para o MONAN 2.0").
!
! The objective of this coordinate is to increase vertical resolution in the
! low atmosphere relative to the pre-existing zeta1 (power-law, alpha = 1.5)
! and zeta2 (piecewise exponential/cubic, controlled by als/alt/zetal) options
! already present in init_atm_case_gfs, while keeping a single smooth cubic
! polynomial in the fractional level index all the way to the model top.
!
! Let F_k = (k-1)/(KI-1), k = 1, ..., KI, where KI = nz1+1 is the number of
! vertical interfaces (F_1 = 0 at the surface, F_KI = 1 at the model top).
! The coordinate is defined as
!
!     zeta3(k) = Ztopo * (A0 + A1*F_k + A2*F_k**2 + A3*F_k**3)
!
! and the four coefficients A0..A3 are obtained by imposing:
!
!   (1) zeta3(1)  = 0                     => A0 = 0
!   (2) zeta3(KI) = Ztopo                 => A1 + A2 + A3 = 1
!   (3) zeta3(2)  = config_bonatti_zeta0  => target thickness of the first
!                                            model layer above the surface
!   (4) zeta3(F_M) = Ztopo/2, with 0.5 <= F_M <= 0.8, a mid-atmosphere blend
!                                            point controlling how quickly the
!                                            profile transitions to its coarser,
!                                            upper-atmosphere spacing
!
! This subroutine is deliberately self-contained (no MPAS pools, no mesh) so
! that it can be exercised by a small standalone test program independently
! of the rest of the init_atmosphere core; see test_bonatti_zeta3.f90.
!
!--------------------------------------------------------------------------------------------------

    implicit none

!   input arguments:
    real (kind=RKIND), intent(in) :: config_ztop            ! model top height, Ztopo (m)
    real (kind=RKIND), intent(in) :: config_bonatti_zeta0    ! target thickness of first layer, zeta0 (m)
    real (kind=RKIND), intent(in) :: config_bonatti_fm       ! mid-blend fractional level, F_M (0.5 <= F_M <= 0.8)
    integer,           intent(in) :: nz1                     ! number of layers, KC = KI - 1 = nVertLevels

!   output arguments:
    real (kind=RKIND), dimension(nz1+1), intent(out) :: zw   ! heights of the KI = nz1+1 zeta interfaces (m)
    real (kind=RKIND), dimension(nz1),   intent(out) :: dzw  ! thicknesses of the KC = nz1 layers (m)

!   local variables:
    real (kind=RKIND) :: zt
    real (kind=RKIND) :: a1, a2, a3
    real (kind=RKIND) :: cf1, cf2, cf3, cf4, cf5, cf6
    real (kind=RKIND) :: f2, fm, z0, zm, factor
    integer :: k

    zt = config_ztop
    z0 = config_bonatti_zeta0

    fm = config_bonatti_fm
    if (fm < 0.5_RKIND .or. fm > 0.8_RKIND) then
       call mpas_log_write('*******************************************************************', &
                            messageType=MPAS_LOG_WARN)
       call mpas_log_write('config_bonatti_fm should satisfy 0.5 <= config_bonatti_fm <= 0.8;', &
                            messageType=MPAS_LOG_WARN)
       call mpas_log_write('resetting config_bonatti_fm to its default value of 0.8', &
                            messageType=MPAS_LOG_WARN)
       call mpas_log_write('*******************************************************************', &
                            messageType=MPAS_LOG_WARN)
       fm = 0.8_RKIND
    end if

    f2 = 1.0_RKIND/real(nz1, RKIND)     ! F_2, the fractional height of the first interface above the surface

    zm = 0.5_RKIND*zt
    if (zm < z0) zm = z0
    if (zm > zt) zm = zt

    cf1 = z0/zt - f2
    cf2 = zm/zt - fm
    cf3 = f2*f2       - f2
    cf4 = fm*fm       - fm
    cf5 = f2*f2*f2    - f2
    cf6 = fm*fm*fm    - fm

    a3 = (cf1*cf4 - cf3*cf2) / (cf5*cf4 - cf3*cf6)
    a2 = (cf2 - cf6*a3) / cf4
    a1 = 1.0_RKIND - a2 - a3

    call mpas_log_write('')
    call mpas_log_write('--- Bonatti zeta3 vertical coordinate (CPTEC/INPE, MONAN 2.0) ---')
    call mpas_log_write('--- config_bonatti_zeta0 (m) = $r', realArgs=(/z0/))
    call mpas_log_write('--- config_bonatti_fm        = $r', realArgs=(/fm/))
    call mpas_log_write('--- A0, A1, A2, A3           = 0.0 $r $r $r', realArgs=(/a1,a2,a3/))
    call mpas_log_write('')

    do k=1,nz1+1
       factor = real(k-1, RKIND)/real(nz1, RKIND)
       zw(k) = zt*(a1*factor + a2*factor**2 + a3*factor**3)
       if (k > 1) dzw(k-1) = zw(k)-zw(k-1)
    end do

 end subroutine bonatti_zeta3_levels

 end module mpas_init_atm_bonatti_zeta
