module modConstants
   !! ## Constants and parameters
   !!
   !! ![](https://i.ibb.co/LNqGy3S/logo-Monan-Color-75x75.png)
   !! ## MONAN
   !!
   !! Author: Rodrigues, L.F. [LFR]
   !!
   !! E-mail: <mailto:luizfrodrigues@protonmail.com>
   !!
   !! Date: 09Fevereiro2023 08:55
   !!
   !! #####Version: 0.1.0
   !!
   !! ---
   !! **Full description**:
   !!
   !! Constants and parameters
   !!
   !! ** History**:
   !!
   !! --- 
   !! ** Licence **:
   !!
   !!  <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
   !!
   !!  This program is free software: you can redistribute it and/or modify
   !!  it under the terms of the GNU General Public License as published by
   !!  the  Free  Software  Foundation, either version 3 of the License, or
   !!  (at your option) any later version.
   !!
   !!  This program is distributed in the hope that it  will be useful, but
   !!  ** WITHOUT  ANY  WARRANTY **;  without  even  the   implied   warranty  of
   !!  **MERCHANTABILITY** or **FITNESS FOR A  PARTICULAR PURPOSE**.  See  the, GNU
   !!  GNU General Public License for more details.
   !!
   !!  You should have received a copy  of the GNU General  Public  License
   !!  along with this program.  If not, see [GNU Public License](https://www.gnu.org/licenses/gpl-3.0.html).
   !!
   
   implicit none
   character(len=*), parameter :: sourceName = 'constants.f90' 
   !! Source code name 
   character(len=*), parameter :: moduleName = 'modConstants' 
   !! module name 

   !integer kinds
   integer, parameter :: i8 = selected_int_kind(8)
   !# 8 byte integer
   integer, parameter :: i4 = selected_int_kind(4)
   !# 4 byte integer
   integer, parameter :: i0 = kind(1)
   !# native integer

   !real kinds
   integer, parameter :: r8 = selected_real_kind(8)
   !# 8 byte real
   integer, parameter :: r4 = selected_real_kind(4)
   !# 4 byte real
   integer, parameter :: r0 = kind(1.0)
   !# native real

   character(len=*), parameter :: c_modelVersion='Rev. 0.1.0'
   character(len=*), parameter :: modelVersion='Rev. 0.1.0'
   character(len=*), parameter :: c_license='[GNU Public License](https://www.gnu.org/licenses/gpl-3.0.html)'
   !# Last version of model

   !# Phys & others constants
   real(kind=r4), parameter :: c_pi       = 3.1415926535897932384626433
   real(kind=r4), parameter :: c_rgas     = 287.
   real(kind=r4), parameter :: c_rgas_rd  = 287.06
   real(kind=r4), parameter :: c_cp       = 1004.
   real(kind=r4), parameter :: c_cp_pd    = 1004.71
   real(kind=r4), parameter :: c_cv       = 717.
   real(kind=r4), parameter :: c_rm       = 461.
   real(kind=r4), parameter :: c_rm_rd    = 461.52
   real(kind=r4), parameter :: c_p00      = 1.e5
   real(kind=r4), parameter :: c_t00      = 273.16
   real(kind=r4), parameter :: c_t01      = 273.155
   real(kind=r4), parameter :: c_t100     = c_t00+100.
   real(kind=r4), parameter :: c_Tice     = 235.16

   !Absolute temperature
   real(kind=r4), parameter :: c_pi180    = c_pi / 180.
   real(kind=r4), parameter :: c_i_pi180  = 1./c_pi180
   real(kind=r4), parameter :: c_pi4      = c_pi * 4.
   real(kind=r4), parameter :: c_spcon    = 111120.
   real(kind=r4), parameter :: c_erad     = 6367000.
   real(kind=r4), parameter :: c_arc      = c_erad*c_pi180
   real(kind=r4), parameter :: c_vonk     = 0.40
   real(kind=r4), parameter :: p_tkmin    = 1.e-5
   !# Minimum TKE [J/kg]
   real(kind=r4), parameter :: c_grav      = 9.80665
   !# Gravity acceleration [m/s]
   real(kind=r4), parameter :: c_alvl     = 2.50e6
   real(kind=r4), parameter :: c_alvi     = 2.834e6
   real(kind=r4), parameter :: c_alli     = 0.334e6
   real(kind=r4), parameter :: c_alvl2    = 6.25e12
   real(kind=r4), parameter :: c_alvi2    = 8.032e12
   real(kind=r4), parameter :: c_solar    = 1.3533e3
   real(kind=r4), parameter :: c_stefan   = 5.6696e-8
   real(kind=r4), parameter :: c_cww      = 4218.
   real(kind=r4), parameter :: c_c0       = 752.55 * 4.18684e4
   real(kind=r4), parameter :: c_viscos   = .15e-4
   real(kind=r4), parameter :: c_rowt     = 1.e3
   real(kind=r4), parameter :: c_dlat     = 111120.
   real(kind=r4), parameter :: c_omega    = 7.292e-5
   real(kind=r4), parameter :: c_rocp     = c_rgas / c_cp
   real(kind=r4), parameter :: c_p00i     = 1. / c_p00
   real(kind=r4), parameter :: c_cpor     = c_cp / c_rgas
   real(kind=r4), parameter :: c_rocv     = c_rgas / c_cv
   real(kind=r4), parameter :: c_cpi      = 1. / c_cp
   real(kind=r4), parameter :: c_cpi4     = 4. * c_cpi
   real(kind=r4), parameter :: c_cp253i   = c_cpi / 253.
   real(kind=r4), parameter :: c_allii    = 1. / c_alli
   real(kind=r4), parameter :: c_aklv     = c_alvl / c_cp
   real(kind=r4), parameter :: c_akiv     = c_alvi / c_cp
   real(kind=r4), parameter :: c_gama     = c_cp / c_cv
   real(kind=r4), parameter :: c_gg       = .5 * c_grav
   real(kind=r4), parameter :: c_ep       = c_rgas / c_rm
   real(kind=r4), parameter :: c_p00k     = 26.870941
   !#  = p00 ** rocp
   real(kind=r4), parameter :: c_p00ki    = 1. / c_p00k

   real(kind=r4), parameter :: c_tcrit    = 258.
   real(kind=r4), parameter :: c_akmin    = 1.0
   real(kind=r4), parameter :: c_ccnclean = 250.
   real(kind=r4), parameter :: c_t_ice    = 235.16
   real(kind=r4), parameter :: c_xlf      = 0.333e6
   !! ! latent heat of freezing (J kg-1)
   real(kind=r4), parameter :: p_max_qsat = 0.5
   real(kind=r4), parameter :: p_smaller_qv = 1.e-16
   real, parameter :: p_mintracer = tiny(1.)
   real, parameter :: p_xmbmaxshal = 0.05
   real, parameter :: p_ccnclean= 250.   ! # cm-3
   real, parameter :: c_temp0 =    298.15
   !! standard temperature [K]
   real, parameter :: c_temp0i= 1./c_temp0
   !! inverse of standard temperature [K]
   real, parameter :: c_rgas_atm = 8.205e-2 
   ! atm M^-1 K^-1 ! 8.314 gas constant [J/(mol*K)]
   real(kind=r4), parameter :: c_hplus = 1.175e-4     
   !!  for cloud water. pH is asuumed to be 3.93: pH=3.93 =>hplus=10**(-pH)
   real(kind=r4), parameter :: c_retv = c_rm_rd/c_rgas_rd - 1.0
   real(kind=r4), parameter :: c_r2es = 611.21*c_rgas_rd/c_rm_rd
   real(kind=r4), parameter :: c_rtt  = 273.16 
   real(kind=r4), parameter :: c_r3les = 17.502
   real(kind=r4), parameter :: c_r3ies = 22.587
   real(kind=r4), parameter :: c_r4les = 32.19
   real(kind=r4), parameter :: c_r4ies = -0.7
   real(kind=r4), parameter :: c_rtwat= c_rtt 
   real(kind=r4), parameter :: c_rtice = c_t00 - 23.
   real(kind=r4), parameter :: c_rticecu = c_t00 - 23.
   real(kind=r4), parameter :: c_rtwat_rtice_r = 1./(c_t00 - c_rtice)
   real(kind=r4), parameter :: c_rtwat_rticecu_r = 1./(c_t00 - c_rticecu)
   real(kind=r4), parameter :: c_r5les = c_r3les*(c_t00 - c_r4les)
   real(kind=r4), parameter :: c_r5ies = c_r3ies*(c_t00 - c_r4ies)
   real(kind=r4), parameter :: c_rlvtt = 2.5008e+6
   real(kind=r4), parameter :: c_rlstt = 2.8345e+6
   real(kind=r4), parameter :: c_rlmlt = c_rlstt - c_rlvtt
   real(kind=r4), parameter :: c_rho_h2o = 1000.
   real(kind=r4), parameter :: c_r5alvcp = c_r5les*c_rlvtt/c_cp_pd
   real(kind=r4), parameter :: c_r5alscp = c_r5ies*c_rlstt/c_cp_pd
   real(kind=r4), parameter :: c_ralvdcp = c_rlvtt/c_cp_pd
   real(kind=r4), parameter :: c_ralsdcp = c_rlstt/c_cp_pd
   real(kind=r4), parameter :: c_ralfdcp = c_rlmlt/c_cp_pd
   real(kind=r4), parameter :: c_rcpv = 4.*c_rm_rd
   real(kind=r4), parameter :: c_rvtmp2 = c_rcpv/c_cp_pd - 1.0
   real(kind=r4), parameter :: c_rtber = c_t00 - 5. 
   real(kind=r4), parameter :: c_rtbercu = c_t00 - 5.0 

   real(kind=r4), parameter :: p_zqmax = 0.5

   real(kind=r4), parameter :: c_onethird  = 1./3.
   !# 1/3
   real(kind=r4), parameter :: c_half  = 1./2.
   !# 1/2

   ! Lower bounds for turbulence-related variables                                        !
   real(kind=r4), parameter :: c_sigwmin     = 1.e-4
   !# Minimum sigma-w                     [m/s]
   real(kind=r4), parameter :: c_abslmomin   = 1.e-4
   !# Minimum abs value of Obukhov length [m]
   real(kind=r4), parameter :: c_ltscalemax  = 1.e5
   !# Maximum Lagrangian timescale        [s]
   real(kind=r4), parameter :: c_abswltlmin  = 1.e-4
   !# Minimum abs value of Theta*         [K m/s]
   real(kind=r4), parameter :: c_lturbmin    = 1.e-3
   !# Minimum abs value of turb. lenght   [m]

   real,parameter :: c_scday=86400.0
   !# Seconds by day

   !Tiny numbers
   real(kind=r4), parameter :: tinyReal = 1.17549435E-38
   !# Tiny real number
   !real(kind=kind_rb), parameter :: tinyDouble = 2.2250738585072014E-308
   !# Tiny double precision number

   !Errors type for dump functions
   integer, parameter :: c_noError=0
   !# No Error, just write
   integer, parameter :: c_notice=1
   !# Notice Error
   integer, parameter :: c_warning=2
   !# warning Error
   integer, parameter :: c_fatal=3
   !# Fatal error
   integer, parameter :: c_kill=4
   !# Kill the model Signal
   integer, parameter :: c_continue=5
   !# Continue run Signal

   logical, parameter :: c_yes=.true.
   !# Yes is the .true. value
   logical, parameter :: c_no=.false.
   !# No is the false value

   character(len=*), parameter :: c_empty=''
   !# empty string

   !Colors strings for print and write commands
   character(len=*), parameter :: c_darkGrey   =achar(27)//'[90m'
   character(len=*), parameter :: c_peach      =achar(27)//'[91m'
   character(len=*), parameter :: c_lightGreen =achar(27)//'[92m'
   character(len=*), parameter :: c_lightYellow=achar(27)//'[93m'
   character(len=*), parameter :: c_lightBlue  =achar(27)//'[94m'
   character(len=*), parameter :: c_pink       =achar(27)//'[95m'
   character(len=*), parameter :: c_lightAqua  =achar(27)//'[96m'
   character(len=*), parameter :: c_pearlWhite =achar(27)//'[97m'
   character(len=*), parameter :: c_black      =achar(27)//'[30m'
   character(len=*), parameter :: c_red        =achar(27)//'[31m'
   character(len=*), parameter :: c_green      =achar(27)//'[32m'
   character(len=*), parameter :: c_yellow     =achar(27)//'[33m'
   character(len=*), parameter :: c_blue       =achar(27)//'[34m'
   character(len=*), parameter :: c_purple     =achar(27)//'[35m'
   character(len=*), parameter :: c_aqua       =achar(27)//'[36m'
   character(len=*), parameter :: c_blink      =achar(27)//'[31;5;95;38;5;214m'
   character(len=*), parameter :: c_noColor    =achar(27)//'[0m'
   character(len=*), parameter :: c_Iblack     =achar(27)//'[40m'
   character(len=*), parameter :: c_Ired       =achar(27)//'[41m'
   character(len=*), parameter :: c_Igreen     =achar(27)//'[42m'
   character(len=*), parameter :: c_Iyellow    =achar(27)//'[43m'
   character(len=*), parameter :: c_Iblue      =achar(27)//'[44m'
   character(len=*), parameter :: c_IMagenta   =achar(27)//'[45m'
   character(len=*), parameter :: c_Icyan      =achar(27)//'[46m'
   character(len=*), parameter :: c_underline  =achar(27)//'[4m'
   character(len=*), parameter :: c_strike     =achar(27)//'[9m'
   character(len=*), parameter :: c_bold       =achar(27)//'[1m'
   character(len=*), parameter :: c_normal     =achar(27)//'[0m'
   character(len=*), parameter :: c_blinking   =achar(27)//'[5m'
   character(len=*), parameter :: c_reverse    =achar(27)//'[7m'
   character(len=*), parameter :: c_inverted   =achar(27)//'[30m'//achar(27)//'[47m'
   character(len=*), parameter :: c_bg_black   =achar(27)//'[40m'
   character(len=*), parameter :: c_bg_red     =achar(27)//'[41m'
   character(len=*), parameter :: c_bg_green   =achar(27)//'[42m'
   character(len=*), parameter :: c_bg_brown   =achar(27)//'[43m'
   character(len=*), parameter :: c_bg_blue    =achar(27)//'[44m'
   character(len=*), parameter :: c_bg_purple  =achar(27)//'[45m'
   character(len=*), parameter :: c_bg_cyan    =achar(27)//'[46m'
   character(len=*), parameter :: c_bg_lgray   =achar(27)//'[47m' 


   integer, parameter :: c_tty=6
   !# Default TTY (terminal) output
!#ifdef splitlog
!   integer, parameter :: logUnit=69
!#else
   integer, parameter :: logUnit=6
!#endif
   !# Number of file to log
   integer :: iErrNumber
   !# integer to use with dumps

   real(kind=r4), parameter :: c_adjust=0.01
   !# to convert from Pascal to mbar

   character(len=3), parameter, dimension(12) :: month_name=(/'jan','feb','mar' &
                                                            , 'apr','may','jun' &
                                                            , 'jul','aug','sep' &
                                                            , 'oct','nov','dec'/)

   character(len=3), parameter, dimension(7)  ::  week_name=(/'sun','mon','tue' &
                                                            , 'wed','thu','fri' &
                                                            , 'sat'/)



contains


end module modConstants

