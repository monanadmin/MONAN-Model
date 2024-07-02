module modGate
   !! ## brief
   !!
   !! ![](https://i.ibb.co/LNqGy3S/logo-Monan-Color-75x75.png)
   !! ## MONAN
   !!
   !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
   !!
   !! E-mail: Saulo Freitas [SRF] e Georg Grell [GAG]
   !!
   !! Date: 2014
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

   character(len=*), parameter :: sourceName = 'module_gate.f90' 
   !! Source code name 
   character(len=*), parameter :: moduleName = 'module_name' 
   !! module name 

   logical :: p_use_gate = .true.
   !! -for BRAMS runs, set use_gate=.false

   !- Here are the place for data related with the GATE soundings
   integer, parameter :: p_gate = 1   
   !! flag to turn on/off : 1/0
   integer, parameter :: p_klon = 161 
   !! number of soundings for gate
   integer, parameter :: p_klev = 41  
   !! number of vertical levels
   integer, parameter :: p_ktrac = 2   
   !! number of chemical tracers
   integer, parameter :: p_levs = p_klev
   integer, parameter :: p_nvar_grads = 300

   type t_cupout_vars
      real, pointer :: varp(:, :)
      character(len=80), allocatable :: varn(:)
   end type t_cupout_vars

   type(t_cupout_vars), allocatable :: cupout(:)

   real, dimension(p_klon, p_klev):: pgeo, ppres, ptemp, pq, pu, pv, pvervel &
                                  ,  zrvten, ztten, zq1, zq2, zqr, zadvt, zadvq

   integer :: jl, klev_sound
   character(len=128) :: runname, runlabel, rundata = "NONE"

end module modGate
