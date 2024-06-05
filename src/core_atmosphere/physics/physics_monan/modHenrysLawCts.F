module modHenrysLawCts
   !! ## Henry's Law Constants
   !!
   !! ![](https://i.ibb.co/LNqGy3S/logo-Monan-Color-75x75.png)
   !! ## MONAN
   !!
   !! Author: Saulo Freitas [SRF]
   !!
   !! E-mail: <mailto:saulo.freitas@inpe.br>
   !!
   !! Date: 09Fevereiro2023 18:01
   !!
   !! #####Version: version
   !!
   !! ---
   !! **Full description**:
   !!
   !! Henry Law Constants
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
   character(len=*), parameter :: sourceName = 'modHenrysLawCts.F90' 
   !! Source code name 
   character(len=*), parameter :: moduleName = 'modHenrysLawCts' 
   !! module name 

   integer, parameter :: p_nspecies_hl = 051
   real, parameter :: p_not_found = -1.

   private
   public :: getHenrysLawCts

   !--- : ak0(ispc), dak(ispc),  hstar(ispc), dhr(ispc)
   !--- corrh=1.+ak0(ispc)*exp(dak(ispc)*tcorr)/hplus
   !--- hplus = 1.175E-4  - for cloud water. pH is assumed to be 3.93: pH=3.93 =>hplus=10**(-pH)
   !--- tcorr = 1./temp - 1./298.15
   !--- fct   = 1.e-3 * rgas * temp
   !--- henry_coef =  hstar(ispc)* exp(dhr(ispc)*tcorr) * fct * corrh

   type t_hcts_vars
      real :: hstar, dhr, ak0, dak
   end type t_hcts_vars
   type(t_hcts_vars) :: hcts(p_nspecies_hl)

   !Name of species
   character(LEN=8), parameter, dimension(p_nspecies_hl) :: spc_name = (/ &
                                                            'O3  ' & !001
                                                          , 'H2O2' & !002
                                                          , 'NO  ' & !003
                                                          , 'NO2 ' & !004
                                                          , 'NO3 ' & !005
                                                          , 'N2O5' & !006
                                                          , 'HONO' & !007
                                                          , 'HNO3' & !008
                                                          , 'HNO4' & !009
                                                          , 'SO2 ' & !010
                                                          , 'SULF' & !011
                                                          , 'CO  ' & !012
                                                          , 'CO2 ' & !013
                                                          , 'N2  ' & !014
                                                          , 'O2  ' & !015
                                                          , 'H2O ' & !016
                                                          , 'H2  ' & !017
                                                          , 'O3P ' & !018
                                                          , 'O1D ' & !019
                                                          , 'HO  ' & !020
                                                          , 'HO2 ' & !021
                                                          , 'CH4 ' & !022
                                                          , 'ETH ' & !023
                                                          , 'ALKA' & !024
                                                          , 'ALKE' & !025
                                                          , 'BIO ' & !026
                                                          , 'ARO ' & !027
                                                          , 'HCHO' & !028
                                                          , 'ALD ' & !029
                                                          , 'KET ' & !030
                                                          , 'CRBO' & !031
                                                          , 'ONIT' & !032
                                                          , 'PAN ' & !033
                                                          , 'OP1 ' & !034
                                                          , 'OP2 ' & !035
                                                          , 'ORA1' & !036
                                                          , 'ORA2' & !037
                                                          , 'MO2 ' & !038
                                                          , 'AKAP' & !039
                                                          , 'AKEP' & !040
                                                          , 'BIOP' & !041
                                                          , 'PHO ' & !042
                                                          , 'ADD ' & !043
                                                          , 'AROP' & !044
                                                          , 'CBOP' & !045
                                                          , 'OLN ' & !046
                                                          , 'XO2 ' & !047
                                                          , 'DMS ' & !048
                                                          , 'NH3 ' & !049
                                                          , 'CFC ' & !050
                                                          , 'N2O ' & !050
                                                          /)

   !Number of each specie
   integer, parameter :: O3 = 001
   integer, parameter :: H2O2 = 002
   integer, parameter :: NO = 003
   integer, parameter :: NO2 = 004
   integer, parameter :: NO3 = 005
   integer, parameter :: N2O5 = 006
   integer, parameter :: HONO = 007
   integer, parameter :: HNO3 = 008
   integer, parameter :: HNO4 = 009
   integer, parameter :: SO2 = 010
   integer, parameter :: SULF = 011
   integer, parameter :: CO = 012
   integer, parameter :: CO2 = 013
   integer, parameter :: N2 = 014
   integer, parameter :: O2 = 015
   integer, parameter :: H2O = 016
   integer, parameter :: H2 = 017
   integer, parameter :: O3P = 018
   integer, parameter :: O1D = 019
   integer, parameter :: HO = 020
   integer, parameter :: HO2 = 021
   integer, parameter :: CH4 = 022
   integer, parameter :: ETH = 023
   integer, parameter :: ALKA = 024
   integer, parameter :: ALKE = 025
   integer, parameter :: BIO = 026
   integer, parameter :: ARO = 027
   integer, parameter :: HCHO = 028
   integer, parameter :: ALD = 029
   integer, parameter :: KET = 030
   integer, parameter :: CRBO = 031
   integer, parameter :: ONIT = 032
   integer, parameter :: PAN = 033
   integer, parameter :: OP1 = 034
   integer, parameter :: OP2 = 035
   integer, parameter :: ORA1 = 036
   integer, parameter :: ORA2 = 037
   integer, parameter :: MO2 = 038
   integer, parameter :: AKAP = 039
   integer, parameter :: AKEP = 040
   integer, parameter :: BIOP = 041
   integer, parameter :: PHO = 042
   integer, parameter :: ADD = 043
   integer, parameter :: AROP = 044
   integer, parameter :: CBOP = 045
   integer, parameter :: OLN = 046
   integer, parameter :: XO2 = 047
   integer, parameter :: DMS = 048
   integer, parameter :: NH3 = 049
   integer, parameter :: CFC = 050
   integer, parameter :: N2O = 051

!     HENRYS LAW COEFFICIENTS
!     Henrys law coefficient
!     [KH298]=mole/(l atm)
!     Referencias em R. Sander (1999)
!     Compilation of Henry Law Constants
!     for Inorganic and Organic Species
!     of Potential Importance in
!     Environmental Chemistry (Version 3)
!     http://www.henrys-law.org
!     * indica artigos nao encontrados nesse endere�o eletronico
   real, parameter, dimension(p_nspecies_hl) :: p_hstar = (/ &
                                              1.10e-2, & ! O3 - 001
                                              8.30e+4, & ! H2O2 - 002
                                              1.90e-3, & ! NO - 003
                                              1.20e-2, & ! NO2 - 004
                                              6.1e-01, & ! NO3 - 005
                                              2.1e+00, & ! N2O5 - 006
                                              5.00e+1, & ! HONO - 007
                                              2.10e+5, & ! HNO3 - 008
                                              1.20e+4, & ! HNO4 - 009
                                              1.40e+0, & ! SO2 - 010
                                              2.10e+5, & ! SULF - 011
                                              9.90e-4, & ! CO - 012
                                              3.6e-02, & ! CO2 - 013
                                              6.1e-04, & ! N2 - 014
                                              1.3e-03, & ! O2 - 015
                                              0.0e+00, & ! H2O - 016
                                              7.8e-04, & ! H2 - 017
                                              0.00e+0, & ! O3P - 018
                                              0.00e+0, & ! O1D - 019
                                              3.00e+1, & ! HO - 020
                                              5.70e+3, & ! HO2 - 021
                                              1.40e-3, & ! CH4 - 022
                                              1.90e-3, & ! ETH - 023
                                              1.00e-3, & ! ALKA - 024
                                              5.00e-3, & ! ALKE - 025
                                              2.80e-2, & ! BIO - 026
                                              1.73e-1, & ! ARO - 027
                                              3.20e+3, & ! HCHO - 028
                                              1.40e+1, & ! ALD - 029
                                              3.00e+1, & ! KET - 030
                                              2.1e+05, & ! CRBO - 031
                                              1.00e+0, & ! ONIT - 032
                                              3.60e+0, & ! PAN - 033
                                              3.10e+2, & ! OP1 - 034
                                              3.40e+2, & ! OP2 - 035
                                              8.90e+3, & ! ORA1 - 036
                                              4.10e+3, & ! ORA2 - 037
                                              2.00e+3, & ! MO2 - 038
                                              0.0e+00, & ! AKAP - 039
                                              0.0e+00, & ! AKEP - 040
                                              0.0e+00, & ! BIOP - 041
                                              0.0e+00, & ! PHO - 042
                                              0.0e+00, & ! ADD - 043
                                              0.0e+00, & ! AROP - 044
                                              1.14e+1, & ! CBOP - 045
                                              0.0e+00, & ! OLN - 046
                                              0.0e+00, & ! XO2 - 047
                                              5.6e-01, & ! DMS - 048
                                              5.9e+01, & ! NH3 - 048
                                              -1., & ! CFC - 048
                                              2.4e-02 & ! N2O - 051
                                              /)

!     -DH/R (for temperature correction)
!     [-DH/R]=K
!     Referencias em R. Sander (1999)
!     Compilation of Henry Law Constants
!     for Inorganic and Organic Species
!     of Potential Importance in
!     Environmental Chemistry (Version 3)
!     http://www.henrys-law.org
   real, parameter, dimension(p_nspecies_hl) :: p_dhr = (/ &
                                              2400., & ! O3 - 001
                                              7400., & ! H2O2 - 002
                                              1400., & ! NO - 003
                                              2500., & ! NO2 - 004
                                              2000., & ! NO3 - 005
                                              3400., & ! N2O5 - 006
                                              4900., & ! HONO - 007
                                              8700., & ! HNO3 - 008
                                              6900., & ! HNO4 - 009
                                              2900., & ! SO2 - 010
                                              0., & ! SULF - 011
                                              1300., & ! CO - 012
                                              2200., & ! CO2 - 013
                                              1300., & ! N2 - 014
                                              1500., & ! O2 - 015
                                              0., & ! H2O - 016
                                              500., & ! H2 - 017
                                              0., & ! O3P - 018
                                              0., & ! O1D - 019
                                              4500., & ! HO - 020
                                              5900., & ! HO2 - 021
                                              1600., & ! CH4 - 022
                                              2300., & ! ETH - 023
                                              2700., & ! ALKA - 024
                                              3000., & ! ALKE - 025
                                              0., & ! BIO - 026
                                              4045., & ! ARO - 027
                                              6800., & ! HCHO - 028
                                              5600., & ! ALD - 029
                                              4600., & ! KET - 030
                                              5300., & ! CRBO - 031
                                              5800., & ! ONIT - 032
                                              6500., & ! PAN - 033
                                              5200., & ! OP1 - 034
                                              6000., & ! OP2 - 035
                                              5700., & ! ORA1 - 036
                                              6300., & ! ORA2 - 037
                                              6600., & ! MO2 - 038
                                              0., & ! AKAP - 039
                                              0., & ! AKEP - 040
                                              0., & ! BIOP - 041
                                              0., & ! PHO - 042
                                              0., & ! ADD - 043
                                              0., & ! AROP - 044
                                              0., & ! CBOP - 045
                                              0., & ! OLN - 046
                                              0., & ! XO2 - 047
                                              3500., & ! DMS - 048
                                              4200., & ! NH3 - 048
                                              -1., & ! CFC - 048
                                              2700. & ! N2O - 048
                                              /)

   real, parameter, dimension(p_nspecies_hl) :: p_weight = (/ &
                                              48., & ! O3 - 001
                                              34., & ! H2O2 - 002
                                              30., & ! NO - 003
                                              46., & ! NO2 - 004
                                              62., & ! NO3 - 005
                                              108., & ! N2O5 - 006
                                              47., & ! HONO - 007
                                              63., & ! HNO3 - 008
                                              79., & ! HNO4 - 009
                                              64., & ! SO2 - 010
                                              98., & ! SULF - 011
                                              28., & ! CO - 012
                                              44., & ! CO2 - 013
                                              28., & ! N2 - 014
                                              32., & ! O2 - 015
                                              18., & ! H2O - 016
                                              2., & ! H2 - 017
                                              16., & ! O3P - 018
                                              16., & ! O1D - 019
                                              17., & ! HO - 020
                                              33., & ! HO2 - 021
                                              16., & ! CH4 - 022
                                              30., & ! ETH - 023
                                              61.6, & ! ALKA - 024
                                              33.0, & ! ALKE - 025
                                              68., & ! BIO - 026
                                              97.9, & ! ARO - 027
                                              30., & ! HCHO - 028
                                              44., & ! ALD - 029
                                              72., & ! KET - 030
                                              68.6, & ! CRBO - 031
                                              119., & ! ONIT - 032
                                              122., & ! PAN - 033
                                              48., & ! OP1 - 034
                                              62., & ! OP2 - 035
                                              46., & ! ORA1 - 036
                                              60., & ! ORA2 - 037
                                              47., & ! MO2 - 038
                                              102., & ! AKAP - 039
                                              88.4, & ! AKEP - 040
                                              117., & ! BIOP - 041
                                              107., & ! PHO - 042
                                              107., & ! ADD - 043
                                              151., & ! AROP - 044
                                              85.4, & ! CBOP - 045
                                              136., & ! OLN - 046
                                              44., & ! XO2 - 047
                                              62.13, & ! DMS - 048
                                              17.03, & ! NH3 - 048
                                              -1., & ! CFC - 048
                                              44. & ! CFC - 048
                                              /)

!    ACID DISSOCIATION CONSTANT AT 298K
!     [mole/liter of liquid water]
!     Referencias: Barth et al. JGR 112, D13310 2007
!     Martell and Smith, 1976, Critical stability
!     vol1-4 Plenum Press New York
   real, parameter, dimension(p_nspecies_hl) :: p_ak0 = (/ &
                                              0.00e+00, & ! O3 - 001
                                              2.20e-12, & ! H2O2 - 002
                                              0.00e+00, & ! NO - 003
                                              0.00e+00, & ! NO2 - 004
                                              0.00e+00, & ! NO3 - 005
                                              0.00e+00, & ! N2O5 - 006
                                              7.10e-04, & ! HONO - 007
                                              1.54e+01, & ! HNO3 - 008
                                              0.00e+00, & ! HNO4 - 009
                                              1.30e-02, & ! SO2 - 010
                                              1.00e-02, & ! SULF - 011
                                              0.00e+00, & ! CO - 012
                                              4.50e-07, & ! CO2 - 013
                                              0.00e+00, & ! N2 - 014
                                              0.00e+00, & ! O2 - 015
                                              0.00e+00, & ! H2O - 016
                                              0.00e+00, & ! H2 - 017
                                              0.00e+00, & ! O3P - 018
                                              0.00e+00, & ! O1D - 019
                                              0.00e+00, & ! HO - 020
                                              3.50e-05, & ! HO2 - 021
                                              0.00e+00, & ! CH4 - 022
                                              0.00e+00, & ! ETH - 023
                                              0.00e+00, & ! ALKA - 024
                                              0.00e+00, & ! ALKE - 025
                                              0.00e+00, & ! BIO - 026
                                              0.00e+00, & ! ARO - 027
                                              0.00e+00, & ! HCHO - 028
                                              0.00e+00, & ! ALD - 029
                                              0.00e+00, & ! KET - 030
                                              0.00e+00, & ! CRBO - 031
                                              0.00e+00, & ! ONIT - 032
                                              0.00e+00, & ! PAN - 033
                                              0.00e+00, & ! OP1 - 034
                                              0.00e+00, & ! OP2 - 035
                                              1.80e-04, & ! ORA1 - 036
                                              1.75e-05, & ! ORA2 - 037
                                              0.00e+00, & ! MO2 - 038
                                              0.00e+00, & ! AKAP - 039
                                              0.00e+00, & ! AKEP - 040
                                              0.00e+00, & ! BIOP - 041
                                              0.00e+00, & ! PHO - 042
                                              0.00e+00, & ! ADD - 043
                                              0.00e+00, & ! AROP - 044
                                              0.00e+00, & ! CBOP - 045
                                              0.00e+00, & ! OLN - 046
                                              0.00e+00, & ! XO2 - 047
                                              0.00e+00, & ! DMS - 048
                                              0.00e+00, & ! NH3 - 049
                                              0.00e+00, & ! NH3 - 049
                                              0.00e+00 & ! CFC - 050
                                              /)

!     Temperature correction factor for
!     acid dissociation constants
!     [K]
!     Referencias: Barth et al. JGR 112, D13310 2007
   real, parameter, dimension(p_nspecies_hl) :: p_dak = (/ &
                                              0., & ! O3 - 001
                                              -3700., & ! H2O2 - 002
                                              0., & ! NO - 003
                                              0., & ! NO2 - 004
                                              0., & ! NO3 - 005
                                              0., & ! N2O5 - 006
                                              0., & ! HONO - 007
                                              0., & ! HNO3 - 008
                                              0., & ! HNO4 - 009
                                              2000., & ! SO2 - 010
                                              0., & ! SULF - 011
                                              0., & ! CO - 012
                                              -1000., & ! CO2 - 013
                                              0., & ! N2 - 014
                                              0., & ! O2 - 015
                                              0., & ! H2O - 016
                                              0., & ! H2 - 017
                                              0., & ! O3P - 018
                                              0., & ! O1D - 019
                                              0., & ! HO - 020
                                              0., & ! HO2 - 021
                                              0., & ! CH4 - 022
                                              0., & ! ETH - 023
                                              0., & ! ALKA - 024
                                              0., & ! ALKE - 025
                                              0., & ! BIO - 026
                                              0., & ! ARO - 027
                                              0., & ! HCHO - 028
                                              0., & ! ALD - 029
                                              0., & ! KET - 030
                                              0., & ! CRBO - 031
                                              0., & ! ONIT - 032
                                              0., & ! PAN - 033
                                              0., & ! OP1 - 034
                                              0., & ! OP2 - 035
                                              -1500., & ! ORA1 - 036
                                              0., & ! ORA2 - 037
                                              0., & ! MO2 - 038
                                              0., & ! AKAP - 039
                                              0., & ! AKEP - 040
                                              0., & ! BIOP - 041
                                              0., & ! PHO - 042
                                              0., & ! ADD - 043
                                              0., & ! AROP - 044
                                              0., & ! CBOP - 045
                                              0., & ! OLN - 046
                                              0., & ! XO2 - 047
                                              0., & ! DMS - 048
                                              0., & ! NH3 - 049
                                              0., & ! NH3 - 049
                                              0. & ! CFC - 050
                                              /)

contains

   ! -----------------------------------------------------------------------
   subroutine getHenrysLawCts(name, c1, c2, c3, c4)
      !! ## Get the constants
      !!
      !! Author: autor
      !!
      !! E-mail: <mailto:email>
      !!
      !! Date: 09Fevereiro2023 18:06
      !!
      !! #####Version: version
      !!
      !! ---
      !! **Full description**:
      !!
      !! Get the constants
      !!
      !! ** History**:
      !!
      !! --- 
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!
   
      implicit none
      !Parameters:
      character(len=*), parameter :: procedureName = 'getHenrysLawCts' 
      !! subroutine name
   
      !Variables (input, output, inout)
      character(len=*), intent(in) :: name

      real, intent(out) :: c1
      real, intent(out) :: c2
      real, intent(out) :: c3
      real, intent(out) :: c4

      !Local variables:
      integer :: l, found
      
      found = 0
      loop2: do l = 1, p_nspecies_hl
         if (trim(spc_name(l)) == trim(name)) then
            c1 = p_hstar(l)
            c2 = p_dhr(l)
            c3 = p_ak0(l)
            c4 = p_dak(l)
            found = 1
            exit loop2
         end if
      end do loop2
      if (found == 0) then
         c1 = p_not_found
         c2 = p_not_found
         c3 = p_not_found
         c4 = p_not_found
      end if

   end subroutine getHenrysLawCts

end module modHenrysLawCts
