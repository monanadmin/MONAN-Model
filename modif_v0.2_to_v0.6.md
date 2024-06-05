git clone https://github.com/saulorfreitas/MONAN_Phys_SRF_v1
git checkout develop



Modifications:
1) included physics_monan directory


2) update Makefile for gfortran
3) update Registry.xml of init_atmosphere
4) update Registry.xml of core_atmosphere
5) update Makefile of core_atmosphere


6) Changes in physics
  a) mpas_atmphys_control.F
  b) mpas_atmphys_driver.F
  c) mpas_atmphys_driver_cloudiness.F
  d) mpas_atmphys_driver_convection.F
  e) mpas_atmphys_init.F
  f) mpas_atmphys_interface.F
  g) mpas_atmphys_packages.F
  h) mpas_atmphys_todynamics.F
  i) mpas_atmphys_vars.F

7) changes in dynamics
  a) mpas_atm_time_integration.F

8) changes in core_atmosphere
  a) mpas_atm_halos.F

