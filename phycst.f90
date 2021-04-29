MODULE phycst
  USE par_oce
  USE in_out_manager
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: phy_cst
  REAL(KIND = wp), PUBLIC :: rpi = 3.141592653589793_wp
  REAL(KIND = wp), PUBLIC :: rad = 3.141592653589793_wp / 180._wp
  REAL(KIND = wp), PUBLIC :: rsmall = 0.5 * EPSILON(1.E0)
  REAL(KIND = wp), PUBLIC :: rday = 24. * 60. * 60.
  REAL(KIND = wp), PUBLIC :: rsiyea
  REAL(KIND = wp), PUBLIC :: rsiday
  REAL(KIND = wp), PUBLIC :: raamo = 12._wp
  REAL(KIND = wp), PUBLIC :: rjjhh = 24._wp
  REAL(KIND = wp), PUBLIC :: rhhmm = 60._wp
  REAL(KIND = wp), PUBLIC :: rmmss = 60._wp
  REAL(KIND = wp), PUBLIC :: omega
  REAL(KIND = wp), PUBLIC :: ra = 6371229._wp
  REAL(KIND = wp), PUBLIC :: grav = 9.80665_wp
  REAL(KIND = wp), PUBLIC :: rt0 = 273.15_wp
  REAL(KIND = wp), PUBLIC :: rau0
  REAL(KIND = wp), PUBLIC :: r1_rau0
  REAL(KIND = wp), PUBLIC :: rcp
  REAL(KIND = wp), PUBLIC :: r1_rcp
  REAL(KIND = wp), PUBLIC :: rau0_rcp
  REAL(KIND = wp), PUBLIC :: r1_rau0_rcp
  REAL(KIND = wp), PUBLIC :: emic = 0.97_wp
  REAL(KIND = wp), PUBLIC :: sice = 6.0_wp
  REAL(KIND = wp), PUBLIC :: soce = 34.7_wp
  REAL(KIND = wp), PUBLIC :: rLevap = 2.5E+6_wp
  REAL(KIND = wp), PUBLIC :: vkarmn = 0.4_wp
  REAL(KIND = wp), PUBLIC :: stefan = 5.67E-8_wp
  REAL(KIND = wp), PUBLIC :: rhos = 330._wp
  REAL(KIND = wp), PUBLIC :: rhoi = 917._wp
  REAL(KIND = wp), PUBLIC :: rhow = 1000._wp
  REAL(KIND = wp), PUBLIC :: rcnd_i = 2.034396_wp
  REAL(KIND = wp), PUBLIC :: rcpi = 2067.0_wp
  REAL(KIND = wp), PUBLIC :: rLsub = 2.834E+6_wp
  REAL(KIND = wp), PUBLIC :: rLfus = 0.334E+6_wp
  REAL(KIND = wp), PUBLIC :: rTmlt = 0.054_wp
  REAL(KIND = wp), PUBLIC :: r1_rhoi
  REAL(KIND = wp), PUBLIC :: r1_rhos
  REAL(KIND = wp), PUBLIC :: r1_rcpi
  CONTAINS
  SUBROUTINE phy_cst
    rsiyea = 365.25_wp * rday * 2._wp * rpi / 6.283076_wp
    rsiday = rday / (1._wp + rday / rsiyea)
    omega = 2._wp * rpi / rsiday
    r1_rhoi = 1._wp / rhoi
    r1_rhos = 1._wp / rhos
    r1_rcpi = 1._wp / rcpi
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'phy_cst : initialization of ocean parameters and constants'
      WRITE(numout, *) '~~~~~~~'
      WRITE(numout, *) '      mathematical constant                 rpi = ', rpi
      WRITE(numout, *) '      day                                rday   = ', rday, ' s'
      WRITE(numout, *) '      sideral year                       rsiyea = ', rsiyea, ' s'
      WRITE(numout, *) '      sideral day                        rsiday = ', rsiday, ' s'
      WRITE(numout, *) '      omega                              omega  = ', omega, ' s^-1'
      WRITE(numout, *)
      WRITE(numout, *) '      nb of months per year               raamo = ', raamo, ' months'
      WRITE(numout, *) '      nb of hours per day                 rjjhh = ', rjjhh, ' hours'
      WRITE(numout, *) '      nb of minutes per hour              rhhmm = ', rhhmm, ' mn'
      WRITE(numout, *) '      nb of seconds per minute            rmmss = ', rmmss, ' s'
      WRITE(numout, *)
      WRITE(numout, *) '      earth radius                         ra   = ', ra, ' m'
      WRITE(numout, *) '      gravity                              grav = ', grav, ' m/s^2'
      WRITE(numout, *)
      WRITE(numout, *) '      freezing point of water              rt0  = ', rt0, ' K'
      WRITE(numout, *)
      WRITE(numout, *) '   reference density and heat capacity now defined in eosbn2.f90'
      WRITE(numout, *)
      WRITE(numout, *) '      thermal conductivity of pure ice          = ', rcnd_i, ' J/s/m/K'
      WRITE(numout, *) '      thermal conductivity of snow is defined in a namelist '
      WRITE(numout, *) '      fresh ice specific heat                   = ', rcpi, ' J/kg/K'
      WRITE(numout, *) '      latent heat of fusion of fresh ice / snow = ', rLfus, ' J/kg'
      WRITE(numout, *) '      latent heat of subl.  of fresh ice / snow = ', rLsub, ' J/kg'
      WRITE(numout, *) '      density of sea ice                        = ', rhoi, ' kg/m^3'
      WRITE(numout, *) '      density of snow                           = ', rhos, ' kg/m^3'
      WRITE(numout, *) '      density of freshwater (in melt ponds)     = ', rhow, ' kg/m^3'
      WRITE(numout, *) '      salinity of ice (for pisces)              = ', sice, ' psu'
      WRITE(numout, *) '      salinity of sea (for pisces and isf)      = ', soce, ' psu'
      WRITE(numout, *) '      latent heat of evaporation (water)        = ', rLevap, ' J/m^3'
      WRITE(numout, *) '      von Karman constant                       = ', vkarmn
      WRITE(numout, *) '      Stefan-Boltzmann constant                 = ', stefan, ' J/s/m^2/K^4'
      WRITE(numout, *)
      WRITE(numout, *) '      conversion: degre ==> radian          rad = ', rad
      WRITE(numout, *)
      WRITE(numout, *) '      smallest real computer value       rsmall = ', rsmall
    END IF
  END SUBROUTINE phy_cst
END MODULE phycst