MODULE icethd_da
  USE par_oce
  USE phycst
  USE ice
  USE ice1D
  USE in_out_manager
  USE lib_mpp
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_thd_da
  PUBLIC :: ice_thd_da_init
  REAL(KIND = wp) :: rn_beta
  REAL(KIND = wp) :: rn_dmin
  CONTAINS
  SUBROUTINE ice_thd_da
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ji
    REAL(KIND = wp) :: zastar, zdfloe, zperi, zwlat, zda
    REAL(KIND = wp), PARAMETER :: zdmax = 300._wp
    REAL(KIND = wp), PARAMETER :: zcs = 0.66_wp
    REAL(KIND = wp), PARAMETER :: zm1 = 3.E-6_wp
    REAL(KIND = wp), PARAMETER :: zm2 = 1.36_wp
    REAL(KIND = wp), DIMENSION(jpij) :: zda_tot
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ice_thd_da', 'r0', 0, 0)
    zastar = 1._wp / (1._wp - (rn_dmin / zdmax) ** (1._wp / rn_beta))
    DO ji = 1, npti
      zdfloe = rn_dmin * (zastar / (zastar - at_i_1d(ji))) ** rn_beta
      zperi = at_i_1d(ji) * rpi / (zcs * zdfloe)
      zwlat = zm1 * (MAX(0._wp, sst_1d(ji) - (t_bo_1d(ji) - rt0))) ** zm2
      zda_tot(ji) = MIN(zwlat * zperi * rdt_ice, at_i_1d(ji))
      IF (a_i_1d(ji) > 0._wp) THEN
        zda = MIN(a_i_1d(ji), zda_tot(ji) * a_i_1d(ji) / at_i_1d(ji))
        sfx_lam_1d(ji) = sfx_lam_1d(ji) + rhoi * h_i_1d(ji) * zda * s_i_1d(ji) * r1_rdtice
        hfx_thd_1d(ji) = hfx_thd_1d(ji) - zda * r1_rdtice * (h_i_1d(ji) * r1_nlay_i * SUM(e_i_1d(ji, 1 : nlay_i)) + h_s_1d(ji) * r1_nlay_s * SUM(e_s_1d(ji, 1 : nlay_s)))
        wfx_lam_1d(ji) = wfx_lam_1d(ji) + zda * r1_rdtice * (rhoi * h_i_1d(ji) + rhos * h_s_1d(ji))
        a_i_1d(ji) = a_i_1d(ji) - zda
        IF (a_i_1d(ji) == 0._wp) THEN
          h_i_1d(ji) = 0._wp
          h_s_1d(ji) = 0._wp
        END IF
      END IF
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ice_thd_da
  SUBROUTINE ice_thd_da_init
    INTEGER :: ios
    NAMELIST /namthd_da/ rn_beta, rn_dmin
    REWIND(UNIT = numnam_ice_ref)
    READ(numnam_ice_ref, namthd_da, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namthd_da in reference namelist', lwp)
    REWIND(UNIT = numnam_ice_cfg)
    READ(numnam_ice_cfg, namthd_da, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namthd_da in configuration namelist', lwp)
    IF (lwm) WRITE(numoni, namthd_da)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'ice_thd_da_init: Ice lateral melting'
      WRITE(numout, *) '~~~~~~~~~~~~~~~'
      WRITE(numout, *) '   Namelist namthd_da:'
      WRITE(numout, *) '      Coef. beta for lateral melting param.               rn_beta = ', rn_beta
      WRITE(numout, *) '      Minimum floe diameter for lateral melting param.    rn_dmin = ', rn_dmin
    END IF
  END SUBROUTINE ice_thd_da_init
END MODULE icethd_da