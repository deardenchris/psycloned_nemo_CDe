MODULE icethd_sal
  USE dom_oce
  USE phycst
  USE ice
  USE ice1D
  USE icevar
  USE in_out_manager
  USE lib_mpp
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_thd_sal
  PUBLIC :: ice_thd_sal_init
  REAL(KIND = wp) :: rn_sal_gd
  REAL(KIND = wp) :: rn_time_gd
  REAL(KIND = wp) :: rn_sal_fl
  REAL(KIND = wp) :: rn_time_fl
  CONTAINS
  SUBROUTINE ice_thd_sal(ld_sal)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    LOGICAL, INTENT(IN) :: ld_sal
    INTEGER :: ji, jk
    REAL(KIND = wp) :: iflush, igravdr
    REAL(KIND = wp) :: zs_sni, zs_i_gd, zs_i_fl, zs_i_si, zs_i_bg
    REAL(KIND = wp) :: z1_time_gd, z1_time_fl
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ice_thd_sal', 'r0', 0, 0)
    SELECT CASE (nn_icesal)
    CASE (2)
      z1_time_gd = 1._wp / rn_time_gd * rdt_ice
      z1_time_fl = 1._wp / rn_time_fl * rdt_ice
      DO ji = 1, npti
        IF (h_i_1d(ji) > 0._wp) THEN
          zs_sni = sss_1d(ji) * (rhoi - rhos) * r1_rhoi
          zs_i_si = (zs_sni - s_i_1d(ji)) * dh_snowice(ji) / h_i_1d(ji)
          zs_i_bg = (s_i_new(ji) - s_i_1d(ji)) * dh_i_bog(ji) / h_i_1d(ji)
          s_i_1d(ji) = s_i_1d(ji) + zs_i_bg + zs_i_si
        END IF
        IF (ld_sal) THEN
          iflush = MAX(0._wp, SIGN(1._wp, t_su_1d(ji) - rt0))
          igravdr = MAX(0._wp, SIGN(1._wp, t_bo_1d(ji) - t_su_1d(ji)))
          zs_i_gd = - igravdr * MAX(s_i_1d(ji) - rn_sal_gd, 0._wp) * z1_time_gd
          zs_i_fl = - iflush * MAX(s_i_1d(ji) - rn_sal_fl, 0._wp) * z1_time_fl
          s_i_1d(ji) = s_i_1d(ji) + zs_i_fl + zs_i_gd
          sfx_bri_1d(ji) = sfx_bri_1d(ji) - rhoi * a_i_1d(ji) * h_i_1d(ji) * (zs_i_fl + zs_i_gd) * r1_rdtice
        END IF
      END DO
      CALL ice_var_salprof1d
    CASE (3)
      CALL ice_var_salprof1d
    END SELECT
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ice_thd_sal
  SUBROUTINE ice_thd_sal_init
    INTEGER :: ios
    NAMELIST /namthd_sal/ nn_icesal, rn_icesal, rn_sal_gd, rn_time_gd, rn_sal_fl, rn_time_fl, rn_simax, rn_simin
    REWIND(UNIT = numnam_ice_ref)
    READ(numnam_ice_ref, namthd_sal, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namthd_sal in reference namelist', lwp)
    REWIND(UNIT = numnam_ice_cfg)
    READ(numnam_ice_cfg, namthd_sal, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namthd_sal in configuration namelist', lwp)
    IF (lwm) WRITE(numoni, namthd_sal)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'ice_thd_sal_init : Ice parameters for salinity '
      WRITE(numout, *) '~~~~~~~~~~~~~~~~'
      WRITE(numout, *) '   Namelist namthd_sal:'
      WRITE(numout, *) '      switch for salinity                                     nn_icesal  = ', nn_icesal
      WRITE(numout, *) '      bulk salinity value if nn_icesal = 1                    rn_icesal  = ', rn_icesal
      WRITE(numout, *) '      restoring salinity for gravity drainage                 rn_sal_gd  = ', rn_sal_gd
      WRITE(numout, *) '      restoring time for for gravity drainage                 rn_time_gd = ', rn_time_gd
      WRITE(numout, *) '      restoring salinity for flushing                         rn_sal_fl  = ', rn_sal_fl
      WRITE(numout, *) '      restoring time for flushing                             rn_time_fl = ', rn_time_fl
      WRITE(numout, *) '      Maximum tolerated ice salinity                          rn_simax   = ', rn_simax
      WRITE(numout, *) '      Minimum tolerated ice salinity                          rn_simin   = ', rn_simin
    END IF
  END SUBROUTINE ice_thd_sal_init
END MODULE icethd_sal