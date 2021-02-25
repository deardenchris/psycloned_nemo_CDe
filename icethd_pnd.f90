MODULE icethd_pnd
  USE phycst
  USE dom_oce
  USE ice
  USE ice1D
  USE icetab
  USE in_out_manager
  USE lib_mpp
  USE lib_fortran
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_thd_pnd_init
  PUBLIC :: ice_thd_pnd
  INTEGER :: nice_pnd
  INTEGER, PARAMETER :: np_pndNO = 0
  INTEGER, PARAMETER :: np_pndCST = 1
  INTEGER, PARAMETER :: np_pndH12 = 2
  CONTAINS
  SUBROUTINE ice_thd_pnd
    SELECT CASE (nice_pnd)
    CASE (np_pndcst)
      CALL pnd_cst
    CASE (np_pndh12)
      CALL pnd_h12
    END SELECT
  END SUBROUTINE ice_thd_pnd
  SUBROUTINE pnd_CST
    INTEGER :: ji
    DO ji = 1, npti
      IF (a_i_1d(ji) > 0._wp .AND. t_su_1d(ji) >= rt0) THEN
        a_ip_frac_1d(ji) = rn_apnd
        h_ip_1d(ji) = rn_hpnd
        a_ip_1d(ji) = a_ip_frac_1d(ji) * a_i_1d(ji)
      ELSE
        a_ip_frac_1d(ji) = 0._wp
        h_ip_1d(ji) = 0._wp
        a_ip_1d(ji) = 0._wp
      END IF
    END DO
  END SUBROUTINE pnd_CST
  SUBROUTINE pnd_H12
    REAL(KIND = wp), PARAMETER :: zrmin = 0.15_wp
    REAL(KIND = wp), PARAMETER :: zrmax = 0.70_wp
    REAL(KIND = wp), PARAMETER :: zpnd_aspect = 0.8_wp
    REAL(KIND = wp), PARAMETER :: zTp = - 2._wp
    REAL(KIND = wp) :: zfr_mlt
    REAL(KIND = wp) :: zdv_mlt
    REAL(KIND = wp) :: z1_Tp
    REAL(KIND = wp) :: z1_rhow
    REAL(KIND = wp) :: z1_zpnd_aspect
    REAL(KIND = wp) :: zfac, zdum
    INTEGER :: ji
    z1_rhow = 1._wp / rhow
    z1_zpnd_aspect = 1._wp / zpnd_aspect
    z1_Tp = 1._wp / zTp
    DO ji = 1, npti
      IF (h_i_1d(ji) < rn_himin) THEN
        a_ip_1d(ji) = 0._wp
        a_ip_frac_1d(ji) = 0._wp
        h_ip_1d(ji) = 0._wp
      ELSE
        v_ip_1d(ji) = h_ip_1d(ji) * a_ip_1d(ji)
        zdv_mlt = - (dh_i_sum(ji) * rhoi + dh_s_mlt(ji) * rhos) * z1_rhow * a_i_1d(ji)
        zfr_mlt = zrmin + (zrmax - zrmin) * a_i_1d(ji)
        v_ip_1d(ji) = MAX(0._wp, v_ip_1d(ji) + zfr_mlt * zdv_mlt)
        IF (zdv_mlt > 0._wp) THEN
          zfac = zfr_mlt * zdv_mlt * rhow * r1_rdtice
          wfx_pnd_1d(ji) = wfx_pnd_1d(ji) - zfac
          zdum = zfac / (wfx_snw_sum_1d(ji) + wfx_sum_1d(ji))
          wfx_snw_sum_1d(ji) = wfx_snw_sum_1d(ji) * (1._wp + zdum)
          wfx_sum_1d(ji) = wfx_sum_1d(ji) * (1._wp + zdum)
        END IF
        v_ip_1d(ji) = v_ip_1d(ji) * EXP(0.01_wp * MAX(zTp + rt0 - t_su_1d(ji), 0._wp) * z1_Tp)
        a_ip_1d(ji) = SQRT(v_ip_1d(ji) * z1_zpnd_aspect * a_i_1d(ji))
        a_ip_frac_1d(ji) = a_ip_1d(ji) / a_i_1d(ji)
        h_ip_1d(ji) = zpnd_aspect * a_ip_frac_1d(ji)
      END IF
    END DO
  END SUBROUTINE pnd_H12
  SUBROUTINE ice_thd_pnd_init
    INTEGER :: ios, ioptio
    NAMELIST /namthd_pnd/ ln_pnd_H12, ln_pnd_CST, rn_apnd, rn_hpnd, ln_pnd_alb
    REWIND(UNIT = numnam_ice_ref)
    READ(numnam_ice_ref, namthd_pnd, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namthd_pnd  in reference namelist', lwp)
    REWIND(UNIT = numnam_ice_cfg)
    READ(numnam_ice_cfg, namthd_pnd, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namthd_pnd in configuration namelist', lwp)
    IF (lwm) WRITE(numoni, namthd_pnd)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'ice_thd_pnd_init: ice parameters for melt ponds'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namicethd_pnd:'
      WRITE(numout, FMT = *) '      Evolutive  melt pond fraction and depth (Holland et al 2012) ln_pnd_H12 = ', ln_pnd_H12
      WRITE(numout, FMT = *) '      Prescribed melt pond fraction and depth                      ln_pnd_CST = ', ln_pnd_CST
      WRITE(numout, FMT = *) '         Prescribed pond fraction                                  rn_apnd    = ', rn_apnd
      WRITE(numout, FMT = *) '         Prescribed pond depth                                     rn_hpnd    = ', rn_hpnd
      WRITE(numout, FMT = *) '      Melt ponds affect albedo or not                              ln_pnd_alb = ', ln_pnd_alb
    END IF
    ioptio = 0
    nice_pnd = np_pndNO
    IF (ln_pnd_cst) THEN
      ioptio = ioptio + 1
      nice_pnd = np_pndcst
    END IF
    IF (ln_pnd_h12) THEN
      ioptio = ioptio + 1
      nice_pnd = np_pndh12
    END IF
    IF (ioptio > 1) CALL ctl_stop('ice_thd_pnd_init: choose one and only one pond scheme (ln_pnd_H12 or ln_pnd_CST)')
    SELECT CASE (nice_pnd)
    CASE (np_pndNO)
      IF (ln_pnd_alb) THEN
        ln_pnd_alb = .FALSE.
        CALL ctl_warn('ln_pnd_alb=false when no ponds')
      END IF
    END SELECT
  END SUBROUTINE ice_thd_pnd_init
END MODULE icethd_pnd