MODULE diatmb
  USE oce
  USE dom_oce
  USE in_out_manager
  USE iom
  USE wet_dry
  IMPLICIT NONE
  PRIVATE
  LOGICAL, PUBLIC :: ln_diatmb
  PUBLIC :: dia_tmb_init
  PUBLIC :: dia_tmb
  CONTAINS
  SUBROUTINE dia_tmb_init
    INTEGER :: ios
    NAMELIST /nam_diatmb/ ln_diatmb
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, nam_diatmb, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'nam_diatmb in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, nam_diatmb, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'nam_diatmb in configuration namelist', lwp)
    IF (lwm) WRITE(numond, nam_diatmb)
    IF (lwp) THEN
      WRITE(numout, *)
      WRITE(numout, *) 'dia_tmb_init : Output Top, Middle, Bottom Diagnostics'
      WRITE(numout, *) '~~~~~~~~~~~~'
      WRITE(numout, *) '   Namelist nam_diatmb : set tmb outputs '
      WRITE(numout, *) '      Switch for TMB diagnostics (T) or not (F)  ln_diatmb  = ', ln_diatmb
    END IF
  END SUBROUTINE dia_tmb_init
  SUBROUTINE dia_calctmb(pfield, ptmb)
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pfield
    REAL(KIND = wp), DIMENSION(jpi, jpj, 3), INTENT(OUT) :: ptmb
    INTEGER :: ji, jj
    INTEGER :: itop, imid, ibot
    REAL(KIND = wp) :: zmdi = 1.E+20_wp
    !$ACC KERNELS
    !$ACC loop independent collapse(2)
    DO jj = 1, jpj
      DO ji = 1, jpi
        itop = mikt(ji, jj)
        ibot = mbkt(ji, jj)
        imid = itop + (ibot - itop + 1) / 2
        ptmb(ji, jj, 1) = pfield(ji, jj, itop) * tmask(ji, jj, itop) + zmdi * (1._wp - tmask(ji, jj, itop))
        ptmb(ji, jj, 2) = pfield(ji, jj, imid) * tmask(ji, jj, imid) + zmdi * (1._wp - tmask(ji, jj, imid))
        ptmb(ji, jj, 3) = pfield(ji, jj, ibot) * tmask(ji, jj, ibot) + zmdi * (1._wp - tmask(ji, jj, ibot))
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE dia_calctmb
  SUBROUTINE dia_tmb
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp) :: zmdi = 1.E+20
    REAL(KIND = wp), DIMENSION(jpi, jpj, 3) :: zwtmb
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('dia_tmb', 'r0', 0, 0)
    CALL dia_calctmb(tsn(:, :, :, jp_tem), zwtmb)
    IF (ll_wd) THEN
      CALL iom_put("sshnmasked", (sshn(:, :) + ssh_ref) * tmask(:, :, 1) + zmdi * (1.0 - tmask(:, :, 1)))
    ELSE
      CALL iom_put("sshnmasked", sshn(:, :) * tmask(:, :, 1) + zmdi * (1.0 - tmask(:, :, 1)))
    END IF
    CALL iom_put("top_temp", zwtmb(:, :, 1))
    CALL iom_put("mid_temp", zwtmb(:, :, 2))
    CALL iom_put("bot_temp", zwtmb(:, :, 3))
    CALL dia_calctmb(tsn(:, :, :, jp_sal), zwtmb)
    CALL iom_put("top_sal", zwtmb(:, :, 1))
    CALL iom_put("mid_sal", zwtmb(:, :, 2))
    CALL iom_put("bot_sal", zwtmb(:, :, 3))
    CALL dia_calctmb(un(:, :, :), zwtmb)
    CALL iom_put("top_u", zwtmb(:, :, 1))
    CALL iom_put("mid_u", zwtmb(:, :, 2))
    CALL iom_put("bot_u", zwtmb(:, :, 3))
    CALL dia_calctmb(vn(:, :, :), zwtmb)
    CALL iom_put("top_v", zwtmb(:, :, 1))
    CALL iom_put("mid_v", zwtmb(:, :, 2))
    CALL iom_put("bot_v", zwtmb(:, :, 3))
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE dia_tmb
END MODULE diatmb