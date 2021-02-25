MODULE zdfiwm
  USE oce
  USE dom_oce
  USE zdf_oce
  USE zdfddm
  USE lbclnk
  USE eosbn2
  USE phycst
  USE prtctl
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: zdf_iwm
  PUBLIC :: zdf_iwm_init
  INTEGER :: nn_zpyc
  LOGICAL :: ln_mevar
  LOGICAL :: ln_tsdiff
  REAL(KIND = wp) :: r1_6 = 1._wp / 6._wp
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: ebot_iwm
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: epyc_iwm
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: ecri_iwm
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: hbot_iwm
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: hcri_iwm
  CONTAINS
  INTEGER FUNCTION zdf_iwm_alloc()
    ALLOCATE(ebot_iwm(jpi, jpj), epyc_iwm(jpi, jpj), ecri_iwm(jpi, jpj), hbot_iwm(jpi, jpj), hcri_iwm(jpi, jpj), STAT = &
&zdf_iwm_alloc)
    CALL mpp_sum('zdfiwm', zdf_iwm_alloc)
    IF (zdf_iwm_alloc /= 0) CALL ctl_stop('STOP', 'zdf_iwm_alloc: failed to allocate arrays')
  END FUNCTION zdf_iwm_alloc
  SUBROUTINE zdf_iwm(kt, p_avm, p_avt, p_avs)
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: p_avm
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: p_avt, p_avs
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zztmp, ztmp1, ztmp2
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zfact
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zhdep
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zwkb
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zweight
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: znu_t
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: znu_w
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zReb
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zemx_iwm
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zav_ratio
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zav_wave
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: z3d
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: z2d
    zemx_iwm(:, :, 1) = 0._wp
    zemx_iwm(:, :, jpk) = 0._wp
    zav_ratio(:, :, 1) = 0._wp
    zav_ratio(:, :, jpk) = 0._wp
    zav_wave(:, :, 1) = 0._wp
    zav_wave(:, :, jpk) = 0._wp
    DO jj = 1, jpj
      DO ji = 1, jpi
        zhdep(ji, jj) = gdepw_0(ji, jj, mbkt(ji, jj) + 1)
        zfact(ji, jj) = rau0 * (1._wp - EXP(- zhdep(ji, jj) / hcri_iwm(ji, jj)))
        IF (zfact(ji, jj) /= 0._wp) zfact(ji, jj) = ecri_iwm(ji, jj) / zfact(ji, jj)
      END DO
    END DO
    DO jk = 2, jpkm1
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (zfact(ji, jj) == 0._wp .OR. wmask(ji, jj, jk) == 0._wp) THEN
            zemx_iwm(ji, jj, jk) = 0._wp
          ELSE
            zemx_iwm(ji, jj, jk) = zfact(ji, jj) * (EXP((gde3w_n(ji, jj, jk) - zhdep(ji, jj)) / hcri_iwm(ji, jj)) - &
&EXP((gde3w_n(ji, jj, jk - 1) - zhdep(ji, jj)) / hcri_iwm(ji, jj))) / (gde3w_n(ji, jj, jk) - gde3w_n(ji, jj, jk - 1))
          END IF
        END DO
      END DO
    END DO
    SELECT CASE (nn_zpyc)
    CASE (1)
      zfact(:, :) = 0._wp
      DO jk = 2, jpkm1
        zfact(:, :) = zfact(:, :) + e3w_n(:, :, jk) * SQRT(MAX(0._wp, rn2(:, :, jk))) * wmask(:, :, jk)
      END DO
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (zfact(ji, jj) /= 0) zfact(ji, jj) = epyc_iwm(ji, jj) / (rau0 * zfact(ji, jj))
        END DO
      END DO
      DO jk = 2, jpkm1
        zemx_iwm(:, :, jk) = zemx_iwm(:, :, jk) + zfact(:, :) * SQRT(MAX(0._wp, rn2(:, :, jk))) * wmask(:, :, jk)
      END DO
    CASE (2)
      zfact(:, :) = 0._wp
      DO jk = 2, jpkm1
        zfact(:, :) = zfact(:, :) + e3w_n(:, :, jk) * MAX(0._wp, rn2(:, :, jk)) * wmask(:, :, jk)
      END DO
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (zfact(ji, jj) /= 0) zfact(ji, jj) = epyc_iwm(ji, jj) / (rau0 * zfact(ji, jj))
        END DO
      END DO
      DO jk = 2, jpkm1
        zemx_iwm(:, :, jk) = zemx_iwm(:, :, jk) + zfact(:, :) * MAX(0._wp, rn2(:, :, jk)) * wmask(:, :, jk)
      END DO
    END SELECT
    zwkb(:, :, :) = 0._wp
    zfact(:, :) = 0._wp
    DO jk = 2, jpkm1
      zfact(:, :) = zfact(:, :) + e3w_n(:, :, jk) * SQRT(MAX(0._wp, rn2(:, :, jk))) * wmask(:, :, jk)
      zwkb(:, :, jk) = zfact(:, :)
    END DO
    DO jk = 2, jpkm1
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (zfact(ji, jj) /= 0) zwkb(ji, jj, jk) = zhdep(ji, jj) * (zfact(ji, jj) - zwkb(ji, jj, jk)) * wmask(ji, jj, jk) / &
&zfact(ji, jj)
        END DO
      END DO
    END DO
    zwkb(:, :, 1) = zhdep(:, :) * wmask(:, :, 1)
    DO jk = 2, jpkm1
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (rn2(ji, jj, jk) <= 0._wp .OR. wmask(ji, jj, jk) == 0._wp) THEN
            zweight(ji, jj, jk) = 0._wp
          ELSE
            zweight(ji, jj, jk) = rn2(ji, jj, jk) * hbot_iwm(ji, jj) * (EXP(- zwkb(ji, jj, jk) / hbot_iwm(ji, jj)) - EXP(- &
&zwkb(ji, jj, jk - 1) / hbot_iwm(ji, jj)))
          END IF
        END DO
      END DO
    END DO
    zfact(:, :) = 0._wp
    DO jk = 2, jpkm1
      zfact(:, :) = zfact(:, :) + zweight(:, :, jk)
    END DO
    DO jj = 1, jpj
      DO ji = 1, jpi
        IF (zfact(ji, jj) /= 0) zfact(ji, jj) = ebot_iwm(ji, jj) / (rau0 * zfact(ji, jj))
      END DO
    END DO
    DO jk = 2, jpkm1
      zemx_iwm(:, :, jk) = zemx_iwm(:, :, jk) + zweight(:, :, jk) * zfact(:, :) * wmask(:, :, jk) / (gde3w_n(:, :, jk) - &
&gde3w_n(:, :, jk - 1))
    END DO
    znu_t(:, :, :) = 1.E-4_wp * (17.91_wp - 0.53810_wp * tsn(:, :, :, jp_tem) + 0.00694_wp * tsn(:, :, :, jp_tem) * tsn(:, :, :, &
&jp_tem) + 0.02305_wp * tsn(:, :, :, jp_sal)) * tmask(:, :, :) * r1_rau0
    DO jk = 2, jpkm1
      znu_w(:, :, jk) = 0.5_wp * (znu_t(:, :, jk - 1) + znu_t(:, :, jk)) * wmask(:, :, jk)
    END DO
    DO jk = 2, jpkm1
      zReb(:, :, jk) = zemx_iwm(:, :, jk) / MAX(1.E-20_wp, znu_w(:, :, jk) * rn2(:, :, jk))
    END DO
    DO jk = 2, jpkm1
      zav_wave(:, :, jk) = znu_w(:, :, jk) * zReb(:, :, jk) * r1_6
    END DO
    IF (ln_mevar) THEN
      DO jk = 2, jpkm1
        DO jj = 1, jpj
          DO ji = 1, jpi
            IF (zReb(ji, jj, jk) > 480.00_wp) THEN
              zav_wave(ji, jj, jk) = 3.6515_wp * znu_w(ji, jj, jk) * SQRT(zReb(ji, jj, jk))
            ELSE IF (zReb(ji, jj, jk) < 10.224_wp) THEN
              zav_wave(ji, jj, jk) = 0.052125_wp * znu_w(ji, jj, jk) * zReb(ji, jj, jk) * SQRT(zReb(ji, jj, jk))
            END IF
          END DO
        END DO
      END DO
    END IF
    DO jk = 2, jpkm1
      zav_wave(:, :, jk) = MIN(MAX(1.4E-7_wp, zav_wave(:, :, jk)), 1.E-2_wp) * wmask(:, :, jk)
    END DO
    IF (kt == nit000) THEN
      zztmp = 0._wp
      DO jk = 2, jpkm1
        DO jj = 1, jpj
          DO ji = 1, jpi
            zztmp = zztmp + e3w_n(ji, jj, jk) * e1e2t(ji, jj) * MAX(0._wp, rn2(ji, jj, jk)) * zav_wave(ji, jj, jk) * wmask(ji, jj, &
&jk) * tmask_i(ji, jj)
          END DO
        END DO
      END DO
      CALL mpp_sum('zdfiwm', zztmp)
      zztmp = rau0 * zztmp
      IF (lwp) THEN
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) 'zdf_iwm : Internal wave-driven mixing (iwm)'
        WRITE(numout, FMT = *) '~~~~~~~ '
        WRITE(numout, FMT = *)
        WRITE(numout, FMT = *) '      Total power consumption by av_wave =  ', zztmp * 1.E-12_wp, 'TW'
      END IF
    END IF
    IF (ln_tsdiff) THEN
      ztmp1 = 0.505_wp + 0.495_wp * TANH(0.92_wp * (LOG10(1.E-20_wp) - 0.60_wp))
      DO jk = 2, jpkm1
        DO jj = 1, jpj
          DO ji = 1, jpi
            ztmp2 = zReb(ji, jj, jk) * 5._wp * r1_6
            IF (ztmp2 > 1.E-20_wp .AND. wmask(ji, jj, jk) == 1._wp) THEN
              zav_ratio(ji, jj, jk) = 0.505_wp + 0.495_wp * TANH(0.92_wp * (LOG10(ztmp2) - 0.60_wp))
            ELSE
              zav_ratio(ji, jj, jk) = ztmp1 * wmask(ji, jj, jk)
            END IF
          END DO
        END DO
      END DO
      CALL iom_put("av_ratio", zav_ratio)
      DO jk = 2, jpkm1
        p_avs(:, :, jk) = p_avs(:, :, jk) + zav_wave(:, :, jk) * zav_ratio(:, :, jk)
        p_avt(:, :, jk) = p_avt(:, :, jk) + zav_wave(:, :, jk)
        p_avm(:, :, jk) = p_avm(:, :, jk) + zav_wave(:, :, jk)
      END DO
    ELSE
      DO jk = 2, jpkm1
        p_avs(:, :, jk) = p_avs(:, :, jk) + zav_wave(:, :, jk)
        p_avt(:, :, jk) = p_avt(:, :, jk) + zav_wave(:, :, jk)
        p_avm(:, :, jk) = p_avm(:, :, jk) + zav_wave(:, :, jk)
      END DO
    END IF
    CALL iom_put("av_wave", zav_wave)
    IF (iom_use("bflx_iwm") .OR. iom_use("pcmap_iwm")) THEN
      ALLOCATE(z2d(jpi, jpj), z3d(jpi, jpj, jpk))
      z3d(:, :, :) = MAX(0._wp, rn2(:, :, :)) * zav_wave(:, :, :)
      z2d(:, :) = 0._wp
      DO jk = 2, jpkm1
        z2d(:, :) = z2d(:, :) + e3w_n(:, :, jk) * z3d(:, :, jk) * wmask(:, :, jk)
      END DO
      z2d(:, :) = rau0 * z2d(:, :)
      CALL iom_put("bflx_iwm", z3d)
      CALL iom_put("pcmap_iwm", z2d)
      DEALLOCATE(z2d, z3d)
    END IF
    CALL iom_put("emix_iwm", zemx_iwm)
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = zav_wave, clinfo1 = ' iwm - av_wave: ', tab3d_2 = avt, clinfo2 = ' avt: ', kdim = jpk)
  END SUBROUTINE zdf_iwm
  SUBROUTINE zdf_iwm_init
    INTEGER :: ji, jj, jk
    INTEGER :: inum
    INTEGER :: ios
    REAL(KIND = wp) :: zbot, zpyc, zcri
    NAMELIST /namzdf_iwm/ nn_zpyc, ln_mevar, ln_tsdiff
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namzdf_iwm, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namzdf_iwm in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namzdf_iwm, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namzdf_iwm in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namzdf_iwm)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'zdf_iwm_init : internal wave-driven mixing'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namzdf_iwm : set wave-driven mixing parameters'
      WRITE(numout, FMT = *) '      Pycnocline-intensified diss. scales as N (=1) or N^2 (=2) = ', nn_zpyc
      WRITE(numout, FMT = *) '      Variable (T) or constant (F) mixing efficiency            = ', ln_mevar
      WRITE(numout, FMT = *) '      Differential internal wave-driven mixing (T) or not (F)   = ', ln_tsdiff
    END IF
    avmb(:) = 1.4E-6_wp
    avtb(:) = 1.E-10_wp
    avtb_2d(:, :) = 1.E0_wp
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) '   Force the background value applied to avm & avt in TKE to be everywhere ', 'the viscous molecular &
&value & a very small diffusive value, resp.'
    END IF
    IF (zdf_iwm_alloc() /= 0) CALL ctl_stop('STOP', 'zdf_iwm_init : unable to allocate iwm arrays')
    CALL iom_open('mixing_power_bot', inum)
    CALL iom_get(inum, jpdom_data, 'field', ebot_iwm, 1)
    CALL iom_close(inum)
    CALL iom_open('mixing_power_pyc', inum)
    CALL iom_get(inum, jpdom_data, 'field', epyc_iwm, 1)
    CALL iom_close(inum)
    CALL iom_open('mixing_power_cri', inum)
    CALL iom_get(inum, jpdom_data, 'field', ecri_iwm, 1)
    CALL iom_close(inum)
    CALL iom_open('decay_scale_bot', inum)
    CALL iom_get(inum, jpdom_data, 'field', hbot_iwm, 1)
    CALL iom_close(inum)
    CALL iom_open('decay_scale_cri', inum)
    CALL iom_get(inum, jpdom_data, 'field', hcri_iwm, 1)
    CALL iom_close(inum)
    ebot_iwm(:, :) = ebot_iwm(:, :) * ssmask(:, :)
    epyc_iwm(:, :) = epyc_iwm(:, :) * ssmask(:, :)
    ecri_iwm(:, :) = ecri_iwm(:, :) * ssmask(:, :)
    zbot = glob_sum('zdfiwm', e1e2t(:, :) * ebot_iwm(:, :))
    zpyc = glob_sum('zdfiwm', e1e2t(:, :) * epyc_iwm(:, :))
    zcri = glob_sum('zdfiwm', e1e2t(:, :) * ecri_iwm(:, :))
    IF (lwp) THEN
      WRITE(numout, FMT = *) '      High-mode wave-breaking energy:             ', zbot * 1.E-12_wp, 'TW'
      WRITE(numout, FMT = *) '      Pycnocline-intensifed wave-breaking energy: ', zpyc * 1.E-12_wp, 'TW'
      WRITE(numout, FMT = *) '      Critical slope wave-breaking energy:        ', zcri * 1.E-12_wp, 'TW'
    END IF
  END SUBROUTINE zdf_iwm_init
END MODULE zdfiwm