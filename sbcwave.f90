MODULE sbcwave
  USE phycst
  USE oce
  USE sbc_oce
  USE zdf_oce, ONLY: ln_zdfswm
  USE bdy_oce
  USE domvvl
  USE iom
  USE in_out_manager
  USE lib_mpp
  USE fldread
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sbc_stokes
  PUBLIC :: sbc_wstress
  PUBLIC :: sbc_wave
  PUBLIC :: sbc_wave_init
  LOGICAL, PUBLIC :: cpl_hsig = .FALSE.
  LOGICAL, PUBLIC :: cpl_phioc = .FALSE.
  LOGICAL, PUBLIC :: cpl_sdrftx = .FALSE.
  LOGICAL, PUBLIC :: cpl_sdrfty = .FALSE.
  LOGICAL, PUBLIC :: cpl_wper = .FALSE.
  LOGICAL, PUBLIC :: cpl_wfreq = .FALSE.
  LOGICAL, PUBLIC :: cpl_wnum = .FALSE.
  LOGICAL, PUBLIC :: cpl_tauwoc = .FALSE.
  LOGICAL, PUBLIC :: cpl_tauw = .FALSE.
  LOGICAL, PUBLIC :: cpl_wdrag = .FALSE.
  INTEGER :: jpfld
  INTEGER :: jp_usd
  INTEGER :: jp_vsd
  INTEGER :: jp_hsw
  INTEGER :: jp_wmp
  INTEGER :: jp_wfr
  TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_cd
  TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_sd
  TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_wn
  TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_tauwoc
  TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_tauw
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: cdn_wave
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: hsw, wmp, wnum
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: wfreq
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: tauoc_wave
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: tauw_x, tauw_y
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: tsd2d
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: div_sd
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: ut0sd, vt0sd
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :, :) :: usd, vsd, wsd
  CONTAINS
  SUBROUTINE sbc_stokes
    INTEGER :: jj, ji, jk
    INTEGER :: ik
    REAL(KIND = wp) :: ztransp, zfac, zsp0
    REAL(KIND = wp) :: zdepth, zsqrt_depth, zexp_depth, z_two_thirds, zsqrtpi
    REAL(KIND = wp) :: zbot_u, zbot_v, zkb_u, zkb_v, zke3_u, zke3_v, zda_u, zda_v
    REAL(KIND = wp) :: zstokes_psi_u_bot, zstokes_psi_v_bot
    REAL(KIND = wp) :: zdep_u, zdep_v, zkh_u, zkh_v
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: zk_t, zk_u, zk_v, zu0_sd, zv0_sd
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: zstokes_psi_u_top, zstokes_psi_v_top
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: ze3divh
    ALLOCATE(ze3divh(jpi, jpj, jpk))
    ALLOCATE(zk_t(jpi, jpj), zk_u(jpi, jpj), zk_v(jpi, jpj), zu0_sd(jpi, jpj), zv0_sd(jpi, jpj))
    IF (ll_st_bv_li) THEN
      zfac = 2.0_wp * rpi / 16.0_wp
      DO jj = 1, jpj
        DO ji = 1, jpi
          ztransp = zfac * hsw(ji, jj) * hsw(ji, jj) / MAX(wmp(ji, jj), 0.0000001_wp)
          tsd2d(ji, jj) = SQRT(ut0sd(ji, jj) * ut0sd(ji, jj) + vt0sd(ji, jj) * vt0sd(ji, jj))
          zk_t(ji, jj) = ABS(tsd2d(ji, jj)) / MAX(ABS(5.97_wp * ztransp), 0.0000001_wp)
        END DO
      END DO
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zk_u(ji, jj) = 0.5_wp * (zk_t(ji, jj) + zk_t(ji + 1, jj))
          zk_v(ji, jj) = 0.5_wp * (zk_t(ji, jj) + zk_t(ji, jj + 1))
          zu0_sd(ji, jj) = 0.5_wp * (ut0sd(ji, jj) + ut0sd(ji + 1, jj))
          zv0_sd(ji, jj) = 0.5_wp * (vt0sd(ji, jj) + vt0sd(ji, jj + 1))
        END DO
      END DO
    ELSE IF (ll_st_peakfr) THEN
      DO jj = 1, jpj
        DO ji = 1, jpi
          zk_t(ji, jj) = (2.0_wp * rpi * wfreq(ji, jj)) * (2.0_wp * rpi * wfreq(ji, jj)) / grav
        END DO
      END DO
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zk_u(ji, jj) = 0.5_wp * (zk_t(ji, jj) + zk_t(ji + 1, jj))
          zk_v(ji, jj) = 0.5_wp * (zk_t(ji, jj) + zk_t(ji, jj + 1))
          zu0_sd(ji, jj) = 0.5_wp * (ut0sd(ji, jj) + ut0sd(ji + 1, jj))
          zv0_sd(ji, jj) = 0.5_wp * (vt0sd(ji, jj) + vt0sd(ji, jj + 1))
        END DO
      END DO
    END IF
    IF (ll_st_bv2014) THEN
      !$OMP parallel default(shared), private(ji,jj,jk,zda_u,zda_v,zdep_u,zdep_v,zkh_u,zkh_v)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zdep_u = 0.5_wp * (gdept_n(ji, jj, jk) + gdept_n(ji + 1, jj, jk))
            zdep_v = 0.5_wp * (gdept_n(ji, jj, jk) + gdept_n(ji, jj + 1, jk))
            zkh_u = zk_u(ji, jj) * zdep_u
            zkh_v = zk_v(ji, jj) * zdep_v
            zda_u = EXP(- 2.0_wp * zkh_u) / (1.0_wp + 8.0_wp * zkh_u)
            zda_v = EXP(- 2.0_wp * zkh_v) / (1.0_wp + 8.0_wp * zkh_v)
            usd(ji, jj, jk) = zda_u * zu0_sd(ji, jj) * umask(ji, jj, jk)
            vsd(ji, jj, jk) = zda_v * zv0_sd(ji, jj) * vmask(ji, jj, jk)
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
    ELSE IF (ll_st_li2017 .OR. ll_st_peakfr) THEN
      ALLOCATE(zstokes_psi_u_top(jpi, jpj), zstokes_psi_v_top(jpi, jpj))
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zstokes_psi_u_top(ji, jj) = 0._wp
          zstokes_psi_v_top(ji, jj) = 0._wp
        END DO
      END DO
      zsqrtpi = SQRT(rpi)
      z_two_thirds = 2.0_wp / 3.0_wp
      !$OMP parallel default(shared),  &
!$omp& private(ji,jj,jk,zbot_u,zbot_v,zda_u,zda_v,zdepth,zexp_depth,zkb_u,zkb_v,zke3_u,zke3_v,zsqrt_depth,zstokes_psi_u_bot, &
!$omp& zstokes_psi_v_bot)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zbot_u = (gdepw_n(ji, jj, jk + 1) + gdepw_n(ji + 1, jj, jk + 1))
            zbot_v = (gdepw_n(ji, jj, jk + 1) + gdepw_n(ji, jj + 1, jk + 1))
            zkb_u = zk_u(ji, jj) * zbot_u
            zkb_v = zk_v(ji, jj) * zbot_v
            zke3_u = MAX(1.E-8_wp, 2.0_wp * zk_u(ji, jj) * e3u_n(ji, jj, jk))
            zke3_v = MAX(1.E-8_wp, 2.0_wp * zk_v(ji, jj) * e3v_n(ji, jj, jk))
            zdepth = zkb_u
            zsqrt_depth = SQRT(zdepth)
            zexp_depth = EXP(- zdepth)
            zstokes_psi_u_bot = 1.0_wp - zexp_depth - z_two_thirds * (zsqrtpi * zsqrt_depth * zdepth * ERFC(zsqrt_depth) + 1.0_wp &
&- (1.0_wp + zdepth) * zexp_depth)
            zda_u = (zstokes_psi_u_bot - zstokes_psi_u_top(ji, jj)) / zke3_u
            zstokes_psi_u_top(ji, jj) = zstokes_psi_u_bot
            zdepth = zkb_v
            zsqrt_depth = SQRT(zdepth)
            zexp_depth = EXP(- zdepth)
            zstokes_psi_v_bot = 1.0_wp - zexp_depth - z_two_thirds * (zsqrtpi * zsqrt_depth * zdepth * ERFC(zsqrt_depth) + 1.0_wp &
&- (1.0_wp + zdepth) * zexp_depth)
            zda_v = (zstokes_psi_v_bot - zstokes_psi_v_top(ji, jj)) / zke3_v
            zstokes_psi_v_top(ji, jj) = zstokes_psi_v_bot
            usd(ji, jj, jk) = zda_u * zu0_sd(ji, jj) * umask(ji, jj, jk)
            vsd(ji, jj, jk) = zda_v * zv0_sd(ji, jj) * vmask(ji, jj, jk)
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
      DEALLOCATE(zstokes_psi_u_top, zstokes_psi_v_top)
    END IF
    CALL lbc_lnk_multi('sbcwave', usd, 'U', - 1., vsd, 'V', - 1.)
    !$OMP parallel default(shared), private(ji,jj,jk)
    !$OMP do schedule(static)
    DO jk = 1, jpkm1
      DO jj = 2, jpj
        DO ji = 2, jpi
          ze3divh(ji, jj, jk) = (e2u(ji, jj) * e3u_n(ji, jj, jk) * usd(ji, jj, jk) - e2u(ji - 1, jj) * e3u_n(ji - 1, jj, jk) * &
&usd(ji - 1, jj, jk) + e1v(ji, jj) * e3v_n(ji, jj, jk) * vsd(ji, jj, jk) - e1v(ji, jj - 1) * e3v_n(ji, jj - 1, jk) * vsd(ji, jj - &
&1, jk)) * r1_e1e2t(ji, jj)
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
    CALL lbc_lnk('sbcwave', ze3divh, 'T', 1.)
    IF (ln_linssh) THEN
      ik = 1
    ELSE
      ik = 2
    END IF
    !$OMP parallel default(shared), private(jk)
    !$OMP do schedule(static)
    DO jk = jpkm1, ik, - 1
      wsd(:, :, jk) = wsd(:, :, jk + 1) - ze3divh(:, :, jk)
    END DO
    !$OMP end do
    !$OMP end parallel
    IF (ln_bdy) THEN
      !$OMP parallel default(shared), private(jk)
      !$OMP do schedule(static)
      DO jk = 1, jpkm1
        wsd(:, :, jk) = wsd(:, :, jk) * bdytmask(:, :)
      END DO
      !$OMP end do
      !$OMP end parallel
    END IF
    div_sd(:, :) = 0._wp
    !$OMP parallel default(shared), private(jk)
    !$OMP do schedule(static)
    DO jk = 1, jpkm1
      div_sd(:, :) = div_sd(:, :) + ze3divh(:, :, jk)
    END DO
    !$OMP end do
    !$OMP end parallel
    CALL iom_put("ustokes", usd)
    CALL iom_put("vstokes", vsd)
    CALL iom_put("wstokes", wsd)
    DEALLOCATE(ze3divh)
    DEALLOCATE(zk_t, zk_u, zk_v, zu0_sd, zv0_sd)
  END SUBROUTINE sbc_stokes
  SUBROUTINE sbc_wstress
    INTEGER :: jj, ji
    IF (ln_tauwoc) THEN
      utau(:, :) = utau(:, :) * tauoc_wave(:, :)
      vtau(:, :) = vtau(:, :) * tauoc_wave(:, :)
      taum(:, :) = taum(:, :) * tauoc_wave(:, :)
    END IF
    IF (ln_tauw) THEN
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          utau(ji, jj) = 0.5_wp * (tauw_x(ji, jj) + tauw_x(ji + 1, jj))
          vtau(ji, jj) = 0.5_wp * (tauw_y(ji, jj) + tauw_y(ji, jj + 1))
          taum(ji, jj) = SQRT(tauw_x(ji, jj) * tauw_x(ji, jj) + tauw_y(ji, jj) * tauw_y(ji, jj))
        END DO
      END DO
      CALL lbc_lnk_multi('sbcwave', utau(:, :), 'U', - 1., vtau(:, :), 'V', - 1., taum(:, :), 'T', - 1.)
    END IF
  END SUBROUTINE sbc_wstress
  SUBROUTINE sbc_wave(kt)
    INTEGER, INTENT(IN) :: kt
    IF (ln_cdgw .AND. .NOT. cpl_wdrag) THEN
      CALL fld_read(kt, nn_fsbc, sf_cd)
      cdn_wave(:, :) = sf_cd(1) % fnow(:, :, 1) * tmask(:, :, 1)
    END IF
    IF (ln_tauwoc .AND. .NOT. cpl_tauwoc) THEN
      CALL fld_read(kt, nn_fsbc, sf_tauwoc)
      tauoc_wave(:, :) = sf_tauwoc(1) % fnow(:, :, 1) * tmask(:, :, 1)
    END IF
    IF (ln_tauw .AND. .NOT. cpl_tauw) THEN
      CALL fld_read(kt, nn_fsbc, sf_tauw)
      tauw_x(:, :) = sf_tauw(1) % fnow(:, :, 1) * tmask(:, :, 1)
      tauw_y(:, :) = sf_tauw(2) % fnow(:, :, 1) * tmask(:, :, 1)
    END IF
    IF (ln_sdw) THEN
      IF (jpfld > 0) THEN
        CALL fld_read(kt, nn_fsbc, sf_sd)
        IF (jp_hsw > 0) hsw(:, :) = sf_sd(jp_hsw) % fnow(:, :, 1) * tmask(:, :, 1)
        IF (jp_wmp > 0) wmp(:, :) = sf_sd(jp_wmp) % fnow(:, :, 1) * tmask(:, :, 1)
        IF (jp_wfr > 0) wfreq(:, :) = sf_sd(jp_wfr) % fnow(:, :, 1) * tmask(:, :, 1)
        IF (jp_usd > 0) ut0sd(:, :) = sf_sd(jp_usd) % fnow(:, :, 1) * tmask(:, :, 1)
        IF (jp_vsd > 0) vt0sd(:, :) = sf_sd(jp_vsd) % fnow(:, :, 1) * tmask(:, :, 1)
      END IF
      IF (ln_zdfswm .AND. .NOT. cpl_wnum) THEN
        CALL fld_read(kt, nn_fsbc, sf_wn)
        wnum(:, :) = sf_wn(1) % fnow(:, :, 1) * tmask(:, :, 1)
      END IF
      IF ((ll_st_bv_li .AND. jp_hsw > 0 .AND. jp_wmp > 0 .AND. jp_usd > 0 .AND. jp_vsd > 0) .OR. (ll_st_peakfr .AND. jp_wfr > 0 &
&.AND. jp_usd > 0 .AND. jp_vsd > 0)) CALL sbc_stokes
    END IF
  END SUBROUTINE sbc_wave
  SUBROUTINE sbc_wave_init
    INTEGER :: ierror, ios
    INTEGER :: ifpr
    CHARACTER(LEN = 100) :: cn_dir
    TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) :: slf_i, slf_j
    TYPE(FLD_N) :: sn_cdg, sn_usd, sn_vsd, sn_hsw, sn_wmp, sn_wfr, sn_wnum, sn_tauwoc, sn_tauwx, sn_tauwy
    NAMELIST /namsbc_wave/ sn_cdg, cn_dir, sn_usd, sn_vsd, sn_hsw, sn_wmp, sn_wfr, sn_wnum, sn_tauwoc, sn_tauwx, sn_tauwy
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namsbc_wave, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namsbc_wave in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namsbc_wave, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namsbc_wave in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namsbc_wave)
    IF (ln_cdgw) THEN
      IF (.NOT. cpl_wdrag) THEN
        ALLOCATE(sf_cd(1), STAT = ierror)
        IF (ierror > 0) CALL ctl_stop('STOP', 'sbc_wave_init: unable to allocate sf_wave structure')
        ALLOCATE(sf_cd(1) % fnow(jpi, jpj, 1))
        IF (sn_cdg % ln_tint) ALLOCATE(sf_cd(1) % fdta(jpi, jpj, 1, 2))
        CALL fld_fill(sf_cd, (/sn_cdg/), cn_dir, 'sbc_wave_init', 'Wave module ', 'namsbc_wave')
      END IF
      ALLOCATE(cdn_wave(jpi, jpj))
    END IF
    IF (ln_tauwoc) THEN
      IF (.NOT. cpl_tauwoc) THEN
        ALLOCATE(sf_tauwoc(1), STAT = ierror)
        IF (ierror > 0) CALL ctl_stop('STOP', 'sbc_wave_init: unable to allocate sf_wave structure')
        ALLOCATE(sf_tauwoc(1) % fnow(jpi, jpj, 1))
        IF (sn_tauwoc % ln_tint) ALLOCATE(sf_tauwoc(1) % fdta(jpi, jpj, 1, 2))
        CALL fld_fill(sf_tauwoc, (/sn_tauwoc/), cn_dir, 'sbc_wave_init', 'Wave module', 'namsbc_wave')
      END IF
      ALLOCATE(tauoc_wave(jpi, jpj))
    END IF
    IF (ln_tauw) THEN
      IF (.NOT. cpl_tauw) THEN
        ALLOCATE(sf_tauw(2), STAT = ierror)
        IF (ierror > 0) CALL ctl_stop('STOP', 'sbc_wave_init: unable to allocate sf_tauw structure')
        ALLOCATE(slf_j(2))
        slf_j(1) = sn_tauwx
        slf_j(2) = sn_tauwy
        ALLOCATE(sf_tauw(1) % fnow(jpi, jpj, 1))
        ALLOCATE(sf_tauw(2) % fnow(jpi, jpj, 1))
        IF (slf_j(1) % ln_tint) ALLOCATE(sf_tauw(1) % fdta(jpi, jpj, 1, 2))
        IF (slf_j(2) % ln_tint) ALLOCATE(sf_tauw(2) % fdta(jpi, jpj, 1, 2))
        CALL fld_fill(sf_tauw, (/slf_j/), cn_dir, 'sbc_wave_init', 'read wave input', 'namsbc_wave')
      END IF
      ALLOCATE(tauw_x(jpi, jpj))
      ALLOCATE(tauw_y(jpi, jpj))
    END IF
    IF (ln_sdw) THEN
      jpfld = 0
      jp_usd = 0
      jp_vsd = 0
      jp_hsw = 0
      jp_wmp = 0
      jp_wfr = 0
      IF (.NOT. cpl_sdrftx) THEN
        jpfld = jpfld + 1
        jp_usd = jpfld
      END IF
      IF (.NOT. cpl_sdrfty) THEN
        jpfld = jpfld + 1
        jp_vsd = jpfld
      END IF
      IF (.NOT. cpl_hsig .AND. ll_st_bv_li) THEN
        jpfld = jpfld + 1
        jp_hsw = jpfld
      END IF
      IF (.NOT. cpl_wper .AND. ll_st_bv_li) THEN
        jpfld = jpfld + 1
        jp_wmp = jpfld
      END IF
      IF (.NOT. cpl_wfreq .AND. ll_st_peakfr) THEN
        jpfld = jpfld + 1
        jp_wfr = jpfld
      END IF
      IF (jpfld > 0) THEN
        ALLOCATE(slf_i(jpfld))
        IF (jp_usd > 0) slf_i(jp_usd) = sn_usd
        IF (jp_vsd > 0) slf_i(jp_vsd) = sn_vsd
        IF (jp_hsw > 0) slf_i(jp_hsw) = sn_hsw
        IF (jp_wmp > 0) slf_i(jp_wmp) = sn_wmp
        IF (jp_wfr > 0) slf_i(jp_wfr) = sn_wfr
        ALLOCATE(sf_sd(jpfld), STAT = ierror)
        IF (ierror > 0) CALL ctl_stop('STOP', 'sbc_wave_init: unable to allocate sf_wave structure')
        DO ifpr = 1, jpfld
          ALLOCATE(sf_sd(ifpr) % fnow(jpi, jpj, 1))
          IF (slf_i(ifpr) % ln_tint) ALLOCATE(sf_sd(ifpr) % fdta(jpi, jpj, 1, 2))
        END DO
        CALL fld_fill(sf_sd, slf_i, cn_dir, 'sbc_wave_init', 'Wave module ', 'namsbc_wave')
      END IF
      ALLOCATE(usd(jpi, jpj, jpk), vsd(jpi, jpj, jpk), wsd(jpi, jpj, jpk))
      ALLOCATE(hsw(jpi, jpj), wmp(jpi, jpj))
      ALLOCATE(wfreq(jpi, jpj))
      ALLOCATE(ut0sd(jpi, jpj), vt0sd(jpi, jpj))
      ALLOCATE(div_sd(jpi, jpj))
      ALLOCATE(tsd2d(jpi, jpj))
      ut0sd(:, :) = 0._wp
      vt0sd(:, :) = 0._wp
      hsw(:, :) = 0._wp
      wmp(:, :) = 0._wp
      usd(:, :, :) = 0._wp
      vsd(:, :, :) = 0._wp
      wsd(:, :, :) = 0._wp
      IF (.NOT. cpl_wnum) THEN
        ALLOCATE(sf_wn(1), STAT = ierror)
        IF (ierror > 0) CALL ctl_stop('STOP', 'sbc_wave_init: unable toallocate sf_wave structure')
        ALLOCATE(sf_wn(1) % fnow(jpi, jpj, 1))
        IF (sn_wnum % ln_tint) ALLOCATE(sf_wn(1) % fdta(jpi, jpj, 1, 2))
        CALL fld_fill(sf_wn, (/sn_wnum/), cn_dir, 'sbc_wave', 'Wave module', 'namsbc_wave')
      END IF
      ALLOCATE(wnum(jpi, jpj))
    END IF
  END SUBROUTINE sbc_wave_init
END MODULE sbcwave