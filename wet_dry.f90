MODULE wet_dry
  USE oce
  USE dom_oce
  USE sbc_oce, ONLY: ln_rnf
  USE sbcrnf
  USE in_out_manager
  USE lbclnk
  USE lib_mpp
  USE timing
  IMPLICIT NONE
  PRIVATE
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: wdmask
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: wdramp, wdrampu, wdrampv
  LOGICAL, PUBLIC :: ln_wd_il
  LOGICAL, PUBLIC :: ln_wd_dl
  REAL(KIND = wp), PUBLIC :: rn_wdmin0
  REAL(KIND = wp), PUBLIC :: rn_wdmin1
  REAL(KIND = wp), PUBLIC :: r_rn_wdmin1
  REAL(KIND = wp), PUBLIC :: rn_wdmin2
  REAL(KIND = wp), PUBLIC :: rn_wd_sbcdep
  REAL(KIND = wp), PUBLIC :: rn_wd_sbcfra
  REAL(KIND = wp), PUBLIC :: rn_wdld
  INTEGER, PUBLIC :: nn_wdit
  LOGICAL, PUBLIC :: ln_wd_dl_bc
  LOGICAL, PUBLIC :: ln_wd_dl_rmp
  REAL(KIND = wp), PUBLIC :: ssh_ref
  LOGICAL, PUBLIC :: ll_wd
  PUBLIC :: wad_init
  PUBLIC :: wad_lmt
  PUBLIC :: wad_lmt_bt
  CONTAINS
  SUBROUTINE wad_init
    INTEGER :: ios, ierr
    NAMELIST /namwad/ ln_wd_il, ln_wd_dl, rn_wdmin0, rn_wdmin1, rn_wdmin2, rn_wdld, nn_wdit, ln_wd_dl_bc, ln_wd_dl_rmp, &
&rn_wd_sbcdep, rn_wd_sbcfra
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namwad, IOSTAT = ios, ERR = 905)
905 IF (ios /= 0) CALL ctl_nam(ios, 'namwad in reference namelist', .TRUE.)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namwad, IOSTAT = ios, ERR = 906)
906 IF (ios > 0) CALL ctl_nam(ios, 'namwad in configuration namelist', .TRUE.)
    IF (lwm) WRITE(numond, namwad)
    IF (rn_wd_sbcfra >= 1) CALL ctl_stop('STOP', 'rn_wd_sbcfra >=1 : must be < 1')
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'wad_init : Wetting and drying initialization through namelist read'
      WRITE(numout, FMT = *) '~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namwad'
      WRITE(numout, FMT = *) '      Logical for Iter Lim wd option   ln_wd_il     = ', ln_wd_il
      WRITE(numout, FMT = *) '      Logical for Dir. Lim wd option   ln_wd_dl     = ', ln_wd_dl
      WRITE(numout, FMT = *) '      Depth at which wet/drying starts rn_wdmin0    = ', rn_wdmin0
      WRITE(numout, FMT = *) '      Minimum wet depth on dried cells rn_wdmin1    = ', rn_wdmin1
      WRITE(numout, FMT = *) '      Tolerance of min wet depth       rn_wdmin2    = ', rn_wdmin2
      WRITE(numout, FMT = *) '      land elevation threshold         rn_wdld      = ', rn_wdld
      WRITE(numout, FMT = *) '      Max iteration for W/D limiter    nn_wdit      = ', nn_wdit
      WRITE(numout, FMT = *) '      T => baroclinic u,v=0 at dry pts: ln_wd_dl_bc = ', ln_wd_dl_bc
      WRITE(numout, FMT = *) '      use a ramp for rwd limiter:  ln_wd_dl_rwd_rmp = ', ln_wd_dl_rmp
      WRITE(numout, FMT = *) '      cut off depth sbc for wd   rn_wd_sbcdep       = ', rn_wd_sbcdep
      WRITE(numout, FMT = *) '      fraction to start sbc wgt rn_wd_sbcfra        = ', rn_wd_sbcfra
    END IF
    IF (.NOT. ln_read_cfg) THEN
      IF (lwp) WRITE(numout, FMT = *) '      No configuration file so seting ssh_ref to zero  '
      ssh_ref = 0._wp
    END IF
    r_rn_wdmin1 = 1 / rn_wdmin1
    ll_wd = .FALSE.
    IF (ln_wd_il .OR. ln_wd_dl) THEN
      ll_wd = .TRUE.
      ALLOCATE(wdmask(jpi, jpj), STAT = ierr)
      ALLOCATE(wdramp(jpi, jpj), wdrampu(jpi, jpj), wdrampv(jpi, jpj), STAT = ierr)
      IF (ierr /= 0) CALL ctl_stop('STOP', 'wad_init : Array allocation error')
    END IF
  END SUBROUTINE wad_init
  SUBROUTINE wad_lmt(sshb1, sshemp, z2dt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: sshb1
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN) :: sshemp
    REAL(KIND = wp), INTENT(IN) :: z2dt
    INTEGER :: ji, jj, jk, jk1
    INTEGER :: jflag
    REAL(KIND = wp) :: zcoef, zdep1, zdep2
    REAL(KIND = wp) :: zzflxp, zzflxn
    REAL(KIND = wp) :: zdepwd
    REAL(KIND = wp) :: ztmp
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zwdlmtu, zwdlmtv
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zflxp, zflxn
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zflxu, zflxv
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zflxu1, zflxv1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    IF (ln_timing) CALL timing_start('wad_lmt')
    DO jk = 1, jpkm1
      !$ACC KERNELS
      un(:, :, jk) = un(:, :, jk) * zwdlmtu(:, :)
      vn(:, :, jk) = vn(:, :, jk) * zwdlmtv(:, :)
      !$ACC END KERNELS
    END DO
    !$ACC KERNELS
    jflag = 0
    zdepwd = 50._wp
    zflxp(:, :) = 0._wp
    zflxn(:, :) = 0._wp
    zflxu(:, :) = 0._wp
    zflxv(:, :) = 0._wp
    zwdlmtu(:, :) = 1._wp
    zwdlmtv(:, :) = 1._wp
    !$ACC END KERNELS
    DO jk = 1, jpkm1
      !$ACC KERNELS
      zflxu(:, :) = zflxu(:, :) + e3u_n(:, :, jk) * un(:, :, jk) * umask(:, :, jk)
      zflxv(:, :) = zflxv(:, :) + e3v_n(:, :, jk) * vn(:, :, jk) * vmask(:, :, jk)
      !$ACC END KERNELS
    END DO
    !$ACC KERNELS
    zflxu(:, :) = zflxu(:, :) * e2u(:, :)
    zflxv(:, :) = zflxv(:, :) * e1v(:, :)
    wdmask(:, :) = 1._wp
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('wad_lmt', 'r0', 0, 0)
    DO jj = 2, jpj
      DO ji = 2, jpi
        IF (tmask(ji, jj, 1) < 0.5_wp) CYCLE
        IF (ht_0(ji, jj) - ssh_ref > zdepwd) CYCLE
        zflxp(ji, jj) = MAX(zflxu(ji, jj), 0._wp) - MIN(zflxu(ji - 1, jj), 0._wp) + MAX(zflxv(ji, jj), 0._wp) - MIN(zflxv(ji, jj - &
&1), 0._wp)
        zflxn(ji, jj) = MIN(zflxu(ji, jj), 0._wp) - MAX(zflxu(ji - 1, jj), 0._wp) + MIN(zflxv(ji, jj), 0._wp) - MAX(zflxv(ji, jj - &
&1), 0._wp)
        zdep2 = ht_0(ji, jj) + sshb1(ji, jj) - rn_wdmin1
        IF (zdep2 <= 0._wp) THEN
          sshb1(ji, jj) = rn_wdmin1 - ht_0(ji, jj)
          IF (zflxu(ji, jj) > 0._wp) zwdlmtu(ji, jj) = 0._wp
          IF (zflxu(ji - 1, jj) < 0._wp) zwdlmtu(ji - 1, jj) = 0._wp
          IF (zflxv(ji, jj) > 0._wp) zwdlmtv(ji, jj) = 0._wp
          IF (zflxv(ji, jj - 1) < 0._wp) zwdlmtv(ji, jj - 1) = 0._wp
          wdmask(ji, jj) = 0._wp
        END IF
      END DO
    END DO
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    wdramp(:, :) = MIN((ht_0(:, :) + sshb1(:, :) - rn_wdmin1) / (rn_wdmin0 - rn_wdmin1), 1.0_wp)
    !$ACC LOOP INDEPENDENT COLLAPSE(2)
    DO jj = 1, jpjm1
      DO ji = 1, jpim1
        wdrampu(ji, jj) = MIN(wdramp(ji, jj), wdramp(ji + 1, jj))
        wdrampv(ji, jj) = MIN(wdramp(ji, jj), wdramp(ji, jj + 1))
      END DO
    END DO
    !$ACC END KERNELS
    DO jk1 = 1, nn_wdit + 1
      !$ACC KERNELS
      zflxu1(:, :) = zflxu(:, :) * zwdlmtu(:, :)
      zflxv1(:, :) = zflxv(:, :) * zwdlmtv(:, :)
      jflag = 0
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('wad_lmt', 'r1', 0, 0)
      DO jj = 2, jpj
        DO ji = 2, jpi
          IF (tmask(ji, jj, 1) < 0.5_wp) CYCLE
          IF (ht_0(ji, jj) > zdepwd) CYCLE
          ztmp = e1e2t(ji, jj)
          zzflxp = MAX(zflxu1(ji, jj), 0._wp) - MIN(zflxu1(ji - 1, jj), 0._wp) + MAX(zflxv1(ji, jj), 0._wp) - MIN(zflxv1(ji, jj - &
&1), 0._wp)
          zzflxn = MIN(zflxu1(ji, jj), 0._wp) - MAX(zflxu1(ji - 1, jj), 0._wp) + MIN(zflxv1(ji, jj), 0._wp) - MAX(zflxv1(ji, jj - &
&1), 0._wp)
          zdep1 = (zzflxp + zzflxn) * z2dt / ztmp
          zdep2 = ht_0(ji, jj) + sshb1(ji, jj) - rn_wdmin1 - z2dt * sshemp(ji, jj)
          IF (zdep1 > zdep2) THEN
            wdmask(ji, jj) = 0._wp
            zcoef = ((zdep2 - rn_wdmin2) * ztmp - zzflxn * z2dt) / (zflxp(ji, jj) * z2dt)
            IF (zcoef > 0._wp) THEN
              jflag = 1
            ELSE
              zcoef = 0._wp
            END IF
            IF (jk1 > nn_wdit) zcoef = 0._wp
            IF (zflxu1(ji, jj) > 0._wp) zwdlmtu(ji, jj) = zcoef
            IF (zflxu1(ji - 1, jj) < 0._wp) zwdlmtu(ji - 1, jj) = zcoef
            IF (zflxv1(ji, jj) > 0._wp) zwdlmtv(ji, jj) = zcoef
            IF (zflxv1(ji, jj - 1) < 0._wp) zwdlmtv(ji, jj - 1) = zcoef
          END IF
        END DO
      END DO
      CALL lbc_lnk_multi('wet_dry', zwdlmtu, 'U', 1., zwdlmtv, 'V', 1.)
      CALL mpp_max('wet_dry', jflag)
      IF (jflag == 0) EXIT
      CALL profile_psy_data1 % PostEnd
    END DO
    DO jk = 1, jpkm1
      !$ACC KERNELS
      un(:, :, jk) = un(:, :, jk) * zwdlmtu(:, :)
      vn(:, :, jk) = vn(:, :, jk) * zwdlmtv(:, :)
      !$ACC END KERNELS
    END DO
    !$ACC KERNELS
    un_b(:, :) = un_b(:, :) * zwdlmtu(:, :)
    vn_b(:, :) = vn_b(:, :) * zwdlmtv(:, :)
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('wad_lmt', 'r2', 0, 0)
    CALL lbc_lnk_multi('wet_dry', un, 'U', - 1., vn, 'V', - 1.)
    CALL lbc_lnk_multi('wet_dry', un_b, 'U', - 1., vn_b, 'V', - 1.)
    IF (jflag == 1 .AND. lwp) WRITE(numout, FMT = *) 'Need more iterations in wad_lmt!!!'
    IF (ln_timing) CALL timing_stop('wad_lmt')
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE wad_lmt
  SUBROUTINE wad_lmt_bt(zflxu, zflxv, sshn_e, zssh_frc, rdtbt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(IN) :: rdtbt
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: zflxu, zflxv, sshn_e, zssh_frc
    INTEGER :: ji, jj, jk, jk1
    INTEGER :: jflag
    REAL(KIND = wp) :: z2dt
    REAL(KIND = wp) :: zcoef, zdep1, zdep2
    REAL(KIND = wp) :: zzflxp, zzflxn
    REAL(KIND = wp) :: zdepwd
    REAL(KIND = wp) :: ztmp
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zwdlmtu, zwdlmtv
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zflxp, zflxn
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zflxu1, zflxv1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    IF (ln_timing) CALL timing_start('wad_lmt_bt')
    !$ACC KERNELS
    jflag = 0
    zdepwd = 50._wp
    z2dt = rdtbt
    zflxp(:, :) = 0._wp
    zflxn(:, :) = 0._wp
    zwdlmtu(:, :) = 1._wp
    zwdlmtv(:, :) = 1._wp
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('wad_lmt_bt', 'r0', 0, 0)
    DO jj = 2, jpj
      DO ji = 2, jpi
        IF (tmask(ji, jj, 1) < 0.5_wp) CYCLE
        IF (ht_0(ji, jj) > zdepwd) CYCLE
        zflxp(ji, jj) = MAX(zflxu(ji, jj), 0._wp) - MIN(zflxu(ji - 1, jj), 0._wp) + MAX(zflxv(ji, jj), 0._wp) - MIN(zflxv(ji, jj - &
&1), 0._wp)
        zflxn(ji, jj) = MIN(zflxu(ji, jj), 0._wp) - MAX(zflxu(ji - 1, jj), 0._wp) + MIN(zflxv(ji, jj), 0._wp) - MAX(zflxv(ji, jj - &
&1), 0._wp)
        zdep2 = ht_0(ji, jj) + sshn_e(ji, jj) - rn_wdmin1
        IF (zdep2 <= 0._wp) THEN
          sshn_e(ji, jj) = rn_wdmin1 - ht_0(ji, jj)
          IF (zflxu(ji, jj) > 0._wp) zwdlmtu(ji, jj) = 0._wp
          IF (zflxu(ji - 1, jj) < 0._wp) zwdlmtu(ji - 1, jj) = 0._wp
          IF (zflxv(ji, jj) > 0._wp) zwdlmtv(ji, jj) = 0._wp
          IF (zflxv(ji, jj - 1) < 0._wp) zwdlmtv(ji, jj - 1) = 0._wp
        END IF
      END DO
    END DO
    CALL profile_psy_data0 % PostEnd
    DO jk1 = 1, nn_wdit + 1
      !$ACC KERNELS
      zflxu1(:, :) = zflxu(:, :) * zwdlmtu(:, :)
      zflxv1(:, :) = zflxv(:, :) * zwdlmtv(:, :)
      jflag = 0
      !$ACC END KERNELS
      CALL profile_psy_data1 % PreStart('wad_lmt_bt', 'r1', 0, 0)
      DO jj = 2, jpj
        DO ji = 2, jpi
          IF (tmask(ji, jj, 1) < 0.5_wp) CYCLE
          IF (ht_0(ji, jj) > zdepwd) CYCLE
          ztmp = e1e2t(ji, jj)
          zzflxp = MAX(zflxu1(ji, jj), 0._wp) - MIN(zflxu1(ji - 1, jj), 0._wp) + MAX(zflxv1(ji, jj), 0._wp) - MIN(zflxv1(ji, jj - &
&1), 0._wp)
          zzflxn = MIN(zflxu1(ji, jj), 0._wp) - MAX(zflxu1(ji - 1, jj), 0._wp) + MIN(zflxv1(ji, jj), 0._wp) - MAX(zflxv1(ji, jj - &
&1), 0._wp)
          zdep1 = (zzflxp + zzflxn) * z2dt / ztmp
          zdep2 = ht_0(ji, jj) + sshn_e(ji, jj) - rn_wdmin1 - z2dt * zssh_frc(ji, jj)
          IF (zdep1 > zdep2) THEN
            zcoef = ((zdep2 - rn_wdmin2) * ztmp - zzflxn * z2dt) / (zflxp(ji, jj) * z2dt)
            IF (zcoef > 0._wp) THEN
              jflag = 1
            ELSE
              zcoef = 0._wp
            END IF
            IF (jk1 > nn_wdit) zcoef = 0._wp
            IF (zflxu1(ji, jj) > 0._wp) zwdlmtu(ji, jj) = zcoef
            IF (zflxu1(ji - 1, jj) < 0._wp) zwdlmtu(ji - 1, jj) = zcoef
            IF (zflxv1(ji, jj) > 0._wp) zwdlmtv(ji, jj) = zcoef
            IF (zflxv1(ji, jj - 1) < 0._wp) zwdlmtv(ji, jj - 1) = zcoef
          END IF
        END DO
      END DO
      CALL lbc_lnk_multi('wet_dry', zwdlmtu, 'U', 1., zwdlmtv, 'V', 1.)
      CALL mpp_max('wet_dry', jflag)
      IF (jflag == 0) EXIT
      CALL profile_psy_data1 % PostEnd
    END DO
    !$ACC KERNELS
    zflxu(:, :) = zflxu(:, :) * zwdlmtu(:, :)
    zflxv(:, :) = zflxv(:, :) * zwdlmtv(:, :)
    !$ACC END KERNELS
    CALL profile_psy_data2 % PreStart('wad_lmt_bt', 'r2', 0, 0)
    CALL lbc_lnk_multi('wet_dry', zflxu, 'U', - 1., zflxv, 'V', - 1.)
    IF (jflag == 1 .AND. lwp) WRITE(numout, FMT = *) 'Need more iterations in wad_lmt_bt!!!'
    IF (ln_timing) CALL timing_stop('wad_lmt_bt')
    CALL profile_psy_data2 % PostEnd
  END SUBROUTINE wad_lmt_bt
END MODULE wet_dry