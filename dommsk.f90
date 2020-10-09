MODULE dommsk
  USE oce
  USE dom_oce
  USE usrdef_fmask
  USE bdy_oce
  USE in_out_manager
  USE iom
  USE lbclnk
  USE lib_mpp
  USE iom
  USE fldread
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dom_msk
  REAL(KIND = wp) :: rn_shlat
  LOGICAL, PUBLIC :: ln_vorlat
  CONTAINS
  SUBROUTINE dom_msk(k_top, k_bot)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, DIMENSION(:, :), INTENT(IN) :: k_top, k_bot
    INTEGER :: ji, jj, jk
    INTEGER :: iif, iil
    INTEGER :: ijf, ijl
    INTEGER :: iktop, ikbot
    INTEGER :: ios, inum
    REAL(KIND = wp) :: zshlat
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zshlat2d
    LOGICAL :: ln_shlat2d
    CHARACTER(LEN = 256) :: cn_shlat2d_file, cn_shlat2d_var
    NAMELIST /namlbc/ rn_shlat, ln_vorlat, ln_shlat2d, cn_shlat2d_file, cn_shlat2d_var
    NAMELIST /nambdy/ ln_bdy, nb_bdy, ln_coords_file, cn_coords_file, ln_mask_file, cn_mask_file, cn_dyn2d, nn_dyn2d_dta, &
&cn_dyn3d, nn_dyn3d_dta, cn_tra, nn_tra_dta, ln_tra_dmp, ln_dyn3d_dmp, rn_time_dmp, rn_time_dmp_out, cn_ice, nn_ice_dta, &
&rn_ice_tem, rn_ice_sal, rn_ice_age, ln_vol, nn_volctl, nn_rimwidth, nb_jpk_bdy
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('dom_msk', 'r0', 0, 0)
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namlbc, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namlbc in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namlbc, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namlbc in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namlbc)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dommsk : ocean mask '
      WRITE(numout, FMT = *) '~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namlbc'
      WRITE(numout, FMT = *) '      lateral momentum boundary cond.    rn_shlat  = ', rn_shlat
      WRITE(numout, FMT = *) '      consistency with analytical form   ln_vorlat = ', ln_vorlat
    END IF
    IF (lwp) WRITE(numout, FMT = *)
    IF (ln_shlat2d) THEN
      IF (lwp) WRITE(numout, FMT = *) '         READ shlat as a 2D coefficient in a file '
      ALLOCATE(zshlat2d(jpi, jpj))
      CALL iom_open(TRIM(cn_shlat2d_file), inum)
      CALL iom_get(inum, jpdom_data, TRIM(cn_shlat2d_var), zshlat2d, 1)
      CALL iom_close(inum)
    ELSE
      IF (rn_shlat == 0.) THEN
        IF (lwp) WRITE(numout, FMT = *) '   ==>>>   ocean lateral  free-slip'
      ELSE IF (rn_shlat == 2.) THEN
        IF (lwp) WRITE(numout, FMT = *) '   ==>>>   ocean lateral  no-slip'
      ELSE IF (0. < rn_shlat .AND. rn_shlat < 2.) THEN
        IF (lwp) WRITE(numout, FMT = *) '   ==>>>   ocean lateral  partial-slip'
      ELSE IF (2. < rn_shlat) THEN
        IF (lwp) WRITE(numout, FMT = *) '   ==>>>   ocean lateral  strong-slip'
      ELSE
        CALL ctl_stop('dom_msk: wrong value for rn_shlat (i.e. a negalive value). We stop.')
      END IF
    END IF
    tmask(:, :, :) = 0._wp
    DO jj = 1, jpj
      DO ji = 1, jpi
        iktop = k_top(ji, jj)
        ikbot = k_bot(ji, jj)
        IF (iktop /= 0) THEN
          tmask(ji, jj, iktop : ikbot) = 1._wp
        END IF
      END DO
    END DO
    CALL lbc_lnk('dommsk', tmask, 'T', 1._wp)
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, nambdy, IOSTAT = ios, ERR = 903)
903 IF (ios /= 0) CALL ctl_nam(ios, 'nambdy in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, nambdy, IOSTAT = ios, ERR = 904)
904 IF (ios > 0) CALL ctl_nam(ios, 'nambdy in configuration namelist', lwp)
    IF (ln_bdy .AND. ln_mask_file) THEN
      CALL iom_open(cn_mask_file, inum)
      CALL iom_get(inum, jpdom_data, 'bdy_msk', bdytmask(:, :))
      CALL iom_close(inum)
      DO jk = 1, jpkm1
        DO jj = 1, jpj
          DO ji = 1, jpi
            tmask(ji, jj, jk) = tmask(ji, jj, jk) * bdytmask(ji, jj)
          END DO
        END DO
      END DO
    END IF
    DO jk = 1, jpk
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          umask(ji, jj, jk) = tmask(ji, jj, jk) * tmask(ji + 1, jj, jk)
          vmask(ji, jj, jk) = tmask(ji, jj, jk) * tmask(ji, jj + 1, jk)
        END DO
        DO ji = 1, jpim1
          fmask(ji, jj, jk) = tmask(ji, jj, jk) * tmask(ji + 1, jj, jk) * tmask(ji, jj + 1, jk) * tmask(ji + 1, jj + 1, jk)
        END DO
      END DO
    END DO
    CALL lbc_lnk_multi('dommsk', umask, 'U', 1., vmask, 'V', 1., fmask, 'F', 1.)
    wmask(:, :, 1) = tmask(:, :, 1)
    wumask(:, :, 1) = umask(:, :, 1)
    wvmask(:, :, 1) = vmask(:, :, 1)
    DO jk = 2, jpk
      wmask(:, :, jk) = tmask(:, :, jk) * tmask(:, :, jk - 1)
      wumask(:, :, jk) = umask(:, :, jk) * umask(:, :, jk - 1)
      wvmask(:, :, jk) = vmask(:, :, jk) * vmask(:, :, jk - 1)
    END DO
    ssmask(:, :) = MAXVAL(tmask(:, :, :), DIM = 3)
    ssumask(:, :) = MAXVAL(umask(:, :, :), DIM = 3)
    ssvmask(:, :) = MAXVAL(vmask(:, :, :), DIM = 3)
    iif = nn_hls
    iil = nlci - nn_hls + 1
    ijf = nn_hls
    ijl = nlcj - nn_hls + 1
    tmask_h(:, :) = 1._wp
    tmask_h(1 : iif, :) = 0._wp
    tmask_h(iil : jpi, :) = 0._wp
    tmask_h(:, 1 : ijf) = 0._wp
    tmask_h(:, ijl : jpj) = 0._wp
    tpol(1 : jpiglo) = 1._wp
    fpol(1 : jpiglo) = 1._wp
    IF (jperio == 3 .OR. jperio == 4) THEN
      tpol(jpiglo / 2 + 1 : jpiglo) = 0._wp
      fpol(1 : jpiglo) = 0._wp
      IF (mjg(nlej) == jpjglo) THEN
        DO ji = iif + 1, iil - 1
          tmask_h(ji, nlej - 1) = tmask_h(ji, nlej - 1) * tpol(mig(ji))
        END DO
      END IF
    END IF
    IF (jperio == 5 .OR. jperio == 6) THEN
      tpol(1 : jpiglo) = 0._wp
      fpol(jpiglo / 2 + 1 : jpiglo) = 0._wp
    END IF
    tmask_i(:, :) = ssmask(:, :) * tmask_h(:, :)
    IF (rn_shlat /= 0 .OR. ln_shlat2d) THEN
      DO jk = 1, jpk
        IF (ln_shlat2d) THEN
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              IF (fmask(ji, jj, jk) == 0._wp) THEN
                fmask(ji, jj, jk) = zshlat2d(ji, jj) * MIN(1._wp, MAX(umask(ji, jj, jk), umask(ji, jj + 1, jk), vmask(ji, jj, jk), &
&vmask(ji + 1, jj, jk)))
              END IF
            END DO
          END DO
        ELSE
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              IF (fmask(ji, jj, jk) == 0._wp) THEN
                fmask(ji, jj, jk) = rn_shlat * MIN(1._wp, MAX(umask(ji, jj, jk), umask(ji, jj + 1, jk), vmask(ji, jj, jk), &
&vmask(ji + 1, jj, jk)))
              END IF
            END DO
          END DO
        END IF
        DO jj = 2, jpjm1
          IF (fmask(1, jj, jk) == 0._wp) THEN
            fmask(1, jj, jk) = rn_shlat * MIN(1._wp, MAX(vmask(2, jj, jk), umask(1, jj + 1, jk), umask(1, jj, jk)))
          END IF
          IF (fmask(jpi, jj, jk) == 0._wp) THEN
            fmask(jpi, jj, jk) = rn_shlat * MIN(1._wp, MAX(umask(jpi, jj + 1, jk), vmask(jpim1, jj, jk), umask(jpi, jj - 1, jk)))
          END IF
        END DO
        DO ji = 2, jpim1
          IF (fmask(ji, 1, jk) == 0._wp) THEN
            fmask(ji, 1, jk) = rn_shlat * MIN(1._wp, MAX(vmask(ji + 1, 1, jk), umask(ji, 2, jk), vmask(ji, 1, jk)))
          END IF
          IF (fmask(ji, jpj, jk) == 0._wp) THEN
            fmask(ji, jpj, jk) = rn_shlat * MIN(1._wp, MAX(vmask(ji + 1, jpj, jk), vmask(ji - 1, jpj, jk), umask(ji, jpjm1, jk)))
          END IF
        END DO
      END DO
      IF (ln_shlat2d) DEALLOCATE(zshlat2d)
      CALL lbc_lnk('dommsk', fmask, 'F', 1._wp)
    END IF
    CALL usr_def_fmask(cn_cfg, nn_cfg, fmask)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE dom_msk
END MODULE dommsk