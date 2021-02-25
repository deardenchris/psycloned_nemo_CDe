MODULE traqsr
  USE oce
  USE phycst
  USE dom_oce
  USE sbc_oce
  USE trc_oce
  USE trd_oce
  USE trdtra
  USE in_out_manager
  USE prtctl
  USE iom
  USE fldread
  USE restart
  USE lib_mpp
  USE lbclnk
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tra_qsr
  PUBLIC :: tra_qsr_init
  LOGICAL, PUBLIC :: ln_traqsr
  LOGICAL, PUBLIC :: ln_qsr_rgb
  LOGICAL, PUBLIC :: ln_qsr_2bd
  LOGICAL, PUBLIC :: ln_qsr_bio
  INTEGER, PUBLIC :: nn_chldta
  REAL(KIND = wp), PUBLIC :: rn_abs
  REAL(KIND = wp), PUBLIC :: rn_si0
  REAL(KIND = wp), PUBLIC :: rn_si1
  INTEGER, PUBLIC :: nksr
  INTEGER, PARAMETER :: np_RGB = 1
  INTEGER, PARAMETER :: np_RGBc = 2
  INTEGER, PARAMETER :: np_2BD = 3
  INTEGER, PARAMETER :: np_BIO = 4
  INTEGER :: nqsr
  REAL(KIND = wp) :: xsi0r
  REAL(KIND = wp) :: xsi1r
  REAL(KIND = wp), DIMENSION(3, 61) :: rkrgb
  TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_chl
  CONTAINS
  SUBROUTINE tra_qsr(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    INTEGER :: irgb
    REAL(KIND = wp) :: zchl, zcoef, z1_2
    REAL(KIND = wp) :: zc0, zc1, zc2, zc3
    REAL(KIND = wp) :: zzc0, zzc1, zzc2, zzc3
    REAL(KIND = wp) :: zz0, zz1
    REAL(KIND = wp) :: zCb, zCmax, zze, zpsi, zpsimax, zdelpsi, zCtot, zCze
    REAL(KIND = wp) :: zlogc, zlogc2, zlogc3
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zekb, zekg, zekr
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: ze0, ze1, ze2, ze3, zea, ztrdt
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: zetot, zchl3d
    IF (ln_timing) CALL timing_start('tra_qsr')
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'tra_qsr : penetration of the surface solar radiation'
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~'
    END IF
    IF (l_trdtra) THEN
      ALLOCATE(ztrdt(jpi, jpj, jpk))
      ztrdt(:, :, :) = tsa(:, :, :, jp_tem)
    END IF
    IF (kt == nit000) THEN
      IF (ln_rstart .AND. iom_varid(numror, 'qsr_hc_b', ldstop = .FALSE.) > 0) THEN
        IF (lwp) WRITE(numout, FMT = *) '          nit000-1 qsr tracer content forcing field read in the restart file'
        z1_2 = 0.5_wp
        CALL iom_get(numror, jpdom_autoglo, 'qsr_hc_b', qsr_hc_b, ldxios = lrxios)
      ELSE
        z1_2 = 1._wp
        qsr_hc_b(:, :, :) = 0._wp
      END IF
    ELSE
      z1_2 = 0.5_wp
      qsr_hc_b(:, :, :) = qsr_hc(:, :, :)
    END IF
    SELECT CASE (nqsr)
    CASE (np_BIO)
      !$OMP parallel default(shared), private(jk)
      !$OMP do schedule(static)
      DO jk = 1, nksr
        qsr_hc(:, :, jk) = r1_rau0_rcp * (etot3(:, :, jk) - etot3(:, :, jk + 1))
      END DO
      !$OMP end do
      !$OMP end parallel
    CASE (np_RGB, np_RGBc)
      ALLOCATE(zekb(jpi, jpj), zekg(jpi, jpj), zekr(jpi, jpj), ze0(jpi, jpj, jpk), ze1(jpi, jpj, jpk), ze2(jpi, jpj, jpk), &
&ze3(jpi, jpj, jpk), zea(jpi, jpj, jpk), zchl3d(jpi, jpj, jpk))
      IF (nqsr == np_RGBc) THEN
        CALL fld_read(kt, 1, sf_chl)
        DO jk = 1, nksr + 1
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              zchl = sf_chl(1) % fnow(ji, jj, 1)
              zCtot = 40.6 * zchl ** 0.459
              zze = 568.2 * zCtot ** (- 0.746)
              IF (zze > 102.) zze = 200.0 * zCtot ** (- 0.293)
              zpsi = gdepw_n(ji, jj, jk) / zze
              zlogc = LOG(zchl)
              zlogc2 = zlogc * zlogc
              zlogc3 = zlogc * zlogc * zlogc
              zCb = 0.768 + 0.087 * zlogc - 0.179 * zlogc2 - 0.025 * zlogc3
              zCmax = 0.299 - 0.289 * zlogc + 0.579 * zlogc2
              zpsimax = 0.6 - 0.640 * zlogc + 0.021 * zlogc2 + 0.115 * zlogc3
              zdelpsi = 0.710 + 0.159 * zlogc + 0.021 * zlogc2
              zCze = 1.12 * (zchl) ** 0.803
              zchl3d(ji, jj, jk) = zCze * (zCb + zCmax * EXP(- ((zpsi - zpsimax) / zdelpsi) ** 2))
            END DO
          END DO
        END DO
      ELSE
        !$OMP parallel default(shared), private(jk)
        !$OMP do schedule(static)
        DO jk = 1, nksr + 1
          zchl3d(:, :, jk) = 0.05
        END DO
        !$OMP end do
        !$OMP end parallel
      END IF
      zcoef = (1. - rn_abs) / 3._wp
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ze0(ji, jj, 1) = rn_abs * qsr(ji, jj)
          ze1(ji, jj, 1) = zcoef * qsr(ji, jj)
          ze2(ji, jj, 1) = zcoef * qsr(ji, jj)
          ze3(ji, jj, 1) = zcoef * qsr(ji, jj)
          zea(ji, jj, 1) = qsr(ji, jj)
        END DO
      END DO
      ! !$OMP parallel default(shared), private(irgb,ji,jj,jk,zc0,zc1,zc2,zc3,zchl)
      ! CDe - race condition on ze0, ze1 etc?
      ! !$OMP do schedule(static)
      DO jk = 2, nksr + 1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zchl = MIN(10., MAX(0.03, zchl3d(ji, jj, jk)))
            irgb = NINT(41 + 20. * LOG10(zchl) + 1.E-15)
            zekb(ji, jj) = rkrgb(1, irgb)
            zekg(ji, jj) = rkrgb(2, irgb)
            zekr(ji, jj) = rkrgb(3, irgb)
          END DO
        END DO
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zc0 = ze0(ji, jj, jk - 1) * EXP(- e3t_n(ji, jj, jk - 1) * xsi0r)
            zc1 = ze1(ji, jj, jk - 1) * EXP(- e3t_n(ji, jj, jk - 1) * zekb(ji, jj))
            zc2 = ze2(ji, jj, jk - 1) * EXP(- e3t_n(ji, jj, jk - 1) * zekg(ji, jj))
            zc3 = ze3(ji, jj, jk - 1) * EXP(- e3t_n(ji, jj, jk - 1) * zekr(ji, jj))
            ze0(ji, jj, jk) = zc0
            ze1(ji, jj, jk) = zc1
            ze2(ji, jj, jk) = zc2
            ze3(ji, jj, jk) = zc3
            zea(ji, jj, jk) = (zc0 + zc1 + zc2 + zc3) * wmask(ji, jj, jk)
          END DO
        END DO
      END DO
      ! !$OMP end do
      ! !$OMP barrier ! CDe
      ! !$OMP do schedule(static)
      DO jk = 1, nksr
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            qsr_hc(ji, jj, jk) = r1_rau0_rcp * (zea(ji, jj, jk) - zea(ji, jj, jk + 1))
          END DO
        END DO
      END DO
      ! !$OMP end do
      ! !$OMP end parallel
      DEALLOCATE(zekb, zekg, zekr, ze0, ze1, ze2, ze3, zea, zchl3d)
    CASE (np_2BD)
      zz0 = rn_abs * r1_rau0_rcp
      zz1 = (1. - rn_abs) * r1_rau0_rcp
      !$OMP parallel default(shared), private(ji,jj,jk,zc0,zc1)
      !$OMP do schedule(static)
      DO jk = 1, nksr
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zc0 = zz0 * EXP(- gdepw_n(ji, jj, jk) * xsi0r) + zz1 * EXP(- gdepw_n(ji, jj, jk) * xsi1r)
            zc1 = zz0 * EXP(- gdepw_n(ji, jj, jk + 1) * xsi0r) + zz1 * EXP(- gdepw_n(ji, jj, jk + 1) * xsi1r)
            qsr_hc(ji, jj, jk) = qsr(ji, jj) * (zc0 * wmask(ji, jj, jk) - zc1 * wmask(ji, jj, jk + 1))
          END DO
        END DO
      END DO
      !$OMP end do
      !$OMP end parallel
    END SELECT
    !$OMP parallel default(shared), private(ji,jj,jk)
    !$OMP do schedule(static)
    DO jk = 1, nksr
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          tsa(ji, jj, jk, jp_tem) = tsa(ji, jj, jk, jp_tem) + z1_2 * (qsr_hc_b(ji, jj, jk) + qsr_hc(ji, jj, jk)) / e3t_n(ji, jj, jk)
        END DO
      END DO
    END DO
    !$OMP end do
    !$OMP end parallel
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        IF (qsr(ji, jj) /= 0._wp) THEN
          fraqsr_1lev(ji, jj) = qsr_hc(ji, jj, 1) / (r1_rau0_rcp * qsr(ji, jj))
        ELSE
          fraqsr_1lev(ji, jj) = 1._wp
        END IF
      END DO
    END DO
    CALL lbc_lnk('traqsr', fraqsr_1lev(:, :), 'T', 1._wp)
    IF (iom_use('qsr3d')) THEN
      ALLOCATE(zetot(jpi, jpj, jpk))
      zetot(:, :, nksr + 1 : jpk) = 0._wp
      ! !$OMP parallel default(shared), private(jk) ! CDe race condition?
      ! !$OMP do schedule(static)
      DO jk = nksr, 1, - 1
        zetot(:, :, jk) = zetot(:, :, jk + 1) + qsr_hc(:, :, jk) * rau0_rcp
      END DO
      ! !$OMP end do
      ! !$OMP end parallel
      CALL iom_put('qsr3d', zetot)
      DEALLOCATE(zetot)
    END IF
    IF (lrst_oce) THEN
      IF (lwxios) CALL iom_swap(cwxios_context)
      CALL iom_rstput(kt, nitrst, numrow, 'qsr_hc_b', qsr_hc, ldxios = lwxios)
      CALL iom_rstput(kt, nitrst, numrow, 'fraqsr_1lev', fraqsr_1lev, ldxios = lwxios)
      IF (lwxios) CALL iom_swap(cxios_context)
    END IF
    IF (l_trdtra) THEN
      ztrdt(:, :, :) = tsa(:, :, :, jp_tem) - ztrdt(:, :, :)
      CALL trd_tra(kt, 'TRA', jp_tem, jptra_qsr, ztrdt)
      DEALLOCATE(ztrdt)
    END IF
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = tsa(:, :, :, jp_tem), clinfo1 = ' qsr  - Ta: ', mask1 = tmask, clinfo3 = 'tra-ta')
    IF (ln_timing) CALL timing_stop('tra_qsr')
  END SUBROUTINE tra_qsr
  SUBROUTINE tra_qsr_init
    INTEGER :: ji, jj, jk
    INTEGER :: ios, irgb, ierror, ioptio
    REAL(KIND = wp) :: zz0, zc0, zc1, zcoef
    REAL(KIND = wp) :: zz1, zc2, zc3, zchl
    CHARACTER(LEN = 100) :: cn_dir
    TYPE(FLD_N) :: sn_chl
    NAMELIST /namtra_qsr/ sn_chl, cn_dir, ln_qsr_rgb, ln_qsr_2bd, ln_qsr_bio, nn_chldta, rn_abs, rn_si0, rn_si1
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namtra_qsr, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namtra_qsr in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namtra_qsr, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namtra_qsr in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namtra_qsr)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'tra_qsr_init : penetration of the surface solar radiation'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namtra_qsr : set the parameter of penetration'
      WRITE(numout, FMT = *) '      RGB (Red-Green-Blue) light penetration       ln_qsr_rgb = ', ln_qsr_rgb
      WRITE(numout, FMT = *) '      2 band               light penetration       ln_qsr_2bd = ', ln_qsr_2bd
      WRITE(numout, FMT = *) '      bio-model            light penetration       ln_qsr_bio = ', ln_qsr_bio
      WRITE(numout, FMT = *) '      RGB : Chl data (=1) or cst value (=0)        nn_chldta  = ', nn_chldta
      WRITE(numout, FMT = *) '      RGB & 2 bands: fraction of light (rn_si1)    rn_abs     = ', rn_abs
      WRITE(numout, FMT = *) '      RGB & 2 bands: shortess depth of extinction  rn_si0     = ', rn_si0
      WRITE(numout, FMT = *) '      2 bands: longest depth of extinction         rn_si1     = ', rn_si1
      WRITE(numout, FMT = *)
    END IF
    ioptio = 0
    IF (ln_qsr_rgb) ioptio = ioptio + 1
    IF (ln_qsr_2bd) ioptio = ioptio + 1
    IF (ln_qsr_bio) ioptio = ioptio + 1
    IF (ioptio /= 1) CALL ctl_stop('Choose ONE type of light penetration in namelist namtra_qsr', ' 2 bands, 3 RGB bands or &
&bio-model light penetration')
    IF (ln_qsr_rgb .AND. nn_chldta == 0) nqsr = np_RGB
    IF (ln_qsr_rgb .AND. nn_chldta == 1) nqsr = np_RGBc
    IF (ln_qsr_2bd) nqsr = np_2BD
    IF (ln_qsr_bio) nqsr = np_BIO
    xsi0r = 1._wp / rn_si0
    xsi1r = 1._wp / rn_si1
    SELECT CASE (nqsr)
    CASE (np_RGB, np_RGBc)
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   R-G-B   light penetration '
      CALL trc_oce_rgb(rkrgb)
      nksr = trc_oce_ext_lev(r_si2, 33._wp)
      IF (lwp) WRITE(numout, FMT = *) '        level of light extinction = ', nksr, ' ref depth = ', gdepw_1d(nksr + 1), ' m'
      IF (nqsr == np_RGBc) THEN
        IF (lwp) WRITE(numout, FMT = *) '   ==>>>   Chlorophyll read in a file'
        ALLOCATE(sf_chl(1), STAT = ierror)
        IF (ierror > 0) THEN
          CALL ctl_stop('tra_qsr_init: unable to allocate sf_chl structure')
          RETURN
        END IF
        ALLOCATE(sf_chl(1) % fnow(jpi, jpj, 1))
        IF (sn_chl % ln_tint) ALLOCATE(sf_chl(1) % fdta(jpi, jpj, 1, 2))
        CALL fld_fill(sf_chl, (/sn_chl/), cn_dir, 'tra_qsr_init', 'Solar penetration function of read chlorophyll', 'namtra_qsr', &
&no_print)
      END IF
      IF (nqsr == np_RGB) THEN
        IF (lwp) WRITE(numout, FMT = *) '   ==>>>   Constant Chlorophyll concentration = 0.05'
      END IF
    CASE (np_2BD)
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   2 bands light penetration'
      nksr = trc_oce_ext_lev(rn_si1, 100._wp)
      IF (lwp) WRITE(numout, FMT = *) '        level of light extinction = ', nksr, ' ref depth = ', gdepw_1d(nksr + 1), ' m'
    CASE (np_BIO)
      IF (lwp) WRITE(numout, FMT = *) '   ==>>>   bio-model light penetration'
      IF (.NOT. lk_top) CALL ctl_stop('No bio model : ln_qsr_bio = true impossible ')
    END SELECT
    qsr_hc(:, :, :) = 0._wp
    IF (iom_varid(numror, 'fraqsr_1lev', ldstop = .FALSE.) > 0) THEN
      CALL iom_get(numror, jpdom_autoglo, 'fraqsr_1lev', fraqsr_1lev, ldxios = lrxios)
    ELSE
      fraqsr_1lev(:, :) = 1._wp
    END IF
    IF (lwxios) THEN
      CALL iom_set_rstw_var_active('qsr_hc_b')
      CALL iom_set_rstw_var_active('fraqsr_1lev')
    END IF
  END SUBROUTINE tra_qsr_init
END MODULE traqsr
